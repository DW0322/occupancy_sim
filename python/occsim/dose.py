"""Building synthetic stiffness maps from a beam shape and a crosslink model.
Optics and math translated over using claude as a TEMPORARY stand in not final
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .geometry import _compartment_centres

VERSION_LINE = "# occupancy-sim stiffness map v1"

BEAMS = ("bessel", "gaussian")
ENVELOPES = ("inverse-r", "flat")

# The scaffold's own stiffness range, kPa (Li 2025 soft-hydrogel regime). The
# floor is the un-illuminated substrate, so a map never reaches zero.
DEFAULT_E_FLOOR_KPA = 0.1
DEFAULT_E_MAX_KPA = 1.5
# Crosslink rate per unit accumulated dose. Chosen so a fully lit compartment
# is part-way up the saturating curve by DIV 25 rather than pinned at E_max,
# which is what leaves any contrast for a coupling to act on.
DEFAULT_KC_PER_DAY = 0.02
DEFAULT_RING_PERIOD_UM = 200.0
DEFAULT_WAIST_UM = 400.0


def _j0(x):
    ax = np.abs(np.asarray(x, dtype=float))
    out = np.empty_like(ax)

    # |x| <= 3: direct polynomial in (x/3)^2
    near = ax <= 3.0
    y = (ax[near] / 3.0) ** 2
    out[near] = (1.0 - 2.2499997 * y + 1.2656208 * y**2 - 0.3163866 * y**3
                 + 0.0444479 * y**4 - 0.0039444 * y**5 + 0.0002100 * y**6)

    # |x| > 3: asymptotic amplitude/phase form
    far = ~near
    z = 3.0 / ax[far]
    f0 = (0.79788456 - 0.00000077 * z - 0.00552740 * z**2 - 0.00009512 * z**3
          + 0.00137237 * z**4 - 0.00072805 * z**5 + 0.00014476 * z**6)
    t0 = (ax[far] - 0.78539816 - 0.04166397 * z - 0.00003954 * z**2
          + 0.00262573 * z**3 - 0.00054125 * z**4 - 0.00029333 * z**5
          + 0.00013558 * z**6)
    out[far] = f0 * np.cos(t0) / np.sqrt(ax[far])
    return out


def intensity(beam: str, grid: int, field_um: float, *,
              ring_period_um: float = DEFAULT_RING_PERIOD_UM,
              waist_um: float = DEFAULT_WAIST_UM,
              envelope: str = "inverse-r"):
    if beam not in BEAMS:
        raise ValueError(f"unknown beam {beam!r}, expected one of {BEAMS}")
    if envelope not in ENVELOPES:
        raise ValueError(
            f"unknown envelope {envelope!r}, expected one of {ENVELOPES}")
    if grid < 1:
        raise ValueError(f"grid must be >= 1, got {grid}")
    if field_um <= 0:
        raise ValueError(f"field_um must be > 0, got {field_um}")

    # same compartment convention as the geometry masks, so a map and a mask
    # built at the same grid describe the same squares
    centres = np.asarray(_compartment_centres(grid, field_um))
    r = np.hypot(centres[None, :], centres[:, None])

    if beam == "gaussian":
        if waist_um <= 0:
            raise ValueError(f"waist_um must be > 0, got {waist_um}")
        return np.exp(-2.0 * (r / waist_um) ** 2)

    if ring_period_um <= 0:
        raise ValueError(
            f"ring_period_um must be > 0, got {ring_period_um}")
    # J0's zeros are asymptotically pi apart, so a radial wavenumber of
    # pi/period puts consecutive rings one period apart in um.
    k = math.pi / ring_period_um
    prof = _j0(k * r) ** 2

    if envelope == "flat":
        # J0(x)^2 decays as 2/(pi x). Multiplying by the reciprocal makes every
        # ring carry comparable power; the clamp at 1 leaves the central lobe
        # alone instead of dividing by an envelope that does not apply there.
        prof = prof * np.maximum(1.0, math.pi * k * r / 2.0)
    peak = prof.max()
    return prof / peak if peak > 0 else prof


def build(beam: str, grid: int, field_um: float, divs, *,
          ring_period_um: float = DEFAULT_RING_PERIOD_UM,
          waist_um: float = DEFAULT_WAIST_UM,
          envelope: str = "inverse-r",
          e_floor_kpa: float = DEFAULT_E_FLOOR_KPA,
          e_max_kpa: float = DEFAULT_E_MAX_KPA,
          kc_per_day: float = DEFAULT_KC_PER_DAY):
    """Build one stiffness map per DIV from a beam and the crosslink model.

    Args:
        beam: One of BEAMS.
        grid: Compartments per side.
        field_um: Field side, um.
        divs: Ascending DIVs, evenly spaced; the cadence becomes `div_step`.
        ring_period_um: Ring spacing for `bessel`, um.
        waist_um: 1/e^2 intensity radius for `gaussian`, um.
        envelope: One of ENVELOPES, `bessel` only.
        e_floor_kpa: Stiffness of un-illuminated substrate, kPa.
        e_max_kpa: Stiffness at full crosslinking, kPa.
        kc_per_day: Crosslink rate per unit accumulated dose.

    Returns:
        A list of (grid, grid) arrays, one per DIV, in kPa.

    Raises:
        ValueError: On a bad beam, envelope, scale, DIV list or stiffness range.
    """
    divs = [float(d) for d in divs]
    if not divs:
        raise ValueError("at least one DIV is required")
    if any(d < 0 for d in divs):
        raise ValueError("DIVs must be >= 0")
    if any(b <= a for a, b in zip(divs, divs[1:])):
        raise ValueError("DIVs must be strictly ascending")
    if not e_max_kpa > e_floor_kpa >= 0:
        raise ValueError(
            f"need e_max > e_floor >= 0, got {e_max_kpa} and {e_floor_kpa}")
    if kc_per_day <= 0:
        raise ValueError(f"kc_per_day must be > 0, got {kc_per_day}")

    prof = intensity(beam, grid, field_um, ring_period_um=ring_period_um,
                     waist_um=waist_um, envelope=envelope)
    # Dose accumulates from DIV 0, so the map at DIV d reflects the whole
    # illumination history rather than an instantaneous intensity. That is what
    # makes the maps change with DIV at all.
    return [e_floor_kpa
            + (e_max_kpa - e_floor_kpa) * (1.0 - np.exp(-kc_per_day * prof * d))
            for d in divs]


def div_step_of(divs) -> int:
    """The declared cadence for a DIV list, as the v1 format requires.

    Args:
        divs: Ascending DIVs.

    Returns:
        The whole-day spacing between consecutive DIVs, or the single DIV's own
        value when there is only one.

    Raises:
        ValueError: If the spacing is not a constant whole number of days, which
            the C++ loader rejects as a gap in the sequence.
    """
    divs = [float(d) for d in divs]
    if len(divs) == 1:
        step = divs[0] if divs[0] >= 1 else 1.0
    else:
        gaps = {round(b - a, 9) for a, b in zip(divs, divs[1:])}
        if len(gaps) != 1:
            raise ValueError(
                f"DIVs must be evenly spaced for a v1 div_step, got gaps "
                f"{sorted(gaps)}")
        step = gaps.pop()
    if step < 1 or abs(step - round(step)) > 1e-9:
        raise ValueError(
            f"div_step must be a whole number of days >= 1, got {step}")
    return int(round(step))


def top_quartile_share(m) -> float:
    #Share of a map's value range sitting above its 75th percentile.

    m = np.asarray(m, dtype=float)
    lo, hi = float(m.min()), float(m.max())
    if hi <= lo:
        return 0.0
    return (hi - float(np.percentile(m, 75))) / (hi - lo)


def write(path: Path, maps, divs, grid: int, field_um: float,
          provenance: str, *, units: str = "kPa") -> Path:
    #Write per-DIV maps as a stiffness CSV the engine's loader accepts.
    
    step = div_step_of(divs)
    lines = [VERSION_LINE,
             f"# {provenance}",
             "# SYNTHETIC - GENERATED by occsim, not measured data.",
             "# E = E_floor + (E_max-E_floor)*alpha, alpha = 1-exp(-k_c*D),",
             "# D the dose accumulated from DIV 0. A fixture for exercising the",
             "# interface, not a photochemical model; the real maps come from",
             "# the optical half of the collaboration.",
             f"grid,{grid},{grid}",
             f"scaffold_um,{field_um:g}",
             f"units,{units}",
             f"div_step,{step}",
             "div,cx,cy,stiffness"]
    for d, m in zip(divs, maps):
        for cy in range(grid):
            for cx in range(grid):
                lines.append(f"{d:g},{cx},{cy},{m[cy][cx]:.6f}")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(lines) + "\n")
    return Path(path)
