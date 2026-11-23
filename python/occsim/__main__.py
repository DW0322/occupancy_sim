"""occsim: one entry point for everything that drives the engine.

    python -m occsim run            one invocation, prints occupancy
    python -m occsim sweep          DIV sweep, seed averaged, absolute and relative
    python -m occsim field          per-cell arbors over the territory fill
    python -m occsim compartments   top-down per-compartment occupancy over DIV
    python -m occsim density        occupancy vs seeding density
    python -m occsim stiffness      carried map + evidence the coupling steers
    python -m occsim regime         which cue dominates, per compartment
    python -m occsim make-geometry  write a scaffold mask from shape + area
    python -m occsim make-stiffness write a synthetic stiffness map from a beam
    python -m occsim all            the whole figure set, one invocation

Every subcommand is a thin wrapper over occsim.engine; nothing here rebuilds
the engine's argument list.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from .engine import (ROOT, DEFAULT_ENGINE, DEFAULT_MORPH, MORPH_DIR, POOL,
                         Run, measure, pool, rel_to_root)
    from . import dose, figures, geometry
except ImportError as exc:  # pragma: no cover, guidance, not logic
    if "relative import" not in str(exc):
        raise
    parent = Path(__file__).resolve().parents[1]
    sys.exit(
        "occsim is a package, not a standalone script, so running this file "
        "directly\nleaves Python without the package context its imports "
        "need.\n\nRun it as a module from the directory ABOVE the package:\n"
        f"    cd {parent}\n"
        "    python -m occsim run --n 15 --div 15\n\n"
        "Or install it once and call it from anywhere:\n"
        f"    pip install -e {parent}\n"
        "    occsim run --n 15 --div 15\n")


# Cell counts for the density curve, in one place so `density` and `all` can't
# drift apart. Stops at 120 because past that it's pinned at 100 %.
DENSITY_COUNTS = [1, 2, 3, 5, 8, 15, 30, 60, 120]
DENSITY_COUNTS_QUICK = [3, 8, 15, 40, 100]

# DIVs the compartments and field figures draw panels for. A CSV still gets
# written per swept DIV, so the data keeps the sweep's resolution while the
# figure stays readable.
PANEL_DIVS = [5, 10, 15, 20, 25]

# --quick trades resolution for a 1-2 minute pass. These four DIVs are the ones
# that carry the curve's shape: the DIV5 baseline everything is normalised to,
# the steep DIV8 rise, DIV15 where the reveal schedule saturates, and DIV25 as
# the culture-lifespan endpoint. Dropping the rest changes no conclusion, only
# how smooth the line looks.
QUICK_DIVS = [5, 8, 15, 25]

# The coupling figures were measured at 12 seeds on a 30x30 grid. Change either
# and the figure stops matching the numbers quoted for it.
COUPLING_SEEDS = 12
COUPLING_COMPARTMENTS = 30

# --stiff-weight values for the gradient-coupling figure. Both signs, since the
# sign is the direction, and 0 as the uncoupled reference.
STIFF_WEIGHTS = [-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9]
STIFF_WEIGHTS_QUICK = [-0.9, -0.45, 0.0, 0.45, 0.9]

# --stiff-speed values for the magnitude-coupling figure. 0..1 only, the
# direction lives in the slope sign, not the weight, and 0 is the reference.
SPEED_WEIGHTS = [0.0, 0.25, 0.5, 0.75, 1.0]
SPEED_WEIGHTS_QUICK = [0.0, 0.5, 1.0]

# The synthetic maps shipped in data/. They demonstrate the format, they aren't
# data, the figures use them to show each mechanism works and label themselves
# as such. A real map with the same grid, units and cadence drops straight in.
#
# One per coupling, because the two fail on opposite map shapes: a smooth
# single-gradient map is what --stiff-weight needs and leaves --stiff-speed
# with nothing to average over, while the ring field is the reverse. Running
# both panels on one map would show one working and one flat for reasons that
# are about the map rather than the mechanism.
DEFAULT_STIFFNESS_MAP = ROOT / "data" / "stiffness_example_synthetic.csv"
DEFAULT_BESSEL_MAP = ROOT / "data" / "stiffness_bessel_synthetic.csv"

# Scaffold area `occsim all` defines when the caller names none. It matches the
# shipped example lattice so the committed results and the worked example in
# docs/RUNNING.md describe the same scaffold. It is a MODELLING INPUT, not a
# measurement: the collaboration supplies the real area, and every output that
# uses it records the number.
DEFAULT_SCAFFOLD_AREA_FRAC = 0.4

# Optical-cue weight for the regime panel. The classification depends on the
# optical cue's DIRECTION, not its magnitude (magnitudes get reported beside the
# map, never folded into it), so all this constant really fixes is the sign
# convention. Positive is soft-ward, the model's default. Flipping it swaps the
# reinforcing and opposing halves of the competitive class. It matches
# --seek-weight's default so the panel reads "if both cues were equally strong".
# An explicit --stiff-weight wins over it, which is why it is applied as the
# regime subcommand's own default rather than read at the call site: 0.0 then
# means the caller asked for 0.0, not that they said nothing.
REGIME_STIFF_WEIGHT = 0.3


def _common(ap: argparse.ArgumentParser) -> None:
    """Add the flags every subcommand shares.

    Args:
        ap: The subcommand's parser, mutated in place.
    """
    ap.add_argument("--morph", type=Path, action="append", default=None,
                    help="morphology; repeat to sample from a pool")
    ap.add_argument("--pool", action="store_true",
                    help="39 distinct cells, one file per parent (clone-free)")
    ap.add_argument("--morph-dir", type=Path, default=None,
                    help="use every morphology in this directory")
    ap.add_argument("--engine", type=Path, default=DEFAULT_ENGINE)
    ap.add_argument("--n", type=int, default=15,
                    help="cells to seed. Not a free choice: occupancy "
                         "saturates above about 95 %% (N~15 at DIV15 with "
                         "migration) and stops discriminating, while the "
                         "clustering readout needs density to mean anything. "
                         "See README.md, Limitations")
    # --field is the clearer name: this is the whole analysis area, while
    # --scaffold-geometry is the scaffold within it. --scaffold stays as an
    # alias so logged invocations still reproduce.
    ap.add_argument("--field", "--scaffold", dest="field", type=float,
                    default=1500.0,
                    help="overall field side, um (--scaffold is a legacy "
                         "alias)")
    ap.add_argument("--territory", type=float, default=17.2)
    ap.add_argument("--res", type=float, default=0.5)
    ap.add_argument("--compartments", type=int, default=30,
                    help="compartments per side of the field; 30 = 50 um "
                         "squares. Any scaffold geometry or stiffness map "
                         "must be on this same grid. 150 (10 um) resolves the "
                         "48 um t-ZnO arms and ring-scale dose fields at "
                         "about 2 %% more runtime")
    # The scaffold area WITHIN the field. Either supply a mask, or state the
    # area and let occsim build one. The area is the stated input, so it gets a
    # flag instead of being a file the caller has to write first.
    ap.add_argument("--scaffold-area-frac", type=float, default=None,
                    help="scaffold area as a fraction of the field, 0-1. "
                         "Builds the geometry mask for the run and makes "
                         "scaffold-restricted occupancy the reported number")
    ap.add_argument("--scaffold-shape", choices=list(geometry.SHAPES),
                    default="lattice",
                    help="shape the scaffold area is arranged in; lattice "
                         "reproduces the printed strand grid (Li 2025), disc "
                         "and square are area-matched controls")
    ap.add_argument("--strand-pitch", type=float,
                    default=geometry.DEFAULT_STRAND_PITCH_UM,
                    help="strand centre-to-centre spacing for "
                         "--scaffold-shape lattice, um")
    ap.add_argument("--scaffold-geometry", type=Path, default=None,
                    help="a ready-made scaffold mask "
                         "(docs/SCAFFOLD_GEOMETRY_FORMAT.md), instead of "
                         "--scaffold-area-frac. Carried and displayed; "
                         "coupled to migration only via --plate-scaffold or "
                         "--seek scaffold")
    ap.add_argument("--stiffness-map", type=Path, default=None,
                    help="stiffness map (docs/STIFFNESS_FORMAT.md); carried "
                         "and displayed; coupled only via --stiff-weight or "
                         "--stiff-speed")
    ap.add_argument("--plate-scaffold", action="store_true",
                    help="COUPLING: seed on and confine to the scaffold "
                         "(needs --scaffold-geometry)")
    ap.add_argument("--stiff-weight", type=float, default=0.0,
                    help="COUPLING: stiffness-gradient bias -1..1, +ve toward "
                         "soft, -ve toward stiff (needs --stiffness-map)")
    ap.add_argument("--stiff-speed", type=float, default=0.0,
                    help="COUPLING: stiffness-magnitude speed modulation 0-1 "
                         "(needs --stiffness-map); the mechanism for "
                         "Bessel/ring dose fields")
    ap.add_argument("--speed-centre", type=float, default=None,
                    help="speed-coupling sigmoid centre, kPa "
                         "(default: auto from the maps)")
    ap.add_argument("--speed-slope", type=float, default=None,
                    help="speed-coupling sigmoid slope, 1/kPa; -ve = slower "
                         "on stiff (default: auto from the maps)")
    ap.add_argument("--seed-frac", type=float, default=1.0,
                    help="fraction of the scaffold side cells are seeded "
                         "across; 1.0 = whole field, ~0.1 = a central cluster")
    ap.add_argument("--mig-speed", type=float, default=0.0)
    ap.add_argument("--mig-persist", type=float, default=50.0)
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--occ-mode", choices=["cum", "inst"], default="cum")
    ap.add_argument("--seek", choices=["none", "scaffold"], default="none",
                    help="decision rule: seek the scaffold or keep the normal "
                         "walk (scaffold needs --scaffold-geometry)")
    ap.add_argument("--seek-weight", type=float, default=0.3)
    ap.add_argument("--p-seek", type=float, default=1.0)
    ap.add_argument("--cohere-weight", type=float, default=0.0,
                    help="COUPLING: 0-1, how hard cells steer toward the "
                         "neighbours their neurites reach. The radius is "
                         "derived from the reveal schedule, not set. Inert at "
                         "sparse N - check with "
                         "'occupancy_engine --cohere-check'")
    ap.add_argument("--cohere-gate", choices=["none", "scaffold"],
                    default="none",
                    help="COUPLING: with 'scaffold', only cells off the "
                         "scaffold cohere - substrate and cell-cell adhesion "
                         "compete (needs a scaffold geometry)")
    ap.add_argument("--out", type=Path, default=ROOT / "results",
                    help="output DIRECTORY, created if missing; each "
                         "subcommand derives its own figure filenames inside "
                         "it, so a filename passed here becomes a directory of "
                         "that name. To name a file, use the subcommand's "
                         "--csv where it has one")


def _resolve_geometry(a) -> Path | None:
    """The scaffold mask this run uses, building one if only an area was given.

    Two ways to say what the scaffold is and they must not both be used: a
    ready-made mask names the compartments directly, an area fraction has
    occsim lay them out. Accepting both would leave which scaffold the run
    actually used ambiguous in the output.

    Args:
        a: The parsed argparse.Namespace; `a.out` must already exist.

    Returns:
        Path to a geometry CSV, or None when the run defines no scaffold.

    Raises:
        SystemExit: If both ways are given, or the area fraction is unusable.
    """
    if a.scaffold_geometry and a.scaffold_area_frac is not None:
        raise SystemExit(
            "--scaffold-geometry and --scaffold-area-frac were both given. "
            "The first supplies a mask and the second builds one, so together "
            "they leave it ambiguous which scaffold the run used. Pass one.")
    if a.scaffold_geometry:
        return a.scaffold_geometry
    if a.scaffold_area_frac is None:
        return None
    try:
        path = geometry.ensure(a.out, a.scaffold_shape, a.compartments,
                               a.field, a.scaffold_area_frac, a.strand_pitch)
    except ValueError as exc:
        raise SystemExit(f"--scaffold-area-frac: {exc}")
    print(f"scaffold: {a.scaffold_shape}, "
          f"{100 * a.scaffold_area_frac:.1f} % of a {a.field:g} um field, "
          f"{a.compartments}x{a.compartments} compartments -> "
          f"{rel_to_root(path)}")
    return path


def _run_from(a) -> Run:
    """Build the base Run from parsed CLI arguments.

    Works out which morphology set to use (--morph, --morph-dir, --pool,
    or the single default) and passes the shared flags through.

    Args:
        a: The parsed argparse.Namespace, with `scaffold_geometry` already
            resolved by _resolve_geometry().

    Returns:
        A Run, ready for per-subcommand fields (div, seed, ...) to be set.

    Raises:
        SystemExit: If the requested morphology directory or pool is empty.
    """
    if getattr(a, "morph_dir", None):
        morph = pool(a.morph_dir)
        if not morph:
            raise SystemExit(f"no morphologies in {a.morph_dir}")
    elif a.pool:
        morph = POOL
        if not morph:
            raise SystemExit(f"no morphologies in {MORPH_DIR}")
    else:
        morph = a.morph if a.morph else DEFAULT_MORPH
    if isinstance(morph, list) and len(morph) == 1:
        morph = morph[0]
    return Run(morph=morph, n=a.n, scaffold=a.field, territory=a.territory,
               res=a.res, seed_frac=a.seed_frac,
               compartments=a.compartments,
               scaffold_geometry=a.scaffold_geometry,
               stiffness_map=a.stiffness_map,
               plate_scaffold=a.plate_scaffold,
               stiff_weight=a.stiff_weight, stiff_speed=a.stiff_speed,
               speed_centre=a.speed_centre, speed_slope=a.speed_slope,
               mig_speed=a.mig_speed, mig_persist=a.mig_persist,
               dt=a.dt, occ_mode=a.occ_mode, seek=a.seek,
               seek_weight=a.seek_weight, p_seek=a.p_seek,
               cohere_weight=a.cohere_weight, cohere_gate=a.cohere_gate)


def _dispatch(argv=None) -> int:
    """Parse the CLI and dispatch to a subcommand.

    Args:
        argv: Arguments to parse, or None for sys.argv.

    Returns:
        Exit code: 0 on success, 1 if a panel in "all" failed.
    """
    ap = argparse.ArgumentParser(prog="occsim", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run", help="single invocation")
    _common(p)
    p.add_argument("--div", type=float, default=15.0)
    p.add_argument("--seed", type=int, default=1)

    p = sub.add_parser("sweep", help="DIV sweep, seed averaged")
    _common(p)
    p.add_argument("--divs", type=float, nargs="+",
                   default=[5, 8, 10, 12, 15, 18, 20, 25])
    p.add_argument("--seeds", type=int, default=20,
                   help="number of seeds to average over")
    p.add_argument("--baseline", type=float, default=5.0,
                   help="DIV used as the denominator for relative occupancy")
    p.add_argument("--csv", type=Path, default=None)
    p.add_argument("--figure", action="store_true")
    p.add_argument("--control", action="store_true",
                   help="also sweep with seek=none and plot it as the "
                        "reference on the SOMA_SCAFFOLD panel. Doubles the "
                        "runtime; `all` does this automatically")

    p = sub.add_parser("field", help="arbors over the territory fill")
    _common(p)
    p.add_argument("--divs", type=float, nargs="+", default=PANEL_DIVS)
    p.add_argument("--seed", type=int, default=1)

    p = sub.add_parser("compartments",
                       help="top-down per-compartment occupancy over DIV")
    _common(p)
    p.add_argument("--divs", type=float, nargs="+", default=PANEL_DIVS,
                   help="DIVs to run; a CSV is written for every one")
    p.add_argument("--panel-divs", type=float, nargs="+", default=None,
                   help="subset of --divs to draw panels for (default: all "
                        "of them). Lets the sweep be finer than the figure")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--grid", type=int, default=None,
                   help="deprecated alias for --compartments")

    p = sub.add_parser("density", help="occupancy vs seeding density")
    _common(p)
    p.add_argument("--div", type=float, default=15.0)
    # dense at low N where it's interesting, then out to saturation
    p.add_argument("--counts", type=int, nargs="+", default=DENSITY_COUNTS)
    p.add_argument("--seeds", type=int, default=8,
                   help="number of seeds to average over; seed noise is "
                        "23-25%% of the mean at N=1-3")
    p.add_argument("--csv", type=Path, default=None)

    p = sub.add_parser("stiffness",
                       help="carried map + evidence the coupling steers")
    _common(p)
    p.add_argument("--div", type=float, default=15.0)
    p.add_argument("--coupling", choices=["gradient", "speed"],
                   default="gradient",
                   help="which coupling to sweep: gradient (--stiff-weight, "
                        "for maps with a large-scale gradient) or speed "
                        "(--stiff-speed, the mechanism for Bessel/ring maps)")
    p.add_argument("--weights", type=float, nargs="+", default=None,
                   help="weight values to sweep; include 0. Defaults to "
                        "both signs for gradient, 0..1 for speed")
    p.add_argument("--seeds", type=int, default=12,
                   help="number of seeds to average over")
    p.add_argument("--grid", type=int, default=None,
                   help="deprecated alias for --compartments")
    p.add_argument("--csv", type=Path, default=None)

    p = sub.add_parser("regime",
                       help="classify the local physical regime per "
                            "compartment: mechanical, optical, competitive")
    _common(p)
    p.add_argument("--div", type=float, default=15.0,
                   help="DIV to classify at; the stiffness map in effect is "
                        "the v1 carry-forward one")
    p.add_argument("--csv", type=Path, default=None)
    # Only the SIGN of this reaches the classification, and 0.0, the shared
    # default from _common, has none, so the panel would report the optical
    # cue as weightless and split reinforcing/opposing the wrong way round.
    p.set_defaults(stiff_weight=REGIME_STIFF_WEIGHT)

    p = sub.add_parser("make-geometry",
                       help="write a scaffold geometry CSV from a shape and "
                            "a target scaffold area fraction")
    p.add_argument("out_path", type=Path,
                   help="where to write the geometry CSV")
    p.add_argument("--shape", choices=list(geometry.SHAPES), default="lattice",
                   help="lattice reproduces the printed strand grid (Li 2025); "
                        "disc and square are area-matched controls")
    p.add_argument("--area-frac", type=float, required=True,
                   help="scaffold area as a fraction of the field, 0-1")
    p.add_argument("--grid", type=int, default=30,
                   help="compartments per side; must match the run's "
                        "--compartments")
    p.add_argument("--field", type=float, default=1500.0,
                   help="overall field side, um")
    p.add_argument("--strand-pitch", type=float,
                   default=geometry.DEFAULT_STRAND_PITCH_UM,
                   help="strand centre-to-centre spacing for --shape lattice, "
                        "um")

    p = sub.add_parser("make-stiffness",
                       help="write a synthetic stiffness CSV from a beam "
                            "shape and the saturating-crosslink model")
    p.add_argument("out_path", type=Path,
                   help="where to write the stiffness CSV")
    p.add_argument("--beam", choices=list(dose.BEAMS), default="bessel",
                   help="gaussian is a single central lobe, the shape "
                        "--stiff-weight suits; bessel is J0^2 rings, the "
                        "shape --stiff-speed suits")
    p.add_argument("--ring-period", type=float,
                   default=dose.DEFAULT_RING_PERIOD_UM,
                   help="distance between consecutive J0 zeros for --beam "
                        "bessel, um; this is the feature scale a cell has to "
                        "cross, so it needs a grid fine enough to resolve it")
    p.add_argument("--waist", type=float, default=dose.DEFAULT_WAIST_UM,
                   help="1/e^2 intensity radius for --beam gaussian, um")
    p.add_argument("--envelope", choices=list(dose.ENVELOPES),
                   default="inverse-r",
                   help="bessel only: inverse-r is physical J0^2, whose 1/r "
                        "decay concentrates the value range centrally; flat "
                        "divides that out so the rings carry comparable power")
    p.add_argument("--grid", type=int, default=30,
                   help="compartments per side; must match the run's "
                        "--compartments")
    p.add_argument("--field", type=float, default=1500.0,
                   help="overall field side, um")
    p.add_argument("--divs", type=float, nargs="+",
                   default=[5, 10, 15, 20, 25],
                   help="DIVs to write, ascending and evenly spaced")
    p.add_argument("--e-floor", type=float, default=dose.DEFAULT_E_FLOOR_KPA,
                   help="stiffness of un-illuminated substrate, kPa")
    p.add_argument("--e-max", type=float, default=dose.DEFAULT_E_MAX_KPA,
                   help="stiffness at full crosslinking, kPa")
    p.add_argument("--kc", type=float, default=dose.DEFAULT_KC_PER_DAY,
                   help="crosslink rate per unit accumulated dose")

    p = sub.add_parser("all", help="the figure set used for results")
    _common(p)
    p.add_argument("--seeds", type=int, default=20,
                   help="seeds for the averaged sweeps")
    p.add_argument("--seed", type=int, default=1,
                   help="seed for the single-run panels (field, seek)")
    p.add_argument("--divs", type=float, nargs="+",
                   default=[5, 8, 10, 12, 15, 18, 20, 25])
    p.add_argument("--baseline", type=float, default=5.0)
    p.add_argument("--quick", action="store_true",
                   help="fewer seeds, coarser sweeps - a few minutes, not an hour")
    p.add_argument("--no-scaffold", action="store_true",
                   help="run over a bare field with no scaffold defined, so "
                        "occupancy is reported over the whole field. "
                        f"Otherwise `all` defaults to a "
                        f"{100 * DEFAULT_SCAFFOLD_AREA_FRAC:.0f} %% lattice")

    a = ap.parse_args(argv)

    # --grid was the old per-subcommand name for the same thing. Kept working,
    # but --compartments is the one flag now, since a stiffness map and a
    # geometry mask both have to sit on that grid and two names for it let a
    # run disagree with its own input files.
    if getattr(a, "grid", None) is not None:
        a.compartments = a.grid

    # make-geometry writes a file and exits: it runs no simulation, so it does
    # not go through _run_from() or need an output directory.
    if a.cmd == "make-geometry":
        provenance = (f"occsim make-geometry --shape {a.shape} "
                      f"--area-frac {a.area_frac:g} --grid {a.grid} "
                      f"--field {a.field:g}"
                      + (f" --strand-pitch {a.strand_pitch:g}"
                         if a.shape == "lattice" else ""))
        try:
            mask = geometry.build(a.shape, a.grid, a.field, a.area_frac,
                                  a.strand_pitch)
            realised = geometry.write(a.out_path, mask, a.grid, a.field,
                                      provenance)
        except ValueError as exc:
            raise SystemExit(f"make-geometry: {exc}")
        # The count is exact by construction; only rounding to a whole
        # compartment separates it from the request, and that gap shrinks with
        # the grid. Report both so the file's own area is never inferred.
        print(f"wrote {a.out_path}")
        print(f"  shape        {a.shape}")
        print(f"  grid         {a.grid}x{a.grid} "
              f"({a.field / a.grid:g} um per compartment)")
        print(f"  field side   {a.field:g} um")
        print(f"  scaffold     {sum(mask)} of {a.grid * a.grid} compartments "
              f"= {100 * realised:.2f} pct (requested "
              f"{100 * a.area_frac:.2f} pct)")
        return 0

    # Same deal as make-geometry: writes a file and exits. The maps it makes
    # are fixtures for exercising the interface, the real ones come from the
    # optical half, so the file labels itself SYNTHETIC and records this
    # command in its own header.
    if a.cmd == "make-stiffness":
        shape_arg = (f"--waist {a.waist:g}" if a.beam == "gaussian"
                     else f"--ring-period {a.ring_period:g} "
                          f"--envelope {a.envelope}")
        provenance = (f"occsim make-stiffness --beam {a.beam} {shape_arg} "
                      f"--grid {a.grid} --field {a.field:g} "
                      f"--divs {' '.join(f'{d:g}' for d in a.divs)} "
                      f"--e-floor {a.e_floor:g} --e-max {a.e_max:g} "
                      f"--kc {a.kc:g}")
        try:
            maps = dose.build(a.beam, a.grid, a.field, a.divs,
                              ring_period_um=a.ring_period,
                              waist_um=a.waist, envelope=a.envelope,
                              e_floor_kpa=a.e_floor, e_max_kpa=a.e_max,
                              kc_per_day=a.kc)
            dose.write(a.out_path, maps, a.divs, a.grid, a.field, provenance)
        except ValueError as exc:
            raise SystemExit(f"make-stiffness: {exc}")
        comp_um = a.field / a.grid
        print(f"wrote {a.out_path}")
        print(f"  beam         {a.beam}"
              + (f", waist {a.waist:g} um" if a.beam == "gaussian"
                 else f", ring period {a.ring_period:g} um, "
                      f"envelope {a.envelope}"))
        print(f"  grid         {a.grid}x{a.grid} ({comp_um:g} um per "
              f"compartment)")
        print(f"  field side   {a.field:g} um")
        print(f"  DIVs         {' '.join(f'{d:g}' for d in a.divs)}")
        print(f"  range        {min(m.min() for m in maps):.4f} to "
              f"{max(m.max() for m in maps):.4f} kPa")
        # The one number that decides whether a ring field couples at all: how
        # much of the value range sits in the top quartile (INTEGRATION_NOTES,
        # "The measured tables to judge against").
        share = dose.top_quartile_share(maps[-1])
        print(f"  top-quartile share of the range at DIV "
              f"{a.divs[-1]:g}: {share:.2f}")
        if a.beam == "bessel" and share > 0.5:
            print("  note: most of the range sits above p75, concentrated "
                  "centrally - a population mean will barely feel it. "
                  "Try --envelope flat.")
        # Nyquist would allow two compartments per period, but that is not
        # enough in practice: 200 um rings at 50 um compartments were measured
        # aliasing into a 337 um feature (INTEGRATION_NOTES), which is what makes
        # the map read as broad and inert. Warn at four or fewer per period.
        if a.beam == "bessel" and comp_um * 4 >= a.ring_period:
            print(f"  note: {comp_um:g} um compartments give only "
                  f"{a.ring_period / comp_um:.1f} per {a.ring_period:g} um "
                  f"ring period, which aliases into a broader feature; "
                  f"raise --grid.")
        return 0

    a.out.mkdir(parents=True, exist_ok=True)

    # `all` produces the committed result set, and the deliverable is occupancy
    # on the scaffold, so it defines one by default rather than running over a
    # bare field. The area fraction is a modelling input with no measured value
    # here; it is stated on every output and overridable both ways.
    if (a.cmd == "all" and not a.no_scaffold
            and a.scaffold_geometry is None and a.scaffold_area_frac is None):
        a.scaffold_area_frac = DEFAULT_SCAFFOLD_AREA_FRAC
    a.scaffold_geometry = _resolve_geometry(a)

    base = _run_from(a)
    if a.cmd == "all":
        seeds = list(range(1, (4 if a.quick else a.seeds) + 1))
        divs = QUICK_DIVS if a.quick else a.divs
        div = 15.0

        # The results figure set: every entry is something the paper needs.
        # Seek panels were dropped since those modes are parked.
        # panels only where the sweep actually has a DIV to draw
        panel_divs = [d for d in PANEL_DIVS if d in divs] or divs[:4]

        stiff_weights = STIFF_WEIGHTS_QUICK if a.quick else STIFF_WEIGHTS
        speed_weights = SPEED_WEIGHTS_QUICK if a.quick else SPEED_WEIGHTS

        # An explicit --stiffness-map applies to both panels; otherwise each
        # gets the shipped map its mechanism suits.
        def _resolve(default):
            path = a.stiffness_map or default
            return path if Path(path).exists() else None

        gradient_map = _resolve(DEFAULT_STIFFNESS_MAP)
        speed_map = _resolve(DEFAULT_BESSEL_MAP)

        # The coupling panels answer "does this mechanism steer", on a shipped
        # synthetic map at that map's own 30x30 grid, and they read stiffness
        # rather than occupancy. Carrying the scaffold in would add a second
        # grid to disagree with the map and change none of the numbers they
        # report, so they run without it.
        #
        # Dropping the geometry means dropping everything that READS it, or the
        # engine refuses the run ("seek=scaffold needs a scaffold geometry"). It
        # also keeps the question clean: a scaffold bias on top would mix up
        # "does stiffness steer" with "does the scaffold steer harder".
        coupling_base = base.replace(scaffold_geometry=None, compartments=30,
                                     seek="none", cohere_gate="none")

        def _on_map_field(path):
            """Pin a coupling run to the field the map itself declares.

            A stiffness map states the field side it was computed for, and the
            engine refuses a run whose field disagrees, correctly, since the
            compartments would not line up. The main operating point is free to
            use whatever field the experiment needs, so these panels follow the
            map rather than the run. Same reasoning as dropping the scaffold
            above: the question is whether the mechanism steers, not what the
            headline field is.
            """
            from .engine import read_stiffness
            try:
                _, _, meta = read_stiffness(path)
                field = float(meta.get("scaffold_um", a.field))
            except Exception:
                return coupling_base
            if field == coupling_base.scaffold:
                return coupling_base
            print(f"  coupling panel: map declares {field:g} um, running it "
                  f"there rather than at the {a.field:g} um operating point")
            return coupling_base.replace(scaffold=field)

        def _skip(which, default):
            # name the path that was actually looked for, which may be the
            # user's --stiffness-map rather than a shipped example
            def skip():
                print(f"  no stiffness map at {a.stiffness_map or default}; "
                      f"skipping the {which} coupling panel")
                return 0
            return skip

        # The committed set covers the two cues the brief names: the scaffold
        # (geometry, coverage on it, the rule that steers to it) and the
        # stiffness map. density and regime were dropped from it on 2026-08-26.
        # density is a methods check on metric saturation, not a deliverable,
        # and regime runs no simulation. Both are still subcommands.
        # The seek rule is the brief's item 8, so the committed set has to show
        # it doing something. It cannot be read off the absolute SOMA_SCAFFOLD
        # level, the central seed cluster starts on the scaffold, so the
        # sweep runs a second time with the rule off, and the figure shows the
        # gap. Doubles the sweep's cost, which is why nothing else gets one.
        control = (base.replace(seek="none")
                   if base.seek != "none" and base.scaffold_geometry else None)

        steps = [
            # headline: occupancy, normalised fold change, rate, and presence
            # on the scaffold against a rule-off control
            ("sweep", lambda: figures.sweep(
                base, divs, seeds, a.baseline, a.engine, a.out,
                a.out / "sweep.csv", True, control=control)),
            # spatial: which parts of the scaffold get covered. CSV per swept
            # DIV, panels every 5 so the figure stays readable.
            ("compartments", lambda: figures.compartments(
                base, divs, a.seed, a.compartments, a.engine, a.out,
                a.stiffness_map, a.scaffold_geometry,
                panel_divs=panel_divs)),
            # methods illustration: the real arbors behind the coverage
            ("field", lambda: figures.field(
                base, panel_divs, a.seed, a.engine, a.out,
                a.scaffold_geometry)),
            # the carried input plus evidence each coupling steers. Skipped
            # rather than failed with no map, since these are the panels whose
            # input doesn't live in the repo's engine. Gradient answers the
            # Gaussian-shaped case, speed the Bessel/ring-shaped one.
            ("stiffness gradient", (lambda: figures.stiffness(
                _on_map_field(gradient_map), gradient_map, stiff_weights,
                seeds[:COUPLING_SEEDS], div,
                COUPLING_COMPARTMENTS,
                a.engine, a.out, a.out / "stiffness_coupling.csv",
                coupling="gradient"))
             if gradient_map else _skip("gradient", DEFAULT_STIFFNESS_MAP)),
            ("stiffness speed", (lambda: figures.stiffness(
                _on_map_field(speed_map), speed_map, speed_weights,
                seeds[:COUPLING_SEEDS], div,
                COUPLING_COMPARTMENTS,
                a.engine, a.out, a.out / "stiffness_speed_coupling.csv",
                coupling="speed"))
             if speed_map else _skip("speed", DEFAULT_BESSEL_MAP)),
        ]

        failed = []
        for i, (name, fn) in enumerate(steps, 1):
            print(f"\n== [{i}/{len(steps)}] {name}", flush=True)
            try:
                if fn() != 0:
                    failed.append(name)
            except Exception as exc:          # one bad panel shouldn't lose the rest
                print(f"!! {name} failed: {exc}")
                failed.append(name)

        print(f"\nfigures in {a.out}")
        if failed:
            print(f"failed: {', '.join(failed)}")
            return 1
        return 0

    if a.cmd == "run":
        r = base.replace(div=a.div, seed=a.seed)
        m = measure(r, a.engine)
        # Both denominators, always, with the reported one named. A single bare
        # percentage was the thing that let whole-field and scaffold numbers be
        # quoted interchangeably.
        if m.scaffold is None:
            print(f"{m.whole_field:.4f} %   (whole field; no scaffold "
                  f"defined - pass --scaffold-area-frac to report occupancy "
                  f"on the scaffold)")
        else:
            print(f"{m.scaffold:.4f} %   on scaffold   "
                  f"[{m.whole_field:.4f} % over the whole field]")
        return 0

    if a.cmd == "sweep":
        seeds = list(range(1, a.seeds + 1))
        control = (base.replace(seek="none")
                   if a.control and base.seek != "none" and base.scaffold_geometry
                   else None)
        return figures.sweep(base, a.divs, seeds, a.baseline, a.engine, a.out,
                             a.csv, a.figure, control=control)

    if a.cmd == "field":
        return figures.field(base, a.divs, a.seed, a.engine, a.out,
                             a.scaffold_geometry)

    if a.cmd == "compartments":
        return figures.compartments(base, a.divs, a.seed, a.compartments,
                                    a.engine,
                                    a.out, a.stiffness_map,
                                    a.scaffold_geometry,
                                    panel_divs=a.panel_divs)

    if a.cmd == "density":
        return figures.density(base, a.counts, a.div, a.engine, a.out,
                               list(range(1, a.seeds + 1)), a.csv)

    if a.cmd == "regime":
        # Needs both cues. Neither is inferred: a missing one is an error
        # rather than a silently one-sided classification.
        if a.scaffold_geometry is None:
            raise SystemExit(
                "regime needs a scaffold geometry for the mechanical cue - "
                "pass --scaffold-area-frac to build one, or "
                "--scaffold-geometry to supply it")
        smap = a.stiffness_map or DEFAULT_STIFFNESS_MAP
        if not Path(smap).exists():
            raise SystemExit(f"regime needs a stiffness map for the optical "
                             f"cue; none at {smap} - pass --stiffness-map")
        return figures.regime(a.scaffold_geometry, smap, a.div,
                              a.compartments, a.field, a.seek_weight,
                              a.stiff_weight, a.out, a.csv)

    if a.cmd == "stiffness":
        smap = a.stiffness_map or DEFAULT_STIFFNESS_MAP
        if not Path(smap).exists():
            raise SystemExit(f"no stiffness map at {smap} - pass one with "
                             "--stiffness-map")
        weights = a.weights
        if weights is None:
            weights = (STIFF_WEIGHTS if a.coupling == "gradient"
                       else SPEED_WEIGHTS)
        return figures.stiffness(base, smap, weights,
                                 list(range(1, a.seeds + 1)), a.div,
                                 a.compartments,
                                 a.engine, a.out, a.csv,
                                 coupling=a.coupling)

    return 1


def main(argv=None) -> int:
    """Entry point: dispatch, and report an engine refusal as an error.

    The engine rejects setups that could never produce what was asked for, a
    stiffness coupling with migration off, a map whose grid disagrees with
    --compartments. Those messages name the problem and the fix, but they are
    raised from inside a worker thread, so unhandled they arrive under a
    concurrent.futures stack that buries the one line worth reading.

    Args:
        argv: Arguments to parse, or None for sys.argv.

    Returns:
        Exit code: 0 on success, 1 on a failed panel or an engine refusal.
    """
    try:
        return _dispatch(argv)
    except RuntimeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
