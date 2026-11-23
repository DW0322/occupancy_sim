""" 
All figures produced with help of claude, including descriptions, comments, and ugly palette choice
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from .engine import (ROOT, Run, measure, occupancy, occupancy_and_soma_stiffness,
                     rel_to_root,
                     seed_mean, seed_stats,
                     read_pgm, read_segments, read_compartments,
                     read_stiffness, select_stiffness_div, read_geometry)

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
           "#17becf", "#e377c2", "#8c564b", "#bcbd22", "#7f7f7f"]
BLUE, RED, DARK = "#2a6f97", "#c1121f", "#333"


def _progress(msg: str = "") -> None:
    """In-place progress on stderr, and only when stderr is a terminal.

    Deliberately not stdout: the printed tables are the committed record and
    have to stay byte-identical when piped or redirected. Skipped entirely off
    a tty so log files don't collect carriage returns.

    Args:
        msg: Text to show; empty clears the line.
    """
    if not sys.stderr.isatty():
        return
    sys.stderr.write(f"\r\033[K{msg}")
    sys.stderr.flush()


def _clean(ax) -> None:
    """Shared panel style: drop the top/right spines, add a light grid.

    Args:
        ax: Axes to style, modified in place.
    """
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(alpha=0.22)


# ---------------------------------------------------------------------------
def _tag(base: Run) -> str:
    """Filename tag for the morphology set, so pool and single runs don't
    overwrite each other's output.

    Args:
        base: The Run whose `morph` field identifies the set.

    Returns:
        "pool{N}" for a list of N, else the file stem with "morph" -> "single".
    """
    m = base.morph
    if isinstance(m, (list, tuple)):
        return f"pool{len(m)}"
    return Path(m).stem.replace("morph", "single")


def _curve(per_seed, divs, baseline_div, prefix: str = ""):
    """Seed statistics and normalised curve for one denominator.

    Baseline and ratio are both taken per seed, so the spread on the relative
    curve is the spread of the ratio, not a ratio of two separately noisy means.
    The rate is the same pairing differenced across the interval.

    Args:
        per_seed: {div: [value per seed]}.
        divs: Ascending DIVs, including baseline_div.
        baseline_div: DIV used as each seed's denominator.
        prefix: Prepended to every key, so two denominators can share a row.

    Returns:
        (rows, ratios): one dict per DIV keyed `{prefix}mean` ... and
        {div: per-seed ratio array}. The first row's rate keys are "".
    """
    baseline = np.array(per_seed[baseline_div], dtype=float)
    rows, ratios = [], {}
    for d in divs:
        v = np.array(per_seed[d], dtype=float)
        ratio = v / baseline
        ratios[d] = ratio
        sd = v.std(ddof=1) if len(v) > 1 else 0.0
        rows.append({
            f"{prefix}mean": v.mean(),
            f"{prefix}sd": sd,
            f"{prefix}sem": sd / len(v) ** 0.5 if len(v) > 1 else 0.0,
            f"{prefix}rel": ratio.mean(),
            f"{prefix}rel_sd": ratio.std(ddof=1) if len(v) > 1 else 0.0,
            f"{prefix}rate": "",
            f"{prefix}rate_sd": "",
        })
    # Rate: d(rel)/dDIV over each interval, differenced per seed then averaged.
    # Stored on the row for the interval's end DIV; the first row stays empty.
    for i in range(1, len(divs)):
        r = (ratios[divs[i]] - ratios[divs[i - 1]]) / (divs[i] - divs[i - 1])
        rows[i][f"{prefix}rate"] = r.mean()
        rows[i][f"{prefix}rate_sd"] = r.std(ddof=1) if len(r) > 1 else 0.0
    return rows, ratios


def _scaffold_area_fraction(path) -> float | None:
    """Scaffold's share of the field, for the chance line on the seek panel."""
    if path is None:
        return None
    try:
        from .engine import read_geometry
        m = read_geometry(path)
        return float(m.sum()) / float(m.size)
    except Exception:
        return None


def sweep(base: Run, divs, seeds, baseline_div, engine, out: Path,
          csv_path: Path | None, make_figure: bool, control: Run | None = None) -> int:
    """DIV sweep, seed-averaged, reported absolute and relative.

    Relative occupancy is fold change over the baseline DIV. Absolute percent
    doesn't discriminate at realistic densities, so the normalised curve is what
    gets reported.

    **Both denominators are always written.** When the run defines a scaffold,
    the reported quantity is coverage over the scaffold, the collaboration's
    "area of occupancy on top of the scaffold", and it is the curve the figure
    draws. The whole-field columns stay in the same unprefixed names they have
    always had, so a column never changes meaning between runs and a scaffold
    run stays comparable with a bare-field one.

    Args:
        base: Base run; `div` is overridden per swept DIV. Carrying a
            scaffold_geometry is what turns the scaffold columns on.
        divs: DIVs to sweep. baseline_div is added if missing.
        seeds: Seeds to average over at each DIV.
        baseline_div: Denominator DIV for relative occupancy.
        engine: Path to the compiled occupancy_engine binary.
        out: Output directory for the CSV/PNG.
        csv_path: CSV output path, or None to derive one from `out`.
        make_figure: Also render the PNG, four panels when a geometry is
            carried, three otherwise.
        control: A second run, swept over the same DIVs and seeds, whose
            SOMA_SCAFFOLD curve is the reference the reported one is read
            against. Pass the same run with `seek="none"` to answer "does the
            decision rule do anything", which the coverage panels cannot: cells
            are seeded as a central cluster that already sits on the scaffold,
            so the absolute percentage is set by the seeding and only the gap
            between the two curves belongs to the rule.

    Returns:
        0. An int for consistency with the other figure functions.
    """
    divs = sorted(set(list(divs) + [baseline_div]))
    print(f"n={base.n}  mig_speed={base.mig_speed}  mode={base.occ_mode}  "
          f"seek={base.seek}  seeds={len(seeds)}")

    per_field: dict[float, list[float]] = {}
    per_scaffold: dict[float, list[float]] = {}
    per_clumping: dict[float, tuple] = {}
    # The table below prints only once every DIV is in, so without this the
    # whole sweep runs silent. Cost is linear in DIV, not per-row, migration
    # takes one step an hour, so DIV25 is five times DIV5, which puts about
    # 40 % of the wall time in the last two rows. Weight the fraction by DIV
    # so it tracks work rather than rows.
    # Somata ON the scaffold, the readout the seek rule is judged on. Occupancy
    # cannot stand in for it: guidance spreads cells ALONG a strand, which
    # sweeps less fresh area than a free walk, so coverage falls exactly when
    # the rule is working. Populated whenever a geometry is carried.
    per_on_scaffold: dict = {}
    per_scaffold_dist: dict = {}
    total_div = sum(divs)
    done_div = 0.0
    for d in divs:
        _progress(f"    sweep: DIV {d:g}, {len(seeds)} seeds "
                  f"({100 * done_div / total_div:.0f} % of DIV-hours done)")
        st = seed_stats(base.replace(div=d), seeds, engine)
        done_div += d
        per_field[d] = st.whole_field[3]
        if st.scaffold is not None:
            per_scaffold[d] = st.scaffold[3]
        if st.clark_evans is not None:
            per_clumping[d] = st.clark_evans
        if st.soma_scaffold is not None:
            per_on_scaffold[d] = st.soma_scaffold
            per_scaffold_dist[d] = st.soma_scaffold_dist
    _progress()

    # The control: same everything, rule off. Only its SOMA_SCAFFOLD curve is
    # kept, its coverage is not comparable, since guidance trades fresh area
    # for time spent on the strand.
    per_control: dict = {}
    if control is not None:
        for i, d in enumerate(divs):
            _progress(f"    control ({control.seek}): DIV {d:g}, "
                      f"{len(seeds)} seeds ({100 * i / len(divs):.0f} %)")
            cst = seed_stats(control.replace(div=d), seeds, engine)
            if cst.soma_scaffold is not None:
                per_control[d] = cst.soma_scaffold
        _progress()

    on_scaffold = len(per_scaffold) == len(divs)
    # Arrangement, alongside coverage. They answer different questions and can
    # move independently: a population that clumps covers less ground but two
    # populations at the same coverage can be arranged quite differently, and
    # at plating density coverage saturates while this does not.
    has_clumping = len(per_clumping) == len(divs)
    has_on_scaffold = len(per_on_scaffold) == len(divs)
    has_control = len(per_control) == len(divs) and has_on_scaffold

    field_rows, _ = _curve(per_field, divs, baseline_div)
    rows = [dict(div=d, **fr) for d, fr in zip(divs, field_rows)]
    if has_clumping:
        for row, d in zip(rows, divs):
            mean, sd, sem, _ = per_clumping[d]
            row.update(clark_evans_mean=mean, clark_evans_sd=sd,
                       clark_evans_sem=sem)
    if on_scaffold:
        # the reported curve is the scaffold one whenever there is a scaffold
        scaf_rows, _ = _curve(per_scaffold, divs, baseline_div, "scaffold_")
        for row, sr in zip(rows, scaf_rows):
            row.update(sr)
    if has_on_scaffold:
        for row, d in zip(rows, divs):
            mean, sd, sem, _ = per_on_scaffold[d]
            dmean, dsd, _, _ = per_scaffold_dist[d]
            row.update(soma_on_scaffold_pct=mean, soma_on_scaffold_sd=sd,
                       soma_on_scaffold_sem=sem,
                       soma_scaffold_dist_um=dmean,
                       soma_scaffold_dist_sd=dsd)
    if has_control:
        for row, d in zip(rows, divs):
            mean, sd, sem, _ = per_control[d]
            row.update(control_soma_on_scaffold_pct=mean,
                       control_soma_on_scaffold_sd=sd,
                       control_soma_on_scaffold_sem=sem)

    # Console: the reported denominator first, the other alongside it.
    if on_scaffold:
        print(f"  {'DIV':>5} {'scaffold %':>11} {'sd':>7} {'sem':>6} "
              f"{'rel':>7} {'rel sd':>7} | {'field %':>8}")
        for row in rows:
            print(f"  {row['div']:>5.0f} {row['scaffold_mean']:>11.3f} "
                  f"{row['scaffold_sd']:>7.3f} {row['scaffold_sem']:>6.3f} "
                  f"{row['scaffold_rel']:>7.3f} {row['scaffold_rel_sd']:>7.3f} "
                  f"| {row['mean']:>8.3f}")
    else:
        print(f"  {'DIV':>5} {'mean %':>9} {'sd':>7} {'sem':>6} "
              f"{'rel':>7} {'rel sd':>7}")
        for row in rows:
            print(f"  {row['div']:>5.0f} {row['mean']:>9.3f} {row['sd']:>7.3f} "
                  f"{row['sem']:>6.3f} {row['rel']:>7.3f} "
                  f"{row['rel_sd']:>7.3f}")

    key = "scaffold_rate" if on_scaffold else "rate"
    print(f"  {'interval':>12} {'rate/DIV':>9} {'rate sd':>8}")
    for i in range(1, len(divs)):
        print(f"  {f'{divs[i - 1]:.0f}->{divs[i]:.0f}':>12} "
              f"{rows[i][key]:>9.4f} {rows[i][key + '_sd']:>8.4f}")

    if base.mig_speed <= 0 and max(divs) > 15:
        print("\n  note: migration is off and the sweep runs past DIV15, where "
              "the reveal\n  schedule saturates. Everything from DIV15 on is "
              "identical by construction.")

    csv_path = csv_path or (out / f"sweep_n{base.n}_{base.occ_mode}_{_tag(base)}.csv")
    with open(csv_path, "w", newline="") as f:
        # Every sweep CSV carries the invocation that made it, one turned up
        # without its conditions once and had to be reproduced to recover them.
        nm = (f"pool of {len(base.morph)} files"
              if isinstance(base.morph, (list, tuple))
              else rel_to_root(base.morph))
        f.write(f"# occsim sweep  n={base.n} mig_speed={base.mig_speed} "
                f"mig_persist={base.mig_persist} dt={base.dt} "
                f"occ_mode={base.occ_mode} seek={base.seek} "
                f"seed_frac={base.seed_frac} "
                f"seek_weight={base.seek_weight} p_seek={base.p_seek} "
                f"territory={base.territory} "
                f"scaffold={base.scaffold} res={base.res} morph={nm} "
                f"seeds={list(seeds)} baseline_div={baseline_div} "
                f"compartments={base.compartments} "
                f"cohere_weight={base.cohere_weight} "
                f"cohere_gate={base.cohere_gate} "
                f"scaffold_geometry="
                f"{rel_to_root(base.scaffold_geometry)}\n")
        if on_scaffold:
            f.write("# REPORTED = the scaffold_* columns: covered area inside "
                    "the scaffold over the area of the scaffold, at "
                    "compartment granularity.\n")
            f.write("# The unprefixed columns are the SAME runs over the whole "
                    "field. They keep their names so no column ever changes "
                    "denominator, and so scaffold runs stay comparable with "
                    "bare-field ones.\n")
        else:
            f.write("# No scaffold was defined, so the unprefixed columns are "
                    "over the whole field. Pass --scaffold-area-frac to get "
                    "scaffold_* columns.\n")
        f.write("# rate = per-seed d(rel)/dDIV over the interval ending at "
                "this row's div; empty on the first row\n")
        if has_on_scaffold:
            f.write("# soma_on_scaffold_pct = percent of SOMATA sitting in a "
                    "scaffold compartment, and soma_scaffold_dist_um their "
                    "mean distance to the nearest one. This is the direct "
                    "analogue of what Li 2025 reports - presence on the "
                    "structure - and the readout the --seek rule is judged "
                    "on.\n")
            f.write("# Do NOT read the coverage columns for that question. "
                    "Guidance spreads cells ALONG a strand, which sweeps less "
                    "fresh area than a free walk, so scaffold coverage falls "
                    "when the rule is working.\n")
            f.write("# Nor read the ABSOLUTE percentage as evidence the rule "
                    "works. Cells are seeded as a central cluster, and at "
                    "--seed-frac 0.1 on the shipped lattice that cluster "
                    "starts almost entirely ON the scaffold, so the run opens "
                    "high whatever the rule does.\n")
        if has_control:
            f.write(f"# control_* = the SAME sweep with seek={control.seek}, "
                    f"the reference the reported curve is read against. The "
                    f"GAP between soma_on_scaffold_pct and "
                    f"control_soma_on_scaffold_pct is what belongs to the "
                    f"decision rule; the level they share is the seeding.\n")
        if has_clumping:
            f.write("# clark_evans_* = arrangement, not coverage: observed "
                    "mean nearest-neighbour distance over what a random "
                    "scatter of the same density would give, edge-corrected. "
                    "<1 clustered, 1 random, >1 dispersed.\n")
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {csv_path}")

    if make_figure:
        pre = "scaffold_" if on_scaffold else ""
        where = "on scaffold" if on_scaffold else "whole field"
        d = [r["div"] for r in rows]
        # The fourth panel is SOMA_SCAFFOLD, presence on the structure, which
        # is what Li 2025 actually reports and what the seek rule is judged on.
        # Without it the figure set shows only coverage, which moves the other
        # way under guidance.
        n_panels = 4 if has_on_scaffold else 3
        fig, ax = plt.subplots(1, n_panels, figsize=(5.2 * n_panels, 4.0),
                               facecolor="white")
        ax[0].errorbar(d, [r[pre + "mean"] for r in rows],
                       yerr=[r[pre + "sem"] for r in rows],
                       fmt="-o", color=BLUE, lw=1.8, markersize=4, capsize=3,
                       label=f"occupancy {where}")
        if on_scaffold:
            # The whole-field curve stays visible as the reference the pinned
            # numbers are in, so the two denominators can be read against each
            # other instead of one quietly standing in for the other.
            ax[0].errorbar(d, [r["mean"] for r in rows],
                           yerr=[r["sem"] for r in rows], fmt="--s",
                           color="#9aa5ad", lw=1.3, markersize=3.5, capsize=2,
                           markerfacecolor="white", label="whole field")
            ax[0].legend(fontsize=8, loc="lower right")
        ax[0].set_xlabel("DIV"); ax[0].set_ylabel("occupancy (%)")
        ax[0].set_ylim(bottom=0)
        ax[1].errorbar(d, [r[pre + "rel"] for r in rows],
                       yerr=[r[pre + "rel_sd"] / len(seeds) ** 0.5
                             for r in rows],
                       fmt="-o", color=RED, lw=1.8, markersize=4, capsize=3)
        ax[1].axhline(1.0, color="#bbb", lw=1, ls=":")
        ax[1].set_xlabel("DIV")
        ax[1].set_ylabel(f"occupancy relative to DIV{baseline_div:.0f}")
        mid = [(divs[i - 1] + divs[i]) / 2 for i in range(1, len(divs))]
        rr = [rows[i][pre + "rate"] for i in range(1, len(rows))]
        rs = [rows[i][pre + "rate_sd"] / len(seeds) ** 0.5
              for i in range(1, len(rows))]
        ax[2].errorbar(mid, rr, yerr=rs, fmt="-o", color=DARK, lw=1.8,
                       markersize=4, capsize=3)
        ax[2].axhline(0.0, color="#bbb", lw=1, ls=":")
        ax[2].set_xlabel("DIV (interval midpoint)")
        ax[2].set_ylabel("rate of occupancy (rel / DIV)")

        if has_on_scaffold:
            ax[3].errorbar(d, [r["soma_on_scaffold_pct"] for r in rows],
                           yerr=[r["soma_on_scaffold_sem"] for r in rows],
                           fmt="-o", color=BLUE, lw=1.8, markersize=4,
                           capsize=3,
                           label=f"seek={base.seek}")
            if has_control:
                ax[3].errorbar(
                    d, [r["control_soma_on_scaffold_pct"] for r in rows],
                    yerr=[r["control_soma_on_scaffold_sem"] for r in rows],
                    fmt="--s", color="#9aa5ad", lw=1.3, markersize=3.5,
                    capsize=2, markerfacecolor="white",
                    label=f"control, seek={control.seek}")
            # Chance: the scaffold's share of the field. A population scattered
            # at random sits here, so the distance above it is enrichment -
            # though at --seed-frac < 1 the seeding supplies most of it, which
            # is what the control curve separates out.
            chance = _scaffold_area_fraction(base.scaffold_geometry)
            if chance is not None:
                ax[3].axhline(100 * chance, color="#c8b0b0", lw=1, ls=":")
                ax[3].annotate(f"chance ({100 * chance:.0f} % of field)",
                               xy=(d[0], 100 * chance), xytext=(0, 4),
                               textcoords="offset points", fontsize=7,
                               color="#a08080")
            ax[3].set_xlabel("DIV")
            ax[3].set_ylabel("somata on scaffold (%)")
            ax[3].set_ylim(0, 100)
            ax[3].legend(fontsize=8, loc="lower right")

        for a in ax:
            _clean(a)
        nm = (f"{len(base.morph)} morphologies sampled per cell"
              if isinstance(base.morph, (list, tuple)) else "one morphology")
        cap = (f"N={base.n}, {len(seeds)} seeds, {nm}, mean +/- sem; "
               f"panels 2-3 follow the \"{where}\" curve.")
        if has_on_scaffold:
            cap += ("  Panel 4 is SOMA_SCAFFOLD - presence on the structure, "
                    "the readout the seek rule is judged on; coverage moves "
                    "the other way under guidance.")
            if has_control:
                cap += (f"  Read the GAP to the control, not the level: cells "
                        f"are seeded as a central cluster that already starts "
                        f"on the scaffold.")
        if on_scaffold:
            cap += (f"  Scaffold {Path(base.scaffold_geometry).name}, "
                    f"{base.compartments}x{base.compartments} compartments.")
        # Four panels make the caption long enough to run off the canvas, so
        # wrap it to the figure width rather than letting it clip.
        import textwrap
        cap = "\n".join(textwrap.wrap(cap, width=42 * n_panels))
        fig.text(0.5, -0.02 - 0.035 * cap.count("\n"), cap, ha="center",
                 va="top", fontsize=8, color="#555")
        p = out / f"sweep_n{base.n}_{base.occ_mode}_{_tag(base)}.png"
        fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
        print(f"wrote {p}")
    return 0


# ---------------------------------------------------------------------------
def field(base: Run, divs, seed, engine, out: Path,
          geometry_path: Path | None = None) -> int:
    """Per-cell arbors drawn over the scaffold they sit on, one panel per DIV.

    The literal picture the collaboration asked for: no biophysics, just the
    shape of the neurons over the scaffold, with the occupancy that shape
    produces printed underneath. The scaffold is drawn as a grey underlay so
    the arbors stay readable on top of it.

    Args:
        base: Base run; `div` is overridden per panel.
        divs: DIVs to render, one panel each.
        seed: Single seed; this is a picture of one realisation, not a
            population statistic.
        engine: Path to the compiled occupancy_engine binary.
        out: Output directory; `_frames/` holds the scratch .pgm/.csv files.
        geometry_path: Scaffold geometry CSV to draw underneath the arbors, or
            None for a bare field (docs/SCAFFOLD_GEOMETRY_FORMAT.md).

    Returns:
        0 on success.
    """
    tmp = out / "_frames"
    tmp.mkdir(parents=True, exist_ok=True)
    L = base.scaffold
    gmask = read_geometry(geometry_path) if geometry_path is not None else None
    fig, axes = plt.subplots(1, len(divs), figsize=(3.7 * len(divs), 4.0),
                             facecolor="white")
    axes = np.atleast_1d(axes)

    for ax, d in zip(axes, divs):
        pgm, seg = tmp / f"f{d:g}.pgm", tmp / f"f{d:g}.csv"
        r = base.replace(div=d, seed=seed, pgm=pgm, segments=seg,
                         scaffold_geometry=geometry_path)
        m = measure(r, engine)
        if gmask is not None:
            # scaffold material, under everything else
            under = np.zeros((*gmask.shape, 4))
            under[..., :3] = 0.42
            under[..., 3] = np.where(gmask, 0.30, 0.0)
            ax.imshow(under, extent=[-L / 2, L / 2, -L / 2, L / 2],
                      origin="lower", interpolation="nearest", zorder=0)
        # The engine writes a binary coverage PGM, 0 covered, 255 bare, so
        # any mid-grey cut separates them. 128 is that midpoint, not a
        # tuned threshold.
        fill = (read_pgm(pgm) < 128).astype(float)
        rgba = np.zeros((*fill.shape, 4))
        rgba[..., 0], rgba[..., 1], rgba[..., 2] = 0.16, 0.43, 0.59
        rgba[..., 3] = fill * 0.10
        ax.imshow(rgba, extent=[-L / 2, L / 2, -L / 2, L / 2], origin="lower",
                  interpolation="nearest", zorder=1)
        for cid, sv in sorted(read_segments(seg).items()):
            c = PALETTE[cid % len(PALETTE)]
            ax.add_collection(LineCollection(sv, colors=c, linewidths=0.35,
                                             alpha=0.85, zorder=2))
        ax.set_xlim(-L / 2, L / 2); ax.set_ylim(-L / 2, L / 2)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_edgecolor("#ccc")
        label = (f"DIV{d:g}   {m.reported:.1f} %"
                 + ("  on scaffold" if m.scaffold is not None else ""))
        ax.set_xlabel(label, fontsize=9, color="#333")

    axes[0].plot([-L / 2 + 90, -L / 2 + 290], [-L / 2 + 100] * 2, color="#222",
                 lw=3, zorder=4)
    axes[0].text(-L / 2 + 90, -L / 2 + 150, "200 um", fontsize=8, color="#222")
    # Same caveat as the compartments panel: one draw, and which draw matters.
    # The spread is quoted from sweep.csv rather than hardcoded, because it
    # depends on N and the operating point has changed once already.
    cap = (f"N={base.n}, seed {seed}, migration {base.mig_speed:g} um/h. "
           f"ONE realisation, not a population mean - see sweep.csv.")
    if gmask is not None:
        cap += (f"  Grey is scaffold ({100 * gmask.mean():.0f} % of the "
                f"field); percentages are coverage of it.")
    fig.text(0.5, -0.01, cap, ha="center", fontsize=8, color="#555")
    p = out / f"field_n{base.n}_{_tag(base)}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"wrote {p}")
    return 0


# ---------------------------------------------------------------------------
def compartments(base: Run, divs, seed, grid, engine, out: Path,
                 stiffness_path: Path | None = None,
                 geometry_path: Path | None = None,
                 panel_divs=None) -> int:
    """Top-down per-compartment occupancy, one panel per DIV, optionally with
    the stiffness map in effect on a second row.

    The scaffold as an n x n grid of squares coloured by their own occupancy, so
    you can see the spatial pattern behind one global percentage. The engine
    already asserts that compartments sum to the global number; the sums get
    re-derived here from the CSV as an independent check on the round-trip.

    Watch out for edge compartments reading low, roughly 58 % at a corner
    against 87 % in the interior on the pinned run, because arbors overhang the
    field and dilation only pulls part of them back. That's geometry, not a
    boundary effect.

    The stiffness row is a carried input and changes nothing: the engine
    validates it on every panel run, so a bad file aborts before anything is
    drawn. Panels use the v1 rule (largest map div <= the panel's DIV), so the
    label may name a different DIV than the panel.

    One seed, like the field view. Colour scales are shared across panels so
    they're comparable.

    Args:
        base: Base run; `div` is overridden per panel.
        divs: DIVs to render.
        seed: Single seed to run.
        grid: Compartments per side of the field.
        engine: Path to the compiled occupancy_engine binary.
        out: Output directory; `compartments/` holds the per-DIV CSVs.
        stiffness_path: Stiffness map CSV to carry and draw as a second row, or
            None for no row (docs/STIFFNESS_FORMAT.md).
        geometry_path: Scaffold geometry CSV to overlay as faded
            non-scaffold compartments, or None (docs/SCAFFOLD_GEOMETRY_FORMAT.md).
        panel_divs: Subset of `divs` to draw panels for; None draws all. A CSV
            is written for every DIV in `divs` either way, so the data stays at
            the sweep's resolution while the figure stays readable.

    Returns:
        0 on success.

    Raises:
        RuntimeError: If the aggregate re-derived here disagrees with the
            engine's by more than 1e-3, meaning the CSV round-trip lost data.
    """
    L = base.scaffold
    divs = list(divs)
    panel_divs = list(panel_divs) if panel_divs is not None else list(divs)
    missing = [d for d in panel_divs if d not in divs]
    if missing:
        raise RuntimeError(
            f"panel_divs {missing} are not in divs {divs} - a panel can only "
            "be drawn for a DIV that was actually run")

    csv_dir = out / "compartments"
    csv_dir.mkdir(parents=True, exist_ok=True)

    # Clear this run's own per-DIV CSVs first. Filenames carry the DIV, so a
    # previous run over a different DIV set would leave files sitting next to
    # the new ones with no way to tell them apart.
    stale = sorted(csv_dir.glob(
        f"compartments_n{base.n}_g{grid}_{_tag(base)}_div*_seed{seed}.csv"))
    for old in stale:
        old.unlink()
    if stale:
        print(f"  cleared {len(stale)} CSV(s) from a previous run")

    gmask = read_geometry(geometry_path) if geometry_path is not None else None

    sdivs, smaps, smeta = (None, None, None)
    if stiffness_path is not None:
        sdivs, smaps, smeta = read_stiffness(stiffness_path)
        svals = np.concatenate([m.ravel() for m in smaps.values()])
        svmin, svmax = float(svals.min()), float(svals.max())
    nrows = 1 + (stiffness_path is not None)
    fig, axes = plt.subplots(nrows, len(panel_divs),
                             figsize=(3.7 * len(panel_divs), 4.4 * nrows),
                             facecolor="white")
    axes = np.array(axes).reshape(nrows, len(panel_divs))

    im = im_s = None
    col = -1
    for d in divs:
        cpath = csv_dir / (f"compartments_n{base.n}_g{grid}_{_tag(base)}"
                           f"_div{d:g}_seed{seed}.csv")
        r = base.replace(div=d, seed=seed, compartments=grid,
                         compartments_csv=cpath, stiffness_map=stiffness_path,
                         scaffold_geometry=geometry_path)
        occ = occupancy(r, engine)   # raises if the engine's checks fail
        g, covered, total, _ = read_compartments(cpath)
        agg = 100.0 * covered.sum() / total.sum()
        # occ comes off stdout at print precision, and the bit-exact check
        # already ran in the engine; this just guards the file round-trip
        if abs(agg - occ) > 1e-3:
            raise RuntimeError(
                f"compartment aggregate {agg} != global {occ} for {cpath}")
        if d not in panel_divs:
            continue          # CSV written, no panel
        col += 1
        ax = axes[0, col]
        im = ax.imshow(g, extent=[-L / 2, L / 2, -L / 2, L / 2],
                       origin="lower", interpolation="nearest",
                       cmap="viridis", vmin=0, vmax=100)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_edgecolor("#ccc")
        label = f"DIV{d:g}   {occ:.1f} % of the field"
        if gmask is not None:
            # white-wash everything that isn't scaffold; the scaffold number
            # leads, since that is the one the run reports
            off = np.ones((*gmask.shape, 4))
            off[..., 3] = np.where(gmask, 0.0, 0.80)
            ax.imshow(off, extent=[-L / 2, L / 2, -L / 2, L / 2],
                      origin="lower", interpolation="nearest", zorder=3)
            restricted = (100.0 * covered[gmask].sum() / total[gmask].sum())
            label = (f"DIV{d:g}   {restricted:.1f} % on scaffold"
                     f"\n({occ:.1f} % of the whole field)")
        ax.set_xlabel(label, fontsize=9, color="#333")

        if stiffness_path is not None:
            axs = axes[1, col]
            sel = select_stiffness_div(sdivs, d)
            if sel is None:
                axs.set_facecolor("#f4f4f4")
                axs.text(0.5, 0.5, f"no map \u2264 DIV{d:g}",
                         transform=axs.transAxes, ha="center", va="center",
                         fontsize=8, color="#777")
                axs.set_xlabel(f"(earliest map DIV{sdivs[0]:g})",
                               fontsize=9, color="#333")
            else:
                im_s = axs.imshow(smaps[sel],
                                  extent=[-L / 2, L / 2, -L / 2, L / 2],
                                  origin="lower", interpolation="nearest",
                                  cmap="magma", vmin=svmin, vmax=svmax)
                axs.set_xlabel(f"stiffness, map DIV{sel:g}",
                               fontsize=9, color="#333")
            axs.set_xticks([]); axs.set_yticks([])
            for sp in axs.spines.values():
                sp.set_edgecolor("#ccc")

    cb = fig.colorbar(im, ax=list(axes[0]), fraction=0.02, pad=0.015)
    cb.set_label("compartment occupancy (%)", fontsize=8, color="#555")
    cb.ax.tick_params(labelsize=7)
    if im_s is not None:
        cbs = fig.colorbar(im_s, ax=list(axes[1]), fraction=0.02, pad=0.015)
        cbs.set_label(f"stiffness ({smeta.get('units', '?')})",
                      fontsize=8, color="#555")
        cbs.ax.tick_params(labelsize=7)

    cap = (f"N={base.n}, seed {seed}, {grid}x{grid} compartments "
           f"({L / grid:.0f} um squares), mig_speed={base.mig_speed:g}. "
           f"ONE realisation, not a population mean - see sweep.csv.")
    if len(panel_divs) != len(divs):
        cap += (f"  Panels at DIV {', '.join(f'{d:g}' for d in panel_divs)}; "
                f"CSVs for all {len(divs)} swept DIVs in {csv_dir.name}/.")
    if stiffness_path is not None:
        weight = getattr(base, "stiff_weight", 0.0)
        if weight:
            # sign is the direction: +ve toward soft, -ve toward stiff
            toward = "soft-ward" if weight > 0 else "stiff-ward"
            cap += (f"\nStiffness {Path(stiffness_path).name}: COUPLED, "
                    f"{toward} bias, weight {weight:g}.")
        else:
            cap += (f"\nStiffness {Path(stiffness_path).name}: carried, "
                    f"not coupled.")
    if geometry_path is not None:
        # Two switches couple a geometry, not one: --plate-scaffold confines
        # the walk to it, --seek scaffold steers toward it and then along it.
        # Naming only the first read "carried, not coupled" on figures whose
        # run was steered by the mask.
        if getattr(base, "plate_scaffold", False):
            coupled = "COUPLED, cells confined to it"
        elif getattr(base, "seek", "none") == "scaffold":
            coupled = (f"COUPLED via --seek scaffold, weight "
                       f"{getattr(base, 'seek_weight', 0):g}")
        else:
            coupled = "carried, not coupled"
        cap += (f"\nScaffold {Path(geometry_path).name}: {coupled}. Faded "
                f"compartments are off-scaffold; the reported number is "
                f"coverage ON the scaffold, whole-field below it.")
    fig.text(0.5, -0.02, cap, ha="center", fontsize=8, color="#555")
    p = out / f"compartments_n{base.n}_g{grid}_{_tag(base)}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"wrote {p}")
    print(f"per-DIV CSVs in {csv_dir}")
    return 0


# ---------------------------------------------------------------------------
def density(base: Run, counts, div, engine, out: Path, seeds=None,
            csv_path: Path | None = None) -> int:
    """Occupancy vs projected areal density, for the dilated territory and the
    bare skeleton.

    The point: at realistic plating density even the undilated skeleton fills
    the plane, so the territory radius R isn't the lever and absolute occupancy
    can't discriminate between runs. What moves it is how many arbors project
    onto the 2D plane.

    Seed-averaged, because one seed is uninterpretable in exactly the region
    this figure is about. At DIV15 with migration the seed sd of the dilated
    curve is 23-25 % of the mean at N=1-3 and only drops under 1 % by N=15. At
    one seed the low-N end was actively misleading: N=3 came out at 78.8 %
    against a 59.4 % eight-seed mean, contradicting `sweep`'s 57.7 % under the
    same conditions.

    Args:
        base: Base run; `n`, `div`, `seed` and `territory` are overridden.
        counts: Cell counts to sweep.
        div: DIV to evaluate at.
        engine: Path to the compiled occupancy_engine binary.
        out: Output directory for the PNG.
        seeds: Seeds to average over. None uses `base.seed` alone, which is only
            good for a quick look, see above.
        csv_path: CSV output path; falsy writes no CSV.

    Returns:
        0 on success.
    """
    seeds = list(seeds) if seeds is not None else [base.seed]
    area_cm2 = (base.scaffold * 1e-4) ** 2
    dens = [n / area_cm2 for n in counts]
    # use this run's R rather than a literal, so --territory reaches here and
    # the panel label always matches what was actually run
    r_um = base.territory
    print(f"  {'N':>6} {'cells/cm2':>10} {f'R={r_um:g}':>9} {'sd':>6} "
          f"{'skeleton':>9} {'sd':>6} {'scaffold':>9}")
    dil, dil_e, skel, skel_e, scaf, scaf_e, rows = [], [], [], [], [], [], []
    # Runtime is linear in N, so the largest counts dominate and they are the
    # last rows printed: at the committed counts the final two are about 74 %
    # of the wall time, which looks like a hang with nothing on the terminal.
    # Announce the count being worked on rather than only reporting it after.
    total_cells = sum(counts) * len(seeds) * 2
    done_cells = 0
    for n in counts:
        _progress(f"    density: N={n}, {len(seeds)} seeds x 2 curves "
                  f"({100 * done_cells / total_cells:.0f} % of cells done)")
        d = seed_stats(base.replace(n=n, div=div), seeds, engine)
        s = seed_stats(base.replace(n=n, div=div, territory=base.res),
                       seeds, engine)
        done_cells += n * len(seeds) * 2
        d_m, d_sd, d_sem, _ = d.whole_field
        s_m, s_sd, s_sem, _ = s.whole_field
        dil.append(d_m); dil_e.append(d_sem)
        skel.append(s_m); skel_e.append(s_sem)
        row = dict(n_cells=n, cells_per_cm2=n / area_cm2,
                   dilated_mean=d_m, dilated_sd=d_sd, dilated_sem=d_sem,
                   skeleton_mean=s_m, skeleton_sd=s_sd, skeleton_sem=s_sem)
        # The saturation this figure is about is a property of the dilated
        # territory, so only that curve gets a scaffold denominator; the
        # skeleton curve is a control for R and stays whole-field.
        if d.scaffold is not None:
            sc_m, sc_sd, sc_sem, _ = d.scaffold
            scaf.append(sc_m); scaf_e.append(sc_sem)
            row.update(dilated_scaffold_mean=sc_m, dilated_scaffold_sd=sc_sd,
                       dilated_scaffold_sem=sc_sem)
        row["n_seeds"] = len(seeds)
        rows.append(row)
        _progress()
        print(f"  {n:>6} {n / area_cm2:>10.0f} {d_m:>9.2f} {d_sd:>6.2f} "
              f"{s_m:>9.2f} {s_sd:>6.2f} "
              + (f"{scaf[-1]:>9.2f}" if scaf else f"{'-':>9}"))
    on_scaffold = len(scaf) == len(counts)

    if csv_path:
        Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
        with open(csv_path, "w", newline="") as f:
            nm = (f"pool of {len(base.morph)} files"
                  if isinstance(base.morph, (list, tuple))
                  else rel_to_root(base.morph))
            f.write(f"# occsim density  DIV={div:g} "
                    f"mig_speed={base.mig_speed:g} territory={r_um:g} "
                    f"res={base.res:g} scaffold={base.scaffold:g} "
                    f"seed_frac={base.seed_frac:g} "
                    f"morph={nm} seeds={seeds} "
                    f"compartments={base.compartments} "
                    f"scaffold_geometry="
                f"{rel_to_root(base.scaffold_geometry)}\n")
            f.write("# dilated = territory radius R; skeleton = R set to one "
                    "pixel, so it moves with --res and curves at different "
                    "resolutions are not comparable.\n")
            # The saturation thresholds and the numbers quoted in the docs were
            # measured on UNIFORM seeding, which is what "occupancy versus cell
            # count" means with no initial condition confusing it. A clustered
            # run answers a different question and its curve sits well below
            # the uniform one at the same N, so say which this is rather than
            # leaving it to be inferred.
            if base.seed_frac < 1.0:
                f.write(f"# CLUSTERED seeding: cells start in the central "
                        f"{base.seed_frac:g} of the field side, so this is not "
                        f"the uniform-seeding density curve the saturation "
                        f"thresholds were measured on. Re-run with "
                        f"--seed-frac 1 for that one.\n")
            else:
                f.write("# Uniform seeding across the whole field.\n")
            if on_scaffold:
                f.write("# REPORTED = dilated_scaffold_*: the dilated curve "
                        "over the scaffold area. dilated_* and skeleton_* are "
                        "the same runs over the whole field.\n")
            w_ = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w_.writeheader()
            w_.writerows(rows)
        print(f"wrote {csv_path}")

    fig, ax = plt.subplots(figsize=(6.0, 4.2), facecolor="white")
    if on_scaffold:
        ax.errorbar(dens, scaf, yerr=scaf_e, fmt="-o", color=BLUE, lw=2,
                    markersize=4, capsize=3,
                    label=f"R = {r_um:g} um, on scaffold")
        ax.errorbar(dens, dil, yerr=dil_e, fmt="--o", color="#9aa5ad", lw=1.4,
                    markersize=3.5, capsize=2, markerfacecolor="white",
                    label=f"R = {r_um:g} um, whole field")
    else:
        ax.errorbar(dens, dil, yerr=dil_e, fmt="-o", color=BLUE, lw=2,
                    markersize=4, capsize=3,
                    label=f"R = {r_um:g} um (dilated territory)")
    ax.errorbar(dens, skel, yerr=skel_e, fmt="-s", color=RED, lw=2,
                markersize=4, capsize=3, markerfacecolor="white",
                label=f"R = {base.res:g} um (raw skeleton, whole field)")
    ax.set_xscale("log")
    ax.set_xlabel("projected areal density (cells/cm$^2$)")
    ax.set_ylabel(f"DIV{div:g} occupancy (%)")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=8.5, loc="lower right")
    _clean(ax)
    # skeleton curve is one pixel of radius, so it moves with --res and curves
    # at different resolutions aren't comparable
    fig.text(0.5, -0.06,
             f"DIV{div:g}, migration {base.mig_speed:g} um/h, "
             f"{len(seeds)} seed{'s' if len(seeds) != 1 else ''}, "
             f"mean +/- sem, res {base.res:g} um/px "
             f"(the skeleton curve moves with resolution).",
             ha="center", fontsize=8, color="#555")
    p = out / f"density_saturation_{_tag(base)}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"wrote {p}")
    return 0


def _seeding_vs_map(smap, field_um: float, seed_frac: float):
    """How much stiffer the seed box is than the map as a whole.

    A clustered run does not sample the map uniformly, and on a centrally
    concentrated field, which is what J0^2 gives, that decides whether the
    population feels the map at all. The docs' couplability tables were
    measured under uniform seeding, so a run at --seed-frac < 1 needs this
    ratio stated beside its result or the two cannot be read together.

    Returns (ratio, box_mean, map_mean), or None when seeding is uniform.
    """
    if seed_frac >= 1.0:
        return None
    import numpy as np
    m = np.asarray(smap)
    n = m.shape[0]
    comp_um = field_um / n
    half = max(1, int(round(0.5 * seed_frac * field_um / comp_um)))
    c = n // 2
    box = m[max(0, c - half):c + half, max(0, c - half):c + half]
    if box.size == 0:
        return None
    box_mean, map_mean = float(box.mean()), float(m.mean())
    if map_mean == 0.0:
        return None
    return box_mean / map_mean, box_mean, map_mean


def stiffness(base: Run, stiffness_path: Path, weights, seeds, div, grid,
              engine, out: Path, csv_path: Path | None,
              coupling: str = "gradient") -> int:
    """The stiffness map, plus evidence that coupling to it actually steers.

    Two things the other figures can't show: what a carried map looks like, and
    whether a coupling weight does anything in the direction it claims.

    Occupancy alone can't answer the second one. It falls under any strong
    bias, since funnelled cells sweep less area than a free walk. So an
    occupancy change proves the switch is live but says nothing about
    direction. The directional readout depends on which coupling is swept:

    - `gradient` sweeps `--stiff-weight` and reads coverage-weighted mean
      stiffness (each compartment's stiffness weighted by its covered pixels).
      Steering down-gradient leaves the network in softer material, so the
      readout falls monotonically with the weight. Only responds when the map
      has a coherent large-scale gradient: on a ring field a biased cell drops
      into the nearest local minimum and the panel is flat. That's the map,
      not the engine.
    - `speed` sweeps `--stiff-speed` and reads mean stiffness AT THE SOMA
      positions (the engine's SOMA_STIFFNESS line). The speed rule collects
      slow-moving somata, and with the default negative slope they linger on
      stiff material, so the readout RISES with the weight. It has to be
      soma-scale: the ~700 um dilated arbor averages over ring wavelengths and
      would smooth the accumulation out of any coverage-based readout. This is
      the panel for Bessel and ring dose maps.

    The shipped map is synthetic. This shows each mechanism works, it is not a
    biological result, and the caption says so.

    Args:
        base: Base run; the coupling weight, `seed`, `div` and the compartment
            fields are overridden.
        stiffness_path: Stiffness map CSV to carry and couple.
        weights: Weight values to sweep. Include 0, which is the uncoupled
            reference. Both signs for gradient, 0..1 for speed, where the
            direction lives in the slope sign instead.
        seeds: Seeds to average over at each weight.
        div: DIV to evaluate at.
        grid: Compartment grid side; must match the map.
        engine: Path to the compiled occupancy_engine binary.
        out: Output directory for the PNG.
        csv_path: CSV output path; falsy writes no CSV.
        coupling: "gradient" or "speed", see above.

    Returns:
        0 on success.

    Raises:
        RuntimeError: If the map's grid doesn't match `grid`, the weights are
            outside the coupling's range, or the sweep produces no rows.
    """
    import shutil
    import tempfile
    from concurrent.futures import ThreadPoolExecutor

    if coupling not in ("gradient", "speed"):
        raise RuntimeError(f"unknown coupling {coupling!r}")
    weights = list(weights)
    if coupling == "speed" and any(w < 0 or w > 1 for w in weights):
        raise RuntimeError(
            "speed-coupling weights must be in [0,1] - the direction is the "
            "slope's sign, not the weight's")

    sdivs, smaps, smeta = read_stiffness(stiffness_path)
    # Careful: this returns the DIV value (what `smaps` is keyed by), while
    # StiffnessMap::select in C++ returns an index. Easy to mix up.
    map_div = select_stiffness_div(sdivs, div)
    if map_div is None:
        raise RuntimeError(
            f"no stiffness map at or before DIV {div:g} (earliest is "
            f"DIV {sdivs[0]:g}) - nothing would be in effect")
    smap = smaps[map_div]
    if smap.shape != (grid, grid):
        raise RuntimeError(
            f"map grid {smap.shape[1]}x{smap.shape[0]} does not match "
            f"--grid {grid}; the engine would reject this pairing too")

    seeds = list(seeds)
    tmp = Path(tempfile.mkdtemp(prefix="occ_stiff_"))

    def one(job):
        # index-named scratch file: weights are arbitrary floats, so naming by
        # a formatted weight collided for values equal to 2 decimals
        index, w, s = job
        if coupling == "gradient":
            cpath = tmp / f"job{index}.csv"
            run = base.replace(div=div, seed=s, stiff_weight=w,
                               stiffness_map=stiffness_path,
                               compartments=grid, compartments_csv=cpath)
            occ = occupancy(run, engine)
            _, covered, _, _ = read_compartments(cpath)
            # coverage-weighted mean stiffness: WHERE the network ended up,
            # not how much it covered
            readout = float((covered * smap).sum() / covered.sum())
        else:
            run = base.replace(div=div, seed=s, stiff_speed=w,
                               stiffness_map=stiffness_path,
                               compartments=grid)
            occ, readout = occupancy_and_soma_stiffness(run, engine)
        return w, occ, readout

    jobs = [(i, w, s) for i, (w, s) in
            enumerate((w, s) for w in weights for s in seeds)]
    try:
        with ThreadPoolExecutor(max_workers=min(len(jobs), 8)) as ex:
            done = list(ex.map(one, jobs))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    readout_name = ("covwt_stiffness" if coupling == "gradient"
                    else "soma_stiffness")
    weight_name = "stiff_weight" if coupling == "gradient" else "stiff_speed"
    rows = []
    print(f"  {'weight':>7} {'occ %':>9} {'sem':>6} "
          f"{readout_name:>15} {'sem':>6}")
    for w in weights:
        occs = np.array([o for ww, o, _ in done if ww == w])
        wtd = np.array([x for ww, _, x in done if ww == w])
        n = len(occs)
        row = {
            weight_name: w,
            "occ_mean": occs.mean(),
            "occ_sem": (occs.std(ddof=1) / n ** 0.5) if n > 1 else 0.0,
            f"{readout_name}_mean": wtd.mean(),
            f"{readout_name}_sem":
                (wtd.std(ddof=1) / n ** 0.5) if n > 1 else 0.0,
            "n_seeds": n}
        rows.append(row)
        print(f"  {w:>7.2f} {row['occ_mean']:>9.3f} {row['occ_sem']:>6.3f} "
              f"{row[f'{readout_name}_mean']:>15.4f} "
              f"{row[f'{readout_name}_sem']:>6.4f}")
    if not rows:
        raise RuntimeError("stiffness sweep produced no rows")

    units = smeta.get("units", "?")
    ref = next((r for r in rows if r[weight_name] == 0), None)

    if csv_path:
        Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
        with open(csv_path, "w", newline="") as f:
            f.write(f"# occsim stiffness --coupling {coupling}  N={base.n} "
                    f"DIV={div:g} mig_speed={base.mig_speed:g} "
                    f"seed_frac={base.seed_frac:g} "
                    f"grid={grid}x{grid} map={Path(stiffness_path).name} "
                    f"map_div={map_div:g} units={units} "
                    f"seeds={seeds}\n")
            f.write("# SYNTHETIC map - demonstrates the coupling mechanism, "
                    "NOT a biological result.\n")
            _svm = _seeding_vs_map(smap, base.scaffold, base.seed_frac)
            if _svm is not None:
                _ratio, _box, _all = _svm
                f.write(f"# CLUSTERED seeding: cells start in the central "
                        f"{base.seed_frac:g} of the field side, where the map "
                        f"averages {_box:.4f} {units} against {_all:.4f} over "
                        f"the whole map - {_ratio:.2f}x. The couplability "
                        f"tables in docs/INTEGRATION_NOTES.md were measured "
                        f"under UNIFORM seeding, so a response here is not "
                        f"evidence the map would couple for a uniformly "
                        f"seeded population, and a flat curve is not evidence "
                        f"it would not.\n")
            if coupling == "gradient":
                f.write("# covwt_stiffness = sum(covered_px * stiffness) / "
                        "sum(covered_px), post-processed from the compartment "
                        "CSV; it cannot affect the simulation.\n")
                f.write("# +ve stiff_weight steers down-gradient (toward "
                        "soft), -ve up-gradient (toward stiff). Direction is "
                        "contested - see docs/INTEGRATION_NOTES.md.\n")
            else:
                f.write("# soma_stiffness = mean map value under the final "
                        "soma positions (engine SOMA_STIFFNESS line); "
                        "soma-scale because coverage smooths over the arbor.\n")
                f.write("# stiff_speed scales the step by the Elayan 2026 "
                        "adhesion sigmoid; the default auto slope is "
                        "negative (slower on stiff), so somata accumulate in "
                        "stiff material and the readout rises.\n")
                # Nothing here scores this map, so the caveat below is a
                # standing note rather than a computed verdict. It is true of
                # every J0^2 envelope and is the one thing a reader of this
                # file alone would otherwise miss.
                f.write("# Before quoting a trend below: on a J0^2 envelope "
                        "most of the map's value range sits in a small "
                        "central area, so a population mean barely moves. "
                        "Judge the map against docs/INTEGRATION_NOTES.md, "
                        "'The measured tables to judge against'.\n")
            w_ = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w_.writeheader()
            w_.writerows(rows)
        print(f"wrote {csv_path}")

    L = base.scaffold
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.3), facecolor="white")
    # the map panel's colorbar eats into the gap before the next axes; without
    # the extra width its label collides with the middle panel's y-label
    fig.subplots_adjust(wspace=0.42)

    im = ax[0].imshow(smap, extent=[-L / 2, L / 2, -L / 2, L / 2],
                      origin="lower", interpolation="nearest", cmap="magma")
    ax[0].set_xticks([]); ax[0].set_yticks([])
    for sp in ax[0].spines.values():
        sp.set_edgecolor("#ccc")
    ax[0].set_xlabel(f"carried map, DIV{map_div:g}", fontsize=9, color="#333")
    cb = fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.02)
    cb.set_label(f"stiffness ({units})", fontsize=8, color="#555")
    cb.ax.tick_params(labelsize=7)

    ws = [r[weight_name] for r in rows]
    ax[1].errorbar(ws, [r[f"{readout_name}_mean"] for r in rows],
                   yerr=[r[f"{readout_name}_sem"] for r in rows],
                   fmt="-o", color=BLUE, lw=1.8, markersize=4, capsize=3)
    if ref:
        ax[1].axhline(ref[f"{readout_name}_mean"], color="#bbb", lw=1, ls=":")
    ax[1].axvline(0.0, color="#bbb", lw=1, ls=":")
    if coupling == "gradient":
        ax[1].set_xlabel(
            "stiff_weight   (-ve toward stiff  ->  +ve toward soft)")
        ax[1].set_ylabel(f"coverage-weighted stiffness ({units})", fontsize=9)
    else:
        ax[1].set_xlabel("stiff_speed   (0 = uncoupled)")
        ax[1].set_ylabel(f"soma-position stiffness ({units})", fontsize=9)

    ax[2].errorbar(ws, [r["occ_mean"] for r in rows],
                   yerr=[r["occ_sem"] for r in rows],
                   fmt="-o", color=RED, lw=1.8, markersize=4, capsize=3)
    if ref:
        ax[2].axhline(ref["occ_mean"], color="#bbb", lw=1, ls=":")
    ax[2].axvline(0.0, color="#bbb", lw=1, ls=":")
    ax[2].set_xlabel(weight_name)
    # Caveat kept next to the panel it is about: occupancy drops under a
    # bias in EITHER direction, because
    # funnelled or slowed cells sweep less area, so this panel shows the
    # switch is live and nothing about which way it points.
    ax[2].set_ylabel("occupancy (%)  - not directional", fontsize=9)
    for a in ax[1:]:
        _clean(a)

    _svm = _seeding_vs_map(smap, base.scaffold, base.seed_frac)
    _seed_line = ""
    if _svm is not None:
        _ratio, _, _ = _svm
        _seed_line = (
            f"\nCells are seeded in the central {base.seed_frac:g} of the "
            f"field, where this map averages {_ratio:.2f}x its whole-map mean, "
            f"so the population does not sample it uniformly - "
            f"docs/INTEGRATION_NOTES.md, 'What makes a dose field usable "
            f"here'.")
    fig.text(0.5, -0.08,
             f"N={base.n}, DIV{div:g}, migration {base.mig_speed:g} um/h, "
             f"{len(seeds)} seeds, mean +/- sem; dotted marks the uncoupled "
             f"run ({weight_name} 0).\n"
             f"Map {Path(stiffness_path).name} is SYNTHETIC: the mechanism, "
             f"not a biological result. Reading the panels: "
             f"docs/RUNNING.md." + _seed_line,
             ha="center", fontsize=8, color="#555")
    suffix = "" if coupling == "gradient" else "_speed"
    p = out / f"stiffness{suffix}_coupling_n{base.n}_{_tag(base)}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"wrote {p}")
    return 0


# ---------------------------------------------------------------------------
def _edt_1d(f):
    """Squared 1D Euclidean distance transform (Felzenszwalb & Huttenlocher).

    A direct port of `distance_transform_1d` in src/distance_transform.cpp, so
    the distance field this module builds is the same one the engine steers
    down rather than an approximation of it.

    Args:
        f: Sequence of squared-distance seeds; 0 at a feature, large elsewhere.

    Returns:
        List of squared distances, same length as `f`.
    """
    n = len(f)
    INF = 1e20
    d = [0.0] * n
    v = [0] * n
    z = [0.0] * (n + 1)
    k = 0
    v[0] = 0
    z[0] = -INF
    z[1] = INF
    for q in range(1, n):
        s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2.0 * q - 2.0 * v[k])
        while s <= z[k]:
            k -= 1
            s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2.0 * q - 2.0 * v[k])
        k += 1
        v[k] = q
        z[k] = s
        z[k + 1] = INF
    k = 0
    for q in range(n):
        while z[k + 1] < q:
            k += 1
        d[q] = (q - v[k]) * (q - v[k]) + f[v[k]]
    return d


def _scaffold_distance(mask):
    """Distance from every compartment to the nearest scaffold one, in
    compartment units.

    Args:
        mask: (ny, nx) bool array, True where the compartment is scaffold.

    Returns:
        (ny, nx) float array of Euclidean distances; 0 on the scaffold itself.
    """
    import numpy as np
    ny, nx = mask.shape
    INF = 1e20
    col = np.zeros((ny, nx))
    for x in range(nx):
        col[:, x] = _edt_1d([0.0 if mask[y, x] else INF for y in range(ny)])
    out = np.zeros((ny, nx))
    for y in range(ny):
        out[y, :] = _edt_1d(list(col[y, :]))
    return np.sqrt(out)


def _unit_bias(field, sign):
    """Per-compartment unit steering vector down (or up) a field's gradient.

    Mirrors `cue_bias()` in src/simulation.cpp: a 3-point central difference,
    normalised to unit length. The engine normalises the magnitude away before
    blending by weight, which is exactly why comparing the two cues needs no
    conversion between their units.

    Args:
        field: (ny, nx) array to differentiate.
        sign: +1 steers down-gradient, -1 up-gradient.

    Returns:
        (bx, by, active): two (ny, nx) arrays of unit-vector components, and a
        bool array that is False where the field is locally flat.
    """
    import numpy as np
    ny, nx = field.shape
    gx = np.zeros((ny, nx))
    gy = np.zeros((ny, nx))
    for y in range(ny):
        for x in range(nx):
            xl, xh = max(x - 1, 0), min(x + 1, nx - 1)
            yl, yh = max(y - 1, 0), min(y + 1, ny - 1)
            gx[y, x] = (field[y, xh] - field[y, xl]) / max(xh - xl, 1)
            gy[y, x] = (field[yh, x] - field[yl, x]) / max(yh - yl, 1)
    mag = np.hypot(gx, gy)
    active = mag > 1e-12
    bx = np.where(active, -sign * gx / np.where(active, mag, 1.0), 0.0)
    by = np.where(active, -sign * gy / np.where(active, mag, 1.0), 0.0)
    return bx, by, active


# Regime labels, in the order they are drawn and reported.
REGIMES = ("mechanical only", "optical only", "competitive", "neither")
# How the two cues meet where both are present. Degrees between the unit
# vectors: under 60 they push the same way, over 120 they fight.
COMPETITION = ("reinforcing", "oblique", "opposing")


def _declares_synthetic(path) -> bool:
    """Whether a carried file says in its own comments that it is synthetic.

    Read from the file, not guessed from its name, so a real map dropped in
    loses the warning and a synthetic one keeps it wherever it came from.

    Args:
        path: A stiffness or geometry CSV.

    Returns:
        True if any comment line mentions SYNTHETIC or GENERATED.
    """
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#"):
            break
        if "SYNTHETIC" in line.upper() or "GENERATED" in line.upper():
            return True
    return False


def regime(geometry_path: Path, stiffness_path: Path, div: float, grid: int,
           field_um: float, seek_weight: float, stiff_weight: float,
           out: Path, csv_path: Path | None = None) -> int:
    """Classify the local physical regime per compartment: mechanical,
    optical, competitive, or neither.

    The collaboration asks for "a compact framework that identifies the local
    physical regime in which the neuronal interface operates ... according to
    the relative strengths and directions of the two spatial cues". This is
    that framework, computed from the two files a run already carries, on the
    grid the engine already forces them to share. It runs no simulation and
    can change no simulated number.

    **No unit conversion happens here, and that is the point.** Both cues reach the
    walk as a UNIT vector scaled by its weight (`cue_bias()` in
    src/simulation.cpp normalises the gradient magnitude away), so kPa per
    compartment never has to be weighed against micrometres per compartment.
    What is spatial is *presence* and *direction*; what is strength is the
    weight ratio, which the caller sets and which is reported alongside rather
    than folded into the classification.

    Presence follows the engine's own rules, so this describes the run rather
    than an idealisation of it: the scaffold cue is silent ON the scaffold
    (distance is zero there and the cell walks normally, which is the decision
    rule), and the stiffness cue is silent wherever the map is locally flat.

    Args:
        geometry_path: Scaffold geometry CSV (docs/SCAFFOLD_GEOMETRY_FORMAT.md).
        stiffness_path: Stiffness map CSV (docs/STIFFNESS_FORMAT.md).
        div: DIV to classify at; the map in effect is the v1 carry-forward one.
        grid: Compartments per side; both files must be on it.
        field_um: Field side, um.
        seek_weight: The run's --seek-weight, reported as the mechanical
            cue's strength.
        stiff_weight: The run's --stiff-weight; its sign sets the optical
            cue's direction and its magnitude that cue's strength. Note that
            --stiff-speed is deliberately not represented: it scales step
            LENGTH by local magnitude and has no direction, so it has no place
            in a classification about directions. A run coupled only through
            --stiff-speed has no directional optical cue to classify.
        out: Output directory for the figure.
        csv_path: Where to write the per-compartment classification, or None.

    Returns:
        0 on success.

    Raises:
        SystemExit: If the two files are not on the same grid as each other or
            as `grid`, or if no stiffness map is in effect at `div`.
    """
    import numpy as np

    mask = read_geometry(geometry_path)
    divs, maps, meta = read_stiffness(stiffness_path)
    sel = select_stiffness_div(divs, div)
    if sel is None:
        raise SystemExit(
            f"no stiffness map at or before DIV {div} in {stiffness_path} "
            f"(earliest is DIV {divs[0]}), so there is no optical cue to "
            f"classify against")
    stiff = maps[sel]
    if mask.shape != (grid, grid) or stiff.shape != (grid, grid):
        raise SystemExit(
            f"regime needs both files on the run's {grid}x{grid} grid; got "
            f"geometry {mask.shape[1]}x{mask.shape[0]} and stiffness "
            f"{stiff.shape[1]}x{stiff.shape[0]}")

    # Mechanical cue: down the distance-to-scaffold field, i.e. toward the
    # nearest scaffold compartment. Silent on the scaffold, as the engine is.
    dist = _scaffold_distance(mask)
    mx, my, m_active = _unit_bias(dist, 1.0)
    m_active = m_active & ~mask

    # Optical cue: along the stiffness gradient, sign as --stiff-weight. The
    # `> 0` test is the engine's own (src/simulation.cpp), so the direction
    # shown here is the direction a run would actually steer.
    ox, oy, o_active = _unit_bias(stiff, 1.0 if stiff_weight > 0 else -1.0)

    both = m_active & o_active
    cos = np.clip(mx * ox + my * oy, -1.0, 1.0)
    angle = np.degrees(np.arccos(np.where(both, cos, 1.0)))

    label = np.full((grid, grid), 3, dtype=int)      # neither
    label[m_active & ~o_active] = 0                  # mechanical only
    label[o_active & ~m_active] = 1                  # optical only
    label[both] = 2                                  # competitive

    comp = np.full((grid, grid), -1, dtype=int)
    comp[both & (angle < 60)] = 0
    comp[both & (angle >= 60) & (angle <= 120)] = 1
    comp[both & (angle > 120)] = 2

    total = grid * grid
    counts = {r: int((label == i).sum()) for i, r in enumerate(REGIMES)}
    sub = {c: int((comp == i).sum()) for i, c in enumerate(COMPETITION)}

    print(f"regime at DIV {div:g} (stiffness map DIV {sel:g}), "
          f"{grid}x{grid} compartments")
    for r in REGIMES:
        print(f"  {r:<16} {counts[r]:5d}  ({100 * counts[r] / total:5.1f} %)")
    if counts["competitive"]:
        print("  where both act, the two cues are:")
        for c in COMPETITION:
            print(f"    {c:<14} {sub[c]:5d}  "
                  f"({100 * sub[c] / counts['competitive']:5.1f} % of competitive)")
        print(f"    mean angle between them {angle[both].mean():.0f} deg")
    print(f"  cue strengths are the weights you set: mechanical "
          f"--seek-weight {seek_weight:g}, optical --stiff-weight "
          f"{stiff_weight:g}")

    if csv_path:
        Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
        with open(csv_path, "w", newline="") as f:
            # comments through f.write, like every other CSV here: csv.writer
            # would quote any line containing a comma
            f.write(f"# occsim regime  DIV={div:g} map_div={sel:g} "
                    f"grid={grid}x{grid} field={field_um:g} "
                    f"geometry={Path(geometry_path).name} "
                    f"stiffness={Path(stiffness_path).name} "
                    f"seek_weight={seek_weight:g} "
                    f"stiff_weight={stiff_weight:g}\n")
            f.write("# regime = which spatial cues are PRESENT in this "
                    "compartment. Both cues reach the walk as unit vectors "
                    "scaled by their weight, so no unit conversion is "
                    "involved and none is implied.\n")
            f.write("# the scaffold cue is silent ON the scaffold by the "
                    "decision rule; the stiffness cue is silent where the "
                    "map is locally flat.\n")
            f.write("# angle_deg = between the two unit cue vectors where "
                    "both act; -1 where they do not.\n")
            for q in (geometry_path, stiffness_path):
                if _declares_synthetic(q):
                    f.write(f"# {Path(q).name} declares itself SYNTHETIC or "
                            f"GENERATED - this demonstrates the "
                            f"classification, not a biological result.\n")
            w = csv.writer(f)
            w.writerow(["cx", "cy", "scaffold", "regime", "competition",
                        "angle_deg", "stiffness"])
            for cy in range(grid):
                for cx in range(grid):
                    w.writerow([
                        cx, cy, int(mask[cy, cx]), REGIMES[label[cy, cx]],
                        COMPETITION[comp[cy, cx]] if comp[cy, cx] >= 0 else "",
                        f"{angle[cy, cx]:.1f}" if both[cy, cx] else "-1",
                        f"{stiff[cy, cx]:g}"])
        print(f"wrote {csv_path}")
    return _regime_figure(label, comp, angle, both, mask, stiff, div, sel,
                          grid, field_um, seek_weight, stiff_weight,
                          geometry_path, stiffness_path, out, counts, total,
                          meta.get("units", "kPa"))


def _regime_figure(label, comp, angle, both, mask, stiff, div, sel, grid,
                   field_um, seek_weight, stiff_weight, geometry_path,
                   stiffness_path, out: Path, counts, total,
                   units: str = "kPa") -> int:
    """Draw the regime map beside the two cue fields it was derived from.

    Args:
        label: (grid, grid) int array indexing REGIMES.
        comp: (grid, grid) int array indexing COMPETITION, -1 where neither.
        angle: (grid, grid) degrees between the cue vectors.
        both: (grid, grid) bool, True where both cues act.
        mask: (grid, grid) bool scaffold mask.
        stiff: (grid, grid) the stiffness map in effect.
        div: DIV classified at.
        sel: DIV of the map actually in effect.
        grid: Compartments per side.
        field_um: Field side, um.
        seek_weight: Mechanical cue strength, for the caption.
        stiff_weight: Optical cue strength and direction, for the caption.
        geometry_path: Named in the caption.
        stiffness_path: Named in the caption.
        out: Output directory.
        counts: Per-regime compartment counts.
        total: grid*grid.
        units: The map's declared units, carried verbatim onto the colourbar
            rather than assumed, a map in Pa or MPa labels itself correctly.

    Returns:
        0 on success.
    """
    import numpy as np
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    ext = [-field_um / 2, field_um / 2, -field_um / 2, field_um / 2]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2))

    # the scaffold, i.e. where the mechanical cue comes from
    axes[0].imshow(mask, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=1.6)
    axes[0].set_title("scaffold (mechanical cue)", fontsize=10)

    im = axes[1].imshow(stiff, origin="lower", extent=ext, cmap="magma")
    axes[1].set_title(f"stiffness, map DIV {sel:g} (optical cue)", fontsize=10)
    fig.colorbar(im, ax=axes[1], fraction=0.046).set_label(
        f"stiffness ({units})", fontsize=9)

    colours = ["#2a6f97", "#c1121f", "#f0a202", "#e9e9e9"]
    cmap = ListedColormap(colours)
    axes[2].imshow(label, origin="lower", extent=ext, cmap=cmap,
                   norm=BoundaryNorm([0, 1, 2, 3, 4], cmap.N))
    axes[2].set_title(f"local regime at DIV {div:g}", fontsize=10)
    axes[2].legend(handles=[
        Patch(facecolor=colours[i],
              label=f"{r}  {100 * counts[r] / total:.0f} %")
        for i, r in enumerate(REGIMES)],
        loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=2,
        fontsize=8, frameon=False)

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    # A PNG travels without its data file, so where the inputs are synthetic
    # the figure has to say so itself.
    synth = [Path(q).name for q in (geometry_path, stiffness_path)
             if _declares_synthetic(q)]
    provenance = ("   SYNTHETIC/GENERATED: " + ", ".join(synth)
                  + " - not a biological result.") if synth else ""

    n_comp = counts["competitive"]
    detail = ""
    if n_comp:
        opp = int((comp == 2).sum())
        detail = (f" Mean angle where both act {angle[both].mean():.0f} deg, "
                  f"{opp} of {n_comp} OPPOSING.")
    fig.text(0.5, 0.015,
             "Cue PRESENCE per compartment; no simulation is run. Strength is "
             f"the weight you set (--seek-weight {seek_weight:g}, "
             f"--stiff-weight {stiff_weight:g})." + detail + "\n"
             f"geometry: {Path(geometry_path).name}   stiffness: "
             f"{Path(stiffness_path).name}" + provenance,
             ha="center", fontsize=7.5, color=DARK)
    # top < 1 leaves room for the panel titles, which clipped against the
    # figure edge before
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    Path(out).mkdir(parents=True, exist_ok=True)
    path = Path(out) / f"regime_g{grid}_div{div:g}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")
    return 0
