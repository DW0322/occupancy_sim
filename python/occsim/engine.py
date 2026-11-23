from __future__ import annotations

import os
import re
import subprocess
import sys
import threading
from dataclasses import dataclass, asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _find_engine() -> Path:
    """Locate the built engine across the layouts the presets can produce.

    Single-config generators (Make and Ninja, the default everywhere but MSVC)
    put the binary straight in build/. Visual Studio is multi-config and puts
    it in build/Release/ or build/Debug/ instead, and names it .exe. Probing
    rather than branching on sys.platform keeps a Ninja build on Windows and
    an MSVC build working through the same path.

    Returns:
        The first binary found. If none is, the plain build/ path, so the
        not-built error in measure() names the location people expect.
    """
    for d in (ROOT / "build", ROOT / "build" / "Release",
              ROOT / "build" / "Debug"):
        for name in ("occupancy_engine", "occupancy_engine.exe"):
            if (d / name).is_file():
                return d / name
    return ROOT / "build" / "occupancy_engine"


DEFAULT_ENGINE = _find_engine()
DEFAULT_MORPH = ROOT / "data" / "morph.asc"
# The morphology pool, sorted: the sampled index maps into it in order, so the
# order has to be reproducible. morphologies_balanced/ is 39 files, one per
# distinct cell, Blue Brain ships 1-7 rescaled "clones" per parent, and
# sampling over all of them would weight parents by clone count, a sampling
# artefact. data/morphologies/ keeps that full set so the selection stays
# reproducible; nothing here reads it (the tests use data/morph.asc). Reach it
# with --morph-dir.
MORPH_DIR = ROOT / "data" / "morphologies_balanced"


# Engine warnings go to stderr, which subprocess captures and would otherwise
# swallow. Forwarded here, deduplicated: a 20-seed sweep hits the same condition
# 20 times and one copy says everything.
_warned: set = set()
_warn_lock = threading.Lock()


def rel_to_root(p) -> str:
    if p is None:
        return "None"
    q = Path(p).resolve()
    try:
        return str(q.relative_to(ROOT))
    except ValueError:
        return str(q)


def _forward_warnings(stderr: str) -> None:
    #Re-emit distinct engine warnings on my own stderr.
    for line in stderr.splitlines():
        line = line.strip()
        if not line.startswith(("warning:", "note:")):
            continue
        key = re.sub(r"[0-9]+(?:\.[0-9]+)?", "#", line)
        with _warn_lock:
            if key in _warned:
                continue
            _warned.add(key)
        print(f"engine {line}", file=sys.stderr)


def pool(directory: Path = MORPH_DIR) -> list:
    """Every .asc/.swc/.h5 directly in a directory, sorted, not recursive."""
    return sorted(p for p in Path(directory).iterdir()
                  if p.suffix.lower() in (".asc", ".swc", ".h5"))


POOL = pool() if MORPH_DIR.is_dir() else []

_OCC = re.compile(r"=\s*([0-9.]+)\s*%")
# the soma-scale stiffness readout, printed whenever a map is carried
_SOMA_STIFF = re.compile(r"^SOMA_STIFFNESS mean ([0-9.eE+-]+) ", re.MULTILINE)
# Coverage restricted to scaffold compartments, the brief's "area of occupancy
# on top of the scaffold". Printed whenever a geometry is carried, on a line
# kept off the "= <num> %" pattern _OCC matches, so the two denominators cannot
# be confused.
_OCC_SCAFFOLD = re.compile(
    r"^OCCUPANCY_SCAFFOLD_RESTRICTED ([0-9.eE+-]+) pct", re.MULTILINE)
# Where the somata ended up relative to the scaffold: the closest thing to Li
# 2025's presence ON the structure, and the readout the seek rule is judged on.
# Coverage moves the OTHER way under guidance, since cells spreading along a
# strand sweep less fresh area than free walkers.
_SOMA_SCAFFOLD = re.compile(
    r"^SOMA_SCAFFOLD on_scaffold ([0-9.eE+-]+) pct \([0-9]+ of [0-9]+ "
    r"somata\), mean distance to scaffold ([0-9.eE+-]+) um",
    re.MULTILINE)
_CLUSTERING = re.compile(
    r"^CLUSTERING clark_evans ([0-9.eE+-]+) \(mean_nn ([0-9.eE+-]+) um",
    re.MULTILINE)


@dataclass
class Run:
    """A single engine invocation. Field names match the --flags."""

    # one path or a list of them; with several, each cell samples from the pool
    morph: object = DEFAULT_MORPH
    n: int = 15
    scaffold: float = 1500.0
    territory: float = 17.2
    res: float = 0.5
    div: float = 15.0
    seed: int = 1
    seed_frac: float = 1.0
    mig_speed: float = 0.0
    mig_persist: float = 50.0
    dt: float = 1.0
    occ_mode: str = "cum"
    seek: str = "none"
    seek_weight: float = 0.3
    p_seek: float = 1.0
    pgm: Path | None = None
    segments: Path | None = None
    compartments: int | None = None
    compartments_csv: Path | None = None
    # Carried inputs, validated by the engine. Each only couples through the
    # switches below: stiffness_map via stiff_weight/stiff_speed, geometry via
    # plate_scaffold or seek="scaffold".
    stiffness_map: Path | None = None
    scaffold_geometry: Path | None = None
    # the coupling switches, default off. See docs/INTEGRATION_NOTES.md.
    plate_scaffold: bool = False
    stiff_weight: float = 0.0
    # magnitude (speed) coupling; centre/slope None = engine auto-calibrates
    # from the maps' own spread
    stiff_speed: float = 0.0
    speed_centre: float | None = None
    speed_slope: float | None = None
    # cell-cell cohesion; the radius is derived from the reveal, so this is the
    # only knob. Whether it can act is a density question, --cohere-check.
    cohere_weight: float = 0.0
    # "scaffold" makes cohesion the complement of seek: only cells off it
    # cohere. Off by default, so carrying a geometry stays inert.
    cohere_gate: str = "none"

    def argv(self, engine: Path = DEFAULT_ENGINE) -> list[str]:
        """The engine command line for this run, ready for subprocess.run."""
        morphs = self.morph if isinstance(self.morph, (list, tuple)) else [self.morph]
        a = [rel_to_root(engine)]
        for m in morphs:
            a += ["--morph", rel_to_root(m)]
        for flag, val in (
            ("--n", self.n), ("--scaffold", self.scaffold),
            ("--territory", self.territory), ("--res", self.res),
            ("--div", self.div), ("--seed", self.seed),
            ("--seed-frac", self.seed_frac), ("--mig-speed", self.mig_speed),
            ("--mig-persist", self.mig_persist), ("--dt", self.dt),
            ("--occ-mode", self.occ_mode), ("--seek", self.seek),
            ("--seek-weight", self.seek_weight), ("--p-seek", self.p_seek),
        ):
            a += [flag, str(val)]
        if self.pgm:
            a += ["--pgm", rel_to_root(self.pgm)]
        if self.segments:
            a += ["--segments", rel_to_root(self.segments)]
        if self.compartments is not None:
            a += ["--compartments", str(self.compartments)]
        if self.compartments_csv:
            a += ["--compartments-csv", rel_to_root(self.compartments_csv)]
        if self.stiffness_map:
            a += ["--stiffness-map", rel_to_root(self.stiffness_map)]
        if self.plate_scaffold:
            a += ["--plate-scaffold"]
        if self.stiff_weight:
            a += ["--stiff-weight", str(self.stiff_weight)]
        if self.stiff_speed:
            a += ["--stiff-speed", str(self.stiff_speed)]
        if self.cohere_weight:
            a += ["--cohere-weight", str(self.cohere_weight)]
        if self.cohere_gate != "none":
            a += ["--cohere-gate", self.cohere_gate]
        if self.speed_centre is not None:
            a += ["--speed-centre", str(self.speed_centre)]
        if self.speed_slope is not None:
            a += ["--speed-slope", str(self.speed_slope)]
        if self.scaffold_geometry:
            a += ["--scaffold-geometry", rel_to_root(self.scaffold_geometry)]
        return a

    def replace(self, **kw) -> Run:
        """A copy with some fields overridden; the original is untouched."""
        d = asdict(self)
        d.update(kw)
        return Run(**d)


@dataclass(frozen=True)
class Measurement:
    """Everything one engine invocation reports back.

    Two occupancies, never one: `whole_field` is always there, `scaffold` is
    None when the run carried no geometry. Both are kept and `reported` picks,
    because every pinned value is whole-field and a run that quietly swapped
    denominators would make old and new outputs incomparable.
    """

    whole_field: float
    scaffold: float | None = None
    soma_stiffness: float | None = None
    # Clark-Evans over the final somata: <1 clustered, 1 random, >1 dispersed.
    # None when the run had under two cells, so there was no neighbour to
    # measure against.
    clark_evans: float | None = None
    mean_nn_um: float | None = None
    # Percent of somata sitting in a scaffold compartment, and their mean
    # distance to the nearest one. None when the run carried no geometry.
    soma_scaffold_pct: float | None = None
    soma_scaffold_dist_um: float | None = None

    @property
    def reported(self) -> float:
        """The headline occupancy: scaffold-restricted when a scaffold was
        defined, whole-field otherwise. Decided in one place so the CSVs,
        figures and CLI cannot drift apart over the denominator.
        """
        return self.whole_field if self.scaffold is None else self.scaffold

    @property
    def denominator(self) -> str:
        """Name of the area `reported` is divided by, for labels and headers."""
        return "whole field" if self.scaffold is None else "scaffold"


@dataclass(frozen=True)
class SeedStats:
    """Seed-averaged occupancy, per denominator.

    Each field is (mean, sd, sem, values) or None. `scaffold` is None when the
    runs carried no geometry. Both come out of one pass over the seeds, so
    asking for the scaffold number costs no extra engine runs.
    """

    whole_field: tuple
    scaffold: tuple | None = None
    # Clark-Evans, averaged the same way; None under two cells. Carried here so
    # arrangement comes out of the same pass over seeds as occupancy.
    clark_evans: tuple | None = None
    # Somata on the scaffold, percent; None with no geometry. The readout the
    # seek rule is judged on, from the same pass rather than a second run set.
    soma_scaffold: tuple | None = None
    soma_scaffold_dist: tuple | None = None

    @property
    def reported(self) -> tuple:
        """The (mean, sd, sem, values) block matching Measurement.reported."""
        return self.whole_field if self.scaffold is None else self.scaffold

    @property
    def denominator(self) -> str:
        """Name of the area `reported` is divided by, for labels and headers."""
        return "whole field" if self.scaffold is None else "scaffold"


def _env() -> dict:
    """Environment for the engine subprocess.

    Returns:
        os.environ with /usr/local/lib prepended to LD_LIBRARY_PATH, which a
        system-linked MorphIO needs and the static-linked default ignores.
    """
    return dict(os.environ, LD_LIBRARY_PATH="/usr/local/lib:"
                + os.environ.get("LD_LIBRARY_PATH", ""))


def measure(run: Run, engine: Path = DEFAULT_ENGINE) -> Measurement:
    """Run the engine once and parse everything it reports.

    The single place engine stdout becomes numbers, so a reporting line can
    only be parsed one way.

    Returns a Measurement; `scaffold` and `soma_stiffness` are filled in when
    the run carried a geometry or a map.

    Raises:
        FileNotFoundError: If `engine` isn't built.
        RuntimeError: If the engine fails, prints no occupancy line, or carried
            a geometry without reporting the restricted number.
    """
    if not Path(engine).exists():
        raise FileNotFoundError(
            f"engine not built at {engine}. Run: cmake --preset default "
            "&& cmake --build --preset default")
    # cwd=ROOT so the ROOT-relative paths in argv resolve, whatever
    # directory the caller is in (see rel_to_root).
    p = subprocess.run(run.argv(engine), capture_output=True, text=True,
                       env=_env(), cwd=ROOT)
    if p.returncode != 0:
        # the engine already prefixes its own diagnostics with "error: ", and
        # so does the CLI that prints this, so drop one of the two
        msg = p.stderr.strip().removeprefix("error: ")
        raise RuntimeError(f"engine failed ({p.returncode}): {msg}")
    _forward_warnings(p.stderr)
    m = _OCC.search(p.stdout)
    if not m:
        raise RuntimeError(f"no occupancy in engine output: {p.stdout!r}")
    scaffold = _OCC_SCAFFOLD.search(p.stdout)
    # A run that asked for a scaffold and got no restricted number back would
    # fall through to the whole-field denominator, which is the exact swap this
    # reporting path exists to prevent. Fail instead.
    if run.scaffold_geometry and not scaffold:
        raise RuntimeError(
            f"carried scaffold geometry {run.scaffold_geometry} but the engine "
            f"reported no OCCUPANCY_SCAFFOLD_RESTRICTED line, so there is no "
            f"scaffold-restricted number to report: {p.stdout!r}")
    soma = _SOMA_STIFF.search(p.stdout)
    clust = _CLUSTERING.search(p.stdout)
    on_scaf = _SOMA_SCAFFOLD.search(p.stdout)
    # The engine prints SOMA_SCAFFOLD whenever a geometry is loaded, coupled or
    # not. A run that carried one and reported none did not do what it was
    # asked. Same check as the restricted-occupancy one above.
    if run.scaffold_geometry and not on_scaf:
        raise RuntimeError(
            f"carried scaffold geometry {run.scaffold_geometry} but the engine "
            f"reported no SOMA_SCAFFOLD line: {p.stdout!r}")
    return Measurement(
        whole_field=float(m.group(1)),
        scaffold=float(scaffold.group(1)) if scaffold else None,
        soma_stiffness=float(soma.group(1)) if soma else None,
        clark_evans=float(clust.group(1)) if clust else None,
        mean_nn_um=float(clust.group(2)) if clust else None,
        soma_scaffold_pct=float(on_scaf.group(1)) if on_scaf else None,
        soma_scaffold_dist_um=float(on_scaf.group(2)) if on_scaf else None)


def occupancy(run: Run, engine: Path = DEFAULT_ENGINE) -> float:
    """Whole-field occupancy in percent, whatever the run carries, the
    denominator every pinned number uses. Use measure() for the
    scaffold-restricted number or both.
    """
    return measure(run, engine).whole_field


def occupancy_and_soma_stiffness(run: Run, engine: Path = DEFAULT_ENGINE):
    """Run the engine once; return (occupancy %, mean soma stiffness).

    The engine's SOMA_STIFFNESS line: mean map value under the final somata.
    This is what resolves the speed coupling on ring fields, where the ~700 um
    arbor smooths coverage-weighted stiffness flat.

    Returns (occupancy percent, soma stiffness in the map's units); `run` must
    carry a stiffness_map.

    Raises:
        RuntimeError: If the engine fails, or prints no SOMA_STIFFNESS line
            (no map carried, or no map in effect at this DIV).
    """
    m = measure(run, engine)
    if m.soma_stiffness is None:
        raise RuntimeError(
            "expected a SOMA_STIFFNESS line - the run carried no stiffness "
            f"map, or no map was in effect at DIV {run.div}")
    return m.whole_field, m.soma_stiffness


def measure_many(runs, engine: Path = DEFAULT_ENGINE, jobs: int | None = None):
    """Measurements for several runs, in parallel.

    Each invocation is its own process with its own seed, so scheduling cannot
    change a number: serial and parallel are bit-identical.

    `jobs` defaults to os.cpu_count(); 1 goes serial. Returns Measurements in
    the order of `runs`.
    """
    runs = list(runs)
    if jobs is None:
        jobs = os.cpu_count() or 1
    jobs = max(1, min(jobs, len(runs) or 1))
    if jobs == 1 or len(runs) < 2:
        return [measure(r, engine) for r in runs]
    from concurrent.futures import ThreadPoolExecutor
    # Threads, not processes: each task blocks in subprocess.run, so the GIL is
    # released and the work happens in the children anyway.
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        return list(ex.map(lambda r: measure(r, engine), runs))


def _stats(vals):
    """(mean, sd, sem, values); sd and sem are 0.0 for fewer than 2 values."""
    vals = list(vals)
    n = len(vals)
    mean = sum(vals) / n
    if n < 2:
        return mean, 0.0, 0.0, vals
    sd = (sum((v - mean) ** 2 for v in vals) / (n - 1)) ** 0.5
    return mean, sd, sd / n ** 0.5, vals


def seed_stats(run: Run, seeds, engine: Path = DEFAULT_ENGINE,
               jobs: int | None = None) -> SeedStats:
    """Both occupancies, averaged over seeds, from one pass.

    `run.seed` is overridden per seed. Returns a SeedStats, whose `scaffold`
    is None when the run carried no geometry.
    """
    seeds = list(seeds)
    ms = measure_many([run.replace(seed=s) for s in seeds], engine, jobs)
    scaffold = ([m.scaffold for m in ms] if ms and ms[0].scaffold is not None
                else None)
    clumping = ([m.clark_evans for m in ms]
                if ms and ms[0].clark_evans is not None else None)
    on_scaf = ([m.soma_scaffold_pct for m in ms]
               if ms and ms[0].soma_scaffold_pct is not None else None)
    scaf_dist = ([m.soma_scaffold_dist_um for m in ms]
                 if ms and ms[0].soma_scaffold_dist_um is not None else None)
    return SeedStats(whole_field=_stats([m.whole_field for m in ms]),
                     scaffold=_stats(scaffold) if scaffold else None,
                     clark_evans=_stats(clumping) if clumping else None,
                     soma_scaffold=_stats(on_scaf) if on_scaf else None,
                     soma_scaffold_dist=(_stats(scaf_dist) if scaf_dist
                                         else None))


def seed_mean(run: Run, seeds, engine: Path = DEFAULT_ENGINE,
              jobs: int | None = None):
    """Whole-field occupancy averaged over seeds: (mean, sd, sem, values)."""
    return seed_stats(run, seeds, engine, jobs).whole_field


def read_pgm(path):
    """The engine's binary PGM (P5) coverage image, as (height, width) uint8."""
    import numpy as np
    with open(path, "rb") as f:
        assert f.readline().strip() == b"P5"
        w, h = map(int, f.readline().split())
        f.readline()
        return np.frombuffer(f.read(), dtype=np.uint8).reshape(h, w)


def read_compartments(path):
    #Read a per-compartment occupancy grid from the engine's CSV.

    import numpy as np
    header, rows = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                header.append(line[1:].strip())
            elif line and not line.startswith("cx,"):
                rows.append(line.split(","))
    if not rows:
        raise RuntimeError(f"no compartment rows in {path}")
    nx = max(int(r[0]) for r in rows) + 1
    ny = max(int(r[1]) for r in rows) + 1
    grid = np.zeros((ny, nx))
    covered = np.zeros((ny, nx), dtype=np.int64)
    total = np.zeros((ny, nx), dtype=np.int64)
    for cx, cy, cov, tot, occ in rows:
        grid[int(cy), int(cx)] = float(occ)
        covered[int(cy), int(cx)] = int(cov)
        total[int(cy), int(cx)] = int(tot)
    return grid, covered, total, header


def read_stiffness(path):
    #Read stiffness maps for DISPLAY only.

    import numpy as np
    meta, rows = {}, []
    nx = ny = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            t = line.split(",")
            if t[0] == "grid":
                nx, ny = int(t[1]), int(t[2])
            elif t[0] == "scaffold_um":
                meta["scaffold_um"] = float(t[1])
            elif t[0] == "units":
                meta["units"] = t[1]
            elif t[0] == "div_step":
                meta["div_step"] = int(t[1])
            elif t[0] == "div" and t[1] == "cx":
                continue
            else:
                rows.append((float(t[0]), int(t[1]), int(t[2]), float(t[3])))
    maps = {}
    for d, cx, cy, v in rows:
        maps.setdefault(d, np.zeros((ny, nx)))[cy, cx] = v
    return sorted(maps), maps, meta


def read_geometry(path):
    #Read a scaffold geometry mask for DISPLAY only.

    import numpy as np
    nx = ny = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            t = line.split(",")
            if t[0] == "grid":
                nx, ny = int(t[1]), int(t[2])
            elif t[0] == "scaffold_um" or (t[0] == "cx" and t[1] == "cy"):
                continue
            else:
                rows.append((int(t[0]), int(t[1]), int(t[2])))
    m = np.zeros((ny, nx), dtype=bool)
    for cx, cy, v in rows:
        m[cy, cx] = bool(v)
    return m


def select_stiffness_div(divs, run_div):
    hits = [d for d in divs if d <= run_div]
    return hits[-1] if hits else None


def read_segments(path) -> dict:
    import numpy as np
    cells: dict[int, list] = {}
    with open(path) as f:
        next(f)
        for line in f:
            c, x0, y0, x1, y1 = line.split(",")
            cells.setdefault(int(c), []).append(
                [[float(x0), float(y0)], [float(x1), float(y1)]])
    return {c: np.array(v) for c, v in cells.items()}


if __name__ == "__main__": 
    sys.exit(
        "engine.py only defines the engine wrapper - there is nothing to run "
        "here.\n\nUse the command line entry point instead:\n"
        f"    cd {Path(__file__).resolve().parents[1]}\n"
        "    python -m occsim run --n 15 --div 15\n"
        "    python -m occsim --help\n")
