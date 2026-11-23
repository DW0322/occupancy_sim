# Building and running occupancy-sim

Every command and number below is tested against this checkout.

## Quickstart

```bash
./bootstrap.sh          # deps, build, tests, python tooling
occsim all --pool --n 25 --field 3000 --seed-frac 0.1 \
    --cohere-weight 0.15 --mig-persist 100 --mig-speed 10 --res 1.0 \
    --seek scaffold --quick
xdg-open results/field_n25_pool39.png
```

Those flags are **the committed operating point**, section 2 says what each one
is for and why the defaults are not a substitute. `--quick` is not part of it:
it runs 4 seeds over 4 DIVs instead of 20 over 8, so drop it to reproduce the
committed figures.

## 1. Build

```bash
./bootstrap.sh                  # one-time; SUDO=sudo ./bootstrap.sh unless you are root
# or, once dependencies are present:
cmake --preset default && cmake --build --preset default && ctest --preset default
```

- Needs cmake, a C++17 compiler, git, and the HDF5 **development** package.
  The full per-platform list, Windows included, is
  [README.md, Requirements](../README.md#requirements); that is the source of
  truth and this section is the short form.
- `bootstrap.sh` installs those, builds, tests and pip-installs the tooling, but
  never calls `sudo` itself. It is Linux/macOS only; it dispatches over
  apt/dnf/pacman/zypper/apk/brew and exits 1 on Windows. Build there with the
  three `cmake`/`ctest` commands above, adding `--config Release` / `-C Release`
  under Visual Studio.
- **Expect 6 tests passing in ~16 s.** Fewer means a partial build.
- There is also a `debug` preset, built into `build-debug/` with ASan and UBSan on.
  GCC/Clang only: its `-fsanitize` flags are not MSVC syntax.
- MorphIO is fetched and statically linked, so no `LD_LIBRARY_PATH` is needed.
- **[MEASURED]** cold build from a fresh tree: 83 s to fetch MorphIO v3.4.0 and
  build everything, then 5/5 tests. 2026-08-26, i5-7300U (2 cores), GCC on Arch.
  The suite was five at that date; `clustering` was added after, hence six now.

### If the MorphIO fetch cannot run

MorphIO is pinned to tag `v3.4.0` and upstream is archived, so the pin cannot
move. Behind a proxy or air-gapped, point the build at a local clone:

```bash
git clone --depth 1 --branch v3.4.0 https://github.com/BlueBrain/MorphIO.git /path/to/morphio
cmake --preset default -DFETCHCONTENT_SOURCE_DIR_MORPHIO=/path/to/morphio
```

## 2. The operating point

**The defaults are not a useful configuration.** At `--n 15 --mig-speed 0` on the
1500 um field, occupancy saturates near 100 % and stops responding to anything -
the engine warns when it does. The committed set instead seeds a central cluster
in a large field, which is what the collaboration reports doing in vitro:

| flag | why |
|---|---|
| `--pool` | 39 distinct morphologies rather than one cell repeated |
| `--n 25` | enough local density for the clustering readout to mean something |
| `--field 3000` | keeps global coverage far from the ceiling at that N |
| `--seed-frac 0.1` | the central cluster; **this is the flag that matters**, and see the warning below |
| `--mig-speed 10 --mig-persist 100` | only moving somata sweep fresh scaffold |
| `--seek scaffold` | the decision rule under test |
| `--cohere-weight 0.15` | cells steer toward the neighbours their arbors reach |
| `--res 1.0` | 1 um raster; the cost knob |

### What the clustered seeding does to the two arrangement readouts

`--seed-frac 0.1` puts every cell inside a 300 um box at the centre of a 3000 um
field. That is the initial condition the collaboration describes, and it is what
keeps occupancy off its ceiling, but it also sets the starting value of both
readouts that answer *where* cells are, so neither can be read as an absolute.

**`SOMA_SCAFFOLD` starts high because of where cells are seeded, not because of
the rule.** On the shipped 40 % lattice the central 2x2 compartments are 100 %
scaffold and the central 4x4 is 75 %, so the cluster lands on strands. The run
opens near 70 % on-scaffold whether or not `--seek scaffold` is set.

> Read the **gap** to the `seek=none` control on panel 4 of the sweep figure,
> never the level. `occsim all` runs that control automatically; `occsim sweep`
> needs `--control`.

**Clark-Evans is scored over the convex hull of the somata**, not the field.
Over the field, a 300 um cluster in a 3000 um field has a denominator 100x
too large, which forces the statistic below 1 no matter how the cells are
arranged, the committed run read 0.64, "clumped", where the same positions over
their own hull read 1.17. The engine names the region it used on the
`CLUSTERING` line (`soma hull`, `scaffold` under `--plate-scaffold`, or `field`
when a hull is degenerate), so a number can always be traced to its denominator.

**A hull Clark-Evans has its null above 1, and you have to read it against
that.** The hull is drawn *by* the points, so it undercovers the region they
came from, and the shortfall grows as N falls. Measured over 2000 uniform-random
trials:

| N | 3 | 6 | 10 | 15 | 25 | 40 | 60 | 120 |
|---|---|---|---|---|---|---|---|---|
| median CE of a random scatter | 3.39 | 1.95 | 1.54 | 1.35 | **1.23** | 1.15 | 1.10 | 1.06 |

So at the committed N=25 the null is 1.23, and the run's 1.17 is a random
arrangement, not a spread one. Compare runs at equal N, or against this row -
never against 1.0, and never across different N.

**The same applies to a carried stiffness map.** A clustered population does not
sample the map uniformly, and every coupling CSV and figure prints its seed
box's mean against the whole map's for that reason. See
[INTEGRATION_NOTES.md](INTEGRATION_NOTES.md), "What makes a dose field usable
here".


## 3. Engine flags

```bash
./build/occupancy_engine --morph data/morph.asc --n 15 --div 15
# = 81.7524 %      reveal only, no migration, pinned in the regression test
```

`occsim --help` lists everything; `./build/morph_stats <files...>` prints per-morphology
geometry.

| flag                                                          | meaning                                                                 | default    |
| ------------------------------------------------------------- | ----------------------------------------------------------------------- | ---------- |
| `--n` / `--div`                                               | cell count / days in vitro                                              | 15 / 15    |
| `--field`                                                     | overall field side, um (`--scaffold` is a legacy alias)                 | 1500       |
| `--territory`                                                 | dilation radius R, um                                                   | 17.2       |
| `--res`                                                       | raster resolution, um                                                   | 0.5        |
| `--compartments`                                              | compartment grid side                                                   | 30         |
| `--morph-dir`                                                 | sample from a pool of morphologies                                      | none       |
| `--mig-speed` / `--mig-persist`                               | soma migration, um/h and persistence                                    | **0** / 50 |
| `--seed-frac`                                                 | fraction of the field side cells are seeded across                      | 1.0        |
| `--occ-mode`                                                  | `cum` (swept) / `inst` (footprint)                                      | cum        |
| `--compartments-csv`                                          | write per-compartment occupancy CSV                                     | none       |
| **inputs**                                                    |                                                                         |            |
| `--stiffness-map`                                             | carry a map ([format](STIFFNESS_FORMAT.md))                             | none       |
| `--scaffold-geometry`                                         | which compartments are scaffold ([format](SCAFFOLD_GEOMETRY_FORMAT.md)) | none       |
| `--stiffness-check` / `--geometry-check`                      | validate a file and exit                                                | none       |
| **couplings**, all default off; each changes results once on |                                                                         |            |
| `--seek`                                                      | `none` / `scaffold` (needs a geometry)                                  | none       |
| `--seek-weight` / `--p-seek`                                  | bias strength / fraction of cells that seek                             | 0.3 / 1.0  |
| `--plate-scaffold`                                            | seed on and confine to the scaffold                                     | off        |
| `--stiff-weight`                                              | stiffness **gradient** bias, -1..1                                      | 0          |
| `--stiff-speed`                                               | stiffness **magnitude** scales step length, 0..1                        | 0          |
| `--speed-centre` / `--speed-slope`                            | speed sigmoid; slope sign is the direction                              | auto       |
| `--cohere-weight`                                             | steer toward neighbours the arbor reaches, 0..1                         | 0          |
| `--cohere-gate`                                               | `none` / `scaffold`; with `scaffold`, only cells OFF it cohere          | none       |
| `--cohere-check`                                              | report reach vs spacing (honours `--seed-frac`), then exit              | off        |

## 4. Python tooling

```bash
pip install -e python          # add --break-system-packages on a PEP 668 system
```

Needs numpy >= 2.0 and matplotlib >= 3.9; `bootstrap.sh` does this step for you.

| command | what you get |
|---|---|
| `all` | the committed figure set, one invocation |
| `sweep --figure` | occupancy + fold-change + rate vs DIV. **The headline figure** |
| `compartments` | per-compartment occupancy heatmap over DIV (+ the map, if given) |
| `field` | arbors drawn over the scaffold. **Look at this first** |
| `stiffness --coupling gradient\|speed` | the carried map, and evidence the coupling steers |
| `density` | occupancy vs cell count, a methods check on metric saturation |
| `regime` | which cue dominates per compartment |
| `make-geometry` / `make-stiffness` | write a scaffold mask / a synthetic map |

## 5. Reading the two analysis panels

### `occsim stiffness`

Three panels: the map in effect, the directional readout, and occupancy.

| `--coupling` | sweeps | middle panel |
|---|---|---|
| `gradient` (default) | `--stiff-weight`, both signs | coverage-weighted stiffness |
| `speed` | `--stiff-speed`, 0..1 | mean stiffness under the somata |

- **Read the middle panel for direction, not the right one.** Occupancy falls
  under any strong bias, funnelled or slowed cells sweep less area, so it shows
  the switch is live, not which way it points.
- **A flat middle panel is about the map, not the engine.** The two couplings
  fail on opposite map shapes; [INTEGRATION_NOTES.md](INTEGRATION_NOTES.md) has
  the decision table for picking one.
- **The map is a drop-in.** Both shipped maps are synthetic and labelled so on
  the figure. A real map on the same grid re-runs the figure unchanged; the CSV
  header records which file and which DIV block was in effect.

### `occsim regime`

Classifies every compartment as mechanically dominated, optically dominated,
competitive or neither. It needs both cues and refuses rather than classifying
one-sidedly. It runs no simulation, so it is fast and cannot move a number.

```bash
occsim regime --scaffold-area-frac 0.4 \
    --stiffness-map data/stiffness_example_synthetic.csv \
    --seek-weight 0.3 --stiff-weight 0.8
#   mechanical only      0  (  0.0 %)
#   optical only       357  ( 39.7 %)
#   competitive        540  ( 60.0 %)
```

- **No unit conversion is involved.** Both cues reach the walk as unit vectors
  scaled by their weight, so kPa is never weighed against micrometres, and
  rescaling a map changes nothing. Strength is just the weight you set, reported
  beside the classification rather than folded into it. Only the *sign* of
  `--stiff-weight` reaches the classification; it defaults to +0.3 here so the
  panel reads "if both cues were equally strong".
- **"Optical only" on the scaffold is correct.** The seek rule goes silent once a
  cell has arrived, so on-scaffold compartments feel only the stiffness cue.
- `--csv PATH` writes the per-compartment classification, with the angle between
  the two cue vectors wherever both act.

## 6. Tested and working examples

Everything referenced already exists in `data/`.

```bash
# the three stated inputs: field, scaffold area within it, cell count
occsim run --field 1500 --scaffold-area-frac 0.4 --n 3 --mig-speed 10
# 68.2266 %   on scaffold   [65.1814 % over the whole field]

# the decision rule: off-scaffold cells steer toward it, on-scaffold cells follow the strand
occsim run --n 15 --div 15 --mig-speed 10 --compartments 30 \
    --scaffold-geometry data/scaffold_geometry_example_lattice.csv \
    --seek scaffold --seek-weight 0.6
# 97.1828 %   on scaffold   [98.4753 % over the whole field]

# magnitude coupling, for a Bessel/ring dose map
occsim run --n 15 --div 15 --mig-speed 10 --compartments 30 \
    --stiffness-map data/stiffness_bessel_synthetic.csv --stiff-speed 0.8
# 99.7865 %   (whole field; no scaffold defined)

# both couplings at once
occsim run --n 15 --div 15 --mig-speed 10 --compartments 30 \
    --scaffold-geometry data/scaffold_geometry_example_lattice.csv --plate-scaffold \
    --stiffness-map data/stiffness_example_synthetic.csv --stiff-weight 0.5
# 86.5047 %   on scaffold   [84.3375 % over the whole field]

# author your own inputs
occsim make-geometry data/my_scaffold.csv \
    --shape lattice --area-frac 0.4 --grid 30 --field 1500
```

Validating an input file without running anything, which `occsim` does not expose:

```bash
./build/occupancy_engine --stiffness-check data/stiffness_example_synthetic.csv
./build/occupancy_engine --geometry-check data/scaffold_geometry_example_lattice.csv
# bare, these are format checks against what the file declares, so any
# well-formed file passes at its own size. Add --compartments / --field to
# check it against a particular run instead.
```


## 7. Validation against independent morphologies

The templates were compared against an independent population (Allen Cell Types,
113 cells) through the same `morph_stats` the simulation uses. Committed to
`results/validation/`. Typical on inter-neurite spacing, the metric that sets
`--territory` R, and larger on extent:

| metric | population median (IQR) | ours | in IQR |
|---|---|---|---|
| mean neurite spacing | 27 (18-43) | 28 | 82 % |
| max path distance | 578 (375-718) | 738 | 44 % |
| 2D arbor span | 550 (433-698) | 867 | 18 % |
| total projected length | 9078 (4293-14473) | 20654 | 23 % |
