# occupancy_sim

A 2D, time-resolved geometric simulation of cortical neuron growth on a microstructured alginate scaffold. It answers one question: **how much of a scaffold's surface gets covered by growing neurons, over days in culture?**

## Pipeline

1. Load a real 3D neuron reconstruction and project it to 2D, a top-down view onto the scaffold's print surface.
2. Seed N cell bodies across the field, each with a randomly rotated morphology sampled from a pool.
3. Over simulated days in vitro (DIV), reveal more of each cell's arbor outward from its soma on a fixed schedule. Somata can migrate by a persistent random walk, biased toward the scaffold or modulated by a loaded stiffness map.
4. Dilate the revealed centrelines by a territory radius to turn them into area, then accumulate coverage across all cells and DIVs.
5. Report occupancy as a percentage. The denominator is the scaffold area when a scaffold is defined and the whole field when one is not. Both numbers are always produced and the reported one is named.

The engine is a C++ library (`occupancy::core`) plus a thin CLI (`occupancy_engine`). Everything except argument parsing lives in the library, so the geometry is testable without going through `argv`. A Python package (`occsim`) drives the engine across many runs and renders the figures.

## Requirements

| | version | why |
|---|---|---|
| CMake | >= 3.16 (4.x fine) | build system |
| C++17 compiler | GCC 9+, Clang 10+, or MSVC 19.2x (VS 2019+) | the engine |
| Git | any | MorphIO is fetched at configure time, submodules included |
| HDF5 | the **development** package, headers, not just the runtime | MorphIO reads `.h5` through HighFive and will not configure without them |
| Python | >= 3.9 | the `occsim` driver and the figures |
| numpy / matplotlib | >= 2.0 / >= 3.9 | pinned in `python/pyproject.toml`; `pip` installs them for you |

MorphIO is **not** something you install first. The build fetches tag `v3.4.0`
and links it statically, so it needs github.com reachable at configure time and
a git that can fetch submodules. Upstream is archived read-only, so that pin can
never move; the offline route is a local clone plus
`-DFETCHCONTENT_SOURCE_DIR_MORPHIO=`, in [docs/RUNNING.md](docs/RUNNING.md) section 1.

### Linux / macOS

```bash
sudo apt install cmake g++ git libhdf5-dev      # Debian/Ubuntu
sudo dnf install cmake gcc-c++ git hdf5-devel   # Fedora/RHEL
sudo pacman -S cmake gcc git hdf5               # Arch
brew install cmake git hdf5                     # macOS
```

`./bootstrap.sh` does exactly that and then builds, tests and pip-installs the
tooling, `SUDO=sudo ./bootstrap.sh` unless you are root. It is a convenience for
these two platforms only; the manual route below is equivalent and works
everywhere.

### Windows

Install, in order:

1. **Visual Studio 2019+** with the *Desktop development with C++* workload, which
   brings MSVC, CMake and Ninja. Standalone CMake plus Build Tools works too.
2. **Git for Windows**.
3. **Python 3.9+** from python.org or the Store, with `pip` on `PATH`.
4. **HDF5 development libraries**, most easily through vcpkg:

```powershell
git clone https://github.com/microsoft/vcpkg
.\vcpkg\bootstrap-vcpkg.bat
.\vcpkg\vcpkg install hdf5:x64-windows
```

Then configure with the vcpkg toolchain so `find_package(HDF5)` resolves:

```powershell
cmake --preset default -DCMAKE_TOOLCHAIN_FILE=C:\path\to\vcpkg\scripts\buildsystems\vcpkg.cmake
cmake --build --preset default --config Release
ctest --preset default -C Release
```

### Python tooling (all platforms)

```bash
pip install -e python
```

That pulls numpy and matplotlib and puts the `occsim` command on `PATH`, so it
runs from any directory. On a distro whose Python is externally managed, either
use a virtualenv or add `--break-system-packages`. Without the install you can
still run it in place as `cd python && python -m occsim ...`.

## Quickstart

```bash
cmake --preset default
cmake --build --preset default
ctest --preset default          # expect 6 tests passing
```

```bash
occsim all --pool --n 25 --field 3000 --seed-frac 0.1 \
    --cohere-weight 0.15 --mig-persist 100 --mig-speed 10 --res 1.0 \
    --seek scaffold --quick
```

That writes five figures into `results/` plus their CSVs. `--quick` cuts the sweep to 4 seeds and 4 DIVs; drop it to reproduce the committed figures, which are 20 seeds over 8 DIVs. Look at `field_n25_pool39.png` first: a cluster of arbors spreading and drifting apart over days is the quickest way to catch a broken run.

```bash
./build/occupancy_engine --morph data/morph.asc --n 15 --div 15
#  = 81.7524 %      smoke test; the number is pinned
```

The offline route and every flag are in [docs/RUNNING.md](docs/RUNNING.md).

## Inputs

A run is specified by three things: the overall field, the scaffold area within it and the initial cell count.

| input | flag |
|---|---|
| overall area | `--field 1500` (side in um) |
| scaffold area within it | `--scaffold-area-frac 0.4` |
| initial cell count | `--n 3` |

`--scaffold-area-frac` builds the compartment mask and writes it next to the figures, so the geometry a run used is recoverable from its output.
`--scaffold-geometry` supplies a ready-made mask instead. Passing both is
refused as ambiguous.

The scaffold denominator is quantised to whole compartments, so at the 30x30 default its boundary is a 50 um staircase. `--compartments 150` costs about 2% more runtime.

## Couplings

Six optional couplings change how cells move. All default off, and a run with a file merely loaded is bit-identical to one without it.

| flag | what it does |
|---|---|
| `--seek scaffold` | the decision rule: a cell off the scaffold steers toward it, a cell on it follows the structure |
| `--plate-scaffold` | seed cells onto the scaffold and confine them to it |
| `--stiff-weight` | durotaxis along a stiffness **gradient** |
| `--stiff-speed` | local stiffness **magnitude** scales the step length |
| `--cohere-weight` | cells steer toward the neighbours their neurites reach |
| `--cohere-gate scaffold` | only cells off the scaffold cohere, so substrate and cell-cell adhesion compete |

The two stiffness couplings are not interchangeable: they suit opposite map shapes, so pick by the map's structure rather than trying both. [docs/INTEGRATION_NOTES.md](docs/INTEGRATION_NOTES.md) has the decision table and the measurements behind it.

Coverage cannot answer where cells sit relative to each other, so runs also report a Clark-Evans clustering statistic and `SOMA_SCAFFOLD`, the fraction of somata on the structure. Both are scored against where the cells actually are, not against the whole field: Clark-Evans uses the convex hull of the somata, and `SOMA_SCAFFOLD` is read against a `seek=none` control on the sweep figure, because a clustered seeding sets the absolute level of either. Which readout answers which question is [docs/RUNNING.md](docs/RUNNING.md) section 2.

## Usage

```bash
# the engine directly
./build/occupancy_engine --morph data/morph.asc --n 15 --div 15 --mig-speed 10

# the Python driver
occsim run     --n 15 --div 15
occsim sweep   --n 3 --mig-speed 10 --seeds 20 --figure
```

Every flag, every subcommand and a set of tested example commands are in [docs/RUNNING.md](docs/RUNNING.md), which I recommend you look at next.

## Layout

```
include/occupancy/   public headers, the documented API
src/                 library implementation
apps/                CLI front end, argument parsing only
tests/               core regression suite, run with ctest
python/occsim/       one package, one entry point (python -m occsim)
cmake/               find or fetch MorphIO
data/                morphologies, plus synthetic stiffness and geometry examples
                     (build more with `occsim make-stiffness` / `make-geometry`)
results/             CSVs and figures from the current run
docs/                see Documentation below
```

## Documentation

| doc | what is in it |
|---|---|
| [docs/RUNNING.md](docs/RUNNING.md) | Every flag, every subcommand, tested example commands and common traps. |
| [docs/INTEGRATION_NOTES.md](docs/INTEGRATION_NOTES.md) | The contract for whoever merges this with the code producing its input files. |
| [docs/STIFFNESS_FORMAT.md](docs/STIFFNESS_FORMAT.md) | Stiffness map file format. |
| [docs/SCAFFOLD_GEOMETRY_FORMAT.md](docs/SCAFFOLD_GEOMETRY_FORMAT.md) | Scaffold geometry file format. |

## Testing

```bash
ctest --preset default
```

Six tests, guarding the numbers directly: `regression`, `params`,
`distance_transform`, `compartments`, `stiffness` and `clustering`. A
`regression` failure means the output changed and is worth investigating before
anything else. Under a multi-config generator (Visual Studio) this is
`ctest --preset default -C Release`.

## Limitations

- **Growth is a reveal, not generative.** The arbor is a finished reconstruction uncovered over time. It does not elongate or branch, and a migrating soma carries its arbor rigidly rather than trailing a process behind it.
- **Migration on this scaffold is an extrapolation.** The scaffold paper reports adhesion and neurite outgrowth but no soma migration and no cell tracks. The cue's direction comes from the paper. The trajectory is this model's.
- **2D at 50 um compartments bounds what the scaffold can express.** A top-down mask cannot represent an internal channel that cells enter.

## References

- Markram et al. 2015, *Cell* 163:456-492. Ramaswamy et al. 2015, *Front
  Neural Circuits* 9:44. Reimann et al. 2015, *Front Comput Neurosci* 9:120.
  Blue Brain NMC Portal reconstructions, redistributed for non-commercial use
  under the portal's terms. Citation is a condition of use.
- Li et al. 2025, *Trends in Biotechnology* 43(2):447-458. The M-Alg/t-ZnO
  scaffold, its Day 5/8/15 milestones and the adhesion differential behind the
  scaffold-seek rule.
- Elayan, Mullet, Mohr, Scherer, Barros & Jornet 2026, *NANOCOM '26*. The
  optical and mechanical channel model that produces the stiffness maps and the
  source of the adhesion-gain sigmoid behind `--stiff-speed`.
- Ballesteros-Esteban et al. 2023, *Chaos, Solitons & Fractals* 173:113764.
  Cultured neurons migrate under tension along their own neurites and pull
  themselves into aggregates. The basis for `--cohere-weight`.
- Felzenszwalb & Huttenlocher 2012, *Theory of Computing* 8(19):415-428. The
  exact distance transform behind the dilation and the scaffold distance fields.
- Schnitzer 1993, *Physical Review E* 48(4):2553-2568. The steady-state density
  of a walker with space-dependent speed goes as 1/v. The basis for
  `--stiff-speed`.
- Donnelly 1978. The edge correction on the Clark-Evans statistic; see
  `include/occupancy/clustering.hpp` for how it is applied.
- Allen Cell Types Database, celltypes.brain-map.org. The independent comparison
  population for the morphology validation. Downloaded rather than
  redistributed.

## License

MIT, see [LICENSE](LICENSE). That covers the software only. The morphology data
under `data/` is third-party and keeps its own terms. LICENSE spells both out.
