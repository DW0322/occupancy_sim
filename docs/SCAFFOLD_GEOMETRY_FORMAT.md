# Scaffold geometry format

Which compartments of the field are scaffold material and which are not.

An optional input. The engine loads it, validates it against the run and
displays it over the occupancy heatmap. On its own it changes no number: two
default-off switches make it act, `--plate-scaffold` and `--seek scaffold`.

## Generating one

The stated inputs for a run are the overall field area, the scaffold area
within it, and the cell count. The first and third are flags (`--field`,
`--n`); this file is the second, and `occsim make-geometry` writes it:

```bash
occsim make-geometry data/my_scaffold.csv \
    --shape lattice --area-frac 0.4 --grid 30 --field 1500
```

| shape | what it is |
|---|---|
| `lattice` | square grid of printed strands with open pores, the M-Alg scaffold (Li 2025, Fig 2D-G). Strand spacing is `--strand-pitch`, default 600 um |
| `disc` | centred disc; an area-matched control |
| `square` | centred square; an area-matched control |

The covered compartment count is `round(area_frac * grid^2)` exactly, whatever
the shape or grid, so `--area-frac` means what it says and the three shapes are
directly comparable at equal area. The command prints the realised fraction.

Companion to [STIFFNESS_FORMAT.md](STIFFNESS_FORMAT.md). That file says how
stiff each compartment is; this one says which compartments are scaffold at
all. They are separate files because geometry is fixed for a given scaffold
while stiffness varies with DIV.

## Example

```
# occupancy-sim scaffold geometry v1
grid,3,3
scaffold_um,1500
cx,cy,scaffold
0,0,1
1,0,1
2,0,0
0,1,1
1,1,1
2,1,0
0,2,0
1,2,0
2,2,0
```

## Metadata

Two required lines, either order, before the data header.

| line | meaning |
|---|---|
| `grid,NX,NY` | compartments across and down |
| `scaffold_um,L` | side of the whole field, um |

## Data rows

`cx,cy,scaffold`, one row per compartment.

| column | meaning | valid values |
|---|---|---|
| `cx` | compartment column, 0-indexed | `0` to `NX-1` |
| `cy` | compartment row, 0-indexed | `0` to `NY-1` |
| `scaffold` | is this compartment scaffold | `1` yes, `0` no |

Only the third column is restricted to 0/1.

## Orientation

Top-down view. The field is a square of side `scaffold_um` centred on the
origin, so x and y each run `-L/2` to `+L/2`.

- `(cx=0, cy=0)` is the corner at `(-L/2, -L/2)`.
- `cx` increases toward `+x`, `cy` increases toward `+y`.
- In this project's figures (y up), `cy=0` is the **bottom** row. Source data
  with row 0 at the top, as in a typical image raster, needs the row index
  flipped when writing `cy`.

## Checking a file

```
./build/occupancy_engine --geometry-check FILE
./build/occupancy_engine --geometry-check FILE --compartments 30 --field 1500
```

## What the engine does with it

The file is validated against the run's `--compartments` grid and `--field`
side, summarised on stdout, and drawn as an overlay by `occsim compartments`,
which fades every non-scaffold compartment under a white wash so the scaffold
region stands out.

Two occupancy numbers are then reported: a scaffold-restricted percentage
measured only over scaffold compartments, and the whole-field percentage.
Neither replaces the other, the engine's `= <num> %` headline and the
unprefixed CSV columns stay whole-field, so every pinned result stays
comparable, but the scaffold number is the one `occsim` reports and plots,
because coverage of the scaffold is what the collaboration asked for.

Two default-off switches couple it:

| flag | effect |
|---|---|
| `--plate-scaffold` | cells are seeded onto the scaffold and a migration step that would leave it is cancelled |
| `--seek scaffold` | a cell **off** the scaffold steers toward the nearest scaffold compartment; **on** it, it keeps its normal walk |

`--seek scaffold` is the decision rule between seeking the scaffold and
continuing normal migration. `--seek-weight` sets how hard the bias competes
with the walk and `--p-seek` what fraction of cells apply it at all, so the
rule is a swept hypothesis rather than a hard switch. The steering direction
comes from an exact Euclidean distance transform of this mask, so it follows
whatever geometry is in the file rather than an idealised lattice.
