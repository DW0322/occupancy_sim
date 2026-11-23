# Stiffness map format (v1)

How stiff the scaffold is, per compartment, over DIV.

An optional input. The engine loads it, validates it against the run and
displays it beside the occupancy heatmap. On its own it changes no number: two
default-off switches make it act on migration, `--stiff-weight` (gradient) and
`--stiff-speed` (magnitude).

## Example

```
# occupancy-sim stiffness map v1
grid,3,3
scaffold_um,1500
units,kPa
div_step,5
div,cx,cy,stiffness
5,0,0,1.20
5,1,0,1.20
5,2,0,1.15
5,0,1,0.90
5,1,1,0.45
5,2,1,0.95
5,0,2,1.10
5,1,2,1.05
5,2,2,1.20
10,0,0,1.20
10,1,0,1.10
10,2,0,1.00
10,0,1,0.80
10,1,1,0.30
10,2,1,0.85
10,0,2,1.05
10,1,2,0.95
10,2,2,1.15
```

## Metadata

Four required lines, any order, before the data header.

| line | meaning |
|---|---|
| `grid,NX,NY` | compartments per side |
| `scaffold_um,L` | side of the whole field, um |
| `units,U` | stiffness units, carried verbatim into figures (expected: `kPa`) |
| `div_step,N` | DIV cadence: days between consecutive maps, any integer >= 1 |

## Data rows

`div,cx,cy,stiffness`, one row per compartment per DIV.

| column | meaning | valid values |
|---|---|---|
| `div` | days in vitro this map applies from | finite, `>= 0` |
| `cx` | compartment column, 0-indexed | `0` to `NX-1` |
| `cy` | compartment row, 0-indexed | `0` to `NY-1` |
| `stiffness` | the value, in the declared units | finite, `>= 0` |

## Orientation

Top-down view. The field is a square of side `scaffold_um` centred on the
origin, so x and y each run `-L/2` to `+L/2`.

- `(cx=0, cy=0)` is the corner at `(-L/2, -L/2)`.
- `cx` increases toward `+x`, `cy` increases toward `+y`.
- In this project's figures (y up), `cy=0` is the **bottom** row. Source data
  with row 0 at the top, as in a typical image raster, needs the row index
  flipped when writing `cy`.

## Rules the loader enforces
- guide for AI generation of geometries

1. First non-blank line is exactly `# occupancy-sim stiffness map v1`.
2. All four metadata lines appear before the data header.
3. `div_step` is an integer >= 1. Any cadence the producing pipeline uses is
   accepted; it only has to be declared so rule 5 has a ruler.
4. Every DIV present carries the complete grid: all `NX*NY` compartments, each
   exactly once.
5. The DIV sequence matches `div_step` with no gaps. Under `div_step,5` from
   DIV 5, the file must contain 5, 10, 15 and so on up to wherever it ends. A
   missing DIV is rejected rather than carried forward from the previous one,
   which would otherwise be invisible.
6. Fields are taken verbatim: no spaces or tabs around a field, nothing
   trailing a number. `grid, 30, 30`, `div, cx, cy, stiffness` and
   `5,0,0,0.5abc` are all rejected. Comment lines are exempt.

## What the engine does with it

At run DIV `d` the map in effect is the one with the **largest `div <= d`**:
piecewise-constant carry-forward, no interpolation. Before the first map, no
map is in effect.

The file is validated against the run's `--compartments` grid and `--scaffold`
side, summarised on stdout, and drawn beside the occupancy heatmap by
`occsim compartments`. A run with a map loaded is bit-identical to the same run
without one.

Two flags couple it, and which one suits a given map is a property of the map:

- `--stiff-weight`, `-1..1`, biases migration along the map **gradient**; sign
  sets the direction, positive toward soft.
- `--stiff-speed`, `0..1`, lets local stiffness **magnitude** scale the step
  length instead.

They suit opposite map shapes. [INTEGRATION_NOTES.md](INTEGRATION_NOTES.md) has
the decision table and the measurements behind it; the implementation is
`Params::stiff_weight` and `Params::stiff_speed_weight` in
`include/occupancy/simulation.hpp`.

## Checking a file

```
./build/occupancy_engine --stiffness-check FILE
./build/occupancy_engine --stiffness-check FILE --compartments 30 --scaffold 1500
```

## See also

[SCAFFOLD_GEOMETRY_FORMAT.md](SCAFFOLD_GEOMETRY_FORMAT.md), which says which
compartments are scaffold at all rather than how stiff each one is.
