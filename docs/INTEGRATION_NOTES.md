# Integration notes

For whoever is wrapping the optical and photochemical model around this
simulation. This document is the contract between the two halves.

You produce per-compartment stiffness over DIV as CSV. This side consumes it and
reports how much of a scaffold gets covered by growing neurons; it never models
optics or photochemistry. **The interface between us is the file format, not the
physics.** Everything below is either that format, or something measured here
that changes what a usable field looks like.

## What we need from you

Two things, and only one of them is a file.

| what | form | required? |
|---|---|---|
| stiffness over DIV | a CSV in the v1 format below, on an N x N compartment grid | **yes**, to couple anything mechanical |
| the scaffold area within the field | a **number**, the fraction of the field the scaffold covers | **yes**, it is one of the three stated run inputs |
| the scaffold's actual shape | a geometry CSV | **only if you have measured geometry** |

The middle row is the one that surprises people, so to be explicit: **you do
not have to author a scaffold geometry file.** The brief specifies the scaffold
*area* as an input, not a mask, and this side turns the number into the mask:

```bash
occsim run --field 1500 --scaffold-area-frac 0.4 --n 3 --mig-speed 10
```

That builds a printed-strand lattice covering exactly 40 % of the field, writes
it next to the run's output so the geometry a result used is recoverable, and
reports occupancy over it. `--scaffold-shape disc|square` gives area-matched
controls, `--strand-pitch` sets the strand spacing.

Send a geometry CSV instead (`--scaffold-geometry`) when you have real printed
geometry to represent, traced from a print, say, rather than a plausible
lattice at a stated area. Passing both is refused as ambiguous. The format is
in [SCAFFOLD_GEOMETRY_FORMAT.md](SCAFFOLD_GEOMETRY_FORMAT.md) and the same C++
loader validates it as any other input.

Nothing shipped in `data/` is measurement. Every file there labels itself
`SYNTHETIC` or `GENERATED` in its own header, and the figures repeat the label,
so a demonstration fixture cannot quietly become a result.

## The interface

This codebase reads two kinds of CSV on the same compartment grid with the same
orientation convention, and *carrying* that data (loading, validating,
reporting, displaying) is entirely separate from *coupling* it (letting it
change the simulation's output), which happens only behind an explicit flag.

| file | what it says | how often | format spec |
|---|---|---|---|
| stiffness map | how stiff each compartment is | per DIV | [STIFFNESS_FORMAT.md](STIFFNESS_FORMAT.md) |
| scaffold geometry | which compartments are scaffold at all | once per scaffold | [SCAFFOLD_GEOMETRY_FORMAT.md](SCAFFOLD_GEOMETRY_FORMAT.md) |

## Which number comes back

Supplying a scaffold geometry changes which occupancy is *reported*, never what
an existing name means. Both denominators are always produced:

| surface | whole field | scaffold |
|---|---|---|
| engine stdout | the `= <num> %` headline | `OCCUPANCY_SCAFFOLD_RESTRICTED` line |
| CSV columns | unprefixed (`mean`, `sd`, `rel`, ...) | the same names with a `scaffold_` prefix |
| figures | grey reference curve | the plotted curve |

A merge has to preserve that separation: no column may change denominator
between runs, or a scaffold run stops being comparable with a bare-field one.
The scaffold denominator is a count of whole compartments, not an area
integral, which is the same quantisation that keeps fractional coverage out of
geometry v1.

## Shape of a stiffness file

```
# occupancy-sim stiffness map v1      <- exact first line
grid,30,30                            <- compartments per side
scaffold_um,1500                      <- side of the whole field, um
units,kPa                             <- carried verbatim
div_step,5                            <- any integer >= 1; declared, not inferred
div,cx,cy,stiffness
5,0,0,1.20
...
```

Scaffold geometry is the same shape with a single 0/1 `scaffold` column instead
of a per-DIV value. The rules the loader enforces, orientation, completeness,
no gaps in the DIV sequence, are in the format specs; they are not repeated
here.

## Before coupling: check WHICH coupling the map suits

There are two stiffness couplings and they fail on opposite map shapes, so a
valid file can be inert under one and fine under the other.

| coupling | flag | inert when |
|---|---|---|
| gradient | `--stiff-weight` | the map has no large-scale gradient, a cell reaches the nearest local minimum immediately. Ring and interference fields. |
| magnitude | `--stiff-speed` | the map is near-uniform, or its features are so broad that a cell crosses fewer than about two during the culture. |

```
./build/occupancy_engine --stiffness-check FILE   # format only, says nothing about coupling
occsim stiffness --coupling gradient --stiffness-map FILE --mig-speed 10
occsim stiffness --coupling speed    --stiffness-map FILE --mig-speed 10
```

Pick from the table rather than sweeping both and keeping whichever moved.
`--mig-speed` defaults to 0 and both couplings act on migration, so without it
the engine refuses the run: a coupling that could never move a cell is an error
here, not a silent null.

## What makes a dose field usable here

This is the part most likely to save you a wasted campaign, and every number
below was measured rather than assumed.

**Grid resolution is not the lever it looks like.** A map's grid must equal
`--compartments`, and a 200 um ring period does need roughly 10 um compartments
(150x150) to be represented at all. But regenerating a J0^2 field at 150x150,
where its rings are fully resolved, changed no screen verdict: the feature
period moved only 337 um to 323 um, because J0^2's central lobe is far wider
than its ring spacing and dominates the measured scale. Sending a finer map, on
its own, will not make a ring field couple.

**The envelope is the lever.** J0^2 decays as 1/r, so 91-93 % of the map's value
range sits above the p75, concentrated in a central spot that a population mean
never feels. Flatten that envelope, which is also closer to an ideal
non-diffracting Bessel beam, whose rings carry comparable power, and the
top-quartile share drops to about 0.13. Only then does ring period bind, and
the screen passes at 150 um or finer.

**A passing screen still is not a large effect.** On a flat-envelope 120 um
field at 150x150, which the screen rates "coupling should respond", the
soma-stiffness readout shifts 1.3 pooled sem at N=3 and 1.6 pooled sem at N=15,
both over 12 seeds. Roughly 42 seeds at N=15 would be needed for a 3-sem
result. Budget for that before designing a comparison around it.

**So what to aim for:** contrast spread across the field rather than
concentrated centrally, features fine enough that a cell crosses at least two
during the culture, and seed counts several times what the committed figures
use. Check a candidate against those before generating the real thing.

**And state the seeding with any result.** All of the above assumes the
population samples the map uniformly. It does not at the committed operating
point, and a clustered population can respond strongly to a field that would be
inert for a spread one, or the reverse. `occsim stiffness` prints the seed
box's mean against the whole map's on every CSV and figure; quote that ratio
whenever you quote a coupling result.

**There is no tooled couplability screen.** Judging a field is a hand
comparison against the measured tables below, and the answer of record is a
run: `occsim stiffness --coupling speed --stiffness-map
your_field.csv --mig-speed 10`.

## The measured tables to judge against

These are the calibration data behind the specifications above, kept here so
this document stands on its own. All measured, none assumed.

**Gradient coupling: basin scale.** The mean distance a cell travels downhill
before hitting a local minimum, against one migration step (`--mig-speed` x
`--dt`, 10 um at the operating point). Of three candidate statistics this was
the one that tracked measured steering (correlation +0.995; gradient coherence
+0.94; local-minimum count -0.31), and it called 7/7 maps correctly:

| map | measured span | basin | verdict |
|---|---|---|---|
| single soft disc | 0.2085 | 561 um | works |
| the same, mirrored | 0.2176 | 559 um | works |
| the same, shuffled | 0.0020 | 66 um | inert |
| bessel fixture 30x30 | 0.0023 | 74 um | inert |
| J0^2 physical 30x30 | 0.0007 | 80 um | inert |
| J0^2 physical 300x300 | 0.0004 | 7 um | inert |
| J0^2 demo 300x300 | 0.0025 | 52 um | inert |

Calibration band: strong at 56 steps, null at <= 8 steps, untested between.

**The two statistics invert.** This is what lets you pick a coupling rather
than sweeping both and keeping whichever moved:

| map | basin (gradient) | crossings (speed) |
|---|---|---|
| one coherent soft region | 344 um -> works | 0.4 -> inert |
| concentric rings | 65 um -> inert | 2.4 -> responds |

**Speed coupling: ring period, flat envelope, 150x150.** Once the envelope is
flattened, ring period binds and finer passes:

| ring period | measured feature | crossings | top share | screen |
|---|---|---|---|---|
| 200 um | 328 um | 1.8 | 0.13 | LIKELY INERT |
| 150 um | 242 um | 2.5 | 0.12 | passes |
| 120 um | 197 um | 3.0 | 0.12 | passes |
| 80 um | 132 um | 4.5 | 0.12 | passes |

**The sign of `--stiff-weight` is a convention, not a settled fact.** Magnitude
blends the bias against the random walk; sign picks the direction, positive
down-gradient (toward soft). Soft-ward is the default following the durotaxis
review, but cortical outgrowth is elsewhere reported stiffness-insensitive, and
single axons run faster on stiff substrate in the same assay that supports
soft-ward steering. If your optical model implies a direction, set the sign
deliberately rather than inheriting ours.

Occupancy alone cannot settle it. At N=3, DIV15, `--mig-speed 10`, 39-parent
pool, 12 seeds, on the example map:

| stiff_weight | occupancy % | coverage-weighted stiffness (kPa) |
|---|---|---|
| -0.9 | 34.56 +- 2.08 | 1.0758 +- 0.0153 |
| -0.6 | 36.00 +- 2.05 | 1.0775 +- 0.0154 |
| -0.3 | 48.69 +- 2.90 | 1.0447 +- 0.0136 |
| 0.0 | 57.08 +- 3.90 | 0.9903 +- 0.0150 |
| +0.3 | 63.34 +- 2.96 | 0.9419 +- 0.0123 |
| +0.6 | 48.25 +- 2.78 | 0.8730 +- 0.0125 |
| +0.9 | 46.50 +- 2.85 | 0.8694 +- 0.0128 |

Occupancy is non-monotonic and falls for *both* signs, because any strong bias
funnels cells into one region where they sweep less area than a free walk. The
directional readout is **coverage-weighted mean stiffness**, which falls
monotonically 1.076 -> 0.869 kPa with the uncoupled run at 0.990 between the
halves. Read direction off that column, never off occupancy.

**That table is uniform seeding, and the committed figure is not.** At the
committed operating point (N=25, `--seed-frac 0.1`) the negative half of the
same sweep is flat and unordered: 1.0116, 1.0109, 1.0143 kPa against sems of
about 0.0007, with the uncoupled run at 1.0132, at the top rather than between
the halves. The cause is where the cells start: the seed box averages 0.328 kPa
against 1.015 over the whole map, **0.32x**, i.e. on the floor of the basin
(map minimum 0.311), so "steer toward stiff" has the entire map uphill and
nowhere in particular to go. Only the positive half is directional there. **Read this
table for what the coupling does to a uniformly seeded population, and the
committed figure for what it does at the operating point; they are different
questions and only the first one generalises.**

Three caveats on those tables, all of which matter to you:

- **The speed half is PROVISIONAL.** Basin scale was scored against seven maps
  measured end to end; the crossing count against three. Treat crossings as a
  pre-flight hint, not a verdict, and run the coupling for the answer.
- **The screen is necessary, not sufficient.** A map can pass and still produce
  no resolvable effect; see the 1.6-sem result above.
- **The shipped Bessel fixture fails both, under uniform seeding.**
  `data/stiffness_bessel_synthetic.csv` is a 30x30 map, so its 200 um rings
  alias into a 337 um feature, and 91 % of its value range sits above the p75
  in a central spot a uniformly spread population never feels. That is a
  property of the fixture, not of the mechanism. Any Gaussian-vs-Bessel claim
  needs a 150x150 map, which you can build with `occsim make-stiffness` below.

  **But the committed figure on that fixture does respond**, monotonically,
  0.138 -> 0.180 kPa over `--stiff-speed` 0..1. There is no contradiction: the
  committed operating point seeds a central cluster, and the seed box averages
  0.382 kPa against 0.133 over the whole map, **2.88x**, so the population is
  seeded straight into the concentration this section warns about. So the figure shows the mechanism
  steering, not the fixture passing. Every table on this page was measured
  under **uniform** seeding; a clustered run samples the map somewhere else
  and its result does not transfer. Each coupling CSV and figure now states
  its own seed-box ratio for exactly this reason.

## Fixtures to develop against

You need something to point the loader at before your own pipeline emits
anything, and the claims above need to stay reproducible, so `occsim
make-stiffness` builds maps from a beam shape and the same saturating-crosslink
relation the shipped Bessel fixture declares in its header:

```bash
# a single central lobe: the shape the GRADIENT coupling suits
occsim make-stiffness gauss.csv --beam gaussian --waist 400

# physical J0^2 rings, and the envelope problem this page describes
occsim make-stiffness bessel.csv --beam bessel --ring-period 200

# the flat-envelope, finely-resolved case that the measurements above used
occsim make-stiffness flat.csv --beam bessel --envelope flat \
    --ring-period 120 --grid 150
```

Each run reports the top-quartile share of its own value range, which is the
statistic that decides whether a ring field couples, and warns when a ring
period is too coarse for the grid. The three commands above give 0.87, 0.91 and
0.11 respectively, the middle one is the concentrated-centrally failure this
page warns about, the last is what a usable field looks like.

This is a fixture generator, not an optical model: no photochemistry beyond
`E = E_floor + (E_max-E_floor)*(1-exp(-k_c*D))`, and no substitute for the real
dose fields.
