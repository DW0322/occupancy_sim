# Morphology pool (archive)

The full clone-weighted pool, reachable via `--morph-dir` and
kept so previously logged pool122 results stay reproducible. The default pool
is `../morphologies_balanced/`.

Drop reconstructions in this folder and there is nothing else to do: `pool()`
picks up every `.asc`, `.swc` and `.h5` here, sorted by filename so the seed to
template mapping stays reproducible. Templates are parsed lazily, so a large
pool costs nothing until sampled.

Sources and formats:
- NeuroMorpho.org: `.swc` (standardised version)
- Allen Cell Types: `.swc`
- Blue Brain NMC portal: `.asc` (all 122 current files)

Keep the original filenames: they encode the cell ID needed for citation.
Blue Brain files require the Markram 2015 / Ramaswamy 2015 citations; log any
new source's citation requirement in LICENSE and the README's references when
you add it.

Prefer distinct cells over clones when choosing what to download: clones of a
parent reconstruction are not independent samples.

`rp101229_L5-1_idD_Clone_0.asc` here is byte-identical to `../morph.asc`, which
is the frozen regression reference the tests pin by its short path. Both copies
are intended; edit one and you must edit the other.

## Composition (2026-08-18 audit)

122 files from 39 distinct parent cells, 1-7 clones per parent. Three
byte-identical duplicates were removed, keeping the short-form filenames because
the tests pin them. Cell-type identity of the additions is unverified. Measure
the pool with `./build/morph_stats data/morphologies/*.asc`.
