# Clone-free morphology pool (DEFAULT)

39 files, one per distinct reconstructed cell. This is what `--pool` uses.

`../morphologies/` is the clone-inclusive set these were drawn from; the
MORPH_DIR comment in `python/occsim/engine.py` says why sampling it would
skew the draw.
Selection here is deterministic (one clone per parent, fixed rule).

Measured: switching between this pool and the full 122 does not move the
normalised occupancy curve - max 0.021 fold divergence against per-seed sd
0.43-0.78 at N=3 with migration, 8 seeds (measured 2026-08-19).

Citations required for any use: Markram et al. 2015; Ramaswamy et al. 2015.

Grouping parents from filenames requires stripping BOTH `_-_Clone_N` and
`_Clone_N` - the pool mixes the two conventions.
