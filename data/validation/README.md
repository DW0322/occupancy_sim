# Validation morphologies

Third-party reconstructions used as the independent comparison population for
the morphology validation summarised in docs/RUNNING.md section 7 and recorded
in `results/validation/`. They are **downloaded, not
redistributed**: `data/validation/` is gitignored, so what is here is whatever
was fetched on this machine, and a fresh clone starts empty.

**The fetch and comparison scripts are not shipped** - the comparison was a
one-time check and `results/validation/` is its record. To repopulate this
directory you need the Allen Cell Types Database directly:

```bash
pip install allensdk morphio     # then fetch SWC reconstructions into allen/
```

Both sides of the comparison are measured with `build/morph_stats`, which is
shipped, so the check can be recreated from these files without those
scripts - that is how `results/validation/` was regenerated on 2026-08-26.

`allen/` receives SWC files from the Allen Cell Types Database plus a
MANIFEST.json recording specimen ids, dendrite type, structure and layer.

Every download is parsed through MorphIO, the same library the engine loads
through, before it is kept, so the directory only holds cells that reach the
statistics. Anything rejected is recorded under `excluded_unreadable` in
MANIFEST.json with its specimen id and reason rather than skipped in silence.
At `--max 120`, seven cells currently fail a section-type rule, leaving 113.
`--keep-unreadable` fetches them anyway.

To check what is actually present rather than assuming:

```bash
ls data/validation/allen/*.swc | wc -l
```

See docs/RUNNING.md section 7 for what the comparison does and does not
establish. In short, the Allen set is acute-slice tissue rather than
dissociated culture, so passing means "our templates are typical cortical
interneurons" and not "this is what a neuron on a scaffold looks like".
