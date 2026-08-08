# bench/

Not part of the package (`^bench$` in `.Rbuildignore`), and not tests — nothing
here asserts a correct answer. `regress.R` prints a table of ratios between two
trees and exits non-zero if anything crossed the threshold.

```sh
Rscript bench/regress.R                          # every dataset, every case
Rscript bench/regress.R --data=synthetic         # one dataset
Rscript bench/regress.R --data=atera
Rscript bench/regress.R --cases=PCA              # one slice of cases (regex)
Rscript bench/regress.R --cases='featStats|cellStats' --reps=5
Rscript bench/regress.R --ref=HEAD~1 --workers=1 --threshold=1.1
Rscript bench/regress.R --ref=A --now=B          # pin both ends
Rscript bench/regress.R --report                 # re-print, no work
```

`--ref` compares against your working tree, which is what you want while
developing. `--now` replaces the working-tree side with a second worktree, so the
comparison no longer depends on where the working tree happens to sit — **use it
for any number you record somewhere durable**, such as an issue, since otherwise
the same command produces a different comparison once the branch moves on.

Three files: `regress.R` drives, `_runner.R` measures one tree, `_cases.R` is
the list. Datasets and cases are independently selectable and everything runs by
default. `atera` needs `GD_BENCH_H5` pointing at a `cell_feature_matrix.h5`;
without it the default run does synthetic only and says so.

Parameters are aligned with the original Atera harness — `ncp = 30`, 2000 HVGs,
`batch_cells = 25000`, 8 workers, cap 300 s — so the overlapping cases are
comparable to its recorded `profile.csv` rather than being a fourth
configuration. HVGs come from HVF output, not `rownames(pe)[1:n]`: feature
selection picks the densest rows on purpose (adr/0008), so an arbitrary prefix
reads a sparser band than any real workflow and makes PCA look faster than it is.

## Why it is shaped this way

**Ratios, not baselines.** Absolute times are not comparable across machines,
arrow versions, or background load, so a recorded-baselines file goes stale
immediately. What is stable is the ratio between two trees measured back to back
on one machine, which is why the ref is an argument and there is nothing to
commit but the script.

**Verbs, at explicit chain states.** The regression this was built for — 22x on
`featStats` — would have been invisible to a micro-benchmark of the accumulator
it lives in. The accumulator was fine; the *verb* was slow because the store's op
placement changed which execution path ran. So the case list names the chain
state, and `featStats [raw]` sits next to `featStats [norm]` precisely because
they can diverge. Confirmed: against `e120c82` the raw case reads 0.99x and the
normalized one 0.37x.

**Public API only.** Each side runs against a different tree, which may not have
the internals — or the exported names — the other does. An earlier ad-hoc version
died on `Giotto::reduceParam` not existing on the older ref.

**One runner, two trees.** The working tree's copy of the script always does the
measuring; only the package under test differs. This holds under `--now` too —
both sides are worktrees, but `_runner.R` and `_cases.R` still come from the
working tree. Never run the ref's own copy, or you are comparing two different
measurements.

## Reading the output

A ratio is `ref / now`, so **below 1 means slower than the ref**. `~0` marks a
case under the timer's useful resolution — `processData logNorm` belongs there,
since it appends a record and does no data pass, and it becoming *measurable*
would itself be the regression.

A `SLOWER` flag is a prompt for an explanation, not a verdict. Some are correct:
`featStats` and `cellStats` are slower than `e120c82` because they now apply the
store's op chain, which that version skipped entirely — they were reporting
statistics on unnormalized values. The flag is doing its job by making you say so.

## Verbs and one pipeline

Most cases time a single verb against a fixture built in setup, which is what
localizes a regression to the step that caused it. But verb-level measurement
misses a slowdown spread too thinly across steps to trip any single threshold,
so there is one end-to-end case: `filter -> norm -> HVF -> PCA` from a raw
store, with the HVGs taken from HVF output rather than faked. It is the only
case exercising the HVF -> PCA handoff, and it runs at `reps = 1` because its
job is a trend line, not precision.

Both PCA methods are timed. `random` (Halko) and `gram` are separate
implementations with their own chain-demotion and fallback paths, and the
feature ratio here (1000 of 4000) trips the transient bake, so that is the path
being measured rather than the un-baked one.

## Memory

Peak footprint is sampled system-wide (`vm_stat` / `MemAvailable`) and reported
as the largest drop in free memory, because arrow allocates in C++ where R's
`gc()` accounting cannot see it. Two consequences:

- anything below 50 MB reads as `-`, since background activity dominates there
- the number is attributed to whichever case was running, so a concurrent
  allocation elsewhere on the machine lands on that case. Treat single-sample
  memory differences with more suspicion than timings, and re-run with
  `--reps=3` before acting on one.

Worker startup is a specific trap: eight mirai daemons each attach arrow and
data.table, and although they are warmed outside any timed case, arrow's
per-worker pools can still grow during the first parallel case.

## Real data

`GD_BENCH_H5` points at a 10x/Xenium `cell_feature_matrix.h5` and the run starts
from ingest, for
confirming that synthetic conclusions survive real skew — feature density
concentrated in the HVGs, a wide library-size spread. On Atera the whole sweep is
a couple of minutes, so this is cheap enough to run when a result looks
surprising. Ingest is a measured case there rather than a skipped one — it is the
largest single memory consumer in the run (~16.7 GB on Atera), and it is the
stage that has historically trailed BPCells and scstream.

A store serialized from an older tree is deliberately not supported as an input:
it will not deserialize against a changed class definition, and going through the
Input path is what yields real ids.

Heavier work — cross-toolbox comparisons, scale studies — does not belong here.
It needs its own dependency set and does not have to track this package's API
commit by commit.
