# bench/

Not part of the package (`^bench$` in `.Rbuildignore`), and not tests — nothing
here asserts a correct answer. `regress.R` prints a table of ratios between two
trees and exits non-zero if anything crossed the threshold.

```sh
Rscript bench/regress.R                    # vs upstream/dev, synthetic data
Rscript bench/regress.R --ref=HEAD~1
Rscript bench/regress.R --ref=e120c82 --reps=5 --threshold=1.1
GD_BENCH_STORE=/path/to/store Rscript bench/regress.R
```

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
measuring; only the package under test differs. Never run the ref's own copy, or
you are comparing two different measurements.

## Reading the output

A ratio is `ref / now`, so **below 1 means slower than the ref**. `~0` marks a
case under the timer's useful resolution — `processData logNorm` belongs there,
since it appends a record and does no data pass, and it becoming *measurable*
would itself be the regression.

A `SLOWER` flag is a prompt for an explanation, not a verdict. Some are correct:
`featStats` and `cellStats` are slower than `e120c82` because they now apply the
store's op chain, which that version skipped entirely — they were reporting
statistics on unnormalized values. The flag is doing its job by making you say so.

## Real data

`GD_BENCH_STORE` points at an existing store instead of generating one, for
confirming that synthetic conclusions survive real skew — feature density
concentrated in the HVGs, a wide library-size spread. On Atera the whole sweep is
a couple of minutes, so this is cheap enough to run when a result looks
surprising. `storeWrite` is skipped in that mode, having nothing to ingest.

Heavier work — cross-toolbox comparisons, scale studies — does not belong here.
It needs its own dependency set and does not have to track this package's API
commit by commit.
