# 0019. Parallel windows run as kernel code on a per-call pool that never loads GiottoDisk

- **Status:** Accepted
- **Date:** 2026-09-29
- **Supersedes:** —
- **Superseded by:** —

## Context

The windowed parquet expression write (adr/0018) was single-core in practice.
Its windows are independent, but running them in parallel was blocked by what a
worker would need. Fork is not available on Windows and fails under Positron.
A socket worker (mirai, or future over mirai) that loads GiottoDisk pays about
7 s and 1.65 GB of reported RSS against 0.14 s and 101 MB for arrow alone. It
also keeps the namespace, and whatever Arrow memory it touched, until it exits.
A user's `future` plan shares one long-lived pool across every call, so one call
that loads the package leaves the pool heavy for the rest of the session.

The alternatives to loading the package were weak. Arrow R cannot serialize an
`Expression`, so a built query cannot be shipped, and the installed arrow has
Substrait compiled out. Rebuilding the scan in the worker by hand duplicates
storeRead's composition. A prototype instead copied the functions the write
reaches out of the namespace, re-parented onto a plain environment. On Atera
(170,044 cells x 18,028 genes, library-normalized and logged, 306,980,503
values, 16 cores) it wrote identical files to the serial path:

| path | time |
|---|---|
| serial, default budget (2 windows) | 22.0 s |
| serial, 16 windows | 27.4 s |
| 16 windows on 4 / 8 / 12 arrow-only daemons | 6.6 / 4.5 / 4.6 s |

The daemons started in 0.5 s. Each held about 1.35 GB at its sort peak, so 8
workers peak close to the serial default. Per-window time rose from 1.6 s to
2.1 s to 3.0 s at 4 / 8 / 12 workers, which is where the gain stops.

Through the shipped `storeWrite()`, at the default budget, the same store
wrote in 7.3 s on 8 workers against 22.0-22.6 s serial, and in 10.1-10.6 s on
4 workers. The shipped call pays more than the prototype's loop: 0.7 s to start
and stop the pool, a 1.0 s density scan (this store has no cached marginals),
0.4 s of stats, and 13-14 smaller windows where the serial default uses 2. With
the windows pinned equal, the loop took 6.4-6.6 s on 8 workers against
26.2-27.0 s serial, and the files matched.

A crawl from the interpreter entry points reached 22 functions in 5 files, none
of them S4 dispatch except a bare `unique` that GiottoDisk masks. So the code a
worker needs was already separable from the S4 layer around it.

## Decision

Give that code a boundary. `R/kernel-*.R` holds functions that take plain data
and a lazy query, call other packages only as `pkg::fn`, and reach no S4
generic; a test enforces it. `.kernel_bundle()` copies a kernel entry and what
it reaches into a plain environment. `.isolated_map()` runs it on a mirai pool
started for the call under its own compute profile, and shuts the pool down on
exit. The windowed write fans out through it when `.par_workers() > 1` and
mirai is installed, and falls back to serial for a window it cannot lower.

## Consequences

- The windowed write is about 3x faster end to end at 8 workers on Atera, with
  the same files. Its window budget is shared between workers (`bytes_per_nz = 96 * n`),
  which gives more, smaller windows.
- Kernel code is constrained. It cannot take a store, dispatch a generic, or
  call a shell helper, and a bare call to a base name GiottoDisk masks is a bug
  there even when it happens to work.
- The pool costs about 0.7 s to start and stop per call, so this pays only for
  work of several seconds. A short pass should not fan out.
- The user's `future` plan still sets the worker count but not the workers.
  While the pool runs, the plan's idle daemons stay alive.
- Worker results are row counts, not objects, so a consumer that needs data
  back has to design its return shape.
- **Revisit if** arrow gains a serializable plan, or if more consumers want the
  kernel than this package does. Moving `R/kernel-*.R` into its own light
  package is then a file move.

## Alternatives considered

- **Fork (`mclapply`)** — free and fastest where it works, but not on Windows
  or under Positron, which is the reason for this ADR.
- **Workers that load GiottoDisk** (`future.packages`, what the ingest path
  does) — 7 s per worker, and the load is permanent for the pool.
- **Hand-written scan in the worker** — the prototype's first form. It
  duplicates storeRead's composition, and only a parity test would keep the two
  in step.
- **Kernel as `inst/` scripts sourced by workers** — loses codetools checking,
  roxygen and coverage, and the package would need a second load path to use
  its own code.
- **A persistent named pool** — skips the 0.7 s, but retains Arrow memory
  across unrelated calls, the problem the per-call pool avoids.

## References

- `R/utils-isolate.R` (`.kernel_bundle`, `.isolated_map`); `R/kernel-*.R`;
  `.pestore_write_windowed()` in `R/methods-parquetExprStore.R`.
- `tests/testthat/test-kernel-isolation.R`,
  `tests/testthat/test-parquetExprStore-write-layout.R`.
- adr/0018 (the windowed write); adr/0004 (op machinery roles).
