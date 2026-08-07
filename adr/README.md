# Architecture Decision Records

One file per architectural decision, dated and immutable. An ADR records *why a
choice was made at a point in time*, including what was rejected and what the
choice costs. It is history, not documentation of current behaviour.

## Why these exist alongside the other docs

GiottoDisk already documents architecture in three places. ADRs answer a
question none of them do, so keep the boundary sharp:

| Doc | Answers | Tense |
|---|---|---|
| `AGENTS.md` | What invariants hold right now, and where the code is | present, terse |
| `vignettes/articles/design.Rmd` | How the architecture works and hangs together | present, narrative |
| `vignettes/articles/roadmap.Rmd` | What we intend to change | future |
| `adr/` (here) | Why we chose this over the alternatives, and when | past, immutable |

Practical test for where something belongs:

- "`nrow()` returns a double" → AGENTS.md (invariant a code change must respect).
- "The lazy op system composes ops in order at read time" → design.Rmd.
- "We plan to move the norm ops back onto `@ops`" → roadmap.Rmd.
- "We rejected a positional row index because a 2B-row index costs ~16GB" → ADR.

The overlap is intentional and one-directional: AGENTS.md and design.Rmd state
the *outcome* of an ADR without rehearsing the argument; the ADR is where the
argument and the discarded options live. When they disagree, AGENTS.md wins for
current behaviour and the ADR wins for intent — and the disagreement is itself a
signal that a superseding ADR is owed.

## Writing one

1. Copy `0000-template.md` to `NNNN-short-kebab-title.md`, taking the next free
   number. Numbers are record order, not decision order — an ADR backfilled
   today for a 2025 decision still takes the next number and carries the older
   date.
2. Fill it in. Keep it to a page; if it needs more, the extra belongs in
   design.Rmd and the ADR should link to it.
3. Add a row to the index below.

## Amending one

An accepted ADR is not edited, with two exceptions: its `Status` line, and links
added to it (`Superseded by`). Everything else changes by writing a new ADR that
supersedes it. A wrong ADR that got reversed is more useful than no record of
the reversal.

Statuses: **Proposed** · **Accepted** · **Superseded by NNNN** · **Reversed**
(tried, undone, nothing replaced it) · **Deprecated** (still true, no longer
load-bearing).

## Scope — what earns an ADR

Something an outsider (or you in six months) would otherwise change by accident.
Roughly: a decision that constrains future code, was contested or non-obvious,
or has a cost worth remembering.

Not: bug fixes, refactors that preserve behaviour, or naming conventions
(those go in AGENTS.md "Conventions").

## Index

| # | Title | Status | Date |
|---|---|---|---|
| [0001](0001-no-positional-row-indexing.md) | No positional row indexing | Accepted | 2026-04-17 |
| [0002](0002-post-ops-phase-split.md) | Two-phase `@ops` / `@post_ops` split on `parquetExprStore` | Superseded by 0005 | 2026-07-22 |
| [0003](0003-op-payloads-keyed-by-on-disk-id.md) | Op payloads keyed by on-disk id; `[` needs no subset slicing | Accepted | 2026-07-22 |
| [0004](0004-op-machinery-roles.md) | Op machinery roles, and the invariants that connect them | Accepted | 2026-08-07 |
| [0005](0005-ops-is-the-pre-materialization-prefix.md) | `@ops` is the pre-materialization prefix, not the set of lowerable ops | Accepted | 2026-08-07 |
| [0006](0006-view-state-is-not-chain-state.md) | View state is not chain state; window-dependent ops bake at push time | Accepted | 2026-08-07 |

## Backfill candidates

Decisions already argued out elsewhere in the repo or in commit messages, not
yet written up. Not a queue — write one when it next comes up in conversation,
so the ADR captures the argument while it is fresh.

- **GeoParquet CRS is written as `crs: null` + `edges: planar`.** sedonadb
  hardcodes `null` → `OGC:CRS84`, so the WKT side must emit
  `ST_GeomFromText(wkt, 4326)` for predicates to match. Constrains anything
  touching `.arrow_meta_add_geoparquet()`.
- **`spat_relate` engine precedence** (per-call > option > auto sedona >
  duckdb > terra) and why the narrow path caches ids as an `id_filter` op
  rather than using a tuple-IN subquery (DataFusion does not support it).
- **Manifest concurrency: compare-and-swap on a single JSON**, with the
  append-only chain deferred. Trigger for revisiting is concurrent users, not
  artifact count.
- **`overlapToMatrix` joins string ids to integers before aggregating.**
  Aggregating strings first leaves dangling `utf8_view` buffers at scale — an
  Arrow constraint, not a preference.
- **`parquetExprStore` schema is 4-column with `source_id` load-bearing**, and
  union `storeRead` uses one composite `source_id`-aware filter rather than
  wrapping per substore.
- **Tabular growth is column-add only** (row set fixed at ingest;
  re-segmentation is a new dataset), which is why sidecar + compact was chosen
  over Iceberg/DuckLake.
