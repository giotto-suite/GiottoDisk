# 0012. Expression-store scan modifications are written once, in dplyr, for every carrier

- **Status:** Accepted
- **Date:** 2026-08-24
- **Supersedes:** —
- **Superseded by:** —

*(Takes 0012 rather than 0011: `refactor/pe-window-seam` already holds 0011
upstream.)*

## Context

`parquetStore` supports `output = "query"`, `"duckdb"` and `"sedona"`, and
rebuilds the scan from store state for each. It does that by compiling `@ops`
into SQL text: `.pstore_sql_inner()` walks the op list a second time, and
`.r_expr_to_sql()` walks the R expressions `subset()` captured.

That is a second implementation of the op registry, and it has drifted from the
first. The SQL side drops `join`, `tail` and `sample` with a warning, and it
collapses a sequential chain into one flat `WHERE` with a single trailing
`LIMIT`, so `head(5) |> subset(...)` means something different depending on the
output. The outputs are translatable in shape, not in semantics.

`parquetExprStore` had no native path at all. `output = "duckdb"` registered the
Arrow scanner via `duckdb_register_arrow()`, so Arrow read the parquet and
DuckDB drained the result — DuckDB's memory limit never saw Arrow's allocation,
and no row-group pruning was DuckDB's to do.

The obvious way to wire it up was to copy the parquetStore approach and write a
`.pe_*_sql()` emitter. `vignettes/expression_windows.Rmd` and ADR 0011, both on
`refactor/pe-window-seam`, assume exactly that, and record the cost as the
reason the path did not exist: it "would require rendering the axis predicates
as SQL."

## Decision

Do not write a SQL emitter. Swap the **carrier** underneath the existing dplyr
code.

`.pestore_to_duckdb()` builds a `tbl_dbi` over DuckDB's `read_parquet` and then
calls the same `.pe_apply_axis_pred()` and `.pe_apply_ops()` that the Arrow path
calls. dplyr targets Acero and dbplyr alike, so one body of code applies the
scan modifications to both engines.

The premise the earlier docs recorded is false. dbplyr already renders all three
ADR-0008 predicate shapes, verbatim, from the `bquote()` objects
`.pe_axis_pred_exprs()` emits:

| shape | SQL |
|---|---|
| gapless | `row_id >= 10 AND row_id <= 20` |
| gaps, few dropped | `(col_id >= 1 AND col_id <= 2) AND (NOT((col_id IN (2))))` |
| gaps, few kept | `row_id IN (1, 3, 5)` |

Nothing had to be re-rendered. `.pe_axis_pred()` stays the single classifier;
what it decides is engine-independent because dplyr, not SQL, is the lowering
layer.

Two carrier-specific escapes survive, and both are narrow:

- `.pe_payload_carrier()` — an Arrow `Table` cannot be `left_join`ed into a
  `tbl_dbi`, so `multiply`'s payload is registered on the connection when the
  carrier is one. The join, the multiply and the column drop are unchanged.
- `.pe_apply_axis_pred()` — dbplyr inlines `%in%` into the query *text*, so
  above a threshold the membership half becomes a registered semi/anti join.
  The classification behind it is untouched.

## Consequences

Equivalence between `"query"` and `"duckdb"` is structural. There is no second
op registry to drift, so the failure mode parquetStore has — an op silently
skipped on one output — cannot arise here without deleting the shared call.

Adding an op means one dplyr branch in `.pe_apply_op`, serving both carriers. A
branch per engine is how the two outputs come apart, so the extension protocol
in `R/utils-pestore-ops.R` now says to reach for a carrier test only where an
engine cannot accept the other's data.

The constraint this buys is that every op must have a dplyr form. `log1p` did
not survive it: DuckDB has no such function and dbplyr does not translate it, so
`.op_transform_log` computes `log(value + 1)`. That cost nothing — measured over
50M doubles it is 7.1% *faster* single-threaded and indistinguishable at the
default thread count, where the transform is memory-bandwidth bound — but the
next op may not get off so lightly, and the right response is to change the
formula rather than to branch.

This does **not** generalize back to `parquetStore`. Its ops genuinely have no
common form — `spat_relate` needs `ST_*`, `id_filter` needs `EXISTS`,
`crop`/`window` are AABB literals, and `subset()` captures arbitrary R
expressions. `.pstore_sql_inner` is the right shape for that store. The
distinction is whether the op vocabulary is closed and dplyr-expressible, not
which store is newer.

`vignettes/expression_windows.Rmd` and ADR 0011 described the expression path as
Arrow-backed and named SQL rendering as the blocker. Both were on
`refactor/pe-window-seam` and could not be corrected here; both were corrected
when that branch merged, and now describe the carrier as a choice available
above the scan.

## Alternatives considered

- **A `.pe_*_sql()` emitter mirroring `.pstore_sql_inner`.** The approach the
  existing docs assume. Rejected: it buys a second op registry and the drift
  that comes with it, to re-derive semantics dplyr already lowers correctly.
- **Keep the Arrow bridge and register the scanner.** What the code did.
  Rejected: Arrow owns the scan, so the memory limit is unenforceable and row
  groups are pruned by Acero rather than by the engine the caller asked for.
- **Keep the bridge as a fallback when an op will not lower.** Rejected: silent
  path-switching is hard to debug, and the tabular path's precedent — dropping
  unsupported ops with a warning — is the failure mode being avoided.
- **`COALESCE(w, 1)` on the multiply join.** Rejected on parity: an unmatched
  key yields `NA` on all three carriers today, and defaulting the factor to 1
  would return an entry's raw value dressed as a normalized one.

## References

- `R/methods-storeRead.R` — `.pestore_to_duckdb()`, `.pstore_to_duckdb()`
- `R/methods-parquetExprStore.R` — `.pe_axis_pred()`, `.pe_apply_axis_pred()`
- `R/utils-pestore-ops.R` — `.pe_apply_op()`, `.pe_payload_carrier()`
- ADR 0004 (executors keyed by record type and carrier) — the rule this applies
- ADR 0008 (axis predicate shapes) — the classification shown to be portable
- Measurements: dbplyr 2.5.1 / duckdb 1.4.4 / arrow 23.0.1.2
