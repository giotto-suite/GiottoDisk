# 0016. Transpose is view state, not a chain step

- **Status:** Accepted
- **Date:** 2026-09-18
- **Supersedes:** —
- **Superseded by:** —

## Context

`spatValues()` in GiottoClass reads an expression matrix by slicing it to the
requested features, transposing, and coercing to `dgCMatrix`. A
`parquetExprBase` store could answer none of that: neither the transpose nor
the coercion had a method for it. The coercion is a straightforward wrapper
over `storeRead(output = "dgcmatrix")`. The transpose is the interesting one,
because where it lives decides whether it can interact with the op chain.

An eager `t()` — materialize, `Matrix::t()` — was rejected. It would satisfy
the signature while discarding the property the class exists for, and once a
caller works, the lazy version never gets written.

So the question is whether a lazy transpose is a **recorded step** on `@ops`,
composing in sequence with everything else, or **view state** in a slot.

### What the store's layout actually is

The store is already orientation-indirected, which is easy to miss:

| | |
|---|---|
| logical `dim()` | 337 genes x 462 cells (`nrow` = genes, Bioconductor convention) |
| on-disk `row_id` | 462 distinct values — the **cell** axis |
| on-disk `col_id` | 337 distinct values — the **gene** axis |

So `nrow()` counts genes while `row_id` *means* cell. Every axis-aware site
already resolves a semantic axis to a physical column rather than assuming
i/j — `key <- if (identical(axis, "feat")) "col_id" else "row_id"`.

## Decision

**Orientation is a `@transposed` logical slot on `parquetExprBase`, flipped by
`t()`. It is never a record on `@ops`.**

The flag is consulted **only** where a *logical position* resolves to a
semantic axis:

- `dim()` / `nrow()` / `ncol()` / `dimnames()`
- `[`, which indexes the logical matrix
- the reshape in `storeRead(output = "dgcmatrix")`
- the `max_rows` / `max_cols` materialization cap, which is a promise about
  the matrix the caller receives

Nothing in the op machinery reads it. The semantic-axis → physical-column map
is **invariant** under the flip: `axis = "cell"` resolves to `row_id` either
way, because a cell *is* `row_id` on disk. That is the whole reason this is
safe, and it is the thing to preserve.

## Why a slot rather than an op

Three independent reasons, in increasing order of how much they constrain
future work:

1. **It is an involution.** `t(t(x))` is `x`. That is last-write-wins, which
   is slot semantics; `@ops` means "and then also".

2. **It is the same category as `@cell_idx` / `@gene_idx`,** which are already
   slots. Those are safe because op payloads are keyed by on-disk id rather
   than view position (`adr/0003`), so narrowing a view cannot invalidate a
   queued op. Orientation rides that identical guarantee. If `[` had to be a
   slot, so does this — the proof is shared, and `adr/0006` already states the
   general rule that view state is not chain state.

3. **Every op is orientation-blind.** The expression chain holds `multiply`
   and `add` (semantic-axis-keyed) and `log` (elementwise). Transpose commutes
   with both kinds, so sequencing cannot matter. Verified in both directions:
   flipping before and after a queued `multiply` yields identical matrices,
   and the op record is carried across untouched.

## The invariant this depends on, and where it ends

> The expression chain admits only elementwise or semantic-axis-keyed ops.
> That is what makes orientation view state rather than a step.

Transpose would **not** commute with a *positional* op. The tabular chain
already has three — `head`, `tail`, `sample`. If one of those ever migrates to
the expression chain, orientation stops being view state and this ADR needs
revisiting rather than patching. The two commutation tests in
`test-parquetExprStore-subset.R` are the tripwire.

## Consequences

**Refused, because bytes would disagree with the view:**

- `storeWrite()` on a transposed store. A write bakes the *chain* into on-disk
  values and returns a store with empty chains; the flip is not part of the
  chain and would not be baked, so the result would report an orientation its
  bytes do not have.
- `unionParquetExprStore()` over transposed substores. `feat_ids` / `cell_ids`
  would mean opposite axes across substores, and the existing ordering check
  would then compare a gene vector against a cell vector — passing or failing
  for the wrong reason. Combine upright, then transpose the union.

**Not changed, deliberately:** the file stays sorted `(row_id, col_id)`, i.e.
cell-major. A flip moves no bytes, so gene-axis access remains the scan-heavy
direction however the matrix is presented. A materializing repartition is the
natural companion and is deliberately a *different operation*, not a mode of
this one.

**Message wording:** subset errors now name the axis (`"logical gene axis"`)
rather than the position (`"logical row"`), because position is
orientation-dependent and the axis is not.
