---
name: giottodisk-method
description: Rules for adding a method, verb, or lazy op to GiottoDisk. Use when writing or reviewing any new setMethod/setGeneric in GiottoDisk, adding an op type to a store's op chain, extending storeRead/storeWrite, or adding a store class. Enforces attaching to existing machinery instead of building a parallel path.
---

# Adding to GiottoDisk

GiottoDisk already has machinery for almost everything a new method needs:
lazy op recording, op lowering to Arrow and to SQL, field/special-column
injection, extent resolution, materialization, tiled streaming, chain
editing. **The failure mode this skill exists to prevent is not a bug — it is
a second implementation of a path that already exists**, diverging from the
first the moment either changes.

So the governing rule:

> Before writing a code path, find the seam it attaches to. If you cannot
> name the existing function you are extending, you are probably about to
> duplicate one.

## Orientation

Read in this order. Do not skip to code.

1. `AGENTS.md` — invariants, class hierarchy, where code lives.
2. The **seam table** below — the attachment point for what you are adding.
3. For op work: the header of `R/utils-pestore-ops.R` (chain model + op
   registry) and `adr/0004-op-machinery-roles.md` (roles + invariants).
4. `adr/` before changing any decision that looks arbitrary.

**Code is authoritative over prose.** AGENTS.md and the ADRs state outcomes
and can drift on names; the `switch()` arms and the registry header are the
real list. ADR 0004 says this explicitly of the op registry. When a doc and
the code disagree on a function name or op type, trust the code and fix the
doc in the same change.

## The four questions, before any code

**1. Does the verb already exist?**

```sh
grep -n "exportMethods\|^export(" NAMESPACE      # what GiottoDisk already dispatches
cat R/AllGenerics.R                              # GiottoDisk's own generics (short list)
grep -n "importFrom" NAMESPACE                   # generics borrowed from elsewhere
```

Generics come from one of three places, and picking wrong creates a parallel
verb that shadows a working one:

| Source | Rule |
|---|---|
| GiottoClass / Giotto (`processData`, `filterData`, `reduceData`, `analyzeData`, `calculateOverlap`, `spatRelate`, `affine`, …) | `importFrom()` in `R/pkg_imports.R`, then `setMethod`. **Never `setGeneric` again** — that masks the upstream generic and breaks every other backend's dispatch. |
| terra / base (`crop`, `ext`, `window`, `rasterize`, `centroids`, `subset`, `unique`, …) | Same: `importFrom`, then `setMethod`. |
| Genuinely new to GiottoDisk | Bare `setGeneric` in `R/AllGenerics.R` under the `# definitions ####` banner. No roxygen there — docs live with the method. |

The suite's process verbs are **distinct generics**, not synonyms:
`processData` / `filterData` / `reduceData` / `clusterData` / `analyzeData`.
Route to the correct one; a method on the wrong generic silently never
dispatches.

**2. Is this a method, or an op on an existing method's chain?**

Most "add a capability" requests are the latter. A new *op type* is cheaper
and composes automatically. Ask: does the user call a new verb, or does an
existing verb need to record a new kind of deferred work?

**3. Is it an op, or a slot?**

Ops compose in recorded order; slots are overwritten and have priority
semantics. If your feature needs last-write-wins or override semantics, it is
a slot; if it needs "and then also", it is an op.

Both chains already draw this line, in the same place for different reasons:

- **Geom** — `@crop` and `@window` are slots on `parquetGeomBase`, deliberately
  **not** ops, because window must be able to *override* crop rather than
  stack with it. `@tile_filter` likewise (overwritten, not accumulated, on
  each `crop()`).
- **Expression** — `@cell_idx` / `@gene_idx` are slots, not ops: `[` narrows a
  *view*, and the narrowing must not become a step in the chain. That is what
  lets op payloads stay keyed by on-disk id and therefore invariant under `[`
  (ADR 0003). An op that records "which rows are visible" would have to be
  re-keyed on every subsequent `[`.

Corollary for both: **view state is not chain state** (ADR 0006). If your
feature answers "which subset am I looking at," it is a slot. `[` narrows the
window and never touches the chain; the chain never consults the window at
read time.

**4. Which class tier dispatches it?**

```
dataStore (V) → fileStore → queryableStore → parquetStore → parquetGeomStore → parquetGeomTileStore
                                           → parquetExprStore
                                           → parquetEdgeStore
parquetBase (V)      @fields @ops @post_ops   — shared by single + union
parquetGeomBase (V)  @window @crop @geomtype
parquetExprBase (V)  — shared by parquetExprStore + unionParquetExprStore
```

Dispatch on the **highest virtual that can express the operation**. A method
written on `parquetStore` when `parquetBase` would do is the most common
source of "works on a single store, silently missing on a union" — and then
of a copy-pasted union method.

**Hard rule:** methods on `parquetBase` must be self-contained and must never
call `callNextMethod()` — `parquetBase` does not traverse the `fileStore`
chain. This is why `storeRead` is deliberately *not* defined on `parquetBase`.

## Seam table

For each capability, the one place to attach, and the duplicate it prevents.

| To add | Attach to | Do not |
|---|---|---|
| a lazy row/column op on tabular stores | append `list(type=…, …)` to `@ops`; add an arm to `.ptabular_apply_op` (`R/methods-ops.R`) **and** to `.pstore_sql_inner` (`R/methods-storeRead.R`) | write a bespoke filter inside your method, or post-filter after `collect()` |
| an op that needs columns not requested by the user | add an arm to `.pstore_op_referenced_cols` | widen `fields` at the call site, or re-read the file |
| an R-side step with no Arrow form (tabular/geom) | `@post_ops` + `.apply_post_ops` (`R/utils-spatial.R`) | `storeRead(output="tibble")` then transform — that strips laziness for every downstream caller |
| an op on the expression chain | `.pe_push_op(pe, list(type=…), phase=)`; executor per `(type, carrier)`: `.pe_apply_op` / `.op_*` for Arrow, `.pe_apply_post_op_df` for the collected triplet frame | edit `@ops`/`@post_ops` directly, or reach back and rewrite an earlier record |
| moving where an existing op runs | chain editors: `.pe_demote_ops`, `.pe_chain_none` | take `output="query"` and re-lower the records by hand (ADR 0004, invariant 5) |
| special-column plumbing on a read path | `.pstore_fields_requested` + `.pstore_lazy_fields`, at the **topmost** method before the first `callNextMethod` | inject columns at a lower tier — `output` is rewritten to `"query"` on the way down, so it is too late |
| spatial extent narrowing | `.pstore_active_extent` (window-over-crop) / `.pgeom_ext_intrinsic` (intrinsic scan) / `.pgeom_ext_estimate` (no scan) | call `terra::ext()` on materialized data, or scan when `exact = FALSE` |
| a spatial predicate | the `spat_relate` op + engine dispatch (`R/methods-spatRelate.R`); narrow via `.spat_relate_narrow` | add a fourth engine branch inside your own method |
| shared read post-processing | `.pbase_storeread_processing(atab, store, …)` — the single op-fold + projection + output switch | re-implement the op loop; `parquetStore` and `unionParquetStore` share this on purpose |
| **any statistic over expression values** | `.pe_accum_raw()` (`R/methods-analyzeData.R`) — pick from the `sum` / `sumsq` / `nnz` / `sum_det` accumulators, `axis` for the margin, `by_cell` for a grouped key | write a bespoke `storeRead` + `summarise`, or materialize to compute a mean. If your statistic is a sum, a count, or anything derivable from them, it is already an accumulator |
| per-tile streaming | `tilework::tileApply` (see `R/methods-aggregate.R`) | loop over tile directories yourself |
| a new materialized output format | a new arm in the `storeRead` `output` switch + a `.pstore_to_*` / `.p*_to_*` builder | return a different class from an existing arm |
| writing to disk | `storeWrite` methods; GeoParquet metadata via `.arrow_meta_add_geoparquet()` | hand-roll `arrow::write_dataset` — you lose special cols, hive layout, and `geo` metadata |
| registering an artifact in a project dir | `sourceWrite` (allocate uid → write → register); adoption via `sourceAdopt` | write into the vault path directly and skip the manifest |

## Two op chains. Never cross them.

This is the distinction most easily missed, because both use slots named
`@ops` and `@post_ops`.

**Tabular / geom chain** (`parquetBase`, `parquetStore`, `parquetGeomBase`)
- Op types: `filter`, `head`, `tail`, `sample`, `distinct`, `join`,
  `spat_relate`, plus internal `id_filter`.
- Carriers: lazy Arrow query (`.ptabular_apply_op`) and **SQL**
  (`.pstore_sql_inner`, serving both duckdb and sedona).
- Fold: `.pbase_storeread_processing`.
- `@post_ops` here currently holds one type: `"transform"` (affine2d),
  applied by `.apply_post_ops`.
- No push helper, no monotonic guard — producers append with
  `x@ops <- c(x@ops, list(...))` directly.

**Expression chain** (`parquetExprBase`)
- Op types: `multiply`, `log`, `add` (stub — recorded and refused).
  Authoritative list: the header of `R/utils-pestore-ops.R`.
- Carriers: two **lazy** ones — an Arrow query and a DuckDB `tbl_dbi`, both
  served by the same `.pe_apply_op` → `.op_*` — and the **collected triplet
  `data.table`** (`.pe_apply_post_op_df` → `.pe_apply_post_op_*_df`).
  **There is still no SQL compile on this chain** — do not go looking for a
  `.pe_*_sql_inner`, and do not add one. The duckdb carrier is reached by
  swapping what sits *under* the dplyr executors (`.pestore_to_duckdb` builds a
  `tbl_dbi` over `read_parquet`, then calls the same fold), not by emitting SQL
  text the way `.pstore_sql_inner` does. That is why the two lazy carriers
  share one executor instead of drifting apart. See ADR 0012.
- Folds: `.pe_apply_ops` (lazy), `.pe_apply_post_ops_df` (post).
- Chain editors: `.pe_push_op`, `.pe_demote_ops`, `.pe_chain_none` — the only
  sanctioned way to move where a step runs.

An executor from one chain never runs an op from the other. If you find
yourself wanting that, you have mis-picked the class tier in question 4.

### Which invariants travel

ADR 0004 states its invariants in general terms, but it was written about the
expression chain and **two of them do not hold on the tabular/geom side**.
Check which chain you are on before applying them.

| Invariant | Tabular / geom | Expression |
|---|---|---|
| Records are pure data, no closures | yes | yes |
| Executors keyed by (type, carrier), not by phase | yes | yes |
| Payloads keyed by on-disk id, per source (ADR 0003) | n/a — no axis payloads | yes |
| **Producers append and never revisit** | **no** — see below | yes, enforced by convention |
| **Monotonic phase rule** (nothing on `@ops` once `@post_ops` is non-empty) | **no** — nothing enforces it | yes, enforced by `.pe_push_op` |
| Chain editors are the only way to move where a step runs | n/a — no chain editors exist here | yes |
| Window-dependent ops bake at push time (ADR 0006) | yes — `crop()` freezes half-plane coefficients into an injected `filter` op | yes — `libraryNormParam` freezes its factors |

The revisit exception is deliberate, not an oversight: `.pgeom_set_transform()`
in `R/methods-transforms.R` filters the existing `"transform"` record out of
`@post_ops` and re-adds the composed one, so the chain holds at most one. That
is sound *there* because affine composition is associative and a store has
exactly one pending transform — collapsing loses nothing. It is not sound on
the expression chain, where
arbitrary steps may sit between two records of the same type, which is
precisely the bug ADR 0004 records.

So, both directions:

- Do **not** "fix" the geom transform collapse to satisfy invariant 2. It is
  correct and AGENTS.md documents the auto-collapse as intended behaviour.
- Do **not** copy the collapse idiom into an expression op. Reaching back to
  find and rewrite an earlier record of your own type is the specific mistake
  that `norm_libsize` made; re-running a producer must compose instead.

A new op type on either chain that is **not** associative-and-unique like the
transform gets append-only semantics regardless of chain.

**Slot means position, not capability** (ADR 0005). `@ops` is the prefix that
runs *before* materialization; `@post_ops` is the suffix that runs *after*.
A perfectly lowerable op can legitimately sit on `@post_ops` because something
ahead of it forced materialization. Never infer what an op can do from which
slot holds it, and never "fix" a lowerable op you find in `@post_ops`.

## Checklist A — new op type on the tabular chain

1. Name the **primitive**, not the intent (ADR 0004). `multiply`, not
   `norm_libsize`. Intent belongs to the param class or the verb.
2. Record it as pure data: `list(type = "…", …params)`. **No closures** — the
   record must survive `saveRDS`/`load` and travel to parallel workers.
3. Append directly — there is no push helper on this chain:
   ```r
   x@ops <- c(x@ops, list(list(type = "myop", ...)))
   ```
   Default to append-only: re-running a producer should compose, not rewrite
   an earlier record, because the chain cannot promise nothing was inserted in
   between. Collapsing to a single record (as `"transform"` does) is only
   sound when the operation is associative *and* the store can hold at most
   one — if you are not sure both hold, append. See *Which invariants travel*.
4. Add the Arrow arm to `.ptabular_apply_op`.
5. **Add the SQL arm to `.pstore_sql_inner`.** If you skip this, `output =
   "duckdb"` and `output = "sedona"` warn-and-skip your op and return *wrong
   rows*, not an error. If the op genuinely cannot be lowered to SQL, add it to
   the explicit warn-skip group next to `tail`/`sample`/`join` so the gap is
   declared rather than falling into the unknown-type fallback.
6. If the op references columns, add an arm to `.pstore_op_referenced_cols` so
   the upstream projection keeps them, and to `.pstore_effective_schema` if the
   op widens the visible schema (as `join` does).
7. Check the blockers: ops that make a store unsafe to `rbind2`/`storeWrite`
   must be refused there (`join` already is).
8. Update AGENTS.md's op list and `specialCols` notes if either changed.

## Checklist B — new op type on the expression chain

Follow ADR 0004's checklist; it is the authority. Condensed:

1. Pick a primitive name; add it to the registry header in
   `R/utils-pestore-ops.R` with its params and phase.
2. Decide which **carriers** can express it. Three on this chain — Arrow
   query, DuckDB `tbl_dbi`, collected triplet `data.table` — but unlike
   Checklist A there is no SQL arm to write: the two lazy carriers are both
   dplyr, so **one branch in `.pe_apply_op` written in plain dplyr serves
   both**. Write a separate executor only where a carrier cannot be shared:
   an elementwise expression on `value` shares one (`.op_transform_log`), and
   anything carrying per-axis state needs the collected-frame executor,
   because Arrow cannot index an R vector from inside a plan.

   Two traps this creates, both already paid for once:
   - **Write dplyr that lowers to every carrier, not just the one you
     tested.** `log1p()` is native to Arrow and to data.table, but DuckDB has
     no such function and dbplyr does not translate it, so it reached the
     engine verbatim and failed at `collect()`. `.op_transform_log` uses
     `log(value + 1)` for that reason. A branch that works on Arrow is not
     thereby correct.
   - **A payload has to live on the engine it joins into.** Arrow cannot read
     a `tbl_dbi` and DuckDB cannot read an Arrow `Table`;
     `.pe_payload_carrier()` is the one sanctioned place that branches on
     carrier. Reach for it only for payload residence — a branch per engine
     anywhere else is how `"query"` and `"duckdb"` start returning different
     values.
3. Key payloads **by on-disk id, per source** (ADR 0003) — a named list of
   per-substore vectors, not one stacked vector. This is what makes payloads
   invariant under `[`. Note that `cell_idx`/`gene_idx` are **not** guaranteed
   sorted (`[` preserves caller order), so derive bounds with `min`/`max`, not
   first/last.
4. Push with `.pe_push_op(pe, record, phase = "lazy" | "post")`. It enforces
   the monotonic rule: once `@post_ops` is non-empty, nothing may go on
   `@ops`. That rule is **sequencing, not a defect** — do not remove it.
5. **If the op's meaning depends on the open window, bake it at push time**
   (ADR 0006). An op defined by a population statistic over the current view —
   as `multiply` from `libraryNormParam` is, its factors being
   `scalefactor / colSums` — must run its aggregate in the producer and freeze
   the answer into the record. From then on the record is a pure function of a
   value and an axis id, so `[` can narrow freely without invalidating it.
   Consequences to accept deliberately: the producer is an eager full pass, not
   a queued recipe, and the frozen payload rides on the store through
   `saveRDS` and every subset.

   An op that *cannot* bake — one that would have to re-derive its statistic
   against whatever window is open at read time — **does not belong on the
   chain.** The escape is materialization: `storeWrite` bakes the chain into
   on-disk values and returns a store with empty chains, after which the
   statistic can be taken over the new population by re-running the verb.
   Re-running the producer is how a caller asks for a statistic over a new
   population; subsetting is not.
6. Have the producing **verb** append it (`processData(·, someParam)`), and
   check which consumer verbs will encounter it: `reduceData`, `analyzeData`,
   `filterData`, `storeRead`, `storeWrite`. Every consumer reaches the chain
   through `storeRead`.
7. A store's `@uid` must match the on-disk `source_id` partition. Minting a
   fresh store from a path via `new()`/`initialize()` gives a new uid, and
   `source_id`-keyed payload filters then silently return zeros. Assert
   non-empty results against independent ground truth — identical-looking
   output is not proof of correctness.

## Checklist C — new plain method

1. Pick the tier (question 4). Prefer the virtual.
2. Return a **store**, not data, unless the verb's contract is
   materialization. GiottoDisk's readers expose only a gobject or a store —
   never a transient `data.table`/arrow object, which would drop the file
   handle.
3. Never materialize to answer a metadata question. Counts go through
   `COUNT(*)` (`nrow()` returns a double for 2^53 headroom, uncached);
   extents through the `ext()` helpers; IDs through an Arrow `DISTINCT`.
4. Do not mutate the input. Record on a copy and return it; tests assert the
   original's `@ops` is untouched.
5. Roxygen pattern — a shared doc block terminated by `NULL`, then one
   `@rdname` + `@export` per method:
   ```r
   #' @name myVerb
   #' @title …
   #' @description …
   #' @param x `parquetBase`-inheriting store object
   #' @returns …
   NULL

   #' @rdname myVerb
   #' @export
   setMethod("myVerb", signature("parquetBase"), function(x, ...) { … })
   ```
6. Internal helpers get the matching prefix and are never exported:
   `.pbase_*` (parquetBase), `.pstore_*` (parquetStore), `.pgeom_*`,
   `.pe_*` (parquetExprStore), `.arrow_*` / `.dplyr_*`, `.gdsrc_*` (gsource).
7. New file? Add it to **`Collate:` in DESCRIPTION** and put
   `#' @include class-….R` + `NULL` at the top for class dependencies. Collate
   is explicit and hand-maintained; a missing entry breaks the build.
8. Re-document, then confirm NAMESPACE gained the expected
   `exportMethods()` — not `export()`. A method appearing as `export()` means
   the roxygen block is on a function, not a `setMethod`.

## Checklist D — new store class

1. Extend the existing hierarchy; do not start a parallel root. Slots
   `@fields`/`@ops`/`@post_ops` come from `parquetBase`; spatial slots from
   `parquetGeomBase`.
2. Define `specialCols()` for it — that is the single declaration of which
   columns are enforced on disk and hidden from `colnames()`.
3. Implement `storeRead` by delegating the op fold to
   `.pbase_storeread_processing`, splitting only at how `atab` is obtained.
4. Implement `show` via the `.show_info` pattern (raw values +
   `GiottoUtils::print_list`), and `storeExists`/`storePaths`/`storeUID`.
5. Wire the gsource verbs: `sourceContains` / `sourceAdopt` dispatch, and
   `depends` if the class references other artifacts — `sourcePrune` protects
   transitively by BFS over `depends`.
6. `exportClasses` only if users construct it directly.

## Anti-patterns

Each of these has produced a parallel implementation in this codebase or a
near miss recorded in an ADR.

- **Re-lowering ops by hand.** Taking `output = "query"` and re-applying the
  records at the call site. Nothing keeps the two copies in step. Use a chain
  editor. (ADR 0004, invariant 5.)
- **Collect-then-transform.** `storeRead(output = "tibble")` inside a method
  that returns a store. Every downstream caller loses laziness.
- **Implementing the Arrow arm only.** On the tabular chain the duckdb/sedona
  compile warns and skips unknown ops — results silently differ by output
  format. On the expression chain the shape differs: one dplyr branch serves
  both lazy carriers, so the risk is not a missing arm but a branch written in
  dplyr that only *lowers* to Arrow (`log1p` being the case that bit). Test
  every `output` you claim either way.
- **Giving `storeRead` engine-side state.** A parameter that names, addresses,
  or reuses a database view or temp table across calls. `storeRead` derives a
  scan from store state and returns a handle; it is not a session manager, and
  the engine is not the orchestrator. `conn` is the sole exception because the
  caller already owns it. A consumer needing a prepared table it modifies
  iteratively owns that itself, as a contained optimization. See AGENTS.md,
  "Output Formats".
- **A union copy of a single-store method.** If the logic is the same, it
  belongs on `parquetBase`/`parquetExprBase`.
- **`setGeneric` on a generic that already exists upstream.** Masks it.
- **Naming an op after its intent.** Lets a producer claim ownership of its
  own output and rewrite it — which is only sound if nothing was inserted in
  between, and the chain cannot promise that.
- **Positional row indexing.** Not supported, by decision (ADR 0001) — a
  2B-row index costs ~16GB. Use `subset()` for value filters, `rowSample()`
  for downsampling.
- **Aggregating string columns in Arrow at scale.** Join string IDs to
  integers *first*, then aggregate the integers — aggregating strings leaves
  dangling `utf8_view` buffers. An Arrow constraint, not a preference.
- **Materializing a range index.** Index with `a:b`, not `seq.int(a, b)` —
  ALTREP compact seq becomes one hyperslab; a materialized vector becomes a
  point selection (measured 4.4x time, 3.9x RSS).
- **Scope creep.** Change what was asked. A targeted op addition is not an
  invitation to refactor the fold.

## Tests

`tests/testthat/test-<area>.R`, one behaviour per `test_that`.

- Assert **laziness**: `expect_length(store2@ops, 1L)` and that the original
  is unchanged (`expect_length(store@ops, 0L)`).
- Assert `output = "query"` returns an `ArrowObject`, and that the `"tibble"`
  round-trip matches the in-memory equivalent.
- If you added an op, test it on **both** a single store and a union, and on
  every `output` you claimed to support — that is what catches a missing SQL
  arm.
- `skip_if_not_installed()` for anything in Suggests: `sedonadb`, `duckdb`,
  `Giotto`, `GiottoData`, `BPCells`, `HDF5Array`.
- Verify against independent ground truth, not against another run of the
  same path.

## Doc obligations

Per `adr/README.md`, the boundary is sharp:

| Change | Where it goes |
|---|---|
| A new invariant future code must respect | `AGENTS.md` |
| How the mechanism works and hangs together | `vignettes/articles/design.Rmd` |
| Something intended but not built | `vignettes/articles/roadmap.Rmd` |
| Why this over the alternatives, and when | a new `adr/NNNN-*.md` |
| A naming convention | `AGENTS.md` "Conventions" — not an ADR |

An ADR is owed when a decision **constrains future code, was contested or
non-obvious, or has a cost worth remembering** — not for bug fixes,
behaviour-preserving refactors, or naming. Accepted ADRs are immutable except
for `Status` and added links; change one by writing a superseding ADR and
updating the index table. `adr/README.md` also carries a **backfill
candidates** list — if your work touches one of those decisions, that is the
moment to write it up.

Branch target: GiottoDisk work PRs to **`dev`**. (Other suite packages differ
— on-disk work there targets `@gsource`.)
