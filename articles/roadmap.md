# Roadmap

High-level direction for GiottoDisk. Specific features land in the
[NEWS](https://giotto-suite.github.io/news/index.md) when implemented.

## In active design

- **Disk-backed mutable metadata** — sidecar-based column groups for
  cell/feature metadata at transcript scale (~50M+ rows).
- **Cross-store query dispatch** — ephemeral SQL engines (DuckDB,
  sedonadb) for joins, transforms, and spatial queries spanning multiple
  on-disk artifacts, exposed through a read-only `giottoView`
  abstraction on the GiottoClass side.
- **`bpcMatrixStore` as a wrapping store** — see below.

### Lowering peak memory in the pairwise marker tail

`analyzeData(parquetExprBase, scranMarkersParam)` builds one statistic
frame per ordered pair of groups — three dense vectors of `n_genes` —
and hands each host’s set to `scran::combineMarkers()`. Since the frames
are built and released per host, the **live** set is `O(genes × G)`
rather than `O(genes × G²)`; cell count does not enter at all.

That does not translate into lower peak RSS on its own. Measured at 5000
genes × 40 groups, one shape per process:

| shape | ΔRSS | time |
|----|----|----|
| all-at-once, `G(G−1)` frames live | 1388 MB | 10.5 s |
| per-host, `G−1` frames live | 1365 MB | 11.4 s |
| per-host + [`gc()`](https://rdrr.io/r/base/gc.html) each host | 688 MB | 28.2 s |

R grows its heap to cumulative allocation and does not return freed
frames without a collection, so bounding the live set is necessary but
not sufficient. Forcing [`gc()`](https://rdrr.io/r/base/gc.html) per
host halves peak at ~2.7× the runtime — a real trade, and the reason it
is not the default. **Do not add an unconditional
[`gc()`](https://rdrr.io/r/base/gc.html) to that loop**; if it ever
becomes worth exposing, it should be an option with these numbers
attached.

Cheaper first move for anyone actually memory- or time-bound at high
cluster counts: `comparison = "one_vs_rest"`, which runs `G` pooled
two-group tests instead of `G(G−1)` pairwise ones.

### A moments seam in scran, to retire our transcription

`R/stream-markers.R` hand-rolls roughly 60 lines of scran’s pairwise
tail — `.pe_welch` (Welch standard error and Satterthwaite d.f.),
`.pe_run_t` (two one-sided [`pt()`](https://rdrr.io/r/stats/TDist.html)
calls), `.pe_choose_lr`, `.pe_logBH`, `.pe_full_stats`. Only
`scran::combineMarkers()` is reused, because it is the one part that
never touches expression values.

That exists because `scran::findMarkers(test.type = "t")` needs a
materialized matrix. It makes one C++ pass for per-(feature, group)
moments and then does every pairwise comparison as arithmetic on a
features × groups table — but exposes no way to hand it precomputed
moments, so a backend that cannot produce a matrix has to reimplement
that arithmetic.

**The risk is drift, and it is guarded but not enforced.** The
transcription is pinned two ways: `.pe_welch` / `.pe_run_t` against
[`stats::t.test`](https://rdrr.io/r/stats/t.test.html), `.pe_logBH`
against [`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html), and
the whole path elementwise against `scran::findMarkers()` at 1e-10. But
the parity test sits behind `skip_if_not_installed("scran")`, so it only
fires where scran is installed. If scran changes its statistic and that
test has not run, the two backends diverge silently.

**The fix belongs upstream.** scran already makes this split internally:
`.test_block_internal()` computes moments via `.compute_mean_var()`,
then hands a `STATFUN` closure to `.pairwise_blocked_template()`.
Exposing that boundary — a `pairwiseTTests`-style entry point taking
per-group `n` / `mean` / `var` instead of `x` — would delete our
transcription outright, and would serve every other out-of-core backend
facing the same wall: BPCells, DelayedArray, dbMatrix, TileDBArray.

Worth recording why this is the ask rather than the obvious one:
teaching scran a parquet schema would mean a C++ tatami backend, and it
would still block-realize rather than hash-aggregate — no win on the
axis that matters. The moments seam is a far smaller change with a far
wider payoff.

Until then: keep the parity test green, ensure scran is installed
wherever the suite runs, and treat any divergence from
`scran::findMarkers()` as a bug in the transcription rather than a
reason to adjust the expectation.

### An `add` op for centred display values

`multiply` is implemented; its counterpart `add` is registered as a
refused stub. Both triplet executors reject it, and no verb emits one.

The asymmetry is real: multiplying preserves sparsity, so it lowers to
Acero and applies over a collected triplet frame alike. Adding does not
— every implicit zero becomes the offset — so it cannot be expressed
over triplets at all.

Recording it is still cheap; only materialising it is not. That is the
bargain the rest of the suite already strikes: `BPCells::add_rows()`
appends a node to the op graph and `ScaledMatrix(center=)` keeps the
sparse seed, and neither expands anything until a chunk is pulled
through. Both are reached via
[`GiottoClass::standardise_flex()`](https://giotto-suite.github.io/GiottoClass/reference/standardise_flex.html),
which is the shape a `parquetExprStore` branch should match — today a
pestore matches neither branch and falls through to arithmetic it cannot
do, which is why `zscoreScaleParam` is refused upstream.

At materialisation the frame densifies **into triplet form**:
`CJ(row_id, col_id)` over the slice, the stored nonzeros joined onto it,
`NA` filled with 0. The offset is then an ordinary mutate, and every op
after it keeps working on triplets in the same phase — no matrix-tier
executor, no chain split. What does change is payload size: from the
`add` onward the frame is dense, so the read window has to be sized
against density 1.0 rather than the store’s actual fill.
`.pe_window_cells()` derives the window per read, so it can account for
that.

**Scope, so this does not get over-built.** The only consumer is
*display* — the `"scaled"` expression slot behind heatmaps and similar.
Analysis never needs it: every path that mathematically requires
centring already folds it into algebra instead of materialising, in
Halko’s rank-1 term, the gram path’s `n·μμᵀ`, and the Pearson residual
zero-block. So the slice is small by nature, and `.pe_check_dgc_dims()`
already refuses anything large. A straightforward `CJ` expansion is
sufficient; it does not need to be fast.

### `bpcMatrixStore` as a wrapping store

Sequenced after the scstream feature-adoption work.

Today `bpcMatrixStore` is `contains = "fileStore"`: a path handle with
no matrix slot, and
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
is a one-liner returning a bare `BPCells::open_matrix_dir(store@path)`.
The store is therefore *discarded on read* — everything downstream holds
a raw `IterableMatrix`. Three problems follow, and all are worked around
rather than solved:

1.  **Identity.** An `IterableMatrix` carries no stable identity slot
    (BPCells assigns uid by manifest registration, not on the handle).
    So `.ss_gdsrc_detect_uid_matrices()` has to keep two branches —
    `direct_uids` for the `parquetExprStore` family, `hash_mats` for
    BPCells / HDF5Array — and `.ss_hash_expr_base()` derives identity by
    *opening the matrix* at each leaf
    (`.hash(BPCells::open_matrix_dir(d))`).
    [`sourceAdopt()`](https://giotto-suite.github.io/GiottoDisk/reference/sourceAdopt.md)
    re-hashes at both the old and new path.

2.  **Dispatch — the `ANY` catches are dead code.** Verified with
    `selectMethod()` (BPCells + Giotto + GiottoDisk all loaded):

    | signature asked for | resolves to |
    |----|----|
    | `IterableMatrix`, `libraryNormParam` | `allMatrix`, `libraryNormParam` |
    | `IterableMatrix`, `logNormParam` | `allMatrix`, `logNormParam` |
    | `IterableMatrix`, `binarizeThreshParam` | `allMatrix`, `binarizeThreshParam` |
    | `IterableMatrix`, `minmaxThreshParam` | `allMatrix`, `minmaxThreshParam` |
    | `IterableMatrix`, `filterParam` (`filterData`) | `IterableMatrix`, `filterParam` |

    `IterableMatrix` is a subclass of Giotto’s `allMatrix` union, and
    `allMatrix` is **more specific than `ANY`**. So the two
    `processData` catches in `methods-giotto.R` — commented “catch for
    ‘suggested’ backends (can’t dispatch on class)” — never fire for an
    `IterableMatrix`, and `BPCells::binarize` / `min_scalar` are
    unreachable. Their premise is wrong whenever an upstream package has
    a union method covering the class.

    Note the last row: direct dispatch on `IterableMatrix` **works**
    (`existsMethod()` is `TRUE`) once BPCells is attached. The
    `no definition for class "IterableMatrix"` warning at install time
    is cosmetic. So the fix is not “avoid dispatching on the class” — it
    is “dispatch on a class specific enough to beat the upstream union”,
    which is exactly what a GiottoDisk-owned wrapper provides.

3.  **No collapse point for the lazy op graph.** Nothing in the workflow
    ever materializes a BPCells matrix, so `svds` walks the whole
    normalization stack on every Lanczos matvec. Measured on Atera (170k
    cells, 2000 HVGs, 30 PCs): **122.3 s** on the lazy stack vs **14.7
    s** reading a matrix baked first, where the bake itself costs **1.1
    s**. Any later consumer of the expression slot re-pays the same
    cost, since nothing collapses it.

    The hooks that *would* do this exist but miss BPCells:

    - `setExpression()` (GiottoClass `slot_accessors.R`) takes
      `write = FALSE` and writes via
      [`GiottoDisk::sourceWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite.md)
      when `inherits(mat, c("matrix", "Matrix")) || isTRUE(write)`.
      **That path is functional and native for BPCells** —
      `sourceWrite(gDirSource, ANY)` branches on
      `inherits(data, "IterableMatrix")` to `store_type = "bpcells"`,
      reaching `storeWrite(bpcMatrixStore, ANY)` →
      `BPCells::write_matrix_dir()` with no coercion (measured: 1.1 s
      for 2000 x 170k). It is simply **never invoked with `TRUE`**:
      there is no `write = TRUE` call site anywhere in Giotto,
      GiottoClass, or downstream analysis code, and `normalizeGiotto`
      reaches the setter through
      `setGiotto(gobject, norm_expr, verbose =, initialize = FALSE)`. An
      `IterableMatrix` also fails the implicit
      `inherits(mat, c("matrix", "Matrix"))` test, so nothing triggers
      it. That predicate is arguably inverted — it writes what is
      already cheap to hold in memory and defers what is expensive to
      defer.
    - `normalizeGiotto` already has an option-gated collapse for the
      *other* lazy backend, immediately above that call:
      `getOption("giotto.dbmatrix_compute")` → `.compute_dbMatrix()`,
      run once for `normalized` and once for `scaled`. There is no
      `IterableMatrix` equivalent.

    BPCells falls between the two: lazy like `dbMatrix` but with no
    compute hook, disk-backed so the `memory_matrix` test skips it.

**Proposal.** Give `bpcMatrixStore` a slot holding the `IterableMatrix`
and keep the wrapper on for the whole workflow. Division of authority:

- the **store** is authoritative for *location* (vault path, which
  adoption may have moved),
- the **held matrix** is authoritative for the *lazy op graph*.

Neither duplicates the other. Workflow verbs dispatch on
`bpcMatrixStore` — a class this package defines — and each method
unwraps (rebasing leaf `@dir` from the store’s path), applies the op
with BPCells semantics, then either re-wraps into the *same* wrapper
when the result is another `IterableMatrix` (a lazy step), or returns
the terminal value directly (`matrix_stats` → `data.table`, `svds` → PCA
list).

**This composes with machinery that already exists**, which is most of
why it is cheap:

- `.im_map_leaves()` / `.im_leaf_dirs()` already walk the op graph —
  down `@matrix`, fanning out over `@matrix_list` — and rewrite
  `slot(leaf, "dir")`. Leaf rebasing is a solved problem, already used
  by
  [`sourceAdopt()`](https://giotto-suite.github.io/GiottoDisk/reference/sourceAdopt.md).
- A lazy `IterableMatrix` round-trips `saveRDS` intact (op graph
  preserved, leaf holds a path, iterators are constructed on demand), so
  the held object survives snapshotting. The externalptr-nulling problem
  noted in `methods-snapshotSave.R` is arrow’s, not BPCells’.
- `.store_nostate()` already establishes that pending transforms are
  **not identity-bearing** — it strips `@ops` / `@post_ops` / subset
  indices so that “two references to the same on-disk thing always hash
  identically regardless of any pending filters”. The BPCells op graph
  is the direct analogue, so `uid` stays **stable across re-wraps**:
  nothing is minted per lazy step, and
  `.store_nostate("bpcMatrixStore")` is one more link in the existing
  chain (reset the held matrix to a bare handle on `@path`).

**What it buys.**

- *Identity* becomes direct `@uid`, collapsing the `direct_uids` /
  `hash_mats` split and removing hashing that currently has to open the
  matrix to compute it.
- *Dispatch* becomes real S4 on a class specific enough to beat the
  `allMatrix` union, so the BPCells fast paths stop being shadowed. It
  also puts a single owned method surface where the
  [`inherits()`](https://rdrr.io/r/base/class.html) branches are today.
- *A place to put a materialization policy.* This is the one the
  workflow needs most and the one nothing else can host: Giotto has
  nowhere to decide “this op stack has grown enough, write it and
  rebase”, because `standardise_flex()` is a GiottoClass primitive with
  no notion of storage and `runPCA` receives an already-stacked matrix.
  A store that owns both the path and the op graph is the natural home —
  the same gated-bake decision `.stream_gram_svd()` /
  `.stream_random_svd()` already make for parquet
  (`giottodisk.pca_bake_max_ratio`).
- *A place to put the thread count.* BPCells has **no global option** —
  `matrix_stats(threads =)` is the only parallelism knob on the
  expression path (`write_matrix_dir`, `svds`, `multiply_cols`,
  `colSums` have none), and it defaults to `0L`, which measures as
  single-threaded (1.52 s vs 0.31 s at `threads = 8`). Every call site
  must pass it explicitly, and none currently do — `standardise_flex()`
  and both `variable_genes.R` sites omit it. An unwrapping store can
  feed it from one setting.

**Open question — compound cardinality.** A compound matrix (`rbind` /
`cbind`) has N `MatrixDir` leaves while a store has one `@path`, so
“rebase from the store’s path” is ambiguous, and adoption currently
registers one manifest artifact *per leaf*. Either the store carries a
per-leaf path map, or compounds become N stores under a parent — the
`unionParquetExprStore` / `@stores` shape, which `.ss_hash_expr_base()`
and `.ss_store_uids()` already handle by recursing into `@stores`.
Following that precedent means neither helper needs a new branch.

**Migration.** Keep a `signature("bpcMatrixStore", "ANY")` catch-all
that unwraps and delegates to the existing path, so ops can move over
one at a time instead of needing full coverage before anything works.
That also keeps Giotto’s
[`inherits()`](https://rdrr.io/r/base/class.html) fast paths useful as a
fallback tier — they sit inside shared LOESS / thresholding logic that
should not be duplicated here.

**Interim fixes, if any of this is wanted before the redesign** (each is
independent and small):

- pass `threads` to `matrix_stats()` at the three call sites — 4.9x on
  that stage, one argument each;
- have `normalizeGiotto` forward `write = TRUE` on its two `setGiotto()`
  calls when the result is a lazy disk-backed matrix. This is the
  smallest fix: everything below the setter already works and is native,
  so it is one argument for ~1.1 s of write against ~107 s of PCA. Only
  covers gsource-managed projects, since the setter’s write path is
  gated on `gobject@source`;
- for the non-gsource case, an `IterableMatrix` branch beside the
  `dbMatrix` one in `normalizeGiotto`, mirroring
  `giotto.dbmatrix_compute`;
- move (or delete) the two shadowed `processData` `ANY` catches in
  `methods-giotto.R`, and check what the `allMatrix` method they lose to
  actually does to an `IterableMatrix`.

**Unrelated but adjacent:** Giotto and GiottoClass both call `BPCells::`
(`svds`, `matrix_stats`, `rowVars`, `add_rows`, `multiply_rows`) without
declaring BPCells in either DESCRIPTION. GiottoDisk is the only package
that declares it. That is an `R CMD check` violation today, independent
of this design; migrating ops here reduces but does not fix it.

## Backend variants

The `gsource` protocol is pluggable. Planned implementations beyond the
current local `gDirSource`:

- **Cloud-native storage** — object storage (S3/GCS/ADLS) for parquet
  artifacts with a remote catalog for metadata.
- **SpatialData interop** — read-only adapter over a published
  SpatialData Zarr store, exposing sdata elements through the gsource
  interface.

## Other items

- Exact polygon masking as a post-materialization op (currently
  AABB-only).
- Lazy virtual partitioning via filesystem hardlinks.
- Direct sedona → parquet write path (blocked upstream on Arrow R
  view-type compute kernels).
