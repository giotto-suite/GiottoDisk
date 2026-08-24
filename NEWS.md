# GiottoDisk 0.0.0.2

## bug fixes
- `createGiottoXeniumObject(backend =)` no longer errors on Xenium-format
  directories that ship no panel json. Feature metadata is generated from the
  expression matrix when the panel is absent.

## changes
- `Giotto (>= 4.2.4)` in `Imports:`, for the `AteraReader` class `R/convenience-atera.R` subclasses. Below that the failure is an S4 inheritance error at load rather than a version message.
- Grouped expression statistics — `analyzeData(x, featStatsParam, groups =)` and
  everything riding it, including scran marker detection — now window the scan
  by cells instead of running one arrow plan over the whole store. The aggregate
  is O(groups) either way, but a grouping puts a join in front of it whose output
  is O(nonzeros), and Acero does not spill; at atlas scale that was the failure.
  Windows are exact rather than approximate because the accumulators are
  additive, they are folded as they arrive so the retained state does not grow
  with window count, and a contiguous cell range prunes row groups (the store is
  sorted cell-major). Sized by `.recommend_chunk_size()` against free RAM, so a
  store the budget already covers is one window and behaves exactly as before.
  Ungrouped statistics are unchanged. Callers batching the **feature** axis to
  work around the old memory cost should stop: gene ids are not the sort key, so
  every batch rescanned the store in full and the cost was linear in batch
  count, not in genes per batch.
- Internal refactor, no user-visible behaviour change: cell-window
  streaming now has one seam — `.pe_windows()` (substores x their cell
  ranges), `.pe_chunk_ranges()` (a sub-range of one substore, for the parallel
  PCA band workers) and `.pe_window_store()`. The walk had been hand-rolled in
  nine places — both statistic accumulators, the `storeWrite` bake, four PCA
  passes and the band split — and the copies had drifted. All now route through
  it, and the seam is recorded in `AGENTS.md` and the `giottodisk-method` seam
  table so the next windowed verb attaches instead of copying. Verified against
  the previous implementation at every site: bitwise for PCA (`u`, `d`, `v`,
  `sdev`, `eigenvalues` and per-column magnitudes — a correlation check cannot
  see a scale change) and for the `storeWrite` bake, and to within 1 ULP for the
  float statistic accumulators, where eager folding reassociates the summation
  (see below). Integer accumulators are exact.
- One behaviour change comes with it: the R-side accumulator
  (`.pe_accum_chunked_dt()`, the path taken when `@post_ops` cannot be lowered)
  now folds each window's partial as it arrives instead of collecting one per
  window and reducing at the end. Held state drops from `O(groups x windows)` to
  `O(groups)`, so tightening the window no longer costs memory. The Acero path
  already did this.

  The one visible consequence: folding on arrival reassociates the summation, so
  a float accumulator can differ from the old reduce-at-the-end result by ~1 ULP
  (measured 1.0-1.2 ULP, max relative 2.7e-16). Counts are unaffected. Results
  are equal to tolerance, not bitwise, and comparisons across different window
  counts should be written that way.
- `analyzeData(x, featStatsParam, groups =)` resolves a grouping by `cell_ID`
  when it is named or factored by one. A per-cell vector is a payload, and
  adr/0003 keys those by on-disk id: keyed by view position it reads the wrong
  entries once `[` has narrowed the store. Cells the grouping does not name now
  drop, matching `.pe_axis_pos_map()`; no overlap at all errors. An unnamed
  vector stays positional against the current view, with a warning.
- Stereo-seq is reachable from the public entry points. `Giotto`'s
  `importStereoSeq()`, `createGiottoStereoSeqObjectBin()` and
  `createGiottoStereoSeqObjectCell()` gained a `backend =` argument that routes
  to `importStereoSeqDisk()`, mirroring what `createGiottoXeniumObject()`
  already did. `StereoSeqDiskReader` and the GEF inputs existed before this but
  had no caller. Requires `Giotto@gsource` at or past the matching commit.
- `binGefInput()` addresses `geneExp/<bin_size>/` with the group key **as
  given** (`bin100`), where it previously stripped the `bin` prefix and looked
  under `geneExp/100/`. Real GEFs carry the prefix — Giotto's in-memory reader
  hardcodes `geneExp/bin1/expression` — so the old key found nothing on an
  actual file. Callers passing a bare `"50"` must now pass `"bin50"`.
- `importStereoSeqDisk()` defaults now match `Giotto::importStereoSeq()`:
  `bin_size` is `"bin100"` (was `"bin50"`) and `gef_type` for `type = "cell"`
  is `"adjusted_cellbin"` (was `"cellbin"`). Adding `backend =` to a working
  call no longer changes which file is read.
- `StereoSeqDiskReader`'s `create_gobject()` took `gef_path` and `mask_path`
  through recursive default argument references (`gef_path = gef_path`), so
  neither auto-detected path ever reached the loader and ingest failed with
  "no .gef path provided". Both now resolve.
- `StereoSeqDiskReader`'s `load_expression()` returns `list(exprObj)`, matching
  the in-memory `StereoSeqReader`. It previously returned a bare `exprObj`.
  Both work through `setGiotto()`, so assembled objects were unaffected, but a
  reader driven piecewise — as the Stereo-seq importer vignette does — has to
  be substitutable with `backend =` set or unset.
- upstream (`Giotto@gsource`): cellBorder polygons from
  `.stereoseq_build_polygons_from_border()` now populate `unique_ID_cache`, as
  every other polygon constructor does. Left at the prototype `NA_character_`,
  the IDs get recomputed downstream with `unique(<spatVector>$poly_ID)`, which
  fails on a backend-managed giotto because `setGiotto()` has by then swapped
  the `SpatVector` for a `parquetGeomStore`. This is what made
  `createGiottoStereoSeqObjectCell(load_polygons = TRUE, backend = ...)` — the
  vignette's recommended default — error with "unique() applies only to
  vectors".
- both GEF inputs now sum records for duplicate gene names that sit in
  different chunks. A gene table is ordered by geneID, so two rows sharing a
  geneName scatter arbitrarily (a mouse `tissue.gef` has 16 such names, up to
  25400 rows apart), and `.gef_safe_chunks()` only keeps *consecutive* runs
  together. Their records were aggregated per chunk and written separately,
  so the store held two rows for one `(cell, gene)` pair — 753 of them on that
  file, inflating nnz and every marginal derived from it. Records for
  duplicated columns are now held back and flushed as one aggregated batch at
  end of stream, bounded by the duplicated genes rather than the matrix.
- for `type = "bin"`, spatial locations are built from the `(x, y) -> bin_ID`
  map accumulated during the expression stream rather than by re-reading the
  gef. Bin coordinates live inside the expression records, so the inherited
  in-memory closure had to pull the whole `geneExp/<bin>/expression` dataset
  into memory — the exact read the disk backend exists to avoid. Cellbin is
  unchanged; its coordinates come from the small `cellBin/cell` table.
- `analyzeData(parquetExprBase, varParam)` now evaluates the same Pearson
  residual as `Giotto`'s in-memory path: negative-binomial denominator
  `sqrt(mu + mu^2/theta)` with `theta = 100` (was Poisson, `sqrt(mu)`) and
  clipping to `±sqrt(n)` (was unclipped). Both omissions inflated the
  variances, so **the features selected by `calculateHVF(method =
  "var_p_resid")` change**: on a Stereo-seq cellbin sample the streaming and
  in-memory results now agree exactly, where before they shared 3% of their
  selections. Verified against an independent dense reference at
  `theta = 100`, `10` and `1e6` to 6e-15. `theta` is settable through
  `analyzeParam("var", theta = )`.
- the result gains a `mean_expr` column, for the mean-versus-variance
  diagnostic in `calculateHVF()`'s plot. Free: the gene totals are already
  computed.
- with a finite `theta` the all-zero block no longer collapses to a per-gene
  scalar (`sum_j z^2 = g_i` held only for Poisson), so it is summed over cells
  explicitly. To keep that affordable the cell totals are collapsed to their
  **unique values with multiplicities**, since the term depends on a cell only
  through `mu_ij = g_i c_j / T`. Measured on `C04687E314.tissue.gef`: bin1 has
  60 unique totals across 5,043,144 bins, turning a 1.3e11 product into 1.6e6.
  The gene axis is chunked so the intermediate stays bounded when a dataset
  has many distinct totals (cellbin: 2,789 across 7,527 cells).
- `DESCRIPTION` now carries a `Remotes:` field pinning the upstream development
  branches this package builds against (`GiottoClass@gsource`, `Giotto@gsource`,
  `GiottoUtils@dev`, `drieslab/tilework`), so `remotes::install_github()` and
  `pak` resolve them without a manual install order. Note that Giotto's general
  integration branch is `suite_dev`, not `dev`, and does not satisfy this
  package. Why each pin exists, and what has to land upstream before it can be
  dropped, is in *Upstream branch pins* in `AGENTS.md`.
- `parquetExprStore` no longer has a `@chunk_size` slot. The streaming window is
  derived per read instead of stored, because it depends on free RAM — a
  property of the machine doing the reading, so a value baked in at write time
  was stale as soon as the store moved. Two options steer it:
  `giottodisk.chunk_ram_frac` scales the RAM budget every streaming pass works
  from (default 0.25), and `giottodisk.chunk_size` pins an absolute window.
- new `@stats` slot on `parquetExprStore`: per-axis marginal nonzero counts,
  keyed by on-disk id and filled by `storeWrite()`. Invariant under `[`, so a
  view's nonzero count is exact on either axis without touching the data.
  `parquetExprStore(scan_stats = TRUE)` fills them for a handle attached to
  Parquet written elsewhere; without them, consumers fall back to counting.
- `analyzeData(featStatsParam)` and `analyzeData(cellStatsParam)` now apply the
  store's op chain. They previously read the stored Parquet directly and
  ignored any queued normalization, so **results change for a normalized
  store** — they were reporting statistics on unnormalized values.
- on those same verbs, `detection_threshold` no longer reduces `total_expr` or
  `mean_expr`. It gates `nr_cells` / `nr_feats` and `mean_expr_det` only,
  matching Giotto, where the threshold selects which entries count as detected
  and never modifies a magnitude that participates.
- `processData(x, libraryNormParam(...))` appends a scaling record rather than
  rewriting an earlier one. Re-running now composes: the same `scalefactor` is a
  no-op, and a new one applies the ratio. Previously it replaced in place, which
  rewrote the earlier record underneath any intervening step.
- op records renamed: `norm_libsize` is now `multiply`, carrying an `axis`
  (`"cell"` / `"feat"` / `"all"`). The old name described the verb that produced
  it rather than the operation. `norm_libsize_log` was earlier split into
  independent `norm_libsize` and `log` records — either can be used alone, in
  either order.
- `sc_recommend_chunk()` is no longer exported. Its value had no user-facing
  destination once `@chunk_size` was removed; use `storeChunkInfo()` to see what
  a store's windows come out to.

## new
- `storeChunkInfo()` reports how a store's streaming windows are chosen — the
  view's shape and density, whether marginals are cached, detected free RAM, and
  the resulting window across RAM budgets for both read shapes (a chunk landing
  in a sparse matrix versus a collected triplet frame, which differ ~4x in
  bytes per stored value).
- Halko PCA accepts `scale = TRUE`, applied without densifying.
- an `add` op is registered as a refused stub for a future centred-display path;
  both executors reject it, and nothing emits one.
- `adr/` — architecture decision records, for why a choice was made and what was
  rejected. See `adr/README.md`.

## bug fixes
- library normalization derived its scale factors from raw column sums even when
  a `log` preceded the norm in the chain, so the factors were for the wrong
  quantity and normalized columns did not sum to the scale factor. They are now
  taken from the values the record actually multiplies.
- the gram-eigen to Halko fallback passed `ncp` where `k` was expected and
  errored for every input, so the fallback never worked.
- `analyzeData(featStatsParam / cellStatsParam)` on a `unionParquetExprStore`
  ignored the parent's op chain entirely, since union substores carry no ops by
  constraint.
- stores produced by `calculateOverlap()`'s COO path were returned without
  cached marginals, because they write Parquet directly rather than through
  `storeWrite()`.

## enhancements
- the per-axis statistic verbs (QC feature/cell stats, HVF) share one grouped
  accumulator pass instead of three near-duplicate aggregates, and a union runs
  as a single Acero plan rather than one plan per substore. Measured 2.1x at
  four substores and 5.1x at sixteen; 11.3x end to end on union `cov_loess`.
- the R-side statistics path reads in cell windows rather than collecting a
  whole store, so it is bounded by the window instead of by the data.
- `filterData()` routes through the same accumulator, picking up the union
  speedup and dropping its serial per-substore loop.
- Halko and gram-eigen PCA no longer require a normalization recipe or
  `feats_to_use`; both run on whatever the store holds.
