# GiottoDisk 0.0.0.4

## breaking changes
- `calculateOverlap()` on disk stores picks its engine the way `spatRelate()`
  does: `engine = NULL` (the new default) reads
  `options(giottodisk.spatial_query_engine)`, and unset or `"auto"` takes the
  first installed of sedona, duckdb and terra. It was always terra before. A
  call that supplies a terra tiling argument (`threshold`, `tiles`, `pad_y`,
  `poly_buf_factor`, `tile_idx`, `prune_tiles`) still gets terra under auto.
- `calculateOverlap()` on disk stores returns its overlap in a
  `queryableStore` instead of a `parquetStore`. The files no longer carry a
  `row_index` column; a row is keyed by `(pt_tile_index, pt_row_index,
  poly_ID)`. Projects saved with an older GiottoDisk drop their disk overlaps
  on load, with a warning; rerun `calculateOverlap()` to rebuild them.

## new
- `calculateOverlap(engine = "sedona")` (sedonadb >= 0.4): one spatial join
  over both stores. On a 624M-point, 170k-cell Atera sample it takes 39 s,
  against about 12 minutes for `engine = "terra"`.

## changes
- `calculateOverlap(engine = "duckdb")` reads both stores through
  `storeRead()`, so filters pending on either store now apply (it used to
  read the raw files and count filtered-out cells and features). It no
  longer numbers the rows in a single pass, which ran the whole join on one
  thread: about 4 minutes on the Atera sample instead of about 12.
- `overlapToMatrix()` to a `parquetExprStore` runs as one DuckDB query when
  duckdb is installed: 17.5 s for the Atera sample's 300M values, against
  80 s through Arrow, which remains the path without duckdb.
- `overlapToMatrix()` drops overlap rows outside the feature and cell ID
  universes instead of writing them with missing keys, so narrowing an
  overlap's `@spat_ids` / `@feat_ids` narrows the matrix.

## bug fixes
- `overlapToMatrix()` on an `overlapPointDisk` read the overlap with the
  input stores' ID column names. The overlap files always name them
  `poly_ID` and `feat_ID`, so non-default `poly_id_col` / `feat_id_col`
  failed.

# GiottoDisk 0.0.0.3

## new
- `importVisiumHDDisk()`: disk-backed Visium HD reader (binned and segmented
  outputs), also reached through `backend =` on `Giotto::importVisiumHD()`
  and `createGiottoVisiumHDObjectBin()` / `createGiottoVisiumHDObjectCell()`. It
  gives the same barcodes, features, spatial locations and counts as the
  in-memory reader; the 2 um bin points stay an in-memory `giottoBinPoints`.
  `create_gobject()` takes several `bin`s, a `gobject` to add into, and
  `load_bin_mapping = TRUE` to record each unit's parent units (from
  `barcode_mappings.parquet`) as cell metadata. Only extracted output
  directories are read.
- Gram-eigen PCA (`gramEigenPcaParam`, and `method = "auto"` when it resolves
  to it) runs both of its passes in the compiled `GiottoKernels` package when
  that is installed (optional, in `Suggests`): one scan of the HVF store each,
  with threads inside the call rather than forked R processes, so PCA no
  longer forks and runs in Positron and on Windows. On a 169,420-cell x
  2,000-feature store, 50 components, the whole call takes 7.3 s on one
  thread and 6.2 s on eight, against 16.9 s serial and 8.0 s on eight forked
  workers before. Threads follow
  `giottodisk.par_workers` / the future plan when above one, else
  `options(gkernels.n_threads)`. Without the package, or with
  `options(giottodisk.use_kernels = FALSE)`, the R path runs as before.
- `Compare` methods (`>`, `>=`, `<`, `<=`, `==`, `!=`) for `parquetExprBase`
  against a numeric scalar. `x >= t` queues a lazy indicator on the op chain
  -- stored entries that pass read back as 1, the rest are dropped -- so
  `rowSums(x >= 1)` and `colSums(x > 0)` stream through the existing margin
  methods, and code written for an in-memory matrix (`rowSums_flex(x >= t)`)
  runs unchanged on a backed one. The comparison sees values after anything
  already queued. A comparison that is `TRUE` at 0 (`x >= 0`, `x < 1`) would
  be dense and is an error; so is comparing two stores.
- `storeRead(<parquetEdgeStore>, output = "arrowstream")` returns a
  `nanoarrow_array_stream` over the same query `output = "arrow"` builds,
  with `@ops` already applied. It exists for readers outside R: the batches
  can be pulled across the Arrow C Data Interface, so a consumer holding
  only the stream still observes a pending subset. The alternative a caller
  otherwise reaches for -- taking `@path` and opening the parquet directly --
  reads the files as they sit on disk, which a pending subset is not, and
  also assumes one file per subdirectory. Respects `minimal`, because a
  stream cannot be reshaped once handed over. The stream reads once;
  narrow the store before asking for it. {nanoarrow} moves to Imports.

## changes
- `storeRead(output = "duckdb")` and `storeRead(output = "sedona")` read
  `source_id` and `tile_index` from hive partition discovery on the store
  root instead of reading each tile directory separately and injecting the
  values as SQL literals. Tile selection (`tile_idx`, `@tile_filter`) becomes
  a partition predicate that both engines prune files on. The partition
  columns keep their types (`source_id` string, `tile_index` int32). The
  per-tile reads crashed sedona on stores with a few hundred tiles or more.
- sedonadb >= 0.4.0 is required, for its hive partition discovery. An older
  install is passed over by `"auto"` spatial engine selection, and
  `storeRead(output = "sedona")` errors on it.

## bug fixes
- The Stereo-seq GEF readers now write cell-major stores, like every other
  input. `cellbinGefInput` reads the cell-major `cellExp` copy; `binGefInput`
  reorders through a temporary spill and positions bins in grid order, keeping
  their first-appearance `bin_<id>` names. Before, each part file held one
  gene batch across all cells, so cell-windowed passes (markers, grouped
  `featStats`, PCA) rescanned the store once per window. Stores written earlier
  keep the old layout until re-imported.
- Parquet writes no longer store data.table's `sorted` / `index` attributes.
  arrow restores R attributes on read; restored onto a store read back from
  several files or after a filter they were stale, and keyed or indexed
  subsets of the collected table could return wrong rows.
- `storeWrite()` of a `parquetExprStore` or `unionParquetExprStore` into a
  `parquetExprStore`, the last expression-store writer that was not
  cell-major, now writes that layout: each file sorted by cell then feature,
  files covering disjoint cell ranges. The sort
  was computed and then discarded by the parallel dataset writer, so the
  layout differed on every run. The write is also windowed by cell, so memory
  is bounded by the window rather than by the whole output (17.9 GB peak ->
  10.8 GB writing a 300M-value store). Stores already written keep their
  layout until rewritten.
- Tile stores written from a `queryableStore` or `parquetStore` (the disk
  readers' transcript and polygon paths) no longer nest a second
  `tile_index=000/` level inside every tile directory. Each tile took the flat
  geom store's `tile_idx = 0L` default on top of its own `tile_index=<NNN>/`.
  That broke `storeRead(output = "duckdb")` ("No files found") and, with
  sedonadb >= 0.4, which discovers hive partitions, `output = "sedona"`
  (ambiguous `tile_index`). `.write_parquet()` now errors instead of nesting a
  `tile_index` level. Stores already written keep the nested layout until they
  are rewritten (#74).
- `gDirSource()` given a relative path for a directory that did not exist yet
  kept the relative path, as did every artifact written to it. Reading those
  artifacts from a different working directory (a knitted document, a
  parallel worker, after `setwd()`) then pointed at the wrong place. The path
  is now made absolute once the directory has been created. Objects created
  before this fix keep the relative paths they were saved with.

## changes
- `cellStatsParam` and `featStatsParam` are imported from GiottoClass, where
  they moved from Giotto. Requires GiottoClass >= 0.7.4. No change in behaviour.
- **View resolution follows GiottoClass 0.7.3's `resolveRecipe()`** (the
  former `materialize()` / `resolveSubobject()`), which now requires
  GiottoClass >= 0.7.3. A view is evaluated once per resolve op into the
  surviving cell set by a `resolveKeep()` method, and each backed subobject
  queues that set as one `id_filter` without being handed the gobject.
  - **The set is a lazy arrow plan**, not a collected table: a backed filter
    owner contributes `subset(owner, pred)` as a query, crops contribute
    their cell_IDs, and the steps are chained with semi-joins, so nothing is
    read until a target is collected. Backed getters (`getCellMetadata(g,
    view = )` and friends) stay lazy this way; previously they took a
    separate per-target pushdown path, and `resolveRecipe()` collected an ID
    table up front. The in-memory `vector` form is collected only if an
    in-memory subobject inside a backed object asks for it.
  - A filter whose owner is a backed store now works: the old path called
    `store[, key, drop = FALSE]`, which no store method accepts.
  - A polygon store's own attributes can be filtered on (`region ==
    "tumor"`), with `poly_ID` answering for the cell axis, as `spatValues()`
    allows in memory.
  - Backed transcript points skip filter steps, as the in-memory leaf does;
    a crop still clips their own geometry. A crop drawn in a named space now
    clips in that frame: the region projection dropped translation (it put
    it in the affine's row 3, which `affine()` ignores), so a shifted space
    clipped the wrong area. It now uses GiottoClass's `project_region()`.
  - Feature metadata and dimension reductions inherit the in-memory leaf.
    Backed feature metadata used to be given a cell_ID `id_filter` it has no
    column for.
- Dropped the `prepareIds` import. GiottoClass removed the generic in 0.7.3:
  it was an exported identity transform with no call sites here or anywhere
  else, and `parquetCoordinator`'s own methods already promote an ID set to
  the form each store wants.

# GiottoDisk 0.0.0.2

## new
- `analyzeData(parquetExprBase, scranMarkersParam)` accepts
  `comparison = "nodes"`, the streaming half of Giotto's `findNodeMarkers()`.
  A `sets` list names the clusters on each side of every branch point of a
  cluster tree; the method takes its **one** grouped moment pass and folds each
  node out of it with the existing `.pe_pool_moments()`, exactly as
  `"one_vs_rest"` already does. The accumulators are additive, so every node
  after the first costs arithmetic rather than another scan -- 35 nodes on
  169,528 cells in 5.2 s, against ~35 s for one `findMarkers()` call per node,
  and sub-linear in node count.
- `parquetCoordinator`: view + space recipe resolution for gobjects whose
  subobjects are `parquetStore`-backed. Recipe steps are pushed onto each
  store's lazy-op queue via the existing `subset()` / `crop()` /
  `spatRelate()` dispatch, so no I/O happens at coordinator time and the
  engine (arrow / duckdb / sedona) is chosen per call.
  Inherits `dataTableCoordinator`, so in-memory subobjects inside an
  otherwise-backed gobject fall through to the in-memory path. Selected
  automatically for any gobject whose `@source` inherits `gsource`.
  Cross-storage narrowing is covered: a predicate whose columns live on a
  different subobject than the target resolves through a lazy `[`-join
  when the owner is backed, and an eager `id_filter` when it is in memory.
  Requires GiottoClass >= 0.7.0.
- `snapshotSave(gDirSource, giottoMulti)` and a `snapshotDelete` cascade
  to per-child snapshots.
- `spatIDs()` / `featIDs()` methods for `parquetGeomBase`, so a backed
  geometry answers the ID question through one dispatch rather than a
  column-forcing idiom at each call site.
- `spatRelate()` methods for `parquetGeomBase`, and
  `as.data.table(parquetGeomBase, geom = c("", "wkb", "XY"))` — the `"XY"`
  vertex expansion goes through `wk::wk_coords`, with no terra
  `SpatVector` intermediate.
- Zarr input for the Xenium/Atera disk readers. `importXeniumDisk()` /
  `importAteraDisk()` now work on zarr-only output directories (the only
  format Atera will ship): transcripts, boundaries and cell metadata are
  converted from the `.zarr.zip` archives to 10x-schema parquet in a
  fingerprint-keyed cache (convert once, reuse on later imports), and
  expression streams straight from the zarr into the vault
  `parquetExprStore` with no intermediate file. Requires `Rarr` and `zip`
  (Suggests). Unzipped `.zarr` directory trees are also accepted.
- `tenxZarrInput()`: `exprInput` over `cell_feature_matrix.zarr.zip`.
  Cell-ordered batch iterator over the feature-major CSC arrays; switches
  to bounded per-cell-block rescans when the triplet buffers exceed the
  RAM budget (or nnz exceeds int32, e.g. Atera whole-transcriptome runs).
- `xeniumZarrToParquet()`: standalone zarr -> parquet converter producing
  the 10x-shipped schemas (verified column-exact against 10x parquet,
  including `transcript_id` packing and instrument FOV names), plus the
  expression triplet layout. Transcripts convert in parallel over tile
  groups (`giottodisk.par_workers` / future plan).
- `detectZarrLayout()`: versioned layout detection for zarr output
  directories; unsupported layouts (zarr v3, unknown structure) fail with
  an actionable message.
- `analyzeData(parquetExprBase, pageEnrichParam)` runs PAGE enrichment on a
  disk-backed store without densifying it, so `Giotto::runPAGEEnrich()` works
  on a backed project. It used to fail outright.

  PAGE streams because its statistic is a function of moments rather than of
  the matrix. In memory it builds a dense `geneFold = expr - mean_gene_expr`
  and then only ever asks it for a per-cell mean, a per-cell sd, and a per-cell
  mean over each cell type's marker rows -- all three recoverable from per-cell
  sums of the stored sparse values plus the per-gene reference vector, so the
  dense matrix never has to exist. Three full passes plus one per cell type,
  each of the latter reading that type's marker rows only.

  Scores match the in-memory method to ~1e-12 relative. The per-cell sd is the
  only term not bit-identical: recovering it from raw moments sums in a
  different order than `stats::sd`.

  `p_value = TRUE` is refused with an explanation -- its permutation branch is
  `n_times` x cell-types extra marker-set means, thousands of passes rather
  than one heavier one. `rankEnrichParam` and `hyperEnrichParam` are refused
  too: rank ranks each gene across every cell, so no cell chunk can be scored
  in isolation, and the hypergeometric per-cell quantile shares nothing with
  the accumulators PAGE uses.
- A `expm1` op, the inverse of the existing `log` transform: `value ->
  base^value - 1`. Sparsity-preserving, so it lowers to Acero. Written as
  `base^value - 1` rather than `expm1(value * log(base))` because the
  in-memory backends compute the former and the two are not bit-identical in
  general; arrow's `^` agrees with R's elementwise, which is what makes the
  streamed PAGE reference level exact rather than merely close.

## bug fixes
- Streaming PAGE works on a store that already carries a `@post_ops` record.
  It pushed its `expm1` and `multiply` records with `phase = "lazy"`
  unconditionally, and per adr/0002 a lazy push is refused once `@post_ops` is
  non-empty, so any such store failed with *cannot queue a lazy op after a post
  op*. The records are now appended at the end of the chain, wherever that end
  is; both have executors for the R-side carrier, so either placement computes
  the same thing.
- `createGiottoXeniumObject(backend =)` no longer errors on Xenium-format
  directories that ship no panel json. Feature metadata is generated from the
  expression matrix when the panel is absent.

## changes
- A view step now narrows every cell-keyed subobject by the same
  `cell_ID` set, backed polygon stores included, so a recipe gives the
  same answer whichever slot you read it through. A crop's `geom` picks
  what represents a cell and `engine` picks what evaluates it. Backed
  polygons previously ignored `geom = "poly"`. See `adr/0015`.
- `storeRead(x, output = "duckdb")` on a `parquetExprStore` /
  `unionParquetExprStore` now rebuilds the scan from DuckDB's own
  `read_parquet` rather than registering an Arrow scanner. DuckDB owns the
  scan, so axis ranges prune parquet row groups in the engine and reads honour
  `giottodisk.duckdb_memory_limit`, which the Arrow bridge never saw. The
  subset predicates and the `@ops` chain are applied by the same functions the
  `"query"` path uses — they are dplyr, which lowers to Acero and DuckDB
  alike — so the two outputs return the same values by construction.
  - `storeRead(pe, output = "duckdb")` with no `duckdb_params$conn` now works;
    it previously errored. An ephemeral connection is created and kept alive by
    the returned `tbl_dbi`, matching the tabular stores.
  - `callback` is now applied on this path. A callback written against Arrow
    rather than plain dplyr will error here instead of being ignored.
- `duckdb_params$name` is now honoured by the tabular and geometry stores. It
  was documented but never read by `.pstore_to_duckdb`, which always generated
  its own view name.
- `storeRead(output = "duckdb")` errors when `conn` or `name` is passed
  directly rather than inside `duckdb_params`. Both previously landed in `...`,
  which no duckdb path reads, so the setting was dropped and the caller got a
  valid `tbl_dbi` on a connection they had not chosen. Applies to the tabular
  and geometry stores as well.
  - New option `giottodisk.duckdb_in_subquery_threshold` (default 1000):
    membership predicates larger than this are registered and joined rather
    than inlined by dbplyr as a literal `IN` list.
- `.op_transform_log` computes `log(value + 1)` rather than `log1p(value)`.
  DuckDB has no `log1p` and dbplyr does not translate it, so one expression now
  serves Acero, DuckDB and the data.table executor alike. Values are unchanged
  to within one ulp.
- New vignette, *Cell windows* (`vignette("expression_windows")`): the two options that steer the window, which passes window and when, what forces one, and why it has to be the cell axis. The package's first installed vignette, so `DESCRIPTION` gains `VignetteBuilder: knitr`. Decision recorded in adr/0011; `storeChunkInfo()` carries the options.
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
