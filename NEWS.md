# GiottoDisk 0.0.0.2

## changes
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
