# File Store

S4 class for disk-backed data storage. Defines a file system `path` and
a `read_fun` method for accessing the data. Create with
[`fileStore()`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore.md).

`fileStore` serves dual purposes:

- base class for disk-backed stores

- wild card customizable class allowing ad-hoc compatibility with
  [`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  for non-standard file formats.

## Slots

- `path`:

  `character`. File system path (local or remote URI)

- `uid`:

  `character` automatically generated unique ID for artifact tracking.
  See
  [artifact_uid](https://giotto-suite.github.io/GiottoDisk/reference/artifact_uid.md)

- `params`:

  `list` Additional params such as remote names that are relevant to the
  store type may be stored here as a named list.

- `read_fun`:

  `function`. Function to access data where `read_fun(path)` returns the
  data in a useful format. Can be customized for compression or other
  access requirements. Do not embed credentials in `read_fun()`

## Extending Store Classes

Specific stores with particular **file formats** (parquet, HDF5, etc.)
extend this class. These specialized stores provide preset `read_fun`
implementations in their `initialize()` methods, but the ability to
customize remains available for edge cases.

When extending, register a lightweight `read_fun` to *access* the data,
while specific optimization on reading/filtering should be in the
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
method.

## Customizing Access

The `read_fun` is primarily useful for custom file formats:

    # Wildcard usage - custom formats
    rds_store <- fileStore(path = "data.rds", read_fun = readRDS)
    fst_store <- fileStore(path = "data.fst", read_fun = fst::read_fst)

Specialized stores (parquetStore, h5ArrayStore, etc.) have preset
readers but can be overridden for edge cases if needed.

## `queryableStore`

Subclass of `fileStore` that is accessible via dplyr query semantics.
Flag that the `fileStore` is queryable by coercing to this class

    as(x, "queryableStore")

## See also

Other store types:
[`binGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput-class.md),
[`cellbinGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput-class.md),
[`csvWideInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput-class.md),
[`dataStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/dataStore-class.md),
[`edgeInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/edgeInput-class.md),
[`exprInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/exprInput-class.md),
[`mtxInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput-class.md),
[`parquetEdgeStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore-class.md),
[`parquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md),
[`parquetStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md),
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
