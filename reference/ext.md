# Spatial Extent

Spatially mapped bounds. `exact = TRUE` (default): scans coordinates
with all spatial filter ops applied; pending transform is projected
during the scan. `exact = FALSE`: fast estimate from metadata bounds
(`@crop` \> `disk_extent`) intersected with `@window` and projected
through any pending transform. Axis-aligned transforms are exact;
rotation/shear gives a conservative AABB. Row-level ops (`subset`,
`head`, etc.) are never reflected in either mode.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
ext(x, exact = TRUE, ...)
```

## Arguments

- x:

  object to use

- exact:

  `logical(1)`. If `TRUE` (default), scans for a true extent. If
  `FALSE`, returns a fast metadata-based estimate without scanning.

- ...:

  additional params to pass (not used)
