# Generate a safe color palette for discrete clusters

Returns a vector of visually distinct colors for plotting discrete
clusters. Uses `colors15_cheng` as the base palette. If the number of
clusters exceeds the base palette, additional colors are generated using
[`scales::hue_pal()`](https://scales.r-lib.org/reference/pal_hue.html).

## Usage

``` r
safeColorPalette(
  n_clusters = NULL,
  base_colors = colors15_cheng,
  verbose = TRUE
)
```

## Arguments

- n_clusters:

  Number of unique clusters.

- base_colors:

  A character vector of base colors. Default is `colors15_cheng`.

- verbose:

  Logical. If `TRUE`, displays a message when extended colors are used.
  Default is `TRUE`.

## Value

A character vector of colors of length `n_clusters`.

## Examples

``` r
safeColorPalette(5)
#> [1] "#6495ED" "#BF3EFF" "#FF3030" "#FFD700" "#ADFF2F"
safeColorPalette(20)
#> Note: 15 base colors provided. Generated 5 additional colors using `scales::hue_pal()` (total = 20) to match 20 clusters.
#>  [1] "#6495ED" "#BF3EFF" "#FF3030" "#FFD700" "#ADFF2F" "#00FA9A" "#48D1CC"
#>  [8] "#FFA500" "#FFC0CB" "#CD1076" "#EE82EE" "#FF00FF" "#8B6914" "#00FFFF"
#> [15] "#E5E5E5" "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
```
