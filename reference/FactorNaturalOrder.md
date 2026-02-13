# Create a factor with natural (human-friendly) ordering

Converts a character or numeric vector into a factor where the levels
are ordered naturally (e.g., `a1`, `a2`, ..., `a10` instead of
lexicographically as `a1`, `a10`, `a2`, ...). This is useful for
plotting or labeling grouped data where numeric substrings should follow
numeric order.

## Usage

``` r
factorNaturalOrder(x)
```

## Arguments

- x:

  A character or numeric vector to convert to a factor with natural
  order.

## Value

A factor with levels sorted in natural (human-readable) order.

## Examples

``` r
# Numeric vector
factorNaturalOrder(10:1)
#>  [1] 10 9  8  7  6  5  4  3  2  1 
#> Levels: 1 2 3 4 5 6 7 8 9 10

# Character vector with embedded numbers
factorNaturalOrder(c("a11", "a12", "a1", "a2", "a"))
#> [1] a11 a12 a1  a2  a  
#> Levels: a a1 a2 a11 a12
```
