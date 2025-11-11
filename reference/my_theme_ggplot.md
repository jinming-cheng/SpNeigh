# Custom ggplot2 theme

Provides a consistent and clean visual style for plots generated within
this package. This theme builds on `theme_classic()` and adjusts text
sizes and margins for better readability in figures.

## Usage

``` r
my_theme_ggplot()
```

## Value

A `ggplot2` theme object that can be added to ggplot visualizations.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(mpg, wt)) +
    geom_point() +
    my_theme_ggplot()
```
