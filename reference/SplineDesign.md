# Generate an orthonormal spline-based design matrix

Constructs an orthonormal design matrix from a numeric covariate (e.g.,
spatial distance) using a natural cubic spline basis. The output matrix
can be used in linear modeling to capture smooth, non-linear trends
along continuous variables.

## Usage

``` r
splineDesign(x, df = 3)
```

## Arguments

- x:

  A numeric vector representing a continuous covariate (e.g., distance
  or pseudotime).

- df:

  Integer. Degrees of freedom for the spline basis. Default is 3.

## Value

A numeric matrix with orthonormal columns (same number of rows as `x`).
The columns represent smoothed trends extracted from the spline basis.
The first column is directionally aligned with the input vector (`x`).

## Details

The first column of the resulting matrix is aligned to show a positive
correlation with the input vector and typically captures the main linear
or monotonic trend.

## Examples

``` r
x <- seq(0, 1, length.out = 100)
Z <- splineDesign(x)
cor(Z[, 1], x) # Should be > 0
#> [1] 1

# Use Z in modeling
y <- sin(2 * pi * x) + rnorm(100, sd = 0.2)
fit <- lm(y ~ Z)
summary(fit)
#> 
#> Call:
#> lm(formula = y ~ Z)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -0.47727 -0.12913 -0.01249  0.10844  0.46247 
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -0.01236    0.01838  -0.672    0.503    
#> Z1          -5.63691    0.18383 -30.664   <2e-16 ***
#> Z2           3.43634    0.18383  18.693   <2e-16 ***
#> Z3           2.98165    0.18383  16.220   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 0.1838 on 96 degrees of freedom
#> Multiple R-squared:  0.9418, Adjusted R-squared:   0.94 
#> F-statistic: 517.6 on 3 and 96 DF,  p-value: < 2.2e-16
#> 
```
