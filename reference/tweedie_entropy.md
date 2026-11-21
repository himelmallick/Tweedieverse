# Entropy of a Vector

Compute Shannon entropy for a vector.

## Usage

``` r
tweedie_entropy(target)
```

## Arguments

- target:

  A vector of values.

## Value

A numeric entropy value.

## Details

Missing values are excluded from category counts but included in the
denominator to preserve the package's metadata-filtering convention.
Empty and entirely missing vectors return zero.

## References

Shannon CE (1948). A Mathematical Theory of Communication. Bell System
Technical Journal, 27, 379-423.
[doi:10.1002/j.1538-7305.1948.tb01338.x](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x)
.

## Examples

``` r
tweedie_entropy(c("A", "A", "B", "B", "C"))
#> [1] 1.521928
```
