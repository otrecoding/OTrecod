# transfo_target()

This function prepares the encoding of the target variable before
running an algorithm using optimal transportation theory.

## Usage

``` r
transfo_target(z, levels_order = NULL)
```

## Arguments

- z:

  a factor variable (ordered or not). A variable of another type will
  be, by default, convert to a factor.

- levels_order:

  a vector corresponding to the values of the levels of z. When the
  target is ordinal, the levels can be sorted by ascending order. By
  default, the initial order is remained.

## Value

The list returned is:

- NEW:

  an object of class factor of the same length as z

- LEVELS_NEW:

  the levels (ordered or not) retained for z

## Details

The function `transfo_target` is an intermediate function direcly
implemented in the functions
[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md)
and
[`OT_joint`](https://otrecoding.github.io/OTrecod/reference/OT_joint.md),
two functions dedicated to data fusion (see (1) and (2) for details).
Nevertheless, this function can also be used separately to assist user
in the conversion of a target variable (outcome) according to the
following rules:

- A character variable is converted in factor if the argument
  `levels_order` is set to NULL. In this case, the levels of the factor
  are assigned by order of appearance in the database.

- A character variable is converted in ordered factor if the argument
  `levels_order` differs from NULL. In this case, the levels of the
  factor correspond to those assigned in the argument.

- A factor stays unchanged if the argument `levels_order` is set to
  NULL. Otherwise the factor is converted in ordered factor and the
  levels are ordered according to the argument `levels_order`.

- A numeric variable, discrete or continuous is converted in factor if
  the argument `levels_order` is set to NULL, and the related levels are
  the values assigned in ascending order.

- A numeric variable, discrete or continuous is converted in ordered
  factor if the argument `levels_order` differed from NULL, and the
  related levels correspond to those assigned in the argument.

## References

1.  Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy
    N (2019). On the use of optimal transportation theory to recode
    variables and application to database merging. The International
    Journal of Biostatistics. Volume 16, Issue 1, 20180106, eISSN
    1557-4679. doi:10.1515/ijb-2018-0106

2.  Gares V, Omer J (2020) Regularized optimal transport of covariates
    and outcomes in data recoding. Journal of the American Statistical
    Association.
    [doi:10.1080/01621459.2020.1775615](https://doi.org/10.1080/01621459.2020.1775615)

## See also

[`compare_lists`](https://otrecoding.github.io/OTrecod/reference/compare_lists.md)

## Author

Gregory Guernec

<otrecod.pkg@gmail.com>

## Examples

``` r
y <- rnorm(100, 30, 10)
ynew1 <- transfo_target(y)
#> Your target y was numeric ... By default, it has been converted in factor of integers
#> 39 remaining levels

newlev <- unique(as.integer(y))
ynew2 <- transfo_target(y, levels_order = newlev)
#> Your target y was numeric ... By default, it has been converted in factor of integers
newlev2 <- newlev[-1]
ynew3 <- transfo_target(y, levels_order = newlev2)
#> Your target y was numeric ... By default, it has been converted in factor of integers
#> Inappropriate number or declared labels of levels
#> The default labels have been kept

outco <- c(rep("A", 25), rep("B", 50), rep("C", 25))
outco_new1 <- transfo_target(outco, levels_order = c("B", "C", "A"))
outco_new2 <- transfo_target(outco, levels_order = c("E", "C", "A", "F"))
#> Inappropriate number or declared labels of levels
#> The default levels have been kept
outco_new3 <- transfo_target(outco)

outco2 <- c(rep("A", 25), NA, rep("B", 50), rep("C", 25), NA, NA)
gg <- transfo_target(outco2)
hh <- transfo_target(outco2, levels_order = c("B", "C", "A"))
```
