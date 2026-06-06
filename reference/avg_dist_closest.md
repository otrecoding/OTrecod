# avg_dist_closest()

This function computes average distances between levels of two
categorical variables located in two distinct databases.

## Usage

``` r
avg_dist_closest(proxim, percent_closest = 1)
```

## Arguments

- proxim:

  a `proxim_dist` object

- percent_closest:

  a ratio between 0 and 1 corresponding to the desired part of rows (or
  statistical units, or individuals) that will participate to the
  computation of the average distances between levels of factors or
  between an individual (a row) and levels of only one factor. Indeed,
  target variables are factors and each level of factor is characterized
  by a subset of rows, themselves characterized by their covariate
  profiles. These rows can be ordered according to their distances at
  their factor level. When this ratio is set to 1 (default setting), all
  rows participate to the computation, nevertheless when this ratio is
  less than 1, only rows with the smallest factor level distances will
  be kept for the computation (see 'Details').

## Value

A list of 3 matrices is returned:

- Davg:

  the cost matrix whose number of rows corresponds to \\n_Y\\, the
  number of levels of the target variable \\Y\\ in the database A, and
  whose number of columns corresponds to \\n_Z\\: the number of levels
  of the target variable in B. In this case, the related cost matrix can
  be interpreted as the ability to move from one level of \\Y\\ in A to
  one level of \\Z\\ in B. Davg\[P,Q\] refers to the average distance
  between the modality P of \\Y\\ (only known in A) and modality \\Q\\
  of \\Z\\ (only known in B).

- DindivA:

  a matrix whose number of rows corresponds to the number of rows of the
  first database A and number of columns corresponds to \\n_Z\\, the
  number of levels of the target variable \\Z\\ in the second
  database B. DindivA\[i,Q\] refers to the average distance between the
  \\i^{th}\\ individual (or row) of the first database and a chosen
  proportion of individuals (`percent_closest` set by the user) of the
  second database having the modality \\Q\\ of \\Z\\.

- DindivB:

  a matrix whose number of rows corresponds to the number of rows of the
  second database B and number of columns corresponds to nA, the number
  of levels of the target variable in the first database A.
  DindivB\[k,P\] refers to the average distance between the \\k^{th}\\
  individual (or row) of the second database and a chosen proportion of
  individuals (depending on `percent_closest`) of the first database
  having the modality P of \\Y\\.

## Details

The function `avg_dist_closest` is an intermediate function for the
implementation of original algorithms dedicated to the solving of
recoding problems in data fusion using Optimal Transportation theory
(for more details, consult the corresponding algorithms called
`OUTCOME`, `R_OUTCOME`, `JOINT` and `R_JOINT`, in the reference (2)).
The function `avg_dist_closest` is so directly implemented in the
`OT_outcome` and `OT_joint` functions but can also be used separately.
The function `avg_dist_closest` uses, in particular, the distance matrix
D (that stores distances between rows of A and B) from the function
[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md)
to produce three distinct matrices saved in a list object. Therefore,
the function requires in input, the specific output of the function
[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md)
which is available in the package and so must be used beforehand. In
consequence, do not use this function directly on your database, and do
not hesitate to consult the provided examples provided for a better
understanding.

DEFINITION OF THE COST MATRIX

Assuming that A and B are two databases with a set of shared variables
and that a same information (referred to a same target population) is
stored as a variable \\Y\\ in A and \\Z\\ in B, such that \\Y\\ is
unknown in B and \\Z\\ is unknown in A, whose encoding depends on the
database (\\n_Y\\ levels in A and \\n_Z\\ levels in B). A distance
between one given level y of \\Y\\ and one given level z of \\Z\\ is
estimated by averaging the distances between the two subsets of
individuals (units or rows) assigned to y in A and z in B, characterized
by their vectors of covariates. The distance between two individuals
depends on the variations between the shared covariates, and so depends
on the chosen distance function using the function `proxim_dist`. For
these computations, all the individuals concerned by these two levels
can be taken into account, or only a part of them, depending on the
argument `percent_closest`. When `percent_closest` \< 1, the average
distance between an individual \\i\\ and a given level of factor z only
uses the corresponding part of individuals related to z that are the
closest to \\i\\. Therefore, this choice influences the estimations of
average distances between levels of factors but also permits to reduce
time computation when necessary.

The average distance between each individual of \\Y\\ (resp. \\Z\\) and
each levels of \\Z\\ (resp. \\Y\\) are returned in output, in the object
`DindivA` (`DindivB` respectively). The average distance between each
levels of \\Y\\ and each levels of \\Z\\ are returned in a matrix saved
in output (the object `Davg`). `Davg` returns the computation of the
cost matrix D, whose dimensions (\\n_Y \times n_Z\\) correspond to the
number of levels of \\Y\\ (rows) and \\Z\\ (columns). This matrix can be
seen as the ability for an individual (row) to move from a given level
of the target variable (\\Y\\) in A to a given level of \\Z\\ in the
database B (or vice versa).

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

[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md)

## Author

Gregory Guernec, Valerie Gares, Jeremy Omer

<otrecod.pkg@gmail.com>

## Examples

``` r
data(simu_data)
### The covariates of the data are prepared according to the distance chosen
### using the transfo_dist function

### Example with The Manhattan distance

man1 <- transfo_dist(simu_data,
  quanti = c(3, 8), nominal = c(1, 4:5, 7),
  ordinal = c(2, 6), logic = NULL, prep_choice = "M"
)
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 5 remaining levels
mat_man1 <- proxim_dist(man1, norm = "M")

# proxim_dist() fixes the chosen distance function,
# and defines neighborhoods between profiles and individuals

# The following row uses only 80 percents of individuals of each level
# of factors for the computation of the average distances:

neig_man1 <- avg_dist_closest(mat_man1, percent_closest = 0.80)
```
