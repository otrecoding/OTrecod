# indiv_grp_optimal()

This function assigns individual predictions to the incomplete
information of two integrated datasources by solving a linear
optimization problem.

## Usage

``` r
indiv_grp_optimal(
  proxim,
  jointprobaA,
  jointprobaB,
  percent_closest = 1,
  solvr = "glpk",
  which.DB = "BOTH"
)
```

## Arguments

- proxim:

  a
  [`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md)
  object or an object of similar structure

- jointprobaA:

  a matrix whose number of columns is equal to the number of modalities
  of the target variable \\Y\\ in database A, and whose number of rows
  is equal to the number of modalities of \\Z\\ in database B. It gives
  an estimation of the joint probability \\(Y,Z)\\ in the database A.
  The sum of cells of this matrix must be equal to 1.

- jointprobaB:

  a matrix whose number of columns is equal to the number of modalities
  of the target variable Y in database A, and whose number of rows is
  equal to the number of modalities of \\Z\\ in database B. It gives an
  estimation of the joint probability \\(Y,Z)\\ in the database B. The
  sum of cells of this matrix must be equal to 1.

- percent_closest:

  a value between 0 and 1 (by default) corresponding to the fixed
  `percent closest` of individuals used in the computation of the
  average distances

- solvr:

  a character string that specifies the type of method selected to solve
  the optimization algorithms. The default solver is "glpk".

- which.DB:

  a character string that indicates which individual predictions are
  computed: only the individual predictions of \\Y\\ in B ("B"), only
  those of \\Z\\ in A ("A") or the both ("BOTH" by default).

## Value

A list of two vectors of numeric values:

- YAtrans:

  a vector corresponding to the predicted values of \\Y\\ in database B
  (numeric form) according to the `which.DB` argument

- ZBtrans:

  a vector corresponding to the predicted values of \\Z\\ in database A
  (numeric form) according to the `which.DB` argument

## Details

A. THE RECODING PROBLEM IN DATA FUSION

Assuming that \\Y\\ and \\Z\\ are two target variables which refered to
the same target population in two separate databases A and B
respectively (no overlapping rows), so that \\Y\\ and \\Z\\ are never
jointly observed. Assuming also that A and B share a subset of common
covariates \\X\\ of any types (same encodings in A and B) completed or
not. Merging these two databases often requires to solve a recoding
problem by creating an unique database where the missing information of
\\Y\\ and \\Z\\ is fully completed.

B. DESCRIPTION OF THE FUNCTION

The function `indiv_grp_optimal` is an intermediate function used in the
implementation of an algorithm called `OUTCOME` (and its enrichment
`R-OUTCOME` (2)) dedicated to the solving of recoding problems in data
fusion using Optimal Transportation theory. The model is implemented in
the function
[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md)
which integrates the function `indiv_grp_optimal` in its syntax as a
possible second step of the algorithm. The function `indiv_grp_optimal`
can nevertheless be used separately providing that the argument `proxim`
receives an output object of the function
[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md).
This latter is available in the package and is so directly usable
beforehand.

The function `indiv_grp_optimal` constitutes an alternative method to
the nearest neighbor procedure implemented in the function
[`indiv_grp_closest`](https://otrecoding.github.io/OTrecod/reference/indiv_grp_closest.md).
As for the function
[`indiv_grp_closest`](https://otrecoding.github.io/OTrecod/reference/indiv_grp_closest.md),
assuming that the objective consists in the prediction of \\Z\\ in the
database A, the first step of the algorithm related to `OUTCOME`
provides an estimate of \\\gamma\\, the solution of the optimization
problem, which can be seen, in this case as an estimation of the joint
distribution \\(Y,Z)\\ in A. Rather than using a nearest neighbor
approach to provide individual predictions, the function
`indiv_grp_optimal` solves an optimization problem using the simplex
algorithm which searches for the individual predictions of \\Z\\ that
minimize the computed total distance satisfying the joint probability
distribution estimated in the first part. More details about the theory
related to the solving of this optimization problem is described in the
section 5.3 of (2).

Obviously, this algorithm runs in the same way for the prediction of
\\Y\\ in the database B. The function `indiv_grp_optimal` integrates in
its syntax the function
[`avg_dist_closest`](https://otrecoding.github.io/OTrecod/reference/avg_dist_closest.md)
and the related argument `percent_closest` is identical in the two
functions. Thus, when computing average distances between an individual
i and a subset of individuals assigned to a same level of \\Y\\ or \\Z\\
is required, user can decide if all individuals from the subset of
interest can participate to the computation (`percent_closest = 1`) or
only a fixed part p (\<1) corresponding to the closest neighbors of i
(in this case `percent_closest` = p).

The arguments `jointprobaA` and `jointprobaB` can be seen as estimations
of \\\gamma\\ (sum of cells must be equal to 1) that correspond to
estimations of the joint distributions of \\(Y,Z)\\ in A and B
respectively.

The argument `solvr` permits user to choose the solver of the
optimization algorithm. The default solver is "glpk" that corresponds to
the GNU Linear Programming Kit (see (3) for more details). The solver
"clp" (see (4)) for Coin-or Linear Programming, convenient in linear and
quadratic situations, is also directly integrated in the function.
Moreover, the function actually uses the `R` optimization infrastructure
of the package ROI which offers a wide choice of solver to users by
easily loading the associated plugins of ROI (see (5)).

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

3.  Makhorin A (2011). GNU Linear Programming Kit Reference Manual
    Version 4.47.<http://www.gnu.org/software/glpk/>

4.  Forrest J, de la Nuez D, Lougee-Heimer R (2004). Clp User Guide.
    <https://www.coin-or.org/Clp/userguide/index.html>

5.  Theussl S, Schwendinger F, Hornik K (2020). ROI: An Extensible R
    Optimization Infrastructure.Journal of Statistical Software,94(15),
    1-64.
    [doi:10.18637/jss.v094.i15](https://doi.org/10.18637/jss.v094.i15)

## See also

[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md),
[`avg_dist_closest`](https://otrecoding.github.io/OTrecod/reference/avg_dist_closest.md),
[`indiv_grp_closest`](https://otrecoding.github.io/OTrecod/reference/indiv_grp_closest.md)

## Author

Gregory Guernec, Valerie Gares, Jeremy Omer

<otrecod.pkg@gmail.com>

## Examples

``` r

### Example using The Euclidean distance on a complete database
# For this example we keep only 200 rows:

data(tab_test)
tab_test2 <- tab_test[c(1:80, 5001:5080), ]
dim(tab_test2)
#> [1] 160   6

# Adding NAs in Y1 and Y2
tab_test2[tab_test2$ident == 2, 2] <- NA
tab_test2[tab_test2$ident == 1, 3] <- NA

# Because all covariates are ordered in numeric form,
# the transfo_dist function is not required here

mat_testm <- proxim_dist(tab_test2, norm = "M")

### Y(Y1) and Z(Y2) are a same variable encoded in 2 different forms:
### 4 levels for Y1 and 3 levels for Y2
### ... Stored in two distinct databases, A and B, respectively
### The marginal distribution of Y in B is unknown,
### as the marginal distribution of Z in A ...

# Assuming that the following matrix called transport symbolizes
# an estimation of the joint distribution L(Y,Z) ...
# Note that, in reality this distribution is UNKNOWN and is
# estimated in the OT function by resolving the optimization problem.

# By supposing:

val_trans <- c(0.275, 0.115, 0, 0, 0, 0.085, 0.165, 0, 0, 0, 0.095, 0.265)
mat_trans <- matrix(val_trans, ncol = 3, byrow = FALSE)

# Getting the individual predictions of Z in A (only)
# by computing average distances on 90% of the nearest neighbors of
# each modality of Z in B
predopt_A <- indiv_grp_optimal(mat_testm,
  jointprobaA = mat_trans,
  jointprobaB = mat_trans, percent_closest = 0.90,
  which.DB = "A"
)

# \donttest{
### Example 2 using The Manhattan distance with incomplete covariates
data(simu_data)

man1 <- transfo_dist(simu_data,
  quanti = c(3, 8), nominal = c(1, 4:5, 7),
  ordinal = c(2, 6), logic = NULL, prep_choice = "M"
)
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 5 remaining levels
mat_man1 <- proxim_dist(man1, norm = "M")


### Y and Z are a same variable encoded in 2 different forms:
### (3 levels for Y and 5 levels for Z)
### ... Stored in two distinct databases, A and B, respectively
### The marginal distribution of Y in B is unknown,
### as the marginal distribution of Z in A ...


# By supposing that the following matrix called transport symbolizes
# an estimation of the joint distribution L(Y,Z) ...
# Note that, in reality this distribution is UNKNOWN and is
# estimated in the OT function by resolving an optimisation problem.

mat_trans2 <- matrix(c(0.3625, 0, 0, 0.07083333, 0.05666667,
                       0, 0, 0.0875, 0, 0, 0.1075, 0,
                       0, 0.17166667, 0.1433333),
                     ncol = 5, byrow = FALSE)


# The predicted values of Y in database B and Z in
# database A are stored in the following object:

predopt2 <- indiv_grp_optimal(mat_man1,
  jointprobaA = mat_trans2,
  jointprobaB = mat_trans2,
  percent_closest = 0.90
)
summary(predopt2)
#>         Length Class  Mode   
#> YAtrans 400    -none- numeric
#> ZBtrans 300    -none- numeric
# }
```
