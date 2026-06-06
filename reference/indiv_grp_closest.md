# indiv_grp_closest()

This function sequentially assigns individual predictions using a
nearest neighbors procedure to solve recoding problems of data fusion.

## Usage

``` r
indiv_grp_closest(
  proxim,
  jointprobaA = NULL,
  jointprobaB = NULL,
  percent_closest = 1,
  which.DB = "BOTH"
)
```

## Arguments

- proxim:

  a
  [`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md)
  object or an object of similar structure

- jointprobaA:

  a matrix whose number of columns corresponds to the number of
  modalities of the target variable \\Y\\ in database A, and which
  number of rows corresponds to the number of modalities of Z in
  database B. It gives an estimation of the joint probability of
  \\(Y,Z)\\ in A. The sum of cells of this matrix must be equal to 1

- jointprobaB:

  a matrix whose number of columns equals to the number of modalities of
  the target variable \\Y\\ in database A, and which number of rows
  corresponds to the number of modalities of \\Z\\ in database B. It
  gives an estimation of the joint probability of \\(Y,Z)\\ in B. The
  sum of cells of this matrix must be equal to 1

- percent_closest:

  a value between 0 and 1 (by default) corresponding to the fixed
  `percent closest` of individuals remained in the computation of the
  average distances

- which.DB:

  a character string (with quotes) that indicates which individual
  predictions need to be computed: only the individual predictions of
  \\Y\\ in B ("B"), only those of \\Z\\ in A ("A") or the both ("BOTH"
  by default)

## Value

A list of two vectors of numeric values:

- YAtrans:

  a vector corresponding to the individual predictions of \\Y\\ (numeric
  form) in the database B using the Optimal Transportation algorithm

- ZBtrans:

  a vector corresponding to the individual predictions of \\Z\\ (numeric
  form) in the database A using the Optimal Transportation algorithm

## Details

A. THE RECODING PROBLEM IN DATA FUSION

Assuming that \\Y\\ and \\Z\\ are two variables which refered to the
same target population in two separate databases A and B respectively
(no overlapping rows), so that \\Y\\ and \\Z\\ are never jointly
observed. Assuming also that A and B share a subset of common covariates
\\X\\ of any types (same encodings in A and B) completed or not.
Integrating these two databases often requires to solve the recoding
problem by creating an unique database where the missing information of
\\Y\\ and \\Z\\ is fully completed.

B. DESCRIPTION OF THE FUNCTION

The function `indiv_grp_closest` is an intermediate function used in the
implementation of an algorithm called OUTCOME (and its enrichment
R-OUTCOME, see the reference (2) for more details) dedicated to the
solving of recoding problems in data fusion using Optimal Transportation
theory. The model is implemented in the function
[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md)
which integrates the function `indiv_grp_closest` in its syntax as a
possible second step of the algorithm. The function `indiv_grp_closest`
can also be used separately provided that the argument `proxim` receives
an output object of the function
[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md).
This latter is available in the package and is so directly usable
beforehand.

The algorithms `OUTCOME` (and `R-OUTCOME`) are made of two independent
parts. Assuming that the objective consists in the prediction of \\Z\\
in the database A:

- The first part of the algorithm solves the optimization problem by
  providing a solution called \\\gamma\\ that corresponds here to an
  estimation of the joint distribution \\(Y,Z)\\ in A.

- From the first part, a nearest neighbor procedure is carried out as a
  second part to provide the individual predictions of \\Z\\ in A: this
  procedure is implemented in the function `indiv_group_closest`. In
  other words, this function sequentially assigns to each individual of
  A the modality of \\Z\\ that is closest.

Obviously, this algorithm runs in the same way for the prediction of
\\Y\\ in the database B. The function `indiv_grp_closest` integrates in
its syntax the function
[`avg_dist_closest`](https://otrecoding.github.io/OTrecod/reference/avg_dist_closest.md).
Therefore, the related argument `percent_closest` is identical in the
two functions. Thus, when computing average distances between an
individual \\i\\ and a subset of individuals assigned to a same level of
\\Y\\ or \\Z\\ is required, user can decide if all individuals from the
subset of interest can participate to the computation
(`percent_closest`=1) or only a fixed part p (\<1) corresponding to the
closest neighbors of \\i\\ (in this case `percent_closest` = p).

The arguments `jointprobaA` and `jointprobaB` correspond to the
estimations of \\\gamma\\ (sum of cells must be equal to 1) in A and/or
B respectively, according to the `which.DB` argument. For example,
assuming that \\n\_{Y_1}\\ individuals are assigned to the first
modality of \\Y\\ in A, the objective consists in the individual
predictions of \\Z\\ in A. Then, if `jointprobaA`\[1,2\] = 0.10, the
maximum number of individuals that can be assigned to the second
modality of \\Z\\ in A, can not exceed \\0.10 \times n_A\\. If
\\n\_{Y_1} \leq 0.10 \times n_A\\ then all individuals assigned to the
first modality of \\Y\\ will be assigned to the second modality of
\\Z\\. At the end of the process, each individual with still no
affectation will receive the same modality of \\Z\\ as those of his
nearest neighbor in B.

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

[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md),[`avg_dist_closest`](https://otrecoding.github.io/OTrecod/reference/avg_dist_closest.md),
,[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md)

## Author

Gregory Guernec, Valerie Gares, Jeremy Omer

<otrecod.pkg@gmail.com>

## Examples

``` r
data(simu_data)

### Example with the Manhattan distance

man1 <- transfo_dist(simu_data,
  quanti = c(3, 8), nominal = c(1, 4:5, 7),
  ordinal = c(2, 6), logic = NULL, prep_choice = "M"
)
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 5 remaining levels
mat_man1 <- proxim_dist(man1, norm = "M")

### Y(Yb1) and Z(Yb2) are a same information encoded in 2 different forms:
### (3 levels for Y and 5 levels for Z)
### ... Stored in two distinct databases, A and B, respectively
### The marginal distribution of Y in B is unknown,
### as the marginal distribution of Z in A ...

# Empirical distribution of Y in database A:
freqY <- prop.table(table(man1$Y))
freqY
#> 
#>   [20-40]   [40-60[   [60-80] 
#> 0.4333333 0.4233333 0.1433333 

# Empirical distribution of Z in database B
freqZ <- prop.table(table(man1$Z))
freqZ
#> 
#>      1      2      3      4      5 
#> 0.3625 0.1275 0.0875 0.1075 0.3150 

# By supposing that the following matrix called transport symbolizes
# an estimation of the joint distribution L(Y,Z) ...
# Note that, in reality this distribution is UNKNOWN and is
# estimated in the OT function by resolving an optimisation problem.


transport1 <- matrix(c(0.3625, 0, 0, 0.07083333, 0.05666667,
                      0, 0, 0.0875, 0, 0, 0.1075, 0,
                      0, 0.17166667, 0.1433333),
                     ncol = 5, byrow = FALSE)

# ... So that the marginal distributions of this object corresponds to freqY and freqZ:
apply(transport1, 1, sum) # = freqY
#> [1] 0.4333333 0.4233333 0.1433333
apply(transport1, 2, sum) # = freqZ
#> [1] 0.3625 0.1275 0.0875 0.1075 0.3150

# The affectation of the predicted values of Y in database B and Z in database A
# are stored in the following object:

pred_man1 <- indiv_grp_closest(mat_man1,
  jointprobaA = transport1, jointprobaB = transport1,
  percent_closest = 0.90
)
summary(pred_man1)
#>         Length Class  Mode   
#> YAtrans 400    -none- numeric
#> ZBtrans 300    -none- numeric

# For the prediction of Z in A only, add the corresponding argument:
pred_man1_A <- indiv_grp_closest(mat_man1,
  jointprobaA = transport1, jointprobaB = transport1,
  percent_closest = 0.90, which.DB = "A"
)
```
