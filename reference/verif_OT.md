# verif_OT()

This function proposes post-process verifications after data fusion by
optimal transportation algorithms.

## Usage

``` r
verif_OT(
  ot_out,
  group.class = FALSE,
  ordinal = TRUE,
  stab.prob = FALSE,
  min.neigb = 1
)
```

## Arguments

- ot_out:

  an otres object from
  [`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md)
  or
  [`OT_joint`](https://otrecoding.github.io/OTrecod/reference/OT_joint.md)

- group.class:

  a boolean indicating if the results related to the proximity between
  outcomes by grouping levels are requested in output (`FALSE` by
  default).

- ordinal:

  a boolean that indicates if \\Y\\ and \\Z\\ are ordinal (`TRUE` by
  default) or not. This argument is only useful in the context of groups
  of levels (`group.class`=TRUE).

- stab.prob:

  a boolean indicating if the results related to the stability of the
  algorithm are requested in output (`FALSE` by default).

- min.neigb:

  a value indicating the minimal required number of neighbors to
  consider in the estimation of stability (1 by default).

## Value

A list of 7 objects is returned:

- nb.profil:

  the number of profiles of covariates

- conf.mat:

  the global confusion matrix between \\Y\\ and \\Z\\

- res.prox:

  a summary table related to the association measures between \\Y\\ and
  \\Z\\

- res.grp:

  a summary table related to the study of the proximity of \\Y\\ and
  \\Z\\ using group of levels. Only if the `group.class` argument is set
  to TRUE.

- hell:

  Hellinger distances between observed and predicted distributions

- eff.neig:

  a table which corresponds to a count of conditional probabilities
  according to the number of neighbors used in their computation (only
  the first ten values)

- res.stab:

  a summary table related to the stability of the algorithm

## Details

In a context of data fusion, where information from a same target
population is summarized via two specific variables \\Y\\ and \\Z\\ (two
ordinal or nominal factors with different number of levels \\n_Y\\ and
\\n_Z\\), never jointly observed and respectively stored in two distinct
databases A and B, Optimal Transportation (OT) algorithms (see the
models `OUTCOME`, `R_OUTCOME`, `JOINT`, and `R_JOINT` of the reference
(2) for more details) propose methods for the recoding of \\Y\\ in B
and/or \\Z\\ in A. Outputs from the functions `OT_outcome` and
`OT_joint` so provides the related predictions to \\Y\\ in B and/or
\\Z\\ in A, and from these results, the function `verif_OT` provides a
set of tools (optional or not, depending on the choices done by user in
input) to estimate:

1.  the association between \\Y\\ and \\Z\\ after recoding

2.  the similarities between observed and predicted distributions

3.  the stability of the predictions proposed by the algorithm

A. PAIRWISE ASSOCIATION BETWEEN \\Y\\ AND \\Z\\

The first step uses standard criterions (Cramer's V, and Spearman's rank
correlation coefficient) to evaluate associations between two ordinal
variables in both databases or in only one database. When the argument
`group.class = TRUE`, these informations can be completed by those
provided by the function
[`error_group`](https://otrecoding.github.io/OTrecod/reference/error_group.md),
which is directly integrate in the function `verif_OT`. Assuming that
\\n_Y \> n_Z\\, and that one of the two scales of \\Y\\ or \\Z\\ is
unknown, this function gives additional informations about the potential
link between the levels of the unknown scale. The function proceeds to
this result in two steps. Firsty,
[`error_group`](https://otrecoding.github.io/OTrecod/reference/error_group.md)
groups combinations of modalities of \\Y\\ to build all possible
variables \\Y'\\ verifying \\n\_{Y'} = n_Z\\. Secondly, the function
studies the fluctuations in the association of \\Z\\ with each new
variable \\Y'\\ by using adapted comparisons criterions (see the
documentation of
[`error_group`](https://otrecoding.github.io/OTrecod/reference/error_group.md)
for more details). If grouping successive classes of \\Y\\ leads to an
improvement in the initial association between \\Y\\ and \\Z\\ then it
is possible to conclude in favor of an ordinal coding for \\Y\\ (rather
than nominal) but also to emphasize the consistency in the predictions
proposed by the algorithm of fusion.

B. SIMILARITIES BETWEEN OBSERVED AND PREDICTED DISTRIBUTIONS

When the predictions of \\Y\\ in B and/or \\Z\\ in A are available in
the `datab` argument, the similarities between the observed and
predicted probabilistic distributions of \\Y\\ and/or \\Z\\ are
quantified from the Hellinger distance (see (1)). This measure varies
between 0 and 1: a value of 0 corresponds to a perfect similarity while
a value close to 1 (the maximum) indicates a great dissimilarity. Using
this distance, two distributions will be considered as close as soon as
the observed measure will be less than 0.05.

C. STABILITY OF THE PREDICTIONS

These results are based on the decision rule which defines the stability
of an algorithm in A (or B) as its average ability to assign a same
prediction of \\Z\\ (or \\Y\\) to individuals that have a same given
profile of covariates \\X\\ and a same given level of \\Y\\ (or \\Z\\
respectively).

Assuming that the missing information of \\Z\\ in base A was predicted
from an OT algorithm (the reasoning will be identical with the
prediction of \\Y\\ in B, see (2) and (3) for more details), the
function `verif_OT` uses the conditional probabilities stored in the
object `estimatorZA` (see outputs of the functions
[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md)
and
[`OT_joint`](https://otrecoding.github.io/OTrecod/reference/OT_joint.md))
which contains the estimates of all the conditional probabilities of
\\Z\\ in A, given a profile of covariates \\x\\ and given a level of \\Y
= y\\. Indeed, each individual (or row) from A, is associated with a
conditional probability \\P(Z= z\|Y= y, X= x)\\ and averaging all the
corresponding estimates can provide an indicator of the predictions
stability.

The function
[`OT_joint`](https://otrecoding.github.io/OTrecod/reference/OT_joint.md)
provides the individual predictions for subject \\i\\:
\\\widehat{z}\_i\\, \\i=1,\ldots,n_A\\ according to the the maximum a
posteriori rule: \$\$\widehat{z}\_i= \mbox{argmax}\_{z\in \mathcal{Z}}
P(Z= z\| Y= y_i, X= x_i)\$\$ The function
[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md)
directly deduces the individual predictions from the probablities \\P(Z=
z\|Y= y, X= x)\\ computed in the second part of the algorithm (see (3)).

It is nevertheless common that conditional probabilities are estimated
from too rare covariates profiles to be considered as a reliable
estimate of the reality. In this context, the use of trimmed means and
standard deviances is suggested by removing the corresponding
probabilities from the final computation. In this way, the function
provides in output a table (`eff.neig` object) that provides the
frequency of these critical probabilities that must help the user to
choose. According to this table, a minimal number of profiles can be
imposed for a conditional probability to be part of the final
computation by filling in the `min.neigb` argument.

Notice that these results are optional and available only if the
argument `stab.prob = TRUE`. When the predictions of \\Z\\ in A and
\\Y\\ in B are available, the function `verif_OT` provides in output,
global results and results by database. The `res.stab` table can produce
NA with `OT_outcome` output in presence of incomplete shared variables:
this problem appears when the `prox.dist` argument is set to 0 and can
be simply solved by increasing this value.

## References

1.  Liese F, Miescke K-J. (2008). Statistical Decision Theory:
    Estimation, Testing, and Selection. Springer

2.  Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy
    N (2019). On the use of optimal transportation theory to recode
    variables and application to database merging. The International
    Journal of Biostatistics. Volume 16, Issue 1, 20180106, eISSN
    1557-4679. doi:10.1515/ijb-2018-0106

3.  Gares V, Omer J (2020) Regularized optimal transport of covariates
    and outcomes in data recoding. Journal of the American Statistical
    Association.
    [doi:10.1080/01621459.2020.1775615](https://doi.org/10.1080/01621459.2020.1775615)

## See also

[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md),
[`OT_joint`](https://otrecoding.github.io/OTrecod/reference/OT_joint.md),
[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md),
[`error_group`](https://otrecoding.github.io/OTrecod/reference/error_group.md)

## Author

Gregory Guernec

<otrecod.pkg@gmail.com>

## Examples

``` r

### Example 1
#-----
# - Using the data simu_data
# - Studying the proximity between Y and Z using standard criterions
# - When Y and Z are predicted in B and A respectively
# - Using an outcome model (individual assignment with knn)
#-----
data(simu_data)
outc1 <- OT_outcome(simu_data,
  quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
  dist.choice = "G", percent.knn = 0.90, maxrelax = 0,
  convert.num = 8, convert.class = 3,
  indiv.method = "sequential", which.DB = "BOTH", prox.dist = 0.30
)
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = OUTCOME
#> Distance                 = Gower
#> Percent closest knn      = 90%
#> Relaxation parameter     = NO
#> Relaxation value         = 0
#> Individual pred process  = Sequential
#> DB imputed               = BOTH
#> ---------------------------------------

verif_outc1 <- verif_OT(outc1)
verif_outc1
#> $nb.profil
#> [1] 232
#> 
#> $conf.mat
#>          predZ
#> predY       1   2   3   4   5 Sum
#>   [20-40]   2  83   2   2 217 306
#>   [40-60[ 253   2   2  40   0 297
#>   [60-80]   1   6  59  31   0  97
#>   Sum     256  91  63  73 217 700
#> 
#> $res.prox
#>          N V_cram rank_cor
#> Global 700   0.89   -0.547
#> 1st DB 300   0.88   -0.537
#> 2nd DB 400   0.89   -0.555
#> 
#> $res.grp
#> NULL
#> 
#> $hell
#>                 YA_YB ZA_ZB
#> Hellinger dist. 0.009 0.015
#> 
#> $eff.neig
#> NULL
#> 
#> $res.stab
#> NULL
#> 

# \donttest{

### Example 2
#-----
# - Using the data simu_data
# - Studying the proximity between Y and Z using standard criterions and studying
#   associations by grouping levels of Z
# - When only Y is predicted in B
# - Tolerated distance between a subject and a profile: 0.30 * distance max
# - Using an outcome model (individual assignment with knn)
#-----

data(simu_data)
outc2 <- OT_outcome(simu_data,
  quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
  dist.choice = "G", percent.knn = 0.90, maxrelax = 0, prox.dist = 0.3,
  convert.num = 8, convert.class = 3,
  indiv.method = "sequential", which.DB = "B"
)
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = OUTCOME
#> Distance                 = Gower
#> Percent closest knn      = 90%
#> Relaxation parameter     = NO
#> Relaxation value         = 0
#> Individual pred process  = Sequential
#> DB imputed               = B
#> ---------------------------------------

verif_outc2 <- verif_OT(outc2, group.class = TRUE, ordinal = TRUE)
verif_outc2
#> $nb.profil
#> [1] 232
#> 
#> $conf.mat
#>          Z
#> predY       1   2   3   4   5 Sum
#>   [20-40]   1  47   1   1 126 176
#>   [40-60[ 144   1   1  24   0 170
#>   [60-80]   0   3  33  18   0  54
#>   Sum     145  51  35  43 126 400
#> 
#> $res.prox
#>        N   V_cram rank_cor 
#>  400.000    0.890   -0.555 
#> 
#> $res.grp
#>   grp levels Z to Y error_rate  Kappa Vcramer RankCor
#> 1         1 2 3/4/5       81.8 -0.243    0.57  -0.586
#> 4         1 2/3 4/5       81.8 -0.241    0.73  -0.443
#> 5         1 2/3/4 5       83.2 -0.206    0.66  -0.330
#> 3         1/2/3 4 5       86.8 -0.209    0.63  -0.308
#> 2         1/2 3 4/5       93.2 -0.411    0.75  -0.650
#> 6         1/2 3/4 5       94.8 -0.369    0.66  -0.492
#> 
#> $hell
#>                 YA_YB ZA_ZB
#> Hellinger dist. 0.009    NA
#> 
#> $eff.neig
#> NULL
#> 
#> $res.stab
#> NULL
#> 


### Example 3
#-----
# - Using the data simu_data
# - Studying the proximity between Y and Z using standard criterions and studying
#   associations by grouping levels of Z
# - Studying the stability of the conditional probabilities
# - When Y and Z are predicted in B and A respectively
# - Using an outcome model (individual assignment with knn)
#-----

verif_outc2b <- verif_OT(outc2, group.class = TRUE, ordinal = TRUE, stab.prob = TRUE, min.neigb = 5)
verif_outc2b
#> $nb.profil
#> [1] 232
#> 
#> $conf.mat
#>          Z
#> predY       1   2   3   4   5 Sum
#>   [20-40]   1  47   1   1 126 176
#>   [40-60[ 144   1   1  24   0 170
#>   [60-80]   0   3  33  18   0  54
#>   Sum     145  51  35  43 126 400
#> 
#> $res.prox
#>        N   V_cram rank_cor 
#>  400.000    0.890   -0.555 
#> 
#> $res.grp
#>   grp levels Z to Y error_rate  Kappa Vcramer RankCor
#> 1         1 2 3/4/5       81.8 -0.243    0.57  -0.586
#> 4         1 2/3 4/5       81.8 -0.241    0.73  -0.443
#> 5         1 2/3/4 5       83.2 -0.206    0.66  -0.330
#> 3         1/2/3 4 5       86.8 -0.209    0.63  -0.308
#> 2         1/2 3 4/5       93.2 -0.411    0.75  -0.650
#> 6         1/2 3/4 5       94.8 -0.369    0.66  -0.492
#> 
#> $hell
#>                 YA_YB ZA_ZB
#> Hellinger dist. 0.009    NA
#> 
#> $eff.neig
#>    Nb.neighbor Nb.Prob
#> 1            3       5
#> 2            4       9
#> 3            5      12
#> 4            6      14
#> 5            7      11
#> 6            8      23
#> 7            9      14
#> 8           10      10
#> 9           11      15
#> 10          12      12
#> 
#> $res.stab
#>          N min.N mean    sd
#> 2nd DB 386     5 0.95 0.153
#> 
# }
```
