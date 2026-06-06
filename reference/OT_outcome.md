# OT_outcome()

The function `OT_outcome` integrates two algorithms called (`OUTCOME`)
and (`R-OUTCOME`) dedicated to the solving of recoding problems in data
fusion using optimal transportation (OT) of the joint distribution of
outcomes.

## Usage

``` r
OT_outcome(
  datab,
  index_DB_Y_Z = 1:3,
  quanti = NULL,
  nominal = NULL,
  ordinal = NULL,
  logic = NULL,
  convert.num = NULL,
  convert.class = NULL,
  FAMD.coord = "NO",
  FAMD.perc = 0.8,
  dist.choice = "E",
  percent.knn = 1,
  maxrelax = 0,
  indiv.method = "sequential",
  prox.dist = 0,
  solvR = "glpk",
  which.DB = "BOTH"
)
```

## Arguments

- datab:

  a data.frame made up of two overlayed databases with at least four
  columns sorted in a random order. One column must be a column
  dedicated to the identification of the two databases ranked in
  ascending order (For example: 1 for the top database and 2 for the
  database from below, or more logically here A and B ...But not B and
  A!). One column (\\Y\\ here but other names are allowed) must
  correspond to the target variable related to the information of
  interest to merge with its specific encoding in the database A
  (corresponding encoding should be missing in the database B). In the
  same way, one column (\\Z\\ here) corresponds to the second target
  variable with its specific encoding in the database B (corresponding
  encoding should be missing in the database A). Finally, the input
  database must have at least one shared covariate with same encoding in
  A and B. Please notice that, if your data.frame has only one shared
  covariate (four columns) with missing values (because no imputation is
  desired) then a warning will appear and the algorithm will only run
  with complete cases.

- index_DB_Y_Z:

  a vector of three indexes of variables. The first index must
  correspond to the index of the databases identifier column. The second
  index corresponds to the index of the target variable in the first
  database (A) while the third index corresponds to the column index
  related to the target variable in the second database (B).

- quanti:

  a vector of column indexes of all the quantitative variables (database
  identifier and target variables included if it is the case for them).

- nominal:

  a vector of column indexes of all the nominal (not ordered) variables
  (database identifier and target variables included if it is the case
  for them).

- ordinal:

  a vector of column indexes of all the ordinal variables (database
  identifier and target variables included if it is the case for them).

- logic:

  a vector of column indexes of all the boolean variables of the
  data.frame.

- convert.num:

  indexes of the continuous (quantitative) variables to convert in
  ordered factors if necessary. All declared indexes in this argument
  must have been declared in the argument `quanti` (no conversion by
  default).

- convert.class:

  a vector indicating for each continuous variable to convert, the
  corresponding desired number of levels. If the length of the argument
  `convert_num` exceeds 1 while the length of `convert_class` equals 1
  (only one integer), each discretization will count the same number of
  levels (quantiles).

- FAMD.coord:

  a logical that must be set to TRUE when user decides to work with
  principal components of a factor analysis for mixed data (FAMD)
  instead of the set of raw covariates (FALSE is the default value).

- FAMD.perc:

  a percent (between 0 and 1) linked to the `FAMD.coord` argument (0.8
  is the default value). When this latter equals TRUE, this argument
  corresponds to the minimum part of variability that must be taken into
  account by the principal components of the FAMD method. This option
  fixes the remaining number of principal components for the rest of the
  study.

- dist.choice:

  a character string (with quotes) corresponding to the distance
  function chosen between: the euclidean distance ("E", by default), The
  Manhattan distance ("M"), the Gower distance ("G"), the Hamming
  distance ("H") for binary covariates only, and the Euclidean or
  Manhattan distance computed from principal components of a factor
  analysis of mixed data ("FAMD"). See (1) for details.

- percent.knn:

  the ratio of closest neighbors involved in the computations of the
  cost matrices. 1 is the default value that includes all rows in the
  computation.

- maxrelax:

  the maximum percentage of deviation from expected probability masses.
  It must be equal to 0 (default value) for the `OUTCOME` algorithm, and
  equal to a strictly positive value for the R-OUTCOME algorithm.
  See (2) for details.

- indiv.method:

  a character string indicating the chosen method to get individual
  predictions from the joint probabilities assessed, "sequential" by
  default, or "optimal". See the `details` section and (2) for details.

- prox.dist:

  a probability (between 0 and 1) used to calculate the distance
  threshold below which an individual (a row) is considered as a
  neighbor of a given profile of covariates. When shared variables are
  all factors or categorical, it is suggested to keep this option to 0.

- solvR:

  a character string that specifies the type of method selected to solve
  the optimization algorithms. The default solver is "glpk".

- which.DB:

  a character string indicating the database to complete ("BOTH" by
  default, for the prediction of \\Y\\ and \\Z\\ in the two databases),
  "A" only for the imputation of \\Z\\ in A, "B" only for the imputation
  of \\Y\\ in B.

## Value

A "otres" class object of 9 elements:

- time_exe:

  the running time of the function

- gamma_A:

  a matrix corresponding to an estimation of the joint distribution of
  \\(Y,Z)\\ in A

- gamma_B:

  a matrix corresponding to an estimation of the joint distribution of
  \\(Y,Z)\\ in B

- profile:

  a data.frame that gives all details about the remaining \\P\\ profiles
  of covariates. These informations can be linked to the `estimatorZA`
  and the `estimatorYB` objects for a better interpretation of the
  results.

- res_prox:

  the outputs of the function `proxim_dist`

- estimatorZA:

  an array that corresponds to estimates of the probability distribution
  of \\Z\\ conditional to \\X\\ and \\Y\\ in database A. The number of
  rows of each table corresponds to the total number of profiles of
  covariates. The first dimension of each table (rownames) correspond to
  the profiles of covariates sorted by order of appearance in the merged
  database. The second dimension of the array (columns of the tables)
  corresponds to the levels of \\Y\\ while the third element corresponds
  to the levels of \\Z\\.

- estimatorYB:

  an array that corresponds to estimates of the probability distribution
  of \\Y\\ conditional to \\X\\ and \\Z\\ in database B. The number of
  rows of each table corresponds to the total number of profiles of
  covariates. The first dimension of each table (rownames) correspond to
  the profiles of covariates sorted by order of appearance in the merged
  database. The second dimension of the array (columns of the tables)
  corresponds to the levels of \\Z\\ while the third element corresponds
  to the levels of \\Y\\.

- DATA1_OT:

  the database A with the individual predictions of \\Z\\ using an
  optimal transportation algorithm (`OUTCOME`) or `R-OUTCOME`

- DATA2_OT:

  the database B with the individual predictions of \\Y\\ using an
  optimal transportation algorithm (`OUTCOME`) or `R-OUTCOME`

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

B. INFORMATIONS ABOUT THE ALGORITHM

The algorithm integrated in the function `OT_outcome` provides a
solution to the recoding problem previously described by proposing an
application of optimal transportation which aims is to search for a
bijective mapping between the distributions of of \\Y\\ in A and \\Z\\
in B. Mathematically, the principle of the algorithm is based on the
resolution of an optimization problem which provides an optimal solution
\\\gamma\\ (as called in the related articles) that transfers the
distribution of \\Y\\ in A to the distribution of \\Z\\ in B (or
conversely, according to the sense of the transport)and can be so
interpreted as an estimator of the joint distribution \\(Y,Z)\\ in A (or
B respetively). According to this result, a second step of the algorithm
provides individual predictions of \\Y\\ in B (resp. of \\Z\\ in A, or
both, depending on the choice specified by user in the argument
`which.DB`). Two possible approaches are available depending on the
argument `indiv.method`:

- When `indiv.method = "sequential"`, a nearest neighbor procedure is
  applied. This corresponds to the use of the function
  [`indiv_grp_closest`](https://otrecoding.github.io/OTrecod/reference/indiv_grp_closest.md)
  implemented in the function `OT_outcome`.

- When `indiv.method = "optimal"`, a linear optimization problem is
  solved to determine the individual predictions that minimize the sum
  of the individual distances in A (resp. in B) with the modalities of
  \\Z\\ in B (resp. \\Y\\ in A). This approach is applied via the
  function
  [`indiv_grp_optimal`](https://otrecoding.github.io/OTrecod/reference/indiv_grp_optimal.md)
  implemented in the function `OT_outcome`.

This algorithm supposes the respect of the two following assumptions:

1.  \\Y\\ must follow the same distribution in A and B. In the same way,
    \\Z\\ follows the same distribution in the two databases.

2.  The conditional distribution \\(Y\|X)\\ must be identical in A
    and B. Respectively, \\(Z\|X)\\ is supposed identical in A and B.

Because the first assumption can be too strong in some situations, a
relaxation of the constraints of marginal distribution is possible using
the argument `maxrelax`. When `indiv.method = "sequential"` and
`maxrelax = 0`, the algorithm called `OUTCOME` (see (1) and (2)) is
applied. In all other situations, the algorithm applied corresponds to
an algorithm called `R_OUTCOME` (see (2)). A posteriori estimates of
conditional probabilities \\P\[Y\|X,Z\]\\ and \\P\[Z\|X,Y\]\\ are
available for each profile of covariates (see the output objects
`estimatorYB` and `estimatorZA`). Estimates of \\\gamma\\ are also
available according to the desired direction of the transport (from A to
B and/or conversely. See \\\gamma_A\\ and \\\gamma_B\\).

C. EXPECTED STRUCTURE FOR THE INPUT DATABASE

The input database is a data.frame that must be saved in a specific form
by users:

- Two overlayed databases containing a common column of database
  identifiers (A and B, 1 or 2, by examples, encoded in numeric or
  factor form)

- A column corresponding to the target variable with its specific
  encoding in A (For example a factor \\Y\\ encoded in \\n_Y\\ levels,
  ordered or not, with NAs in the corresponding rows of B)

- A column corresponding to the second target outcome with its specific
  endoded in B (For example a factor \\Z\\ in \\n_Z\\ levels, with NAs
  in rows of A)

- The order of the variables in the database have no importance but the
  column indexes related to the three columns previously described (ie
  ID, \\Y\\ and \\Z\\) must be rigorously specified in the argument
  `index_DB_Y_Z`.

- A set of shared common covariates (at least one but more is
  recommended) of any type, complete or not (provided that the number of
  covariates exceeds 1) is required.

The function
[`merge_dbs`](https://otrecoding.github.io/OTrecod/reference/merge_dbs.md)
is available in this package to assist user in the preparation of their
databases, so please, do not hesitate to use it beforehand if necessary.

Remarks about the target variables:

- A target variable can be of categorical type, but also discrete,
  stored in factor, ordered or not. Nevertheless, notice that, if the
  variable is stored in numeric it will be automatically converted in
  ordered factors.

- If a target outcome is incomplete, the corresponding rows will be
  automatically dropped during the execution of the function.

The type of each variables (including \\ID\\, \\Y\\ and \\Z\\) of the
database must be rigorously specified once, in one of the four arguments
`quanti`,`nominal`, `ordinal` and `logic`.

D. TRANSFORMATIONS OF CONTINUOUS COVARIATES

The function `OT_outcome` integrates in its syntax a process dedicated
to the categorization of continuous covariates. For this, it is
necessary to rigorously fill in the arguments `convert.num` and
`convert.class`. The first one informs about the indexes in database of
the continuous variables to transform in ordered factor while the second
one specifies the corresponding number of desired balanced levels (for
unbalanced levels, users must do transformations by themselves).
Therefore `convert.num` and `convert.class` must be vectors of same
length, but if the length of `convert.num` exceeds 1, while the length
of `convert.class` is 1, then, by default, all the covariates to convert
will have the same number of classes, that corresponds to the value
specified in the argument `convert.class`. Please notice that only
covariates can be transformed (not outcomes) and missing informations
are not taken into account for the transformations. Moreover, all the
indexes informed in the argument `convert.num` must also be informed in
the argument `quanti`.

E. INFORMATIONS ABOUT DISTANCE FUNCTIONS

Each individual (or row) of a given database is here characterized by
their covariates, so the distance between two individuals or groups of
individuals depends on similarities between covariates according to the
distance function chosen by user (via the argument `dist.choice`).
Actually four distance functions are implemented in `OT_outcome` to take
into account the most frequently encountered situation (see (3)):

- the Manhattan distance ("M")

- the Euclidean distance ("E")

- the Gower distance for mixed data (see (4): "G")

- the Hamming distance for binary data ("H")

Moreover, it is also possible to directly apply the first three
distances mentioned on coordinates extracted from a multivariate
analysis (Factor Analysis for Mixed Data, see (5)) applied on raw
covariates using the arguments `FAMD.coord` and `FAMD.perc`. This method
is used (1).

As a decision rule, for a given profile of covariates \\P_j\\, an
individual \\i\\ will be considered as a neighbor of \\P_j\\ if
\\dist(i,P_j) \< \mbox{prox.dist} \times max(dist(i,P_j))\\ where
\\prox.dist\\ must be fixed by user.

F. INFORMATIONS ABOUT THE SOLVER

The argument `solvR` permits user to choose the solver of the
optimization algorithm. The default solver is "glpk" that corresponds to
the GNU Linear Programming Kit (see (6) for more details). Moreover, the
function actually uses the `R` optimization infrastructure of the
package ROI which offers a wide choice of solver to users by easily
loading the associated plugins of ROI (see (7)).

For more details about the algorithms integrated in `OT_outcome`, please
consult (1) and (2).

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

3.  Anderberg, M.R. (1973), Cluster analysis for applications, 359 pp.,
    Academic Press, New York, NY, USA.

4.  Gower J.C. (1971). A general coefficient of similarity and some of
    its properties. Biometrics, 27, 623–637.

5.  Pages J. (2004). Analyse factorielle de donnees mixtes. Revue
    Statistique Appliquee. LII (4). pp. 93-111.

6.  Makhorin A (2011). GNU Linear Programming Kit Reference Manual
    Version 4.47.<http://www.gnu.org/software/glpk/>

7.  Theussl S, Schwendinger F, Hornik K (2020). ROI: An Extensible R
    Optimization Infrastructure.Journal of Statistical Software,94(15),
    1-64.
    [doi:10.18637/jss.v094.i15](https://doi.org/10.18637/jss.v094.i15)

## See also

[`transfo_dist`](https://otrecoding.github.io/OTrecod/reference/transfo_dist.md),[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md),
[`avg_dist_closest`](https://otrecoding.github.io/OTrecod/reference/avg_dist_closest.md),
[`indiv_grp_closest`](https://otrecoding.github.io/OTrecod/reference/indiv_grp_closest.md),
[`indiv_grp_optimal`](https://otrecoding.github.io/OTrecod/reference/indiv_grp_optimal.md)

## Author

Gregory Guernec, Valerie Gares, Jeremy Omer

<otrecod.pkg@gmail.com>

## Examples

``` r

### Using a sample of simu_data dataset
### Y and Z are a same variable encoded in 2 different forms:
### (3 levels for Y and 5 levels for Z)
#--------
data(simu_data)
simu_dat <- simu_data[c(1:200, 301:500), ]

### An example of OUTCOME algorithm that uses:
#-----
# - A nearest neighbor procedure for the estimation of individual predictions
# - The Manhattan distance function
# - 90% of individuals from each modalities to calculate average distances
#   between individuals and modalities
# Predictions are assessed for Y in B and Z in A
#-----

OUTC1 <- OT_outcome(simu_dat,
                    quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
                    dist.choice = "M", maxrelax = 0,
                    indiv.method = "sequential"
)
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = OUTCOME
#> Distance                 = Manhattan
#> Percent closest knn      = 100%
#> Relaxation parameter     = NO
#> Relaxation value         = 0
#> Individual pred process  = Sequential
#> DB imputed               = BOTH
#> ---------------------------------------
head(OUTC1$DATA1_OT) # Part of the completed database A
#>   DB       Y    Z Gender_2 Treatment_2 Treatment_3 Smoking_2 Dosage     Age
#> 1  A [40-60[ <NA>        0           1           0         1      3  3.0342
#> 2  A [20-40] <NA>        1          NA          NA         0      2  0.3472
#> 3  A [40-60[ <NA>        0           0           0         1      2 -0.1796
#> 4  A [40-60[ <NA>        0           0           1        NA      4  1.2620
#> 5  A [40-60[ <NA>        0           1           0         1      4 -1.0325
#> 6  A [20-40] <NA>        0           1           0        NA      1      NA
#>   OTpred
#> 1      4
#> 2      5
#> 3      1
#> 4      1
#> 5      1
#> 6      3
head(OUTC1$DATA2_OT) # Part of the completed database B
#>     DB    Y Z Gender_2 Treatment_2 Treatment_3 Smoking_2 Dosage     Age  OTpred
#> 301  B <NA> 5        0           0           0         1      2 -1.0701 [20-40]
#> 302  B <NA> 1        0           0           1        NA      4  2.9942 [40-60[
#> 303  B <NA> 2        0           0           0         0     NA  0.3189 [20-40]
#> 304  B <NA> 2        0           1           0         0     NA  0.0256 [20-40]
#> 305  B <NA> 1        0           0           1         1      4  2.2648 [40-60[
#> 306  B <NA> 1        0           0           1        NA      3  2.0712 [40-60[

head(OUTC1$estimatorZA[, , 1])
#>     [20-40] [40-60[ [60-80]
#> P_1     0.2     0.0     0.2
#> P_2     0.0     0.2     0.2
#> P_3     0.2     1.0     0.2
#> P_4     0.2     1.0     0.2
#> P_5     0.2     1.0     0.2
#> P_6     0.0     1.0     0.0
# ... Corresponds to P[Z = 1|Y,P1] when P1 corresponds to the 1st profile of covariates (P_1)
# detailed in the 1st row of the profile object:
OUTC1$profile[1, ] # Details of P_1
#>    ID Gender_2 Treatment_2 Treatment_3 Smoking_2 Dosage    Age
#> 1 P_1        0           1           0         1      3 3.0342

# So estimatorZA[1,1,1]= 0.2 corresponds to an estimation of:
# P[Z = 1|Y=[20-40],Gender_2=0,Treatment_2=1,Treatment_3=0,Smoking_2=1,Dosage=3,Age=65.44]
# Thus, we can conclude that all individuals with the P_1 profile of covariates have
# 20% of chance to be affected to the 1st level of Z in database A.
# ... And so on, the reasoning is the same for the estimatorYB object.

# \donttest{

### An example of OUTCOME algorithm with same conditions as the previous example, excepted that;
# - Only the individual predictions of Y in B are required
# - The continuous covariates "age" (related index = 8) will be converted in an ordinal factors
#   of 3 balanced classes (tertiles)
# - The Gower distance is now used
### -----

OUTC2_B <- OT_outcome(simu_dat,
                      quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
                      dist.choice = "G", maxrelax = 0,
                      convert.num = 8, convert.class = 3,
                      indiv.method = "sequential", which.DB = "B"
)
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = OUTCOME
#> Distance                 = Gower
#> Percent closest knn      = 100%
#> Relaxation parameter     = NO
#> Relaxation value         = 0
#> Individual pred process  = Sequential
#> DB imputed               = B
#> ---------------------------------------


### An example of OUTCOME algorithm with same conditions as the first example, excepted that;
# - Only the individual predictions of Z in A are required
# - The continuous covariates "age" (related index = 8) will be converted in an ordinal factors
#   of 3 balanced classes (tertiles)
# - Here, the Hamming distance can be applied because, after conversion, all covariates are factors.
#   Disjunctive tables of each covariates will be automatically used to work with a set of binary
#   variables.
### -----

OUTC3_B <- OT_outcome(simu_data,
                      quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
                      dist.choice = "H", maxrelax = 0,
                      convert.num = 8, convert.class = 3,
                      indiv.method = "sequential", which.DB = "B"
)
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = OUTCOME
#> Distance                 = Hamming
#> Percent closest knn      = 100%
#> Relaxation parameter     = NO
#> Relaxation value         = 0
#> Individual pred process  = Sequential
#> DB imputed               = B
#> ---------------------------------------
#> Hamming 1/3
#> Hamming 2/3
#> Hamming 3/3


### An example of R-OUTCOME algorithm using:
# - An optimization procedure for individual predictions on the 2 databases
# - The Manhattan distance
# - Raw covariates
### -----

R_OUTC1 <- OT_outcome(simu_data,
                      quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
                      dist.choice = "M", maxrelax = 0,
                      indiv.method = "optimal"
)
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = R-OUTCOME
#> Distance                 = Manhattan
#> Percent closest knn      = 100%
#> Relaxation parameter     = NO
#> Relaxation value         = 0
#> Individual pred process  = Optimal
#> DB imputed               = BOTH
#> ---------------------------------------


### An example of R-OUTCOME algorithm with:
# - An optimization procedure for individual predictions on the 2 databases
# - The use of Euclidean distance on coordinates from FAMD
# - Raw covariates
### -----

R_OUTC2 <- OT_outcome(simu_data,
                      quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
                      dist.choice = "E",
                      FAMD.coord = "YES", FAMD.perc = 0.8,
                      indiv.method = "optimal"
)
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = R-OUTCOME
#> Distance                 = Euclidean
#> Percent closest knn      = 100%
#> Relaxation parameter     = NO
#> Relaxation value         = 0
#> Individual pred process  = Optimal
#> DB imputed               = BOTH
#> ---------------------------------------


### An example of R-OUTCOME algorithm with relaxation on marginal distributions and:
# - An optimization procedure for individual predictions on the 2 databases
# - The use of the euclidean distance
# - An arbitrary coefficient of relaxation
# - Raw covariates
#-----

R_OUTC3 <- OT_outcome(simu_data,
                      quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
                      dist.choice = "E", maxrelax = 0.4,
                      indiv.method = "optimal"
)
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = R-OUTCOME
#> Distance                 = Euclidean
#> Percent closest knn      = 100%
#> Relaxation parameter     = YES
#> Relaxation value         = 0.4
#> Individual pred process  = Optimal
#> DB imputed               = BOTH
#> ---------------------------------------
# }
```
