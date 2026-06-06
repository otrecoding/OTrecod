# OT_joint()

The function `OT_joint` integrates two algorithms called (`JOINT`) and
(`R-JOINT`) dedicated to the solving of recoding problems in data fusion
using optimal transportation of the joint distribution of outcomes and
covariates.

## Usage

``` r
OT_joint(
  datab,
  index_DB_Y_Z = 1:3,
  nominal = NULL,
  ordinal = NULL,
  logic = NULL,
  convert.num = NULL,
  convert.class = NULL,
  dist.choice = "E",
  percent.knn = 1,
  maxrelax = 0,
  lambda.reg = 0,
  prox.X = 0.3,
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

  indexes of the continuous (quantitative) variables. They will be
  automatically converted in ordered factors. By default, no continuous
  variables is assumed in the database.

- convert.class:

  a vector indicating for each continuous variable to convert, the
  corresponding desired number of levels. If the length of the argument
  `convert_num` exceeds 1 while the length of `convert_class` equals 1
  (only one integer), each discretization will count the same number of
  levels (quantiles).

- dist.choice:

  a character string (with quotes) corresponding to the distance
  function chosen between: the euclidean distance ("E", by default), the
  Manhattan distance ("M"), the Gower distance ("G"), and the Hamming
  distance ("H") for binary covariates only.

- percent.knn:

  the ratio of closest neighbors involved in the computations of the
  cost matrices. 1 is the default value that includes all rows in the
  computation.

- maxrelax:

  the maximum percentage of deviation from expected probability masses.
  It must be equal to 0 (default value) for the `JOINT` algorithm, and
  equal to a strictly positive value for the R-JOINT algorithm.

- lambda.reg:

  a coefficient measuring the importance of the regularization term. It
  corresponds to the `R-JOINT` algorithm for a value other than 0
  (default value).

- prox.X:

  a probability (betwen 0 and 1) used to calculate the distance
  threshold below which two covariates' profiles are supposed as
  neighbors. If `prox.X = 1`, all profiles are considered as neighbors.

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

  running time of the function

- gamma_A:

  estimate of \\\gamma\\ for the completion of A. A matrix that
  corresponds to the joint distribution of \\(Y,Z,X)\\ in A

- gamma_B:

  estimate of \\\gamma\\ for the completion of B. A matrix that
  corresponds to the joint distribution of \\(Y,Z,X)\\ in B

- profile:

  a data.frame that gives all details about the remaining P profiles of
  covariates. These informations can be linked to the `estimatorZA` and
  the `estimatorYB` objects for a better interpretation of the results.

- res_prox:

  a `proxim_dist` object

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
  optimal transportation algorithm (`JOINT`) or `R-JOINT`

- DATA2_OT:

  the database B with the individual predictions of \\Y\\ using an
  optimal transportation algorithm (`JOINT`) or `R-JOINT`

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

As with the function
[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md),
the function `OT_joint` provides a solution to the recoding problem by
proposing an application of optimal transportation which aims is to
search for a bijective mapping between the joint distributions of
\\(Y,X)\\ and \\(Z,X)\\ in A and B (see (2) for more details). The
principle of the algorithm is also based on the resolution of an
optimization problem, which provides a solution \\\gamma\\ (as called in
(1) and (2)), estimate of the joint distribution of \\(X,Y,Z)\\
according to the database to complete (see the argument `which.DB` for
the choice of the database). While the algorithms `OUTCOME` and
`R_OUTCOME` integrated in the function
[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md)
require post-treatment steps to provide individual predictions, the
algorithm `JOINT` directly uses estimations of the conditional
distributions \\(Y\|Z,X)\\ in B and \\(Z\|Y,X)\\ in A to predict the
corresponding incomplete individuals informations of \\Y\\ and/or \\Z\\
respectively. This algorithm supposes that the conditional distribution
\\(Y\|X)\\ must be identical in A and B. Respectively, \\(Z\|X)\\ is
supposed identical in A and B. Estimations a posteriori of conditional
probabilities \\P\[Y\|X,Z\]\\ and \\P\[Z\|X,Y\]\\ are available for each
profiles of covariates in output (See the objects `estimatorYB` and
`estimatorZA`). Estimations of \\\gamma\\ are also available according
to the chosen transport distributions (See the arguments `gamma_A` and
`gamma_B`).

The algorithm `R-JOINT` gathers enrichments of the algorithm `JOINT` and
is also available via the function `OT_joint`. It allows users to add a
relaxation term in the algorithm to relax distributional assumptions
(`maxrelax`\>0), and (or) add also a positive regularization term
(`lamdba.reg`\>0) expressing that the transportation map does not vary
to quickly with respect of covariates \\X\\. Is is suggested to users to
calibrate these two parameters a posteriori by studying the stability of
the individual predictions in output.

C. EXPECTED STRUCTURE FOR THE INPUT DATABASE

The input database is a data.frame that must satisfy a specific form:

- Two overlayed databases containing a common column of databases
  identifiers (A and B, 1 or 2, by examples, encoded in numeric or
  factor form)

- A column corresponding to the target variable with its specific
  encoding in A (For example a factor \\Y\\ encoded in \\n_Y\\ levels,
  ordered or not, with NAs in the corresponding rows of B)

- A column corresponding to another target outcome summarizing the same
  latent information with its specific encoding in B (By example a
  factor \\Z\\ with \\n_Z\\ levels, with NAs in rows of A)

- The order of the variables in the database have no importance but the
  column indexes related to the three columns previously described (ie
  ID, \\Y\\ and \\Z\\) must be rigorously specified in the argument
  `index_DB_Y_Z`.

- A set of shared common categorical covariates (at least one but more
  is recommended) with or without missing values (provided that the
  number of covariates exceeds 1) is required. On the contrary to the
  function `OT_outcome`, please notice, that the function `OT_joint`
  does not accept continuous covariates therefore these latters will
  have to be categorized beforehand or using the provided input process
  (see `convert.num`).

The function
[`merge_dbs`](https://otrecoding.github.io/OTrecod/reference/merge_dbs.md)
is available in this package to assist user in the preparation of their
databases.

Remarks about the target variables:

- A target variable can be of categorical type, but also discrete,
  stored in factor, ordered or not. Nevertheless, notice that, if the
  variable is stored in numeric it will be automatically converted in
  ordered factors.

- If a target variable is incomplete, the corresponding rows will be
  automatically dropped during the execution of the function.

The type of each variables (including \\ID\\, \\Y\\ and \\Z\\) of the
database must be rigorously specified, in one of the four arguments
`quanti`, `nominal`, `ordinal` and `logic`.

D. TRANSFORMATIONS OF CONTINUOUS COVARIATES

Continuous shared variables (predictors) with infinite numbers of values
have to be categorized before being introduced in the function. To
assist users in this task, the function `OT_joint` integrates in its
syntax a process dedicated to the categorization of continuous
covariates. For this, it is necessary to rigorously fill in the
arguments `quanti` and `convert.class`. The first one informs about the
column indexes of the continuous variables to be transformed in ordered
factor while the second one specifies the corresponding number of
desired balanced levels (for unbalanced levels, users must do
transformations by themselves). Therefore `convert.num` and
`convert.class` must be vectors of same length, but if the length of
`quanti` exceeds 1, while the length of `convert.class` is 1, then, by
default, all the covariates to convert will have the same number of
classes (transformation by quantiles), that corresponds to the value
specified in the argument `convert.class`. Notice that only covariates
can be transformed (not target variables) and that any incomplete
information must have been taken into account beforehand (via the
dedicated functions
[`merge_dbs`](https://otrecoding.github.io/OTrecod/reference/merge_dbs.md)
or
[`imput_cov`](https://otrecoding.github.io/OTrecod/reference/imput_cov.md)
for examples). Moreover, all the indexes informed in the argument
`convert.num` must also be informed in the argument `quanti`. Finally,
it is recommended to declare all discrete covariates as ordinal factors
using the argument `ordinal`.

E. INFORMATIONS ABOUT DISTANCE FUNCTIONS AND RELATED PARAMETERS

Each individual (or row) of a given database is here characterized by a
vector of covariates, so the distance between two individuals or groups
of individuals depends on similarities between covariates according to
the distance function chosen by user (via the argument `dist.choice`).
Actually four distance functions are implemented in `OT_joint` to take
into account the most frequently encountered situation (see (3)):

- the Manhattan distance ("M")

- the Euclidean distance ("E")

- the Gower distance for mixed data (see (4): "G")

- the Hamming distance for binary data ("H")

Finally, two profiles of covariates \\P_1\\ (\\n_1\\ individuals) and
\\P_2\\ (\\n_2\\ individuals) will be considered as neighbors if
\\dist(P_1,P_2) \< prox.X \times max(dist(P_i,P_j))\\ where \\prox.X\\
must be fixed by user (\\i = 1,\dots,n_1\\ and \\j = 1,\dots,n_2\\).
This choice is used in the computation of the `JOINT` and `R_JOINT`
algorithms. The `prox.X` argument influences a lot the running time of
the algorithm. The greater, the more the value will be close to 1, the
more the convergence of the algorithm will be difficult or even
impossible.

Each individual \\i\\ from A or B is here considered as a neighbor of
only one profile of covariates \\P_j\\.

F. INFORMATIONS ABOUT THE SOLVER

The argument `solvR` permits user to choose the solver of the
optimization algorithm. The default solver is "glpk" that corresponds to
the GNU Linear Programming Kit (see (5) for more details). Moreover, the
function actually uses the `R` optimization infrastructure of the
package ROI which offers a wide choice of solver to users by easily
loading the associated plugins of ROI (see (6)).

For more details about the algorithms integrated in `OT_joint`, please
consult (2).

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
    its properties. Biometrics, 27, 623–637

5.  Makhorin A (2011). GNU Linear Programming Kit Reference Manual
    Version 4.47.<http://www.gnu.org/software/glpk/>

6.  Theussl S, Schwendinger F, Hornik K (2020). ROI: An Extensible R
    Optimization Infrastructure.Journal of Statistical Software,94(15),
    1-64.
    [doi:10.18637/jss.v094.i15](https://doi.org/10.18637/jss.v094.i15)

## See also

[`merge_dbs`](https://otrecoding.github.io/OTrecod/reference/merge_dbs.md),
[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md),
[`proxim_dist`](https://otrecoding.github.io/OTrecod/reference/proxim_dist.md),
[`avg_dist_closest`](https://otrecoding.github.io/OTrecod/reference/avg_dist_closest.md)

## Author

Gregory Guernec, Valerie Gares, Jeremy Omer

<otrecod.pkg@gmail.com>

## Examples

``` r

### An example of JOINT algorithm with:
#-----
# - A sample of the database tab_test
# - Y1 and Y2 are a 2 outcomes encoded in 2 different forms in DB 1 and 2:
#   4 levels for Y1 and 3 levels for Y2
# - n1 = n2 = 40
# - 2 discrete covariates X1 and X2 defined as ordinal
# - Distances estimated using the Gower function
# Predictions are assessed for Y1 in B only
#-----

data(tab_test)
tab_test2 <- tab_test[c(1:40, 5001:5040), 1:5]


OUTJ1_B <- OT_joint(tab_test2,
                    nominal = c(1, 4:5), ordinal = c(2, 3),
                    dist.choice = "G", which.DB = "B"
)
#> ---------------------------------------
#> OT JOINT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                  = JOINT
#> Distance              = Gower
#> Percent closest       = 100%
#> Relaxation term       = 0
#> Regularization term   = 0
#> Aggregation tol cov   = 0.3
#> DB imputed            = B
#> ---------------------------------------
#> Your target DB[, "Y"] was numeric ... By default, it has been converted in factor of integers
#> 4 remaining levels
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 3 remaining levels

# \donttest{

### An example of R-JOINT algorithm using the previous database,
### and keeping the same options excepted for:
#-----
# - The distances are estimated using the Gower function
# - Inclusion of an error term in the constraints on
#   the marginals (relaxation term)
# Predictions are assessed for Y1 AND Y2 in A and B respectively
#-----

R_OUTJ1 <- OT_joint(tab_test2,
                    nominal = c(1, 4:5), ordinal = c(2, 3),
                    dist.choice = "G", maxrelax = 0.4,
                    which.DB = "BOTH"
)
#> ---------------------------------------
#> OT JOINT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                  = R-JOINT
#> Distance              = Gower
#> Percent closest       = 100%
#> Relaxation term       = 0.4
#> Regularization term   = 0
#> Aggregation tol cov   = 0.3
#> DB imputed            = BOTH
#> ---------------------------------------
#> Your target DB[, "Y"] was numeric ... By default, it has been converted in factor of integers
#> 4 remaining levels
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 3 remaining levels

### The previous example of R-JOINT algorithm with:
# - Adding a regularization term
# Predictions are assessed for Y1 and Y2 in A and B respectively
#-----

R_OUTJ2 <- OT_joint(tab_test2,
                    nominal = c(1, 4:5), ordinal = c(2, 3),
                    dist.choice = "G", maxrelax = 0.4, lambda.reg = 0.9,
                    which.DB = "BOTH"
)
#> ---------------------------------------
#> OT JOINT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                  = R-JOINT
#> Distance              = Gower
#> Percent closest       = 100%
#> Relaxation term       = 0.4
#> Regularization term   = 0.9
#> Aggregation tol cov   = 0.3
#> DB imputed            = BOTH
#> ---------------------------------------
#> Your target DB[, "Y"] was numeric ... By default, it has been converted in factor of integers
#> 4 remaining levels
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 3 remaining levels



### Another example of JOINT algorithm with:
#-----
# - A sample of the database simu_data
# - Y1 and Y2 are a 2 outcomes encoded in 2 different forms in DB A and B:
#   (3 levels for Y and 5 levels for Z)
# - n1 = n2 = 100
# - 3 covariates: Gender, Smoking and Age in a qualitative form
# - Complete Case study
# - The Hamming distance
# Predictions are assessed for Y1 and Y2 in A and B respectively
#-----

data(simu_data)
simu_data2 <- simu_data[c(1:100, 401:500), c(1:4, 7:8)]
simu_data3 <- simu_data2[!is.na(simu_data2$Age), ]

OUTJ2 <- OT_joint(simu_data3, prox.X = 0.10,
                  convert.num = 6, convert.class = 3,
                  nominal = c(1, 4:5), ordinal = 2:3,
                  dist.choice = "H", which.DB = "B"
)
#> ---------------------------------------
#> OT JOINT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                  = JOINT
#> Distance              = Hamming
#> Percent closest       = 100%
#> Relaxation term       = 0
#> Regularization term   = 0
#> Aggregation tol cov   = 0.1
#> DB imputed            = B
#> ---------------------------------------
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 5 remaining levels
#> Hamming 1/3
#> Hamming 2/3
#> Hamming 3/3
# }
```
