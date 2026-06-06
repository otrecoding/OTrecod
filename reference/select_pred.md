# select_pred()

Selection of a subset of non collinear predictors having relevant
relationships with a given target outcome using a random forest
procedure.

## Usage

``` r
select_pred(
  databa,
  Y = NULL,
  Z = NULL,
  ID = 1,
  OUT = "Y",
  quanti = NULL,
  nominal = NULL,
  ordinal = NULL,
  logic = NULL,
  convert_num = NULL,
  convert_class = NULL,
  thresh_cat = 0.3,
  thresh_num = 0.7,
  thresh_Y = 0.2,
  RF = TRUE,
  RF_ntree = 500,
  RF_condi = FALSE,
  RF_condi_thr = 0.2,
  RF_SEED = sample(1:1e+06, 1)
)
```

## Arguments

- databa:

  a data.frame with a column of identifiers (of row or of database in
  the case of two concatened databases), an outcome, and a set of
  predictors. The number of columns can exceed the number of rows.

- Y:

  the label of a first target variable with quotes

- Z:

  the label of a second target variable with quotes when `databa` is the
  result of two overlayed databases.

- ID:

  the column index of the database identifier (The first column by
  default) in the case of two concatened databases, a row identifier
  otherwise

- OUT:

  a character that indicates the outcome to predict in the context of
  overlayed databases. By default, the outcome declared in the argument
  `Y` is predicted. Another possible outcome to predict can be set with
  the related argument `Z`.

- quanti:

  a vector of integers corresponding to the column indexes of all the
  numeric predictors.

- nominal:

  a vector of integers which corresponds to the column indexes of all
  the categorical nominal predictors.

- ordinal:

  a vector of integers which corresponds to the column indexes of all
  the categorical ordinal predictors.

- logic:

  a vector of integers indicating the indexes of logical predictors. No
  index remained by default

- convert_num:

  a vector of integers indicating the indexes of quantitative variables
  to convert in ordered factors. No index remained by default. Each
  index selected has to be defined as quantitative in the argument
  `quanti`.

- convert_class:

  a vector of integers indicating the number of classes related to each
  transformation of quantitative variable in ordered factor. The length
  of this vector can not exceed the length of the argument
  `convert_num`. Nevertheless, if length(`convert_num`) \> 1 and
  length(`convert_class`) = 1, all quantitative predictors selected for
  discretization will have by default the same number of classes.

- thresh_cat:

  a threshold associated to the Cramer's V coefficient (= 0.30 by
  default)

- thresh_num:

  a threshold associated to the Spearman's coefficient of correlation (=
  0.70 by default)

- thresh_Y:

  a threshold linked to the RF approach, that corresponds to the minimal
  cumulative percent of importance measure required to be kept in the
  final list of predictors.

- RF:

  a boolean sets to TRUE (default) if a random forest procedure must be
  applied to select the best subset of predictors according to the
  outcome.Otherwise, only pairwise associations between predictors are
  used for the selection.

- RF_ntree:

  the number of bootsrap samples required from the row datasource during
  the random forest procedure

- RF_condi:

  a boolean specifying if the conditional importance measures must be
  assessed from the random forest procedure (`TRUE`) rather than the
  standard variable importance measures (`FALSE` by default)

- RF_condi_thr:

  a threshold linked to (1 - pvalue) of an association test between each
  predictor \\X\\ and the other variables, given that a threshold value
  of zero will include all variables in the computation of the
  conditional importance measure of \\X\\ (0.20 is the default value).
  Conversely, a larger threshold will only keeps the subset of variables
  that is strongly correlated to \\X\\ for the computation of the
  variable importance measure of \\X\\.

- RF_SEED:

  an integer used as argument by the set.seed() for offsetting the
  random number generator (random integer by default). This value is
  only used for RF method.

## Value

A list of 14 (if `RF = TRUE`) or 11 objects (Only the first ten objects
if `RF = FALSE`) is returned:

- seed:

  the random number generator related to the study

- outc:

  the identifier of the outcome to predict

- thresh:

  a summarize of the different thresholds fixed for the study

- convert_num:

  the labels of the continuous predictors transformed in categorical
  form

- DB_USED:

  the final database used after potential transformations of predictors

- vcrm_OUTC_cat:

  a table of pairwise associations between the outcome and the
  categorical predictors (Cramer's V)

- cor_OUTC_num:

  a table of pairwise associations between the outcome and the
  continuous predictors (Rank correlation)

- vcrm_X_cat:

  a table of pairwise associations between the categorical predictors
  (Cramer's V)

- cor_X_num:

  a table of pairwise associations between the continuous predictors
  (Cramer's V)

- FG_test:

  the results of the Farrar and Glauber tests, with and without
  approximation form

- collinear_PB:

  a table of predictors with problem of collinearity according to the
  fixed thresholds

- drop_var:

  the labels of predictors to drop after RF process (optional output:
  only if `RF`=TRUE)

- RF_PRED:

  the table of variable importance measurements, conditional or not,
  according to the argument `condi_RF` (optional output: Only if
  `RF`=TRUE)

- RF_best:

  the labels of the best predictors selected (optional output: Only if
  `RF`=TRUE) according to the value of the argument `thresh_Y`

## Details

The `select_pred` function provides several tools to identify, on the
one hand, the relationships between predictors, by detecting especially
potential problems of collinearity, and, on the other hand, proposes a
parcimonious subset of relevant predictors (of the outcome) using
appropriate random forest procedures. The function which can be used as
a preliminary step of prediction in regression areas is particularly
adapted to the context of data fusion by providing relevant subsets of
predictors (the matching variables) to algorithms dedicated to the
solving of recoding problems.

A. REQUIRED STRUCTURE FOR THE DATABASE

The expected input database is a data.frame that especially requires a
specific column of row identifier and a target variable (or outcome)
having a finite number of values or classes (ordinal, nominal or
discrete type). Notice that if the chosen outcome is in numeric form, it
will be automatically converted in ordinal type. The number of
predictors is not a constraint for `select_pred` (even if, with less
than three variables a process of variables selection has no real
sense...), and can exceed the number of rows (no problem of high
dimensionality here). The predictors can be continuous (quantitative),
boolean, nominal or ordinal with or without missing values. In presence
of numeric variables, users can decide to discretize them or a part of
them by themselves beforehand. They can also choose to use the internal
process directly integrated in the function. Indeed, to assist users in
this task, two arguments called `convert_num` and `convert_class`
dedicated to these transformations are available in input of the
function. These options make the function `select_pred` particularly
adapted to the function
[`OT_joint`](https://otrecoding.github.io/OTrecod/reference/OT_joint.md)
which only allows data.frame with categorical covariates. With the
argument `convert_num`, users choose the continuous variables to convert
and the related argument `convert_class` specifies the corresponding
number of classes chosen for each discretization. It is the reason why
these two arguments must be two vectors of indexes of same length.
Nevertheless, an unique exception exists when `convert_class` is
equalled to a scalar \\S\\. In this case, all the continuous predictors
selected for conversion will be discretized with a same number of
classes S. By example, if `convert_class = 4`, all the continuous
variables specified in the `convert_num` argument will be discretized by
quartiles. Moreover, notice that missing values from incomplete
predictors to convert are not taken into account during the conversion,
and that each predictor specified in the argument `convert_num` must be
also specified in the argument `quanti`. In this situation, the label of
the outcome must be entered in the argument `Y`, and the arguments `Z`
and `OUT` must keep their default values. Finally, the order of the
column indexes related to the identifier and the outcome have no
importance.

For a better flexibility, the input database can also be the result of
two overlayed databases. In this case, the structure of the database
must be similar to those observed in the datasets
[`simu_data`](https://otrecoding.github.io/OTrecod/reference/simu_data.md)
and
[`tab_test`](https://otrecoding.github.io/OTrecod/reference/tab_test.md)
available in the package with a column of database identifier, one
target outcome by database (2 columns), and a subset of shared
predictors. Notice that, overlaying two separate databases can also be
done easily using the function
[`merge_dbs`](https://otrecoding.github.io/OTrecod/reference/merge_dbs.md)
beforehand. The labels of the two outcomes will have to be specified in
the arguments `Y` for the top database, and in `Z` for the bottom one.
Notice also that the function `select_pred` deals with only one outcome
at a time that will have to be specified in the argument `OUT` which
must be equalled to "Y" for the study of the top database or "Z" for the
study of the bottom one.

Finally, whatever the structure of the database declared in input, each
column index related to the database variable must be entered once (and
only once) in one of the following four arguments: `quanti`, `nominal`,
`ordinal`, `logic`.

B. PAIRWISE ASSOCIATIONS BETWEEN PREDICTORS

In a first step of process, `select_pred` calculates standard pairwise
associations between predictors according to their types.

1.  Between categorical predictors (ordinal, nominal and logical):
    Cramer's V (and Bias-corrected Cramer's V, see (1) for more details)
    are calculated between categorical predictors and the argument
    `thres_cat` fixed the associated threshold beyond which two
    predictors can be considered as redundant. A similar process is done
    between the target variable and the subset of categorical variables
    which provides in output a first table ranking the top scoring
    predictors. This table summarizes the ability of each variable to
    predict the target outcome.

2.  Between continuous predictors: If the `ordinal` and `logic`
    arguments differ from NULL, all the corresponding predictors are
    beforehand converted in rank values. For numeric (quantitative),
    logical and ordinal predictors, pairwise correlations between ranks
    (Spearman) are calculated and the argument `thresh_num` fixed the
    related threshold beyond which two predictors can be considered as
    redundant. A similar process is done between the outcome and the
    subset of discrete variables which provides in output, a table
    ranking the top scoring predictor variates which summarizes their
    abilities to predict the target. In addition, the result of a Farrar
    and Glauber test is provided. This test is based on the determinant
    of the correlation matrix of covariates and the related null
    hypothesis of the test corresponds to an absence of collinearity
    between them (see (2) for more details about the method). In
    presence of a large number of numeric covariates and/or ordered
    factors, the approximate Farrar-Glauber test, based on the normal
    approximation of the null distribution is more adapted and its
    result is also provided in output. These two tests are highly
    sensitive and, by consequence, it suggested to consider these
    results as simple indicators of collinearity between predictors
    rather than an essential condition of acceptability.

If the initial number of predictors is not too important, these
informations can be sufficient to the user for the visualization of
potential problems of collinearity and for the selection of a subset of
predictors (`RF = FALSE`). It is nevertheless often necessary to
complete this visualization by an automatical process of selection like
the Random Forest approach (see Breiman 2001, for a better understanding
of the method) linked to the function `select_pred` (`RF = TRUE`).

C. RANDOM FOREST PROCEDURE

As a final step of the process, a random forest approach (RF(3)) is here
prefered (to regression models) for two main reasons: RF methods allow
notably the number of variables to exceed the number of rows and remain
applicable whatever the types of covariates considered. The function
`select_pred` integrates in its algorithm the functions
[`cforest`](https://rdrr.io/pkg/party/man/cforest.html) and
[`varimp`](https://rdrr.io/pkg/party/man/varimp.html) of the package
party (Hothorn, 2006) and so gives access to their main arguments.

A RF approach generally provides two types of measures for estimating
the mean variable importance of each covariate in the prediction of an
outcome: the Gini importance and the permutation importance. These
measurements must be used with caution, by taking into account the
following constraints:

1.  The Gini importance criterion can produce bias in favor of
    continuous variables and variables with many categories. To avoid
    this problem, only the permutation criterion is available in the
    function.

2.  The permutation importance criterion can overestimate the importance
    of highly correlated predictors.

The function `select_pred` proposes three different scenarios according
to the types of predictors:

1.  The first one consists in boiling down to a set of categorical
    variables (ordered or not) by discretizing all the continuous
    predictors beforehand, using the internal `convert_num` argument or
    another one, and then works with the conditional importance measures
    (`RF_condi = TRUE`) which give unbiased estimations. In the spirit
    of a partial correlation, the conditional importance measure related
    to a variable \\X\\ for the prediction of an outcome \\Y\\, only
    uses the subset of variables the most correlated to \\X\\ for its
    computation. The argument `RF_condi_thr` that corresponds exactly to
    the argument `threshold` of the function
    [`varimp`](https://rdrr.io/pkg/party/man/varimp.html), fixes a ratio
    below which a variable Z is considered sufficiently correlated to
    \\X\\ to be used as an adjustment variable in the computation of the
    importance measure of \\X\\ (In other words, Z is included in the
    conditioning for the computation, see (4) and (5) for more details).
    A threshold value of zero will include all variables in the
    computation of conditional importance measure of each predictor
    \\X\\, while a threshold \\\< 1\\, will only include a subset of
    variables. Two remarks related to this method: firstly, notice that
    taking into account only subsets of predictors in the computation of
    the variable importance measures could lead to a relevant saving of
    execution time. Secondly, because this approach does not take into
    account incomplete information, the method will only be applied to
    complete data (incomplete rows will be temporarily removed for the
    study).

2.  The second possibility, always in presence of mixed types
    predictors, consists in the execution of two successive RF
    procedures. The first one will be used to select an unique candidate
    in each susbset of correlated predictors (detecting in the 1st
    section), while the second one will extract the permutation measures
    from the remaining subset of uncorrelated predictors
    (`RF_condi = FALSE`, by default). This second possibility has the
    advantage to work in presence of incomplete predictors.

3.  The third scenario consists in running a first time the function
    without RF process (`RF = FALSE`), and according to the presence of
    highly correlated predictors or not, users can choose to extract
    redundant predictors manually and re-runs the function with the
    subset of remaining non-collinear predictors to avoid potential
    biases introduced by the standard permutations measures.

The three scenarios finally lead to a list of uncorrelated predictors of
the outcome sorted in importance order. The argument `thresh_Y`
corresponds to the minimal percent of importance required (and fixed by
user) for a variable to be considered as a reliable predictor of the
outcome. Finally, because all random forest results are subjects to
random variation, users can check whether the same importance ranking is
achieved by varying the random seed parameter (`RF_SEED`) or by
increasing the number of trees (`RF_ntree`).

## References

1.  Bergsma W. (2013). A bias-correction for Cramer's V and
    Tschuprow's T. Journal of the Korean Statistical Society, 42,
    323–328.

2.  Farrar D, and Glauber R. (1968). Multicolinearity in regression
    analysis. Review of Economics and Statistics, 49, 92–107.

3.  Breiman L. (2001). Random Forests. Machine Learning, 45(1), 5–32.

4.  Hothorn T, Buehlmann P, Dudoit S, Molinaro A, Van Der Laan M (2006).
    “Survival Ensembles.” Biostatistics, 7(3), 355–373.

5.  Strobl C, Boulesteix A-L, Kneib T, Augustin T, Zeileis A (2008).
    Conditional Variable Importance for Random Forests. BMC
    Bioinformatics, 9, 307.
    <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-307>

## See also

[`simu_data`](https://otrecoding.github.io/OTrecod/reference/simu_data.md),
[`tab_test`](https://otrecoding.github.io/OTrecod/reference/tab_test.md),
[`OT_joint`](https://otrecoding.github.io/OTrecod/reference/OT_joint.md)

## Author

Gregory Guernec

<gregory.guernec@inserm.fr>

## Examples

``` r

### Example 1
#-----
# - From two overlayed databases: using the table simu_data
# - Searching for the best predictors of "Yb1"
# - Using the row database
# - The RF approaches are not required
#-----

data(simu_data)
sel_ex1 <- select_pred(simu_data,
  Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Y",
  quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
  thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
  RF = FALSE
)
#> The select_pred function is running for outcome= Yb1. Please wait ...
#> The process is now successfully completed
#> ---------
#> For comparison with another outcome from two overlayed tables  : 
#> just adapt the OUT option keeping all the others unchanged in the function
#> ---
#> For comparison with another outcome from two unoverlayed tables:
#> just adapt the arguments from Y to convert_class
#> ---------

### Example 2
#-----
# - With same conditions as example 1
# - Searching for the best predictors of "Yb2"
#-----

sel_ex2 <- select_pred(simu_data,
  Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Z",
  quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
  thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
  RF = FALSE
)
#> The select_pred function is running for outcome= Yb2. Please wait ...
#> The process is now successfully completed
#> ---------
#> For comparison with another outcome from two overlayed tables  : 
#> just adapt the OUT option keeping all the others unchanged in the function
#> ---
#> For comparison with another outcome from two unoverlayed tables:
#> just adapt the arguments from Y to convert_class
#> ---------

# \donttest{
### Example 3
#-----
# - With same conditions as example 1
# - Using a RF approach to estimate the standard variable importance measures
#   and determine the best subset of predictors
# - Here a seed is required
#-----

sel_ex3 <- select_pred(simu_data,
  Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Y",
  quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
  thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
  RF = TRUE, RF_condi = FALSE, RF_SEED = 3023
)
#> The select_pred function is running for outcome= Yb1. Please wait ...
#> Risk of collinearity between predictors detected: Some predictors will be dropped during RF process
#> The process is now successfully completed
#> ---------
#> For comparison with another outcome from two overlayed tables  : 
#> just adapt the OUT option keeping all the others unchanged in the function
#> ---
#> For comparison with another outcome from two unoverlayed tables:
#> just adapt the arguments from Y to convert_class
#> ---------

### Example 4
#-----
# - With same conditions as example 1
# - Using a RF approach to estimate the conditional variable importance measures
#   and determine the best subset of predictors
# - This approach requires to convert the numeric variables: Only "Age" here
#   discretized in 3 levels
#-----

sel_ex4 <- select_pred(simu_data,
  Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Z",
  quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
  convert_num = 8, convert_class = 3,
  thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
  RF = TRUE, RF_condi = TRUE, RF_condi_thr = 0.60, RF_SEED = 3023
)
#> The select_pred function is running for outcome= Yb2. Please wait ...
#> Risk of collinearity between predictors detected: Some predictors will be dropped during RF process
#> The process is now successfully completed
#> ---------
#> For comparison with another outcome from two overlayed tables  : 
#> just adapt the OUT option keeping all the others unchanged in the function
#> ---
#> For comparison with another outcome from two unoverlayed tables:
#> just adapt the arguments from Y to convert_class
#> ---------

### Example 5
#-----
# - Starting with a unique database
# - Same conditions as example 1
#-----
simu_A <- simu_data[simu_data$DB == "A", -3] # Base A

sel_ex5 <- select_pred(simu_A,
  Y = "Yb1",
  quanti = 7, nominal = c(1, 3:4, 6), ordinal = c(2, 5),
  thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
  RF = FALSE
)
#> The select_pred function is running for outcome= Yb1. Please wait ...
#> The process is now successfully completed
#> ---------
#> For comparison with another outcome from two overlayed tables  : 
#> just adapt the OUT option keeping all the others unchanged in the function
#> ---
#> For comparison with another outcome from two unoverlayed tables:
#> just adapt the arguments from Y to convert_class
#> ---------

### Example 6
#-----
# - Starting with an unique database
# - Using a RF approach to estimate the conditional variable importance measures
#   and determine the best subset of predictors
# - This approach requires to convert the numeric variables: Only "Age" here
#   discretized in 3 levels
#-----

simu_B <- simu_data[simu_data$DB == "B", -2] # Base B

sel_ex6 <- select_pred(simu_B,
  Y = "Yb2",
  quanti = 7, nominal = c(1, 3:4, 6), ordinal = c(2, 5),
  convert_num = 7, convert_class = 3,
  thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
  RF = TRUE, RF_condi = TRUE, RF_condi_thr = 0.60, RF_SEED = 3023
)
#> The select_pred function is running for outcome= Yb2. Please wait ...
#> Risk of collinearity between predictors detected: Some predictors will be dropped during RF process
#> The process is now successfully completed
#> ---------
#> For comparison with another outcome from two overlayed tables  : 
#> just adapt the OUT option keeping all the others unchanged in the function
#> ---
#> For comparison with another outcome from two unoverlayed tables:
#> just adapt the arguments from Y to convert_class
#> ---------
# }
```
