# merge_dbs()

Harmonization and merging before data fusion of two databases with
specific outcome variables and shared covariates.

## Usage

``` r
merge_dbs(
  DB1,
  DB2,
  row_ID1 = NULL,
  row_ID2 = NULL,
  NAME_Y,
  NAME_Z,
  order_levels_Y = levels(DB1[, NAME_Y]),
  order_levels_Z = levels(DB2[, NAME_Z]),
  ordinal_DB1 = NULL,
  ordinal_DB2 = NULL,
  impute = "NO",
  R_MICE = 5,
  NCP_FAMD = 3,
  seed_choice = sample(1:1e+06, 1)
)
```

## Arguments

- DB1:

  a data.frame corresponding to the 1st database to merge (top database)

- DB2:

  a data.frame corresponding to the 2nd database to merge (bottom
  database)

- row_ID1:

  the column index of the row identifier of DB1 if it exists (no
  identifier by default)

- row_ID2:

  the column index of the row identifier of DB2 if it exists (no
  identifier by default)

- NAME_Y:

  the name of the outcome (with quotes) in its specific scale/encoding
  from the 1st database (DB1)

- NAME_Z:

  the name of the outcome (with quotes) in its specific scale/encoding
  from the 2nd database (DB2)

- order_levels_Y:

  the levels of \\Y\\ stored in a vector and sorted in ascending order
  in the case of ordered factors. This option permits to reorder the
  levels in the 1st database (DB1) if necessary.

- order_levels_Z:

  the levels of \\Z\\ stored in a vector and sorted in ascending order
  in the case of ordered factors. This option permits to reorder the
  levels in the 2nd database (DB2) if necessary.

- ordinal_DB1:

  a vector of column indexes corresponding to ordinal variables in the
  1st database (no ordinal variable by default)

- ordinal_DB2:

  a vector of column indexes corresponding to ordinal variables in the
  2nd database (no ordinal variable by default)

- impute:

  a character equals to "NO" when missing data on covariates are kept
  (Default option), "CC" for Complete Case by keeping only covariates
  with no missing information , "MICE" for MICE multiple imputation
  approach, "FAMD" for single imputation approach using Factorial
  Analysis for Mixed Data

- R_MICE:

  the chosen number of multiple imputations required for the MICE
  approach (5 by default)

- NCP_FAMD:

  an integer corresponding to the number of components used to predict
  missing values in FAMD imputation (3 by default)

- seed_choice:

  an integer used as argument by the set.seed() for offsetting the
  random number generator (Random integer by default, only useful with
  MICE)

## Value

A list containing 12 elements (13 when `impute` equals "MICE"):

- DB_READY:

  the database matched from the two initial databases with common
  covariates and imputed or not according to the impute option

- ID1_drop:

  the row numbers or row identifiers excluded of the data merging
  because of the presence of missing values in the target variable of
  DB1. NULL otherwise

- ID2_drop:

  the row numbers or row identifiers excluded of the data merging
  because of the presence of missing values in the target variable of
  DB2. NULL otherwise

- Y_LEVELS:

  the remaining levels of the target variable \\Y\\ in the DB1

- Z_LEVELS:

  the remaining Levels of the target variable \\Z\\ in the DB2

- REMOVE1:

  the labels of the deleted covariates because of type incompatibilies
  of type from DB1 to DB2

- REMOVE2:

  the removed factor(s) because of levels incompatibilities from DB1 to
  DB2

- REMAINING_VAR:

  labels of the remained covariates for data fusion

- IMPUTE_TYPE:

  a character with quotes that specify the method eventually chosen to
  handle missing data in covariates

- MICE_DETAILS:

  a list containing the details of the imputed datasets using `MICE`
  when this option is chosen. Raw and imputed databases imputed for DB1
  and DB2 according to the number of multiple imputation selected (Only
  if impute = "MICE")

- DB1_raw:

  a data.frame corresponding to DB1 after merging

- DB2_raw:

  a data.frame corresponding to DB2 after merging

- SEED:

  an integer used as argument by the `set.seed` function for offsetting
  the random number generator (random selection by default)

## Details

Assuming that DB1 and DB2 are two databases (two separate data.frames
with no overlapping rows) to be merged vertically before data fusion,
the function `merge_dbs` performs this merging and checks the
harmonization of the shared variables. Firslty, the two databases
declared as input to the function (via the argument `DB1` and `DB2`)
must have the same specific structure. Each database must contain a
target variable (whose label must be filled in the argument `Y` for DB1
and in `Z` for DB2 respectively, so that the final synthetic database in
output will contain an incomplete variable `Y` whose corresponding
values will be missing in DB2 and another incomplete target `Z` whose
values will be missing in DB1), a subset of shared covariates (by
example, the best predictors of \\Y\\ in DB1, and \\Z\\ in DB2). Each
database can have a row identifier whose label must be assigned in the
argument `row_ID1` for DB1 and `row_ID2` for DB2. Nevertheless, by
default DB1 and DB2 are supposed with no row identifiers. The merging
keeps unchanged the order of rows in the two databases provided that
\\Y\\ and \\Z\\ have no missing values. By building, the first declared
database (in the argument `DB1`) will be placed automatically above the
second one (declared in the argument `DB2`) in the final database.

Firstly, by default, a variable with the same name in the two databases
is abusively considered as shared. This condition is obviously
insufficient to be kept in the final subset of shared variables, and the
function `merge_dbs` so performs checks before merging described below.

A. Discrepancies between shared variables

- Shared variables with discrepancies of types between the two databases
  (for example, a variable with a common name in the two databases but
  stored as numeric in DB1, and stored as character in DB2) will be
  removed from the merging and the variable name will be saved in output
  (`REMOVE1`).

- Shared factors with discrepancies of levels (or number of levels) will
  be also removed from the merging and the variable name will be saved
  in output (`REMOVE2`).

- covariates whose names are specific to each database will be also
  deleted from the merging.

- If some important predictors have been improperly excluded from the
  merging due to the above-mentioned checks, it is possible for user to
  transform these variables a posteriori, and re-run the function.

B. Rules for the two outcomes (target variables)

The types of `Y` and `Z` must be suitable:

- Categorical (ordered or not) factors are allowed.

- Numeric and discrete outcomes with a finite number of values are
  allowed but will be automatically converted as ordered factors using
  the function `transfo_target` integrated in the function `merge_dbs`.

C. The function `merge_dbs` handles incomplete information of shared
variables, by respecting the following rules:

- If `Y` or `Z` have missing values in DB1 or DB2, corresponding rows
  are excluded from the database before merging. Moreover, in the case
  of incomplete outcomes, if A and B have row identifiers, the
  corresponding identifiers are removed and these latters are stored in
  the objects `DB1_ID` and `DB2_ID` of the output.

- Before overlay, the function deals with incomplete covariates
  according to the argument `impute`. Users can decide to work with
  complete case only ("CC"), to keep ("NO") or impute incomplete
  information ("MICE","FAMD").

- The function `imput_cov`, integrated in the syntax of `merge_dbs`
  deals with imputations. Two approaches are actually available: the
  multivariate imputation by chained equation approach (MICE, see (3)
  for more details about the approach or the corresponding package
  mice), and an imputation approach from the package missMDA that uses a
  dimensionality reduction method (here a factor analysis for mixed data
  called FAMD (4)), to provide single imputations. If multiple
  imputation is required (`impute` = "MICE"), the default imputation
  methods are applied according to the type of the variables. The
  average of the plausible values will be kept for a continuous
  variable, while the most frequent candidate will be kept as a
  consensus value for a categorical variable or factor (ordinal or not).

As a finally step, the function checks that all values related to \\Y\\
in B are missing and inversely for \\Z\\ in A.

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

3.  van Buuren, S., Groothuis-Oudshoorn, K. (2011). mice: Multivariate
    Imputation by Chained Equations in R. Journal of Statistical
    Software, 45(3), 1–67. urlhttps://www.jstatsoft.org/v45/i03/

4.  Josse J, Husson F (2016). missMDA: A Package for Handling Missing
    Values in Multivariate Data Analysis. Journal of Statistical
    Software, 70(1), 1–31.
    [doi:10.18637/jss.v070.i01](https://doi.org/10.18637/jss.v070.i01)

## See also

[`imput_cov`](https://otrecoding.github.io/OTrecod/reference/imput_cov.md),
[`transfo_target`](https://otrecoding.github.io/OTrecod/reference/transfo_target.md),
[`select_pred`](https://otrecoding.github.io/OTrecod/reference/select_pred.md)

## Author

Gregory Guernec

<otrecod.pkg@gmail.com>

## Examples

``` r

### Assuming two distinct databases from simu_data: data_A and data_B
### Some transformations will be made beforehand on variables to generate
### heterogeneities between the two bases.
data(simu_data)
data_A <- simu_data[simu_data$DB == "A", c(2, 4:8)]
data_B <- simu_data[simu_data$DB == "B", c(3, 4:8)]

# For the example, a covariate is added (Weight) only in data_A
data_A$Weight <- rnorm(300, 70, 5)

# Be careful: the target variables must be in factor (or ordered) in the 2 databases
# Because it is not the case for Yb2 in data_B, the function will convert it.
data_B$Yb2 <- as.factor(data_B$Yb2)

# Moreover, the Dosage covariate is stored in 3 classes in data_B (instead of 4 classes in data_B)
# to make the encoding of this covariate specific to each database.
data_B$Dosage <- as.character(data_B$Dosage)
data_B$Dosage <- as.factor(ifelse(data_B$Dosage %in% c("Dos 1", "Dos 2"), "D1",
  ifelse(data_B$Dosage == "Dos 3", "D3", "D4")
))

# For more diversity, this covariate iis placed at the last column of the data_B
data_B <- data_B[, c(1:3, 5, 6, 4)]

# Ex 1: The two databases are merged and incomplete covariates are imputed using MICE
merged_ex1 <- merge_dbs(data_A, data_B,
  NAME_Y = "Yb1", NAME_Z = "Yb2",
  ordinal_DB1 = c(1, 4), ordinal_DB2 = c(1, 6),
  impute = "MICE", R_MICE = 2, seed_choice = 3011)
#> DBS MERGING in progress. Please wait ...
#> DBS MERGING OK
#> -----------------------
#> 
#> SUMMARY OF DBS MERGING:
#> Nb of removed subjects due to NA on targets: 0(0%)
#> Nb of removed covariates due to differences between the 2 bases: 1
#> Nb of remained covariates: 4
#> Imputation on incomplete covariates: MICE

summary(merged_ex1$DB_READY)
#>        DB              Y         Z            Age           Gender    Smoking  
#>  Min.   :1.000   [20-40]:130   1  :145   Min.   :34.37   Female:414   NO :365  
#>  1st Qu.:1.000   [40-60[:127   2  : 51   1st Qu.:46.81   Male  :286   YES:335  
#>  Median :2.000   [60-80]: 43   3  : 35   Median :50.09                         
#>  Mean   :1.571   NAs    :400   4  : 43   Mean   :50.07                         
#>  3rd Qu.:2.000                 5  :126   3rd Qu.:53.31                         
#>  Max.   :2.000                 NAs:300   Max.   :65.44                         
#>    Treatment  
#>  Placebo:356  
#>  Trt A  :212  
#>  Trt B  :132  
#>               
#>               
#>               


# Ex 2: The two databases are merged and missing values are kept
merged_ex2 <- merge_dbs(data_A, data_B,
  NAME_Y = "Yb1", NAME_Z = "Yb2",
  ordinal_DB1 = c(1, 4), ordinal_DB2 = c(1, 6),
  impute = "NO", seed_choice = 3011
)
#> DBS MERGING in progress. Please wait ...
#> DBS MERGING OK
#> -----------------------
#> 
#> SUMMARY OF DBS MERGING:
#> Nb of removed subjects due to NA on targets: 0(0%)
#> Nb of removed covariates due to differences between the 2 bases: 1
#> Nb of remained covariates: 4
#> Imputation on incomplete covariates: NO

# Ex 3: The two databases are merged by only keeping the complete cases
merged_ex3 <- merge_dbs(data_A, data_B,
  NAME_Y = "Yb1", NAME_Z = "Yb2",
  ordinal_DB1 = c(1, 4), ordinal_DB2 = c(1, 6),
  impute = "CC", seed_choice = 3011
)
#> DBS MERGING in progress. Please wait ...
#> DBS MERGING OK
#> -----------------------
#> 
#> SUMMARY OF DBS MERGING:
#> Nb of removed subjects due to NA on targets: 0(0%)
#> Nb of removed covariates due to differences between the 2 bases: 1
#> Nb of remained covariates: 4
#> Imputation on incomplete covariates: CC

# Ex 4: The two databases are merged and incomplete covariates are imputed using FAMD
merged_ex4 <- merge_dbs(data_A, data_B,
  NAME_Y = "Yb1", NAME_Z = "Yb2",
  ordinal_DB1 = c(1, 4), ordinal_DB2 = c(1, 6),
  impute = "FAMD", NCP_FAMD = 4, seed_choice = 2096
)
#> DBS MERGING in progress. Please wait ...
#> DBS MERGING OK
#> -----------------------
#> 
#> SUMMARY OF DBS MERGING:
#> Nb of removed subjects due to NA on targets: 0(0%)
#> Nb of removed covariates due to differences between the 2 bases: 1
#> Nb of remained covariates: 4
#> Imputation on incomplete covariates: FAMD

# Conclusion:
# The data fusion is successful in each situation.
# The Dosage and Weight covariates have been normally excluded from the fusion.
# The covariates have been imputed when required.
```
