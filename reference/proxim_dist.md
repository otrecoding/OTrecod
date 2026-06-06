# proxim_dist()

`proxim_dist` computes the pairwise distance matrix of a database and
cross-distance matrix between two databases according to various
distances used in the context of data fusion.

## Usage

``` r
proxim_dist(data_file, indx_DB_Y_Z = 1:3, norm = "E", prox = 0.3)
```

## Arguments

- data_file:

  a data.frame corresponding ideally to an output object of the function
  [`transfo_dist`](https://otrecoding.github.io/OTrecod/reference/transfo_dist.md).
  Otherwise this data.frame is the result of two overlayed databases
  with a column of database identifier ("A" and "B", 1 and 2, for
  example), a target variable (called \\Y\\ by example) only known in
  the first database, a target variable (\\Z\\) only stored in the
  second database, such that \\Y\\ and \\Z\\ summarize a same
  information differently encoded in the two databases and set of common
  covariates (at least one) of any type. The order of the variables in
  the data.frame have no importance. The type of the covariates must be
  in accordance with the chosen distance function in the `norm` option.

- indx_DB_Y_Z:

  a vector of three column indexes corresponding to the database
  identifier, the target variable of the above database and the target
  variable of the below database. The indexes must be declared in this
  specific order.

- norm:

  a character string indicating the choice of the distance function.
  This latest depends on the type of the common covariates: the Hamming
  distance for binary covariates only (`norm` = "H"), the Manhattan
  distance ("M", by default) and the euclidean distance ("E") for
  continuous covariates only, or the Gower distance for mixed covariates
  ("G").

- prox:

  a ratio (betwen 0 and 1) used to calculate the distance threshold
  below which an individual (a row or a given statistical unit) is
  considered as a neighbor of a given profile of covariates. 0.3 is the
  default value.

## Value

A list of 16 elements (the first 16 detailed below) is returned
containing various distance matrices and lists useful for the algorithms
that used Optimal Transportation theory. Two more objects (the last two
of the following list) will be returned if distance matrices contain
NAs.

- FILE_NAME:

  a simple reminder of the name of the raw database

- nA:

  the number of rows of the first database (A)

- nB:

  the number of rows of the second database (B)

- Xobserv:

  the subset of the two overlayed databases composed of the shared
  variables only

- profile:

  the different encountered profiles of covariates according to the
  data.frame

- Yobserv:

  the numeric values of the target variable in the first database

- Zobserv:

  the numeric values of the target variable in the second database

- D:

  a distance matrix corresponding to the computed distances between
  individuals of the two databases

- Y:

  the \\n_Y\\ levels of the target variable in numeric form, in the
  first database

- Z:

  the \\n_Z\\ levels of the target variable in numeric form, in the
  second database

- indY:

  a list of \\n_Y\\ groups of individual (or row) numbers where each
  group corresponds to the individuals indexes related to a given level
  of \\Y\\ in the first database

- indZ:

  a list of \\n_Z\\ groups of individual (or row) numbers where each
  group corresponds to the individuals indexes related to a given level
  of \\Z\\ in the second database

- indXA:

  a list of individual (row) indexes from the first database, sorted by
  profiles of covariates according to their proximities. See the
  `Details` part for more information

- indXB:

  a list of individual (row) indexes from the second database, sorted by
  profiles of covariates according to their proximities. See the
  `Details` part for more information

- DA:

  a distance matrix corresponding to the pairwise distances between
  individuals of the first database

- DB:

  a distance matrix corresponding to the pairwise distances between
  individuals of the second database

- ROWS_TABLE:

  combinations of row numbers of the two databases that generate NAs in
  D

- ROWS_TO_RM:

  number of times a row of the first or second database is involved in
  the NA process of D

## Details

This function is the first step of a family of algorithms that solve
recoding problems of data fusion using optimal transportation theory
(see the details of these corresponding models `OUTCOME`, `R_OUTCOME`,
`JOINT` and `R_JOINT` in (1) and (2)). The function `proxim_dist` is
directly implemented in the functions
[`OT_outcome`](https://otrecoding.github.io/OTrecod/reference/OT_outcome.md)
and
[`OT_joint`](https://otrecoding.github.io/OTrecod/reference/OT_joint.md)
but can also be used separately as long as the input database has as
suitable structure. Nevertheless, its preparation will have to be
rigorously made in two steps detailled in the following sections.

A. EXPECTED STRUCTURE FOR THE INPUT DATABASE

Firsly, the initial database required is a data.frame that must be
prepared in a specific form by users. From two separate databases, the
function
[`merge_dbs`](https://otrecoding.github.io/OTrecod/reference/merge_dbs.md)
available in this package can assist users in this initial merging,
nevertheless notice that this preliminary transformation can also be
made directly by following the imposed structure described below: two
overlayed databases containing a common column of database identifiers
(A and B for examples, encoded in numeric or factor form), a column
corresponding to the target variable with its specific encoding in A
(for example a factor \\Y\\ encoded in \\n_Y\\ levels, ordered or not,
with NAs in the corresponding rows of B), a column corresponding to the
same variable with its specific endoded in B (for example a factor \\Z\\
in \\n_Z\\ levels, with NAs in database A), and a set of shared
covariates (at least one) between the two databases.

The order of these variables in the database have no importance but the
column indexes related to database identifier, \\Y\\ and \\Z\\, must be
specified in the `indx_DB_Y_Z` option. Users can refer to the structure
of the table
[`simu_data`](https://otrecoding.github.io/OTrecod/reference/simu_data.md)
available in the package to adapt their databases to the inital format
required.

Missing values are allowed on covariates only, and are excluded from all
computations involving the rows within which they occur. In the
particular case where only one covariate with NAs is used, we recommend
working with imputed or complete case only to avoid the presence of NA
in the distance matrix that will be computed a posteriori. If the
database counts many covariates and some of them have missing data, user
can keep them or apply beforehand the
[`imput_cov`](https://otrecoding.github.io/OTrecod/reference/imput_cov.md)
function on data.frame to deal with this problem.

B. DISTANCE FUNCTIONS AND TYPES OF COVARIATES

In a second step, the shared variables of the merged database will have
to be encoded according to the choice of the distance function fixed by
user, knowing that it is also frequent that it is the type of the
variables which fixes the distance function to choose. The function
[`transfo_dist`](https://otrecoding.github.io/OTrecod/reference/transfo_dist.md)
is available in the package to assist users in this task but a user can
also decide to make this preparation by themselves. Thus, with the
Euclidean or Manhattan distance ((3), `norm` = "E" or "M"), if all types
of variables are allowed, logical variables are transformed in binary
variables, and categorical variables (factors ordered or not) are
replaced by their related disjunctive tables (the function
[`transfo_quali`](https://otrecoding.github.io/OTrecod/reference/transfo_quali.md)
can make these specific transformations). The Hamming distance (`norm` =
"H") only requires binary variables (all other forms are not allowed).
In this context, continuous variables could have been converted in
factor of k levels (\\k\>2\\) beforehand. The categorical covariates are
then transformed in disjunctive tables (containing the (\\k-1\\)
corresponding binary variables) before use. With this distance,
categorical variables are also transformed in disjunctive tables. Notice
that, using the Hamming distance could be quite long in presence of NAs
on covariates. Finally, the Gower distance ((4), `norm` = "G") uses the
([`gower.dist`](https://rdrr.io/pkg/StatMatch/man/gower.dist.html))
function (5) and so allows logical, categorical and numeric variables
without preliminary transformations.

In conclusion, the structure of the data.frame required in input of the
function `proxim_dist` corresponds to two overlayed databases with two
target outcomes and a set of shared covariates whose encodings depend on
the distance function choosen by user.

If some columns are excluded when computing an Euclidean, Manhattan, or
Hamming distance between two rows, the sum is scaled up proportionally
to the number of columns used in the computation as proposed by the
standard ([`dist`](https://rdrr.io/r/stats/dist.html)) function. If all
pairs are excluded when computing a particular distance, instead of
putting NA in the corresponding cell of the distance matrix, the process
stops and an object listing the problematic rows is proposed in output.
It suggests users to remove these rows before running the process again
or impute NAs related to these rows (see (6) for more details).

C. PROFILES OF COVARIATES AND OUTPUT DETAILS

Whatever the type (mixed or not) and the number of covariates in the
data.frame of interest, the function `proxim_dist` firstly detects all
the possible profiles (or combinations) of covariates from the two
databases, and saves them in the output `profile`. For example, assuming
that a data.frame in input (composed of two overlayed data.frames A and
B) have three shared binary covariates (identically encoded in A and B)
so the sequences `011` and `101` will be considered as two distinct
profiles of covariates. If each covariate is a factor of \\n_1\\,
\\n_2\\ and \\n_3\\ levels respectively, so it exists at most \\n_1
\times n_2 \times n_3\\ possible profiles of covariates. This number is
considered as a maximum here because only the profiles of covariates met
in at least one of the two databases will be kept for the study.

`proxim_dist` classifies individuals from the two databases according to
their proximities to each profile of covariates and saves the
corresponding indexes of rows from A and B in two lists `indXA` and
`indXB` respectively. `indXA` and `indXB` thus contain as many objects
as covariates profiles and the proximity between a given profile and a
given individual is defined as follows. The function also provides in
output the list of all the encountered profiles of covariates. As a
decision rule, for a given profile of covariates \\P_j\\, an individual
\\i\\ will be considered as a neighbor of \\P_j\\ if \\dist(i,P_j) \<
prox \times max(dist(i,P_j))\\ where `prox` will be fixed by user. Set
the value 0 to the `prox` parameter assures that each individual of A
(and B respectively) is exactly the profile of one profile of
covariates. Therefore, it is not recommended in presence of continuous
coavariates. Conversely, assign the value 1 to `prox` is not recommended
because it assumes that each individual is neighbor with all the
encountered profiles of covariates.

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

4.  Gower, J. C. (1971). A general coefficient of similarity and some of
    its properties. Biometrics, 27, 623–637.

5.  D'Orazio M. (2015). Integration and imputation of survey data in R:
    the StatMatch package. Romanian Statistical Review, vol. 63(2)

6.  Borg, I. and Groenen, P. (1997) Modern Multidimensional Scaling.
    Theory and Applications. Springer.

## See also

[`transfo_dist`](https://otrecoding.github.io/OTrecod/reference/transfo_dist.md),
[`imput_cov`](https://otrecoding.github.io/OTrecod/reference/imput_cov.md),
[`merge_dbs`](https://otrecoding.github.io/OTrecod/reference/merge_dbs.md),
[`simu_data`](https://otrecoding.github.io/OTrecod/reference/simu_data.md)

## Author

Gregory Guernec, Valerie Gares, Jeremy Omer

<otrecod.pkg@gmail.com>

## Examples

``` r

data(simu_data)
### The covariates of the data are prepared according to the chosen distance
### using the transfo_dist function

### Ex 1: The Manhattan distance

man1 <- transfo_dist(simu_data,
  quanti = c(3, 8), nominal = c(1, 4:5, 7),
  ordinal = c(2, 6), logic = NULL, prep_choice = "M"
)
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 5 remaining levels
mat_man1 <- proxim_dist(man1, norm = "M") # man1 compatible with norm = "E" for Euclidean


### Ex 2: The Euclidean and Manhattan distance applied on coordinates from FAMD

eucl_famd <- transfo_dist(simu_data,
  quanti = c(3, 8), nominal = c(1, 4:5, 7),
  ordinal = c(2, 6), logic = NULL, prep_choice = "FAMD", info = 0.80
)
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 5 remaining levels
#> Warning: Presence of NA on covariates. You should work with complete or imputed DB before using this method.
#>               By default, only rows with no NA on covariates have been kept.
#> Only 528 rows (75%)are kept here corresponding to complete cases
mat_e_famd <- proxim_dist(eucl_famd, norm = "E")

# \donttest{
mat_m_famd <- proxim_dist(eucl_famd, norm = "M")
# }

### Ex 3: The Gower distance with mixed covariates

gow1 <- transfo_dist(simu_data[c(1:100, 301:400), ],
  quanti = c(3, 8), nominal = c(1, 4:5, 7),
  ordinal = c(2, 6), logic = NULL, prep_choice = "G"
)
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 5 remaining levels
mat_gow1 <- proxim_dist(gow1, norm = "G")

# \donttest{
### Ex 4a: The Hamming distance with binary (but incomplete) covariates only

# categorization of the continuous covariates age by tertiles
ham1 <- transfo_dist(simu_data,
  quanti = c(3, 8), nominal = c(1, 4:5, 7), ordinal = c(2, 6),
  convert_num = 8, convert_class = 3, prep_choice = "H"
)
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 5 remaining levels
mat_ham1 <- proxim_dist(ham1, norm = "H")
#> Hamming 1/3
#> Hamming 2/3
#> Hamming 3/3
# Be patient ... It could take few minutes

### Ex 4b: The Hamming distance with complete cases on nominal and ordinal covariates only
simu_data_CC <- simu_data[(!is.na(simu_data[, 5])) & (!is.na(simu_data[, 6])) &
  (!is.na(simu_data[, 7])), 1:7]
ham2 <- transfo_dist(simu_data_CC,
  quanti = 3, nominal = c(1, 4:5, 7), ordinal = c(2, 6),
  prep_choice = "H"
)
#> Your target DB[, "Z"] was numeric ... By default, it has been converted in factor of integers
#> 5 remaining levels
mat_ham2 <- proxim_dist(ham2, norm = "H")


### Ex 5: PARTICULAR CASE, If only one covariate with no NAs

man2 <- man1[, c(1:3, 7)] # Only Smoking variable
man2_nona <- man2[!is.na(man2[, 4]), ] # Keep complete case
mat_man2_nona <- proxim_dist(man2_nona, norm = "M", prox = 0.10)

mat_man2_nona_H <- proxim_dist(man2_nona, norm = "H") # Hamming


### Ex 6: PARTICULAR CASE, many covariates but NAs in distance matrix

# We generated NAs in the man1 object so that:
# dist(A4,B102) and dist(A122,B102) returns NA whatever the norm chosen:
man1b <- man1
man1b[4, 7:9] <- NA
man1b[122, 6:9] <- NA
man1b[300 + 102, 4:6] <- NA
mat_man3 <- proxim_dist(man1b, norm = "M")
#> Warning: THE PROCESS STOPPED
#> !!! Because of the presence of NAs in distance matrix, the process stopped. Combinations of rows of A and B with NAs cause pbs and have to be removed or imputed. To help you, the indexes of rows are listed in the returned object.
# The process stopped indicates 2 NAs and the corresponding row numbers
# The 2nd output of mat_man3 indicates that removing first the 102th row of the database
# B is enough to solve the pb:
man1c <- man1b[-402, ]
mat_man4 <- proxim_dist(man1c, norm = "M")
# }
```
