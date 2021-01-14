OTrecod package
================

# A package dedicated to data fusion

## Introduction

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/otrecoding/OTrecod.svg?branch=master)](https://travis-ci.org/otrecoding/OTrecod)
[![CRAN
status](https://www.r-pkg.org/badges/version/OTrecod)](https://cran.r-project.org/package=OTrecod)
[![codecov](https://codecov.io/gh/otrecoding/OTrecod/branch/master/graph/badge.svg)](https://codecov.io/gh/otrecoding/OTrecod)
[![Launch
binder](http://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/otrecoding/OTrecod/master)
[![Website](https://img.shields.io/website?url=https%3A%2F%2Fotrecoding.github.io%2FOTrecod%2F)](https://otrecoding.github.io/OTrecod/)
<!-- badges: end -->

The **OTrecod** package gives access to a set of original functions
dedicated to data fusion.

<p align="justify">

From two separate data sources with no overlapping units, sharing only a
set of common variables X and a same target information not jointly
observed in a same encoding from one data source to another (Y in A and
Z in B), the functions **OT\_outcome** and **OT\_joint** aim at
providing users a complete synthetic database where the missing
information is available for every unit.

</p>

<p align="justify">

This recoding problem is solved using the optimal transportation theory
which provides a map that transfers the joint distribution of the first
target variable and X to the joint distribution of the second one and X,
or inversely. Algorithms used in these two functions come from the
references (1) and (2).

</p>

 

## Package installation

If the package **OTrecod** is not installed in their current R versions,
users can install it by following the standard instruction:

``` r
install.packages("OTrecod")
```

Obviously, each time an R session is opened, the **OTrecod** library
must be loaded with:

``` r
library(OTrecod)
```

Moreover, the development version of **OTrecod** can be installed
actually from [GitHub](https://github.com/otrecoding/OTrecod) with:

``` r
# Install development version from GitHub
devtools::install_github("otrecoding/OTrecod")
```

 

## Database examples and expected structure before data fusion

<p align="justify">

The available databases called **tab\_test** and **simu\_data**
correspond to overlayed databases used as examples in the documentation
of all the functions. Their structures can help users understanding the
database structure expected as input argument of the functions
**OT\_outcome** and **OT\_joint**. The first rows of the two overlayed
data sources of **simu\_data** are visualized as follows to inform about
the expected database structure:

</p>

``` r
data(simu_data)
dim(simu_data)
[1] 700   8
simu_data[c(1:5,301:305),]
    DB     Yb1 Yb2 Gender Treatment Dosage Smoking      Age
1    A [40-60[  NA Female     Trt A  Dos 3     YES 65.44273
2    A [20-40]  NA   Male      <NA>  Dos 2      NO 51.78596
3    A [40-60[  NA Female   Placebo  Dos 2     YES 49.10844
4    A [40-60[  NA Female     Trt B  Dos 4    <NA> 56.43524
5    A [40-60[  NA Female     Trt A  Dos 4     YES 44.77365
301  B    <NA>   5 Female   Placebo  Dos 2     YES 44.58233
302  B    <NA>   1 Female     Trt B  Dos 4    <NA> 65.23921
303  B    <NA>   2 Female   Placebo   <NA>      NO 51.64228
304  B    <NA>   2 Female     Trt A   <NA>      NO 50.15125
305  B    <NA>   1 Female     Trt B  Dos 4     YES 61.53242
```

<p align="justify">

The first column called *DB* corresponds here to the database identifier
(two data sources called here 1 and 2 with the data source 1 placed
above the data source 2). The second column called *Yb1* is the target
variable of the data source 1. The values of *Yb1* in the data source 2
are missing and will be predicted using an optimal transportation
algorithm integrated in one of the two functions called **OT\_outcome**
and **OT\_joint**. In the same way, the variable *Yb2* (third column) is
the target variable of the data source 2 whose values in 1 are unknown.
These missing values can also be predicted using **OT\_outcome** and
**OT\_joint**.

</p>

<p align="justify">

The presence of these three variables is essential in any database
dedicated to datafusion in the **OTrecod** package whatevever their
names and whatever their orders in the database. The following columns
correspond to shared variables of any type, complete or not. Note that
continuous variables (like age in years) are not allowed with the
**OT\_joint** function.

</p>

<p align="justify">

Support functions are available in the package (**merge\_dbs**,
**imput\_cov**) to assist user in this preparation.

</p>

<p align="justify">

Finally, the supplementary datasets **api29** and **api35** are simple
datasets extracted from the API program
(<https://www.cde.ca.gov/re/pr/api.asp>) to allow users to practice with
convenient databases.

</p>

 

## Support functions

<p align="justify">

Among the available functions, the **OTrecod** package provides a set of
support functions to assist users in each step of their data fusion
projects.

</p>

### merge\_dbs

<p align="justify">

The **merge\_dbs** function is a pre-process data fusion function
dedicated to the harmonization of two data sources. By default,
variables (not target variables) with same labels are considered as
shared between the two databases. The **merge\_dbs** function detects
potential discrepancies between the variables before merging by:

  - firstly excluding variables with different labels from the first
    database to the second one and inversely.
  - excluding a priori shared variables with different types.
  - excluding a priori shared factors with different levels.

The actual form of the function does not propose automatic
reconciliation actions to reintroduce the problematic variables but
gives user enough information in output to do it by himself if
necessary. The call of the **merge\_dbs** function is
actually:

</p>

``` r
merge_dbs = function(DB1, DB2, row_ID1 = NULL, row_ID2 = NULL, NAME_Y, NAME_Z, order_levels_Y = levels(DB1[, NAME_Y]), order_levels_Z = levels(DB2[, NAME_Z]), ordinal_DB1 = NULL, ordinal_DB2 = NULL,
                     impute = "NO", R_MICE = 5, NCP_FAMD = 3, seed_func = sample(1:1000000, 1))
```

<p align="justify">

The **merge\_dbs** function notably provides in output an unique
database, result of the overlayed of the two initial data sources, in
the structure expected by the **OT\_outcome** and **OT\_joint**
functions.

</p>

 

### select\_pred

<p align="justify">

The **select\_pred** function is a pre-process data fusion function
dedicated to the selection of matching variables. This selection is
essential when the initial set of shared variables is important, but
also because the choice of predictors greatly influences the quality of
the data fusion whatever the optimal transportation algorithms chosen a
posteriori.

The call of the **select\_pred** function is
actually:

</p>

``` r
select_pred = function(databa,Y = NULL, Z = NULL, ID = 1, OUT = "Y", quanti = NULL, nominal = NULL, ordinal = NULL, logic = NULL,
                       convert_num = NULL, convert_clss = NULL, thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
                       RF = TRUE, RF_ntree = 500, RF_condi = FALSE, RF_condi_thr = 0.20, RF_SEED = sample(1:1000000, 1))
```

 

### verif\_OT

<p align="justify">

The **verif\_OT** function is a post-process data fusion function
dedicated to the validation of the fusion. The function provides a set
of tools to assess the quality of the optimal transportation recoding
proposed by the algorithms to predict the missing information of the
target variables in one or both datasources.

</p>

The call of the **verif\_OT** function is
actually:

``` r
verif_OT = function(ot_out, group.clss = FALSE, ordinal = TRUE, stab.prob = FALSE, min.neigb = 1, R = 10, seed.stab = sample(1:1000000, 1))
```

 

## Optimal transportation functions

<p align="justify">

The **OTrecod** package provides two algorithms that use optimal
transportation theory to solve recoding problems in data fusion contexts
(see (1) and (2) for more details). Each algorithm is stored in one
function and each function provides in output a unique and synthetic
database where the two initial data sources are overlayed and the
missing information from only one or both target variables are fully
completed.

Each of the two alogorithms also proposed enrichments by relaxing the
initial distributional constraints and adding regularization terms as
described in (2).

</p>

 

### OT\_outcome

<p align="justify">

The **OT\_outcome** function can provide individual predictions of the
incomplete target variables by considering the recoding problem
involving only optimal transportation of outcomes (see (1) and (2) for
more details).

The call of the **OT\_outcome** function
is:

</p>

``` r
OT_outcome = function(datab, index_DB_Y_Z = 1:3, quanti = NULL, nominal = NULL, ordinal = NULL,logic = NULL,
                      convert.num = NULL, convert.clss = NULL, FAMD.coord = "NO", FAMD.perc = 0.8,
                      dist.choice = "E", percent.knn = 1, maxrelax = 0, indiv.method = "sequential",
                      prox.dist = 0.30, solvR = "glpk", which.DB = "BOTH")
```

 

### OT\_joint

<p align="justify">

The **OT\_joint** function can provide individual predictions of the
incomplete target variables by considering the recoding problem
involving optimal transportation of shared variables and outcomes
(see(2) for more details).

The call of the **OT\_joint** function
is:

</p>

``` r
OT_joint = function(datab, index_DB_Y_Z = 1:3, nominal = NULL, ordinal = NULL,logic = NULL,
                    convert.num = NULL, convert.clss = NULL, dist.choice = "E", percent.knn = 1,
                    maxrelax = 0, lambda.reg = 0.0, prox.X = 0.10, solvR = "glpk", which.DB = "BOTH")
```

 

## References

1)  Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy
    N (2019). On the use of optimal transportation theory to recode
    variables and application to database merging. The International
    Journal of Biostatistics.Volume 16, Issue 1, 20180106, eISSN
    1557-4679.

2)  Gares V, Omer J (2020). Regularized optimal transport of covariates
    and outcomes in data recoding. Journal of the American Statistical
    Association.
