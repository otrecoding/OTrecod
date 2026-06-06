# an-application-of-the-OTrecod-package

## Package installation

If the package **OTrecod** is not installed in their current R versions,
users can install it by following the standard instruction:

``` r

install.packages("OTrecod")
```

Each time an R session is opened, the **OTrecod** library must be loaded
with:

``` r

library(OTrecod)
```

Moreover, the development version of **OTrecod** can be installed
actually from [GitHub](https://github.com/otrecoding/OTrecod) with:

``` r

# Install development version from GitHub
devtools::install_github("otrecoding/OTrecod")
```

 

## I. The context

This vignette illustrates how to use the tools contained in the
**OTrecod** package to solve a variable recoding problem frequently
encountered in the context of data fusion. For more details about the
theory of the algorithms used in the functions **OT_outcome** and
**OT_joint** of the package, user can consult (1),(2) and the
documentation linked to each function.

For this example, we have transformed an available database called
*samp.A* from the **StatMatch** package (see (3) and the help of
*samp.A* for more details). The *samp.A* database provides a limited
number of variables observed at persons levels among those usually
collected in the *European Union Statistics on Income and Living
Conditions Survey* (the EU–SILC survey). In this database, the variable
*c.neti* is a factor that corresponds to the persons net income of
thousand of Euros categorized in 7 ordered classes. The raw distribution
of *c.neti* is done in the following results:

``` r

library(StatMatch)
data(samp.A)
```

``` r

dim(samp.A)
#> [1] 3009   13
head(samp.A)
#>        HH.P.id area5 urb hsize hsize5 age    c.age sex marital edu7 n.income
#> 21384 10149.01    NE   1     1      1  85 (64,104]   2       3    3     1677
#> 35973 17154.02    NE   1     2      2  78 (64,104]   1       2    3    13520
#> 11774  5628.01    NO   2     1      1  48  (44,54]   1       3    3    20000
#> 32127 15319.01     S   1     2      2  78 (64,104]   1       2    1    12428
#> 6301   2973.05     S   2     5    >=5  17  [16,34]   1       1    1        0
#> 12990  6206.02     C   2     2      2  28  [16,34]   2       2    5        0
#>         c.neti        ww
#> 21384   (0,10] 3591.8939
#> 35973  (10,15]  415.1592
#> 11774  (15,20] 2735.4029
#> 32127  (10,15] 1239.5086
#> 6301  (-Inf,0] 5362.7588
#> 12990 (-Inf,0] 2077.7137
table(samp.A$c.neti)   # Repartition of c.neti in the sample
#> 
#>  (-Inf,0]    (0,10]   (10,15]   (15,20]   (20,25]   (25,35] (35, Inf] 
#>       564       671       541       500       307       272       154
```

To construct a standard recoding problem in data fusion, the *samp.A*
database has been transformed as follows:

- From the *c.neti* variable, we built a variable called *c.neti.bis*
  that corresponds to a second encoding of the persons net income now
  coded in 4 ordered levels (1: by grouping the first two levels of
  *c.neti*, 2: by grouping the two subsequent levels of *c.neti*, 3: by
  grouping the last two levels).
- To improve the running time of the example, we only kept the first 350
  subjects (rows) of *samp.A*.
- This new sample has been randomly divided in two distinct subsamples
  called *data1* (200 rows: subjects 1 to 200) and *data2* (150 rows:
  sujects 201 to 350).
- The *c.neti.bis* variable has been excluded from *data1* and
  respectively *c.neti* have been excluded from *data2*.
- Among the other variables, a random subset of them have been assigned
  in each database.

The R code related to these transformations was:

``` r

c.neti            = as.numeric(samp.A$c.neti)
samp.A$c.neti.bis = as.factor(ifelse(c.neti %in% c(1,2),1, 
                              ifelse(c.neti %in% c(3,4),2, 
                              ifelse(c.neti %in% c(5,6),3,4)))) 
data1 = samp.A[1:200,c(2:3,5,7:9,12:13)]
colnames(data1)[4] = "age" 
head(data1)
#>       area5 urb hsize5      age sex marital   c.neti        ww
#> 21384    NE   1      1 (64,104]   2       3   (0,10] 3591.8939
#> 35973    NE   1      2 (64,104]   1       2  (10,15]  415.1592
#> 11774    NO   2      1  (44,54]   1       3  (15,20] 2735.4029
#> 32127     S   1      2 (64,104]   1       2  (10,15] 1239.5086
#> 6301      S   2    >=5  [16,34]   1       1 (-Inf,0] 5362.7588
#> 12990     C   2      2  [16,34]   2       2 (-Inf,0] 2077.7137
data2 = samp.A[201:350,c(3,5:6,8:11,13:14)]
head(data2)
#>       urb hsize5 age sex marital edu7 n.income       ww c.neti.bis
#> 39565   3      2  81   2       2    1     6448 1129.274          1
#> 36490   1    >=5  48   1       3    2    16423 3082.331          2
#> 27529   2      2  66   1       2    3    15600 2433.020          2
#> 201     2      4  26   1       1    3    12876 1869.286          2
#> 31375   2      2  53   1       2    3    25633 2125.361          3
#> 3226    3      1  45   2       3    3     2611 1855.230          1
```

In conclusion, after transformation, we had:

- two separate databases *data1* and *data2* with no common subjet.
- *c.neti* and *c.neti.bis*: two variables that summarize a same
  information (the persons net income) encoded in two different scales
  in *data1* and *data2*.
- a subset of shared variables between the two databases (variables with
  same labels in the 2 databases).
- subsets of specific variables

``` r

table(data1$c.neti)        # 7 levels in data1
#> 
#>  (-Inf,0]    (0,10]   (10,15]   (15,20]   (20,25]   (25,35] (35, Inf] 
#>        37        40        34        49        18        13         9

table(data2$c.neti.bis)    # 4 levels in data2
#> 
#>  1  2  3  4 
#> 60 57 26  7

colnames(data1)
#> [1] "area5"   "urb"     "hsize5"  "age"     "sex"     "marital" "c.neti" 
#> [8] "ww"

colnames(data2)
#> [1] "urb"        "hsize5"     "age"        "sex"        "marital"   
#> [6] "edu7"       "n.income"   "ww"         "c.neti.bis"

intersect(colnames(data1), colnames(data2))   # the susbet of a priori shared variables
#> [1] "urb"     "hsize5"  "age"     "sex"     "marital" "ww"
```

**OBJECTIVE**

**Assuming that the encoding of *c.neti* is unknown for the subjects of
*data2* and that the encoding of *c.neti.bis* is unknown for the
subjects of *data1*, the functions of the OTrecod package solve this
recoding problem by predicting the missing scale of the persons net
income in one or the two databases. This solution allows the user to
fusion his two databases and finally works with a bigger, unique and
synthetic dataset.**

For the rest of the study, *c.neti* and *c.neti.bis* are called the
target variables of *data1* and *data2* respectively. we deliberately
limit this example to the prediction of the variable *c.neti.bis* in
*data1* but note that the proposed approach would be the same for the
prediction of *c.neti* in *data1*.

 

## II. Harmonization of the data sources

Knowing the objective of the study, we first prepare the 2 databases to
data fusion. The two functions dedicated to this data fusion in the
**OTrecod** package expect a specific structure of database as argument.
The **merge_dbs** function assists user in this task by:

- overlaying the two databases
- detecting the false shared variables
- imputing incomplete information of shared variables if desired.

For this, we fill in the arguments of the **merge_dbs** function by
declaring as follows, the column indexes of all the ordinal variables in
*data1* and *data2* (including the target variables indexes). Here there
is no need to handle missing information because all the shared
variables are complete, but note that it exists a specific argument
(*impute*) in the function to take them into account when necessary.

``` r

db_test  = merge_dbs(data1, data2, NAME_Y = "c.neti", NAME_Z = "c.neti.bis",
                     ordinal_DB1 = c(2,3,4,7), ordinal_DB2 = c(1:2,6,9))
#> DBS MERGING in progress. Please wait ...
#> DBS MERGING OK
#> -----------------------
#> 
#> SUMMARY OF DBS MERGING:
#> Nb of removed subjects due to NA on targets: 0(0%)
#> Nb of removed covariates due to differences between the 2 bases: 1
#> Nb of remained covariates: 5
#> Imputation on incomplete covariates: NO

summary(db_test)
#>               Length Class      Mode     
#> DB_READY      8      data.frame list     
#> ID1_drop      0      -none-     character
#> ID2_drop      0      -none-     character
#> Y_LEVELS      7      -none-     character
#> Z_LEVELS      4      -none-     character
#> REMOVE1       1      -none-     character
#> REMOVE2       0      -none-     NULL     
#> REMAINING_VAR 5      -none-     character
#> IMPUTE_TYPE   1      -none-     character
#> DB1_raw       8      data.frame list     
#> DB2_raw       9      data.frame list     
#> SEED          1      -none-     numeric

db_test$REMAINING_VAR
#> [1] "hsize5"  "marital" "sex"     "urb"     "ww"

db_test$REMOVE1
#> [1] "age"

db_test$REMOVE2
#> NULL

db_test$ID1_drop; db_test$ID2_drop
#> character(0)
#> character(0)

db_test$DB_READY[c(1:5,201:205),]   # The 5 1st subjects of the two databases 
#>       DB        Y    Z hsize5 marital sex urb        ww
#> 21384  1   (0,10] <NA>      1       3   2   1 3591.8939
#> 35973  1  (10,15] <NA>      2       2   1   1  415.1592
#> 11774  1  (15,20] <NA>      1       3   1   2 2735.4029
#> 32127  1  (10,15] <NA>      2       2   1   1 1239.5086
#> 6301   1 (-Inf,0] <NA>    >=5       1   1   2 5362.7588
#> 39565  2     <NA>    1      2       2   2   3 1129.2739
#> 36490  2     <NA>    2    >=5       3   1   1 3082.3314
#> 27529  2     <NA>    2      2       2   1   2 2433.0201
#> 201    2     <NA>    2      4       1   1   2 1869.2859
#> 31375  2     <NA>    3      2       2   1   2 2125.3614
```

In output:

- The *REMAINING_VAR* object gives the identity of the shared variables
  kept for the merging of the two databases.
- The *REMOVE1* object gives the identity of potential shared variables
  dropped because of dicrepancies of types: the variable *age* stored as
  factor in *data1* and as numeric in *data2* is in this case. Sometimes
  it is possible to reconciliate the two encodings a posteriori,
  sometimes not. Here it would be possible to categorize the variable
  *age* in *data2* in the same way as *age* in *data1* but we decide for
  the example to discard it for the rest of the study.
- The *REMOVE2* object is null, and this result means that there is no
  factor dropped because of discrepancies of levels.
- The *ID1_drop* and *ID2_drop* objects are null which means that the
  target variables *c.neti* and *c.neti.bis* have no missing values in
  there respective databases.
- Finally the *DB_READY* object provides users a synthetic database that
  overlays *data1* and *data2* (in this order) where the *DB* variable
  is the database identifier.

The variable *Y* now corresponds to the target variable *c.neti* with
its specific encoding in *data1* and missing values in *data2*. In the
same way, the variable *Z* corresponds to *c.neti.bis*. This database
have the structure expected as argument of the **OT_outcome** and
**OT_joint** functions.

 

## III. Selection of the matching variables

Among the set of shared variables kept in the output *DB.READY* database
of the **merge_dbs**, it is important to discard those that never appear
not good predictors of the persons net income whatever the considered
database. The subset of remaining variables will be the matching
variables. This selection of matching variables can be done using the
**select_pred** function of the package.

This function proposes two levels of study to conclude:

- When the sample of shared variables is small, standard correlation
  studies are often enough to select the best set of predictors and get
  rid of potential multicolinearities between the candidates.
- When the conclusion of the first part appears not so clear, this study
  can be completed by a random forest (RF) process of selection (see
  *RF* argument). This part is notably particularly convenient when the
  sample of covariates is large. The function proposes two random forest
  procedures: The standard RF which estimates the permutation importance
  criterion of each shared variable and the conditional importance
  measures from the *proxy* package ((4)).

At this step, the choice of the recoding algorithm (see (2)) to use a
posteriori must be anticipated in the handling of the potential numeric
shared variables (In this example, the *ww* variable is the only one
concerned. Indeed, if the **OT_outcome** function runs whatever the type
considered, this is not the case of the actual version of the
**OT_joint** version which not allows the numeric variable. For this
vignette we will test the two algorithms, so we decide to transform *ww*
as categorial using the quartiles of its distribution as thresholds of
the four classes. This transformation can be dealt by standard R
function (like *cut()*) or directly in the functions using the
*convert_num* and *convert_class* arguments dedicated to this task.

The selection of the best predictors of the persons net income must be
done in the two databases separately but the **selec_pred** functions
allowed overlayed databases as arguments using the following code:

``` r

# for data1
test_DB1 = select_pred(db_test$DB_READY,Y = "Y", Z = "Z", ID = 1, OUT = "Y",
                       quanti = 8, nominal = c(1,5:6,7), ordinal = c(2:4),
                       convert_num = 8, convert_class = 4,
                       thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.10,
                       RF = TRUE, RF_SEED = 3017)
#> The select_pred function is running for outcome= Y. Please wait ...
#> Risk of collinearity between predictors detected: Some predictors will be dropped during RF process
#> The process is now successfully completed
#> ---------
#> For comparison with another outcome from two overlayed tables  : 
#> just adapt the OUT option keeping all the others unchanged in the function
#> ---
#> For comparison with another outcome from two unoverlayed tables:
#> just adapt the arguments from Y to convert_class
#> ---------

# for data2
test_DB2 = select_pred(db_test$DB_READY,Y = "Y", Z = "Z", ID = 1, OUT = "Z",
                       quanti = 8, nominal = c(1,5:6,7), ordinal = c(2:4),
                       convert_num = 8, convert_class = 4,
                       thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.10,
                       RF = TRUE, RF_SEED = 3017)
#> The select_pred function is running for outcome= Z. Please wait ...
#> Risk of collinearity between predictors detected: Some predictors will be dropped during RF process
#> The process is now successfully completed
#> ---------
#> For comparison with another outcome from two overlayed tables  : 
#> just adapt the OUT option keeping all the others unchanged in the function
#> ---
#> For comparison with another outcome from two unoverlayed tables:
#> just adapt the arguments from Y to convert_class
#> ---------
```

As input:

- The *OUT* argument specifies the target variable to predict.
- The *thresh_cat* argument corresponds to the threshold of the V Cramer
  coefficient (for the categorical variables) .
- The *thresh_num* argument corresponds to the threshold of the Spearman
  correlation coefficient (for the continuous and ordinal variables).
  Under these values, two shared variables present acceptable
  colinearities. Under these values, a shared variable is not considered
  as a good predictor of a target variable.
- The *thresh_Y* argument is reserved to random forest procedures. It
  corresponds to an acceptability threshold related to the importance
  mesure criterions. Variables with a cumulative percent of importance
  measures less than this threshold will be dropped from the final list
  of RF predictors.

Notice that it is important to keep the same arguments in the two
selections (*test_DB1* and *test_DB2*) for an optimal comparability.
Here the standard RF process is used, and the *ww* variable is converted
in categorical type before the selection using *convert_num*.

Let see the result for *data1*:

``` r

summary(test_DB1)
#>               Length Class      Mode     
#> seed          1      -none-     numeric  
#> outc          1      -none-     character
#> thresh        3      -none-     numeric  
#> convert_num   1      -none-     character
#> DB_USED       8      data.frame list     
#> vcrm_OUTC_cat 5      data.frame list     
#> cor_OUTC_num  5      data.frame list     
#> vcrm_X_cat    5      data.frame list     
#> cor_X_num     5      data.frame list     
#> FG_test       3      -none-     numeric  
#> collinear_PB  2      -none-     list     
#> drop_var      1      -none-     character
#> RF_PRED       4      -none-     numeric  
#> RF_best       2      -none-     character

test_DB1$vcrm_OUTC_cat
#>   name1   name2 V_Cramer CorrV_Cramer   N
#> 4     Y     sex   0.4384       0.4036 200
#> 3     Y marital   0.2452       0.1740 200
#> 5     Y     urb   0.1703       0.0000 200
#> 2     Y  hsize5   0.1634       0.0000 200
#> 6     Y      ww   0.1465       0.0000 200

test_DB1$collinear_PB
#> $VCRAM
#>     name1   name2 V_Cramer CorrV_Cramer   N
#> 21 hsize5 marital   0.5038       0.4858 200
#> 
#> $SPEARM
#> [1] name1       name2       RANK_COR    pv_COR_test N          
#> <0 rows> (or 0-length row.names)

# Results from RF
test_DB1$drop_var
#> [1] "marital"

test_DB1$RF_PRED
#>     sex  hsize5     urb      ww 
#> 86.9697  8.1313  4.8990  0.0000
```

In the set of shared variables, all variables are now categorical
(factors or ordered factors). According to the *vcrm_OUTC_cat* output
object, the best predictors of *c.neti* are: *sex*, *marital*, *urb*,
*hsize5*, and *ww* in that order. Nevertheless, a risk of collinearity
is detected between *marital* and *hsize5*. The RF process finally
suggests to drop *marital* and to keep only *sex*, *hsize5* and
eventually *urb* as matching variables.

Let see the result for *data2*:

``` r

summary(test_DB2)
#>               Length Class      Mode     
#> seed          1      -none-     numeric  
#> outc          1      -none-     character
#> thresh        3      -none-     numeric  
#> convert_num   1      -none-     character
#> DB_USED       8      data.frame list     
#> vcrm_OUTC_cat 5      data.frame list     
#> cor_OUTC_num  5      data.frame list     
#> vcrm_X_cat    5      data.frame list     
#> cor_X_num     5      data.frame list     
#> FG_test       3      -none-     numeric  
#> collinear_PB  2      -none-     list     
#> drop_var      1      -none-     character
#> RF_PRED       4      -none-     numeric  
#> RF_best       2      -none-     character

test_DB2$vcrm_OUTC_cat
#>   name1   name2 V_Cramer CorrV_Cramer   N
#> 4     Z     sex   0.4783       0.4583 150
#> 2     Z  hsize5   0.2343       0.1693 150
#> 3     Z marital   0.1828       0.1160 150
#> 5     Z     urb   0.1386       0.0000 150
#> 6     Z      ww   0.0949       0.0000 150

test_DB2$collinear_PB
#> $VCRAM
#>     name1   name2 V_Cramer CorrV_Cramer   N
#> 21 hsize5 marital   0.3555       0.3177 150
#> 
#> $SPEARM
#> [1] name1       name2       RANK_COR    pv_COR_test N          
#> <0 rows> (or 0-length row.names)

# Results from RF
test_DB2$drop_var
#> [1] "marital"

test_DB2$RF_PRED
#>     sex  hsize5      ww     urb 
#> 86.1852 13.1555  0.6593  0.0000
```

According to the *vcrm_OUTC_cat* output object, the best predictors of
*c.neti.bis* are: *sex*, *hsize5*, *marital*, *urb*, and *ww* in that
order. A risk of collinearity is also detected between *marital* and
*hsize5*. The RF process here suggests to drop *marital* and to keep
only *sex*, *hsize5* as matching variables.

**In conclusion:**

- The variables *sex* and *hsize5* are the best predictors of the target
  variables in the two databases.
- The variable *marital* is dropped whatever the database because of
  risks of multicollinearity with *hsize5*.
- The variable *ww* is never a good predictor of the target information.
- The variable *urb* is an acceptable predictor in the first database
  only.

The matching variable groups can be:

- *sex*, *hsize5*
- *sex*, *hsize5*, and *urb*

**We finally keep the last group and dropped the other ones for the rest
of the example.**

``` r

match_var = db_test$DB_READY[,-c(5,8)]
match_var[c(1:5,201:205),]
#>       DB        Y    Z hsize5 sex urb
#> 21384  1   (0,10] <NA>      1   2   1
#> 35973  1  (10,15] <NA>      2   1   1
#> 11774  1  (15,20] <NA>      1   1   2
#> 32127  1  (10,15] <NA>      2   1   1
#> 6301   1 (-Inf,0] <NA>    >=5   1   2
#> 39565  2     <NA>    1      2   2   3
#> 36490  2     <NA>    2    >=5   1   1
#> 27529  2     <NA>    2      2   1   2
#> 201    2     <NA>    2      4   1   2
#> 31375  2     <NA>    3      2   1   2
```

**This overlayed database is now ready for recoding.**

 

## IV. Predicting the missing scales in the databases

The **OTrecod** package proposes actually two algorithms using optimal
transportation theory (see (3) for details) to solve the recoding
problem previously introduced. Each algorithm is stored in a unique
function called **OT_outcome** and **OT_joint**. These two functions can
predict the missing values of *c.neti* in *data2*, the missing values of
*c.neti.bis* in *data1* or the both using a same argument called
*which.DB*.

As with the **select_pred** function, it is possible to transform
directly the continuous matching variables if necessary using the
*convert_num* and *convert_class* arguments.

Let see the R approach for the prediction of *c.neti.bis* in *data1*.

 

### A) Transporting target variables to predict the missing scales

The algorithm from **OT_outcome** function solves an optimization
problem to transfer the distributions of the target variables (or
outcome) from one database to another. Using this result, the initial
version of the algorithm (see (2)) transfers the distribution of
*c.neti.bis* in *data2* to the distribution of *c.neti.bis* in *data1*,
and (inversely for *c.neti* if necessary). The algorithm executes in
another distinct step, a nearest neighbor procedure to affect the
indidividual predictions of *c.neti.bis* in *data1*. This version of the
algorithm is actually available by writing *sequential* as argument
*method* of the **OT_outcome** function.

The algorithm assumes the strong hypothesis that the variable
*c.neti.bis* has the same distribution in *data1* and *data2* (and so on
for the variable *c.neti*). If in this example, by construction, this
hypothesis is obviously verified, there are also several contexts where
this latter no longer holds.

Enrichments have been thus proposed (see (2)) to relax this hypothesis
via the *maxrelax* argument of **OT_outcome** and directly provides the
individual predictions without using the nearest neighbor process. This
algorithm is actually available by assigning the *method* argument to
*optimal* in input of the **OT_outcome** function. For our example, the
corresponding R code for the prediction of *c.neti.bis* in *data1* is:

``` r

# sequential algorithm
mod1_seq = OT_outcome(match_var, nominal = c(1,5:6), ordinal = 2:4, dist.choice = "E",
                      maxrelax = 0 , indiv.method = "sequential", which.DB = "A")
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = OUTCOME
#> Distance                 = Euclidean
#> Percent closest knn      = 100%
#> Relaxation parameter     = NO
#> Relaxation value         = 0
#> Individual pred process  = Sequential
#> DB imputed               = A
#> ---------------------------------------
summary(mod1_seq)
#>             Length Class      Mode   
#> time_exe      1    difftime   numeric
#> gamma_A      28    -none-     numeric
#> gamma_B      28    -none-     numeric
#> profile       5    data.frame list   
#> res_prox     16    -none-     list   
#> estimatorZA 812    -none-     numeric
#> estimatorYB   0    -none-     NULL   
#> DATA1_OT      8    data.frame list   
#> DATA2_OT      7    data.frame list

# optimal algorithm with no relaxation term
mod2_opt = OT_outcome(match_var, nominal = c(1,5:6), ordinal = 2:4, dist.choice = "E",
                      maxrelax = 0, indiv.method = "optimal", which.DB = "A")
#> ---------------------------------------
#> OT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                     = R-OUTCOME
#> Distance                 = Euclidean
#> Percent closest knn      = 100%
#> Relaxation parameter     = NO
#> Relaxation value         = 0
#> Individual pred process  = Optimal
#> DB imputed               = A
#> ---------------------------------------
head(mod2_opt$profile)
#>        ID sex_2 urb_2 urb_3 hsize5
#> 21384 P_1     1     0     0      1
#> 35973 P_2     0     0     0      2
#> 11774 P_3     0     1     0      1
#> 6301  P_4     0     1     0      5
#> 12990 P_5     1     1     0      2
#> 39835 P_6     1     1     0      4

dim(mod2_opt$profile)
#> [1] 29  5

mod2_opt$gamma_A
#>               1     2          3          4
#> (-Inf,0]  0.185 0.000 0.00000000 0.00000000
#> (0,10]    0.170 0.030 0.00000000 0.00000000
#> (10,15]   0.000 0.105 0.06500000 0.00000000
#> (15,20]   0.000 0.245 0.00000000 0.00000000
#> (20,25]   0.000 0.000 0.04333333 0.04666667
#> (25,35]   0.000 0.000 0.06500000 0.00000000
#> (35, Inf] 0.045 0.000 0.00000000 0.00000000

head(mod2_opt$DATA1_OT)
#>       DB        Y    Z sex_2 urb_2 urb_3 hsize5 OTpred
#> 21384  1   (0,10] <NA>     1     0     0      1      1
#> 35973  1  (10,15] <NA>     0     0     0      2      3
#> 11774  1  (15,20] <NA>     0     1     0      1      2
#> 32127  1  (10,15] <NA>     0     0     0      2      3
#> 6301   1 (-Inf,0] <NA>     0     1     0      5      1
#> 12990  1 (-Inf,0] <NA>     1     1     0      2      1
```

As input:

- The target *c.neti.bis* is here considered as an ordinal scale but
  sometimes this information can be unknown: In the situation where
  there is no information about the encoding, is recommended to consider
  it as nominal.
- There is many possible distance functions integrated in **OT_outcome**
  to evaluate the distances between individuals: Here, the Euclidean
  distance have been chosen.
- The *percent.knn* argument corresponds to the fixed part (from 0 to 1)
  of individuals that participate to the computations of average
  distances between levels of target variables.

In output, these two models provides lists of same structure:

- *gamma_A* stores the solution of the optimization problem to determine
  the missing scales of *c.neti.bis* in *data1*.
- The *profile* object gives the profiles of covariates detected by the
  algorithm from the two databases: Here we have 29 profiles.
- *estimatorZA* is an array that provides for each profile, the
  probability for *c.neti.bis* to belong to a level given the covariates
  and the corresponding level of *c.neti*.
- The *DATA1_OT* object provides the individual predictions of
  *c.neti.bis* in *data1* (variable *OTpred* of the data.frame).

 

### B) Transporting target and shared variables to predict the missing scales

The objective of the algorithm integrated in the **OT_joint** function
is the same as that of the **OT_outcome** function, but unlike
**OT_outcome** transfers only the distributions of the target variables
to solve the optimization problem, the **OT_joint** function now
transports the joint distribution of outcomes and covariates: This
approach makes it possible to no longer assumes the strong
distributional hypothesis mentionned previously.

Nevertheless, as with the **OT_outcome** function, enrichments have been
proposed to relax distributional constraints (see *maxrelax* argument
and (2)) and eventually adds to the algorithm a regularization term to
smooth the variations of target variables with respect to matching
variables (see *lambda.reg* arguments).

To guarantee the convergence of the algorithm:

- It is suggested to work with a limited number of matching variables
  and so to select them rigorously using *select_pred* or another
  process of selection.
- Choosing a small value of the *prox.X* argument improves the
  convergence and the running time of the function.

For our example, the corresponding R code for the prediction of
*c.neti.bis* in *data1* is:

``` r

# Algorithms with no enrichments
mod3_joint = OT_joint(match_var, nominal = c(1,5), ordinal = c(2:4,6), dist.choice = "E", 
                      prox.X = 0.10, which.DB = "A")
#> ---------------------------------------
#> OT JOINT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                  = JOINT
#> Distance              = Euclidean
#> Percent closest       = 100%
#> Relaxation term       = 0
#> Regularization term   = 0
#> Aggregation tol cov   = 0.1
#> DB imputed            = A
#> ---------------------------------------
summary(mod3_joint)
#>             Length Class      Mode   
#> time_exe      1    difftime   numeric
#> gamma_A      28    -none-     numeric
#> gamma_B       0    -none-     NULL   
#> profile       4    data.frame list   
#> res_prox     16    -none-     list   
#> estimatorZA 812    -none-     numeric
#> estimatorYB   0    -none-     NULL   
#> DATA1_OT      7    data.frame list   
#> DATA2_OT      6    data.frame list
```

For a better understanding, the input arguments and output objects of
the **OT_joint** function have been thought to be similar to those
proposed by the **OT_outcome** function. As previously, the *OTpred*
variable of the *DATA1_OT* object stores the individual predictions of
*c.neti.bis* in *data1*.

 

## V. Validation of the individual predictions

The **verif_OT** function gives access to different tools for assessing
the reliability of the individual predictions proposed by the
**OT_outcome** and **OT_joint** functions.

``` r

# Validation of the mod1_seq model
verif_out1 = verif_OT(mod1_seq, stab.prob = TRUE, min.neigb = 3)
verif_out1$conf.mat
#>            predZ
#> Y             1   2   3   4 Sum
#>   (-Inf,0]   37   0   0   0  37
#>   (0,10]     34   6   0   0  40
#>   (10,15]     1  21  12   0  34
#>   (15,20]     1  48   0   0  49
#>   (20,25]     1   1   9   7  18
#>   (25,35]     1   1  11   0  13
#>   (35, Inf]   9   0   0   0   9
#>   Sum        84  77  32   7 200

verif_out1$res.prox
#>        N   V_cram rank_cor 
#>  200.000    0.740    0.645

verif_out1$res.stab
#>          N min.N mean    sd
#> 1st DB 112     3 0.97 0.119


# Validation of the mod2_seq model
verif_out2 = verif_OT(mod2_opt, stab.prob = TRUE, min.neigb = 3)
verif_out2$conf.mat
#>            predZ
#> Y             1   2   3   4 Sum
#>   (-Inf,0]   37   0   0   0  37
#>   (0,10]     34   6   0   0  40
#>   (10,15]     0  21  13   0  34
#>   (15,20]     0  49   0   0  49
#>   (20,25]     0   0   9   9  18
#>   (25,35]     0   0  13   0  13
#>   (35, Inf]   9   0   0   0   9
#>   Sum        80  76  35   9 200
rate_good_pred = (37+40+31+45+18+13+9)/200
rate_good_pred
#> [1] 0.965

verif_out2$res.prox
#>        N   V_cram rank_cor 
#>  200.000    0.800    0.688

verif_out2$res.stab
#>          N min.N mean sd
#> 1st DB 112     3    1  0


# Validation of the mod3_opt model
verif_jt   = verif_OT(mod3_joint, stab.prob = TRUE, min.neigb = 3)
verif_jt$conf.mat
#>            predZ
#> Y             1   2   3   4 Sum
#>   (-Inf,0]   33   4   0   0  37
#>   (0,10]     32   8   0   0  40
#>   (10,15]     4  29   1   0  34
#>   (15,20]    11  19  19   0  49
#>   (20,25]     0  10   5   3  18
#>   (25,35]     1   2  10   0  13
#>   (35, Inf]   6   1   0   2   9
#>   Sum        87  73  35   5 200

verif_jt$res.prox
#>        N   V_cram rank_cor 
#>  200.000    0.550    0.613

verif_jt$res.stab
#>          N min.N  mean    sd
#> 1st DB 112     3 0.869 0.175
```

As input:

- When the *stab.prob* argument is set to *TRUE*, the function provides
  information about the average stability of the individual predictions.
- The *min.neigb* argument specifies that individuals whose
  probabilities of assignments have been estimated using less than 3
  neighbors are dropped from this average (see the documentation of the
  function for more details).

As output, we have selected the most interesting one to make the
comparison between the three tested models:

- The *conf.mat* object provides the confusion matrix between *c.neti*
  and the individual predictions of *c.neti.bis* in *data1*. Here the
  two target variables *c.neti* and *c.neti.bis* are supposed to
  summarize a same information, so if the two scales are ordered, a
  clear structure in the confusion matrix traduces a good logical in the
  prediction. In our example, by construction, we also know the real
  encoding of *c.neti.bis*: We can so easily verify that the rate of
  good predictions equals 0.965 for the *mod2_opt* model, 0.93 for the
  *mod1_seq* model and 0.635 for the *mod3_jt* model.
- Studying the proximity between the distributions of *c.neti* and
  *c.neti.bis* in *data1* via the V Cramer criterion (see *res.prox*
  object) confirms that the *mod2_opt* model is the more adapted for our
  example.
- Finally, studying the stability of the predictions (see *res.stab*
  object) shows good results whatever the considered model.

**CONCLUSION**

**The two first models are here more adapted when the strong
distributional hypothesis previously introduced holds.** **Nevertheless,
the mod3_jt model could be significantly improved by adding appropriate
relaxation and regularization terms, so now, R user,** **test it by
yourself !**

## References

1.  Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy
    N (2019). On the use of optimal transportation theory to recode
    variables and application to database merging. The International
    Journal of Biostatistics.Volume 16, Issue 1, 20180106, eISSN
    1557-4679.

2.  Gares V, Omer J (2020). Regularized optimal transport of covariates
    and outcomes in data recoding. Journal of the American Statistical
    Association.

3.  D’Orazio, M (2015). Integration and imputation of survey data in R:
    the StatMatch package. Romanian Statistical Review, vol. 63(2),
    pages 57-68

4.  Strobl C, Boulesteix A-L, Kneib T, Augustin T, Zeileis A (2008).
    Conditional Variable Importance for Random Forests. BMC
    Bioinformatics, 9, 307 -
    <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-307>

 
