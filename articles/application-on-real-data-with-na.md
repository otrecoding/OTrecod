# application-on-real-data-with-na

## The NCDS project

### DISCLAIMER

The data owner is the ‘University of London, Institute of Education,
Centre for Longitudinal Studies’ and those safeguarded CLS data are
available to UK Data Service users subject to the End User Licence
Agreement. Clause 4 of the agreement states:

“Not to give access to the Data Collection(s), in whole or in part, or
to any Dataset(s) derived from the Data Collection(s) (including
synthetic Dataset(s)), except to Registered Users for use in accordance
with clause 2 and subject to the EUL Agreement, any relevant Special
Conditions or additional agreements, save where Data Collection(s) or
derived material are supplied for the stated purpose of teaching as
described in the Project Information and shared solely in accordance
with the Access Agreement for Teaching.”

 

The National Child Development Study project also called NCDS project
(<https://cls.ucl.ac.uk/cls-studies/1958-national-child-development-study/>)
is a continuing survey which follows the lives of over 17,000 people
born in England, Scotland and Wales in a same week of the year 1958.

The project is well-known today because its results have greatly
contributed to improve the quality of maternity services in UK.

This survey collects specific information on many distinct fields like
physical and educational development, economic circumstances,
employment, family life, health behaviour, well-being, social
participation and attitudes.

 

### The problem

 

Two possible measurements scales of the social class of the participants
built from profession can have been collected at each updates of data
collection.

The first one corresponds to the **Goldthorp Social Class 90** scale
(**GO90**), a scale with 12 categories, while the second one corresponds
to the **RG Social Class 91** scale (**RG91**), a scale with only 7
categories respectively described in the documentation of the *ncds_14*
and *ncds_5* databases.

If during the first survey waves, the information related to the social
level of each participant is collected using the **GO90** scale, it is
the **RG91** scale that was used exclusively from wave 5. This variation
in data collection could pose a problem related to the exploitation of
this variable, when including new participants in cohort studies,
especially because there is actually no overlaps bethween the two
databases.

**Question: how to work with the greatest number of participants
encountered in the study without losing the information related to the
individual social level?**

Answer: We will see in this vignette how the use of the **OTrecod**
functions can provide a concrete answer to this question.

 

### The solution

 

Functions from the **OTrecod** package can help users to solve the
problem introduced in the previous section by recoding the missing scale
**GO90** or/and **RG91** in the corresponding databases as described in
the following paragraphs dedicated to the imputation of **GO90** in
*ncds_5*. However, it would be obvioulsy possible to impute **RG91** in
*ncds_14* or both in a similar process.

 

### Package installation

 

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

 

### databases installation

 

Participants from waves 1 to 4 were stored in a database called
*ncds_14* provided in the **OTrecod** package. In the same way,
informations related to new participants included in the study from the
wave 5 were stored in the *ncds_5* database. These two base are
available in **OTrecod** and can be simply loaded as follows:

``` r

data(ncds_14); data(ncds_5)
```

Here is the description of all the variables contained in each of the
databases:

``` r

summary(ncds_14); summary(ncds_5)
```

From these descriptive statistics, we can make the following remarks:

1.  Each of the files contains only one scale measuring the social class
    of the participants: the *ncds_14* database measures the social
    class from the **G090** scale only, which currently makes it
    impossible to study this variable on all the waves of data
    collection carried out.
2.  The other variables seem at first sight labeled and coded in a
    similar way from one dataset to another.
3.  Some variables have incomplete informations (including social class
    scales in *ncds_14*) and the missingness structures seems
    independant from the databases.
4.  There are only a few variables contained in the two databases.
5.  Same labels predictors are not sorted in the same order from one
    database to another.

 

### Handling missing information

 

Firsly, **JOINT** algorithms from the **OT_joint** function of the
package requires complete predictors to run, which is not the case for
**OUTCOME** algorithms implemented in the **OT_outcome** function.
Keeping this information in mind, we know that for a comparison of these
two families of algorithm, it will be necessary to work either with
complete cases predictors or to carry out a preliminary imputation on
them before the data fusion.

The **merge_dbs** function can solve problems 2 and 3 introduced in the
previous paragraph:

- by removing from the upcoming data fusion, participants whose social
  class informations are missing from their databases (corresponding to
  806 participants in *ncds_14*).
- by detecting and excluding predictors, labeled and/or coded
  differently from one database to another.
- possibly imputing incomplete predictors using MICE or FAMD procedures
  (see the documentation of **imput_cov** for mored details about the
  methods).
- finally, by superimposing the two datasets by only keeping the shared
  predictors.

Here is the R code related to a MICE imputation (with only 2
replications for running time reasons of the vignette). Using the MICE
procedure directly integrated in **merge_dbs** supposes that the user
accepts the hypothesis that all predictors will participate to the
imputation.

In situation where this hypothesis does not hold but imputation still
appear as necessary, the imputation of missing values will have to be
done through another function (not provided by the package).

``` r

merged_tab <- merge_dbs(ncds_14, ncds_5,
                        row_ID1 = 1, row_ID2 = 1,
                        NAME_Y = "GO90", NAME_Z = "RG91",
                        ordinal_DB1 = 3, ordinal_DB2 = 4,
                        impute = "MICE", R_MICE = 2,
                        seed_choice = 3023)
```

Process details indicate here that 806 participants have been deleted
from the data fusion because of missing information on **GO90** in the
*ncds_14* database.

By definition, **GO90** in *ncds_14* and **RG91** in *ncds_5* does not
seem to correspond clearly to ordered scales and can be considered as
simple factors. This is the case here, because the only variable defined
as an ordered factor in the study is the *health* variable. This
variable correspond to the third column of *ncds_14* and to the fourth
column of *ncds_5*. A seed value is arbitrarily fixed to 3023 here to
guarantee the reproductibility of the results in the vignette.

In the output list, The *DB_READY* object provides the overlayed
database with the shared predictors imputed as described by:

``` r

head(merged_tab$DB_READY); dim(merged_tab$DB_READY)
```

All predictors same labels were finally kept: their coding being visibly
in adequacy from one dataset to another.

 

### Matching predictors evaluation

 

Here, the limited number of predictors at the beginning of the study
makes it unnecessary to use a selection procedure like random forest,
nevertheless the **select_pred** function can provide users interesting
information about the association of predictors, their ability to fit
the target scales in each dataset, and their risks of collinearities,
especially if it is necessary to further reduce the set of sharing
predictors before recoding scales.

Keeping too many predictors for the fusion does not improve its quality
and greatly increases the resolution time of the algorithms, that is
why, as an advice, keeping only 3 predictors for each fusion seems to be
a good compromise between computation time and efficiency.

Let’s use the results of the **select_pred** function to only keep 3
predictors for data fusion (instead of the four currently present).

Here is the related R code taking the database obtained after executing
**merge_dbs** as argument, and the variables **GO90** and **RG91** as
respective outcomes in the two datasets:

``` r

# ncds_14
sel_pred_GO90 <- select_pred(merged_tab$DB_READY,
                             Y = "Y", Z = "Z",
                             ID = 1 , OUT = "Y",
                             nominal = c(1:5,7), ordinal = 6,
                             RF = FALSE)
```

``` r

# ncds_5
sel_pred_RG91 <- select_pred(merged_tab$DB_READY, Y = "Y", Z = "Z",
                             ID = 1, OUT = "Z",
                             nominal = c(1:5,7), ordinal = 6,
                             RF = FALSE)
```

Automatic exits inform the user about the state of the process and
assist him in the implementation of the function.

We first study the existence of potential collinearity between the
predictors. At the acceptability thresholds set by default as input to
the **merge_dbs** function, and in particular, the one associated with
the V-Cramer criterion, we detect collinearity in the two bases, between
the variables *employ* and *gender*, as described below:

``` r

## ncds_14
sel_pred_GO90$collinear_PB$VCRAM

## ncds_5
sel_pred_RG91$collinear_PB$VCRAM
```

This result suggests that we should keep only one of these two
predictors. But which of the two should we keep? *employ* or *gender* ?

To help us decide, we can simply observe the ability of the two
variables to predict the target variables in the two databases:

``` r

## ncds_14
sel_pred_GO90$vcrm_OUTC_cat

## ncds_5
sel_pred_RG91$vcrm_OUTC_cat
```

These results present two variables whose ability to predict the target
variable is equivalent in the two databases. How do these variables
interact with the other predictors?

``` r

## ncds_14
sel_pred_GO90$vcrm_X_cat

## ncds_5
sel_pred_RG91$vcrm_X_cat
```

The stronger association in both datasets of the variable *employ* with
*study*, compared to that observed with the variable *gender*, would
rather encourage us to keep *gender* for the data fusion. Even if, it is
also true that this selection criterion remains very tenuous here.

We therefore decide to keep the *study*, *gender* and *health*
predictors to solve our recoding problem.

``` r

merged_fin = merged_tab$DB_READY[, -4]
head(merged_fin, 3)
```

 

### Imputation of GO90 using optimal transportation theory

The imputation of **GO90** in the *ncds_5* can be the result of two
family of algorithms in **OTrecod**: **OUTCOME** or **JOINT**. Both of
them used optimal transportation theory and can be eventually improved
by relaxing the distributional assumptions.

**OUTCOME** algorithms are implemented in the **OT_outcome** function
while **JOINT** algorithms are implemented in **OT_JOINT**.

It is suggested to fit and compare several models using the output
criteria of the **verif_OT** function, to determine which one provides
the better set of **GO90** predictions in *ncds_5*.

Here is the R commands to fit an **OUTCOME** model and a **JOINT** model
with no relaxation and regularization parameter.

``` r

outc1 = OT_outcome(merged_fin, nominal = c(1:4), ordinal = 5:6,
                   dist.choice = "E", indiv.method = "sequential",
                   which.DB = "B")
```

``` r

outj1 = OT_joint(merged_fin, nominal = 1:4, ordinal = 5:6,
                 dist.choice = "E", which.DB = "B")
```

 

Now let’s see the validity criteria of the two models using the
**verif_OT** function. we simply want to retrieve the stability criteria
and Cramer’s V which compares the association of the two different
(unordered) scales and so complete the arguments according to this
choice:

``` r


### For the model outc1
verif_outc1   = verif_OT(outc1, ordinal = FALSE, 
                         stab.prob = TRUE, min.neigb = 5)

### For the model outj1
verif_outj1   = verif_OT(outj1, ordinal = FALSE, 
                         stab.prob = TRUE, min.neigb = 5)
```

Firstly, the estimates of the Hellinger distance for the two models
(both upper 0.05) confirm that the two scales could respectively share
the same distribution in the two databases. Consequently, there is no
need here to add relaxation parameters in the models to have acceptable
predictions of **GO90** in *ncds_5*.

``` r


### For the model outc1: Hellinger distance
verif_outc1$hell

### For the model outj1: Hellinger distance
verif_outj1$hell
```

The V Cramer criterion informs about the strength of the association
between the two scales. Even if they have specific encodings, **GO90**
and **RG91** summarize a same information: the social class of the
participant. Therefore, we can reasonably think that the stronger will
be the association, the more reliable will be the predictions.

``` r


### For the model outc1: Hellinger distance
verif_outc1$res.prox

### For the model outj1: Hellinger distance
verif_outj1$res.prox
```

The large number of modalities of **GO90** explains the relative
stability observed for the predictions:

``` r


### For the model outc1: Hellinger distance
verif_outc1$res.stab

### For the model outj1: Hellinger distance
verif_outj1$res.stab
```

### Conclusion

The two models provide acceptable predictions of **GO90** in *ncds_5*.
Nevertheless, the **JOINT** algorithm is here prefered for the reasons
previously exposed (no relaxation and regularization parameters
required, better reliability of predictions with this specific
database).
