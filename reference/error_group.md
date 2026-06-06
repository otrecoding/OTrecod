# error_group()

This function studies the association between two categorical
distributions with different numbers of modalities.

## Usage

``` r
error_group(REF, Z, ord = TRUE)
```

## Arguments

- REF:

  a factor with a reference number of levels.

- Z:

  a factor with a number of levels greater than the number of levels of
  the reference.

- ord:

  a boolean. If TRUE, only neighboring levels of \\Z\\ will be grouped
  and tested together.

## Value

A data.frame with five columns:

- combi:

  the first column enumerates all possible groups of modalities of \\Y\\
  to obtain the same number of levels as the reference.

- error_rate:

  the second column gives the corresponding rate error from the
  confusion matrix (ratio of non-diagonal elements)

- Kappa:

  this column indicates the result of the Cohen's kappa coefficient
  related to each combination of \\Y\\

- Vcramer:

  this column indicates the result of the Cramer's V criterion related
  to each combination of \\Y\\

- RankCor:

  this column indicates the result of the Spearman's coefficient of
  correlation related to each combination of \\Y\\

## Details

Assuming that \\Y\\ and \\Z\\ are categorical variables summarizing a
same information, and that one of the two related encodings is unknown
by user because this latter is, for example, the result of predictions
provided by a given model or algorithm, the function `error_group`
searches for potential links between the modalities of \\Y\\ to approach
at best the distribution of \\Z\\.

Assuming that \\Y\\ and \\Z\\ have \\n_Y\\ and \\n_Z\\ modalities
respectively so that \\n_Y \> n_Z\\, in a first step, the function
`error_group` combines modalities of \\Y\\ to build all possible
variables \\Y'\\ verifying \\n\_{Y'} = n_Z\\. In a second step, the
association between \\Z\\ and each new variable \\Y'\\ generated is
measured by studying the ratio of concordant pairs related to the
confusion matrix but also using standard criterions: the Cramer's V (1),
the Cohen's kappa coefficient (2) and the Spearman's rank correlation
coefficient.

According to the type of \\Y\\, different combinations of modalities are
tested:

- If \\Y\\ and \\Z\\ are ordinal (`ord = TRUE`), only consecutive
  modalities of \\Y\\ will be grouped to build the variables \\Y'\\.

- If \\Y\\ and \\Z\\ are nominal (`ord = FALSE`), all combinations of
  modalities of \\Y\\ (consecutive or not) will be grouped to build the
  variables \\Y'\\.

All the associations tested are listed in output as a data.frame object.
The function `error_group` is directly integrated in the function
[`verif_OT`](https://otrecoding.github.io/OTrecod/reference/verif_OT.md)
to evaluate the proximity of two multinomial distributions, when one of
them is estimated from the predictions of an OT algorithm.

Example: Assuming that \\Y = (1,1,2,2,3,3,4,4)\\ and \\Z =
(1,1,1,1,2,2,2,2)\\, so \\n_Y = 4\\ and \\n_Z = 2\\ and the related
coefficient of correlation \\cor(Y,Z)\\ is 0.89. Are there groupings of
modalities of \\Y\\ which contribute to improving the proximity between
\\Y\\ and \\Z\\ ? From \\Y\\, the function `error_group` gives an answer
to this question by successively constructing the variables: \\Y_1 =
(1,1,1,1,2,2,2,2)\\, \\Y_2 = (1,1,2,2,1,1,2,2)\\, \\Y_3 =
(1,1,2,2,2,2,1,1)\\ and tests \\\mbox{cor}(Z,Y_1) = 1\\,
\\\mbox{cor}(Z,Y_2) = 0\\, \\\mbox{cor}(Z,Y_3) = 0\\. Here, the tests
permit to conclude that the difference of encodings between \\Y\\ and
\\Z\\ resulted in fact in a simple grouping of modalities.

## References

1.  Cramér, Harald. (1946). Mathematical Methods of Statistics.
    Princeton: Princeton University Press.

2.  McHugh, Mary L. (2012). Interrater reliability: The kappa statistic.
    Biochemia Medica. 22 (3): 276–282

## Author

Gregory Guernec

<otrecod.pkg@gmail.com>

## Examples

``` r

# Basic examples:
sample1 <- as.factor(sample(1:3, 50, replace = TRUE))
length(sample1)
#> [1] 50
sample2 <- as.factor(sample(1:2, 50, replace = TRUE))
length(sample2)
#> [1] 50
sample3 <- as.factor(sample(c("A", "B", "C", "D"), 50, replace = TRUE))
length(sample3)
#> [1] 50
sample4 <- as.factor(sample(c("A", "B", "C", "D", "E"), 50, replace = TRUE))
length(sample4)
#> [1] 50

# By only grouping consecutive levels of sample1:
error_group(sample1, sample4)
#>       combi error_rate  Kappa Vcramer RankCor
#> 3 A/B/C D E         56  0.156    0.21   0.249
#> 6 A/B C/D E         58  0.116    0.19   0.213
#> 5 A B/C/D E         62  0.086    0.19   0.203
#> 2 A/B C D/E         68 -0.051    0.23   0.209
#> 4 A B/C D/E         72 -0.076    0.26   0.208
#> 1 A B C/D/E         74 -0.077    0.25   0.146
# By only all possible levels of sample1, consecutive or not:
error_group(sample2, sample1, ord = FALSE)
#>   combi error_rate  Kappa Vcramer RankCor
#> 1 1 2/3         46  0.068    0.03   0.071
#> 3 1/2 3         52 -0.022    0.00  -0.025
#> 2 1 3/2         54 -0.090    0.05  -0.092

# \donttest{

### using a sample of the tab_test object (3 complete covariates)
### Y1 and Y2 are a same variable encoded in 2 different forms in DB 1 and 2:
### (4 levels for Y1 and 3 levels for Y2)

data(tab_test)
# Example with n1 = n2 = 70 and only X1 and X2 as covariates
tab_test2 <- tab_test[c(1:70, 5001:5070), 1:5]

### An example of JOINT model (Manhattan distance)
# Suppose we want to impute the missing parts of Y1 in DB2 only ...
try1J <- OT_joint(tab_test2,
  nominal = c(1, 4:5), ordinal = c(2, 3),
  dist.choice = "M", which.DB = "B"
)
#> ---------------------------------------
#> OT JOINT PROCEDURE in progress ...
#> ---------------------------------------
#> Type                  = JOINT
#> Distance              = Manhattan
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

# Error rates between Y2 and the predictions of Y1 in the DB 2
# by grouping the levels of Y1:
error_group(try1J$DATA2_OT$Z, try1J$DATA2_OT$OTpred)
#>     combi error_rate  Kappa Vcramer RankCor
#> 3 1/2/3 4       47.1  0.286    0.68   0.844
#> 1 1 2/3/4       60.0  0.098    0.95   0.664
#> 2 1/2 3/4       82.9 -0.182    0.74   0.316
table(try1J$DATA2_OT$Z, try1J$DATA2_OT$OTpred)
#>    
#>      1  2  3  4
#>   1 12 16  0  0
#>   2  2  0  0 15
#>   3  0  0 25  0
# }
```
