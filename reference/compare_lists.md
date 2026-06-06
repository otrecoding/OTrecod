# compare_lists()

This function compares the elements of two lists of same length.

## Usage

``` r
compare_lists(listA, listB)
```

## Arguments

- listA:

  a first list

- listB:

  a second list

## Value

A boolean vector of same length as the two lists, which ith element is
`TRUE` if the ith element is different between the 2 lists, or `FALSE`
otherwise

## Author

Gregory Guernec

<otrecod.pkg@gmail.com>

## Examples

``` r
data1 <- data.frame(Gender = rep(c("m", "f"), 5), Age = rnorm(5, 20, 4))
data2 <- data.frame(Gender = rep(c("m", "f"), 5), Age = rnorm(5, 21, 5))

list1 <- list(A = 1:4, B = as.factor(c("A", "B", "C")), C = matrix(1:6, ncol = 3))
list2 <- list(A = 1:4, B = as.factor(c("A", "B")), C = matrix(1:6, ncol = 3))
list3 <- list(A = 1:4, B = as.factor(c("A", "B", "C")), C = matrix(c(1:5, 7), ncol = 3))
list4 <- list(A = 1:4, B = as.factor(c("A", "B", "C")), C = matrix(1:6, ncol = 2))
list5 <- list(A = 1:4, B = as.factor(c("A", "B")), C = matrix(1:6, ncol = 2))
list6 <- list(A = 1:4, B = as.factor(c("A", "B")), C = data1)
list7 <- list(A = 1:4, B = as.factor(c("A", "B")), C = data2)

OTrecod::compare_lists(list1, list2)
#> [1] FALSE  TRUE FALSE
OTrecod::compare_lists(list1, list3)
#> [1] FALSE FALSE  TRUE
OTrecod::compare_lists(list1, list4)
#> [1] FALSE FALSE  TRUE
OTrecod::compare_lists(list1, list5)
#> [1] FALSE  TRUE  TRUE
OTrecod::compare_lists(list6, list7)
#> [1] FALSE FALSE  TRUE
```
