# ham()

This function computes a matrix distance using the Hamming distance as
proximity measure.

## Usage

``` r
ham(mat_1, mat_2)
```

## Arguments

- mat_1:

  a vector, a matrix or a data.frame of binary values that may contain
  missing data

- mat_2:

  a vector, a matrix or a data.frame of binary values with the same
  number of columns as `mat_1` that may contain missing data

## Value

A distance matrix

## Details

`ham` returns the pairwise distances between rows (observations) of a
single matrix if `mat_1` equals `mat_2`. Otherwise `ham` returns the
matrix distance between rows of the two matrices `mat_1` and `mat_2` if
this 2 matrices are different in input. Computing the Hamming distance
stays possible despite the presence of missing data by applying the
following formula. Assuming that A and B are 2 matrices such as
`ncol(A) = ncol(B)`. The Hamming distance between the \\i^{th}\\ row of
A and the \\k^{th}\\ row of B equals:

\$\$\mbox{ham}(A_i,B_k) = \frac{\sum_j 1\_{\left\\A\_{ij} \neq
B\_{kj}\right\\}}{\sum_j 1}\times\left(\frac{\sum_j 1}{\sum_j
1\_{\left\\!\mbox{is.na}(A\_{ij}) \\ !\mbox{is.na}(
B\_{kj})\right\\}}\right)\$\$

where: \\i = 1,\dots,\mbox{nrow}(A)\\ and \\k =
1,\dots,\mbox{nrow}(B)\\; And the expression located to the right term
of the multiplication corresponds to a specific weigh applied in
presence of NAs in \\A_i\\ and/or \\B_k\\.

This specificity is not implemented in the `cdist` function and the
Hamming distance can not be computed using the `dist` function either.

The Hamming distance can not be calculated in only two situations:

1.  If a row of A or B has only missing values (ie for each of the
    columns of A or B respectively).

2.  The union of the indexes of the missing values in row i of A with
    the indexes of the missing values in row j of B concerns the indexes
    of all considered columns.

Example: Assuming that \\\mbox{ncol}(A) = \mbox{ncol}(B) = 3\\, if \\A_i
= (1,\mbox{NA},0)\\ and \\B_j = (\mbox{NA},1,\mbox{NA})\\, for each
column, either the information in row i is missing in A, or the
information is missing in B, which induces: \\\mbox{ham}(A_i,B_k) =
\mbox{NA}\\.

If `mat_1` is a vector and `mat_2` is a matrix (or data.frame) or vice
versa, the length of `mat_1` must be equal to the number of columns of
`mat_2`.

## References

Roth R (2006). Introduction to Coding Theory. Cambridge University
Press.

## Author

Gregory Guernec

<otrecod.pkg@gmail.com>

## Examples

``` r
set.seed(3010)
sample_A <- sample(c(0, 1), 12, replace = TRUE)
set.seed(3007)
sample_B <- sample(c(0, 1), 15, replace = TRUE)
A <- matrix(sample_A, ncol = 3)
B <- matrix(sample_B, ncol = 3)

# These 2 matrices have no missing values

# Matrix of pairwise distances with A:
ham(A, A)
#>           [,1]      [,2]      [,3]      [,4]
#> [1,] 0.0000000 0.6666667 0.6666667 1.0000000
#> [2,] 0.6666667 0.0000000 0.6666667 0.3333333
#> [3,] 0.6666667 0.6666667 0.0000000 0.3333333
#> [4,] 1.0000000 0.3333333 0.3333333 0.0000000

# Matrix of distances between the rows of A and the rows of B:
ham(A, B)
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.0000000 0.3333333 0.6666667 0.6666667
#> [2,] 0.3333333 0.6666667 1.0000000 0.6666667 0.6666667
#> [3,] 0.3333333 0.6666667 0.3333333 0.0000000 0.6666667
#> [4,] 0.0000000 1.0000000 0.6666667 0.3333333 0.3333333

# If mat_1 is a vector of binary values:
ham(c(0, 1, 0), B)
#> [1] 0.6666667 0.3333333 0.6666667 0.3333333 1.0000000

# Now by considering A_NA and B_NA two matrices built from A and B respectively,
# where missing values have been manually added:
A_NA <- A
A_NA[3, 1] <- NA
A_NA[2, 2:3] <- rep(NA, 2)

B_NA <- B
B_NA[2, 2] <- NA

ham(A_NA, B_NA)
#>      [,1] [,2]      [,3]      [,4]      [,5]
#> [1,]  1.0    0 0.3333333 0.6666667 0.6666667
#> [2,]  0.0    1 1.0000000 0.0000000 1.0000000
#> [3,]  0.5    0 0.0000000 0.0000000 0.5000000
#> [4,]  0.0    1 0.6666667 0.3333333 0.3333333
```
