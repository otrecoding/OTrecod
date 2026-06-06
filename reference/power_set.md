# power_set()

A function that gives the power set \\P(S)\\ of any non empty set S.

## Usage

``` r
power_set(n, ordinal = FALSE)
```

## Arguments

- n:

  an integer. The cardinal of the set

- ordinal:

  a boolean. If TRUE the power set is only composed of subsets of
  consecutive elements, FALSE (by default) otherwise.

## Value

A list of \\2^n -1\\ subsets (The empty set is excluded)

## References

Devlin, Keith J (1979). Fundamentals of contemporary set theory.
Universitext. Springer-Verlag

## Author

Gregory Guernec

<otrecod.pkg@gmail.com>

## Examples

``` r
# Powerset of set of 4 elements
set1 <- power_set(4)

# Powerset of set of 4 elements by only keeping
# subsets of consecutive elements
set2 <- power_set(4, ordinal = TRUE)
```
