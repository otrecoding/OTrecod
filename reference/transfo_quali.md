# transfo_quali()

A function that transforms a factor of n(\>1) levels in (n-1) binary
variables.

## Usage

``` r
transfo_quali(x, labx = NULL)
```

## Arguments

- x:

  a factor

- labx:

  a new label for the generated binary variables (By default the name of
  the factor is conserved)

## Value

A matrix of (n-1) binary variables

## Author

Gregory Guernec

<otrecod.pkg@gmail.com>

## Examples

``` r

treat <- as.factor(c(rep("A", 10), rep("B", 15), rep("C", 12)))
treat_bin <- transfo_quali(treat, "trt")
```
