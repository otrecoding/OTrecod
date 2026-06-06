# A simulated dataset to test the functions of the OTrecod package

The first 300 rows belong to the database A, while the next 400 rows
belong to the database B. Five covariates: `Gender`, `Treatment`,
`Dosage`, `Smoking` and `Age` are common to both databases (same
encodings). `Gender` is the only complete covariate. The variables `Yb1`
and `Yb2` are the target variables of A and B respectively, summarizing
a same information encoded in two different scales. that summarize a
same information saved in two distinct encodings, that is why, `Yb1` is
missing in the database B and `Yb2` is missing in the database A.

## Usage

``` r
simu_data
```

## Format

A data.frame made of 2 overlayed databases (A and B) with 700
observations on the following 8 variables.

- DB:

  the database identifier, a character with 2 possible classes: `A` or
  `B`

- Yb1:

  the target variable of the database A, stored as factor and encoded in
  3 ordered levels: `[20-40]`, `[40-60[`,`[60-80]` (the values related
  to the database B are missing)

- Yb2:

  the target variable of the database B, stored as integer (an unknown
  scale from 1 to 5) in the database B (the values related to A are
  missing)

- Gender:

  a factor with 2 levels (`Female` or `Male`) and no missing values

- Treatment:

  a covariate of 3 classes stored as a character with 2% of missing
  values: `Placebo`, `Trt A`, `Trt B`

- Dosage:

  a factor with 4 levels and 5% of missing values: from `Dos 1` to
  `dos 4`

- Smoking:

  a covariate of 2 classes stored as a character and 10% of missing
  values: `NO` for non smoker, `YES` otherwise

- Age:

  a numeric corresponding to the age of participants in years. This
  variable counts 5% of missing values

## Source

randomly generated

## Details

The purpose of the functions contained in this package is to predict the
missing information on `Yb1` and `Yb2` in database A and database B
using the Optimal Transportation Theory.

Missing information has been simulated to some covariates following a
simple MCAR process.
