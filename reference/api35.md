# Student performance in California schools: the results of the county 35

This database is a sample of the API program
<https://www.cde.ca.gov/re/pr/api.asp> that ended in 2018. The sample is
extracted from the data [`api`](https://rdrr.io/pkg/survey/man/api.html)
of the package survey, related to the results of the county 35 (San
Benito). The database contains information for the 362 schools of this
county having at least 100 students. Missing information has been
randomly (and voluntary) added to the `awards` and `ell` variables (4%
and 7% respectively). Several variables have been voluntary categorized
from their initial types.

## Usage

``` r
api35
```

## Format

A data.frame with 362 schools (rows) and 12 variables

- cds:

  the school identifier

- apicl_1999:

  the API score in 1999 classed in 4 ordered levels: `G1`,`G2`,`G3`,
  `G4`

- stype:

  the school type in a 3 ordered levels factor: `Elementary`, `Middle`
  or `High School`

- awards:

  the school eligible for awards program ? Two possible answers: `No` or
  `Yes`. This variable counts 4% of missing information.

- acs.core:

  the number of core academic courses in the school

- api.stu:

  the number of students tested in the school

- acs.k3.20:

  the average class size years K-3 in the school. This variable is
  stored in a 3-levels factor: `Unknown`, `<=20`, `>20`.

- grad.sch:

  the percentage of parents with postgraduate education stored in a 3
  ordered levels factor of percents: `0`, `1-10`, `>10`

- ell:

  the percentage of English language learners stored in a 4 ordered
  levels factor: `[0-10]`,`(10-30]`,`(30-50]`,`(50-100]`. This variable
  counts 7% of missing information.

- mobility:

  the percentage of students for whom this is the first year at the
  school, stored in 2 levels: `1` and `2`

- meals:

  the percentage of students eligible for subsidized meals stored in a 4
  balanced levels factor (By quartiles): `[0-25]`, `(25-50]`, `(50-75]`,
  `(75-100]`

- full:

  the percentage of fully qualified teachers stored in a 2-levels
  factor: `1`: For strictly less than 90%, `2` otherwise

## Source

This database is a sample of the data
[`api`](https://rdrr.io/pkg/survey/man/api.html) from the package
survey.
