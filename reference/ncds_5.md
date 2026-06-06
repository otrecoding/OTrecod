# National Child Development Study: a sample of the fifth wave of data collection

This database is a sample of the fifth wave of data collection of the
National Child Development Study (NCDS) started in 1958
(<https://cls.ucl.ac.uk/cls-studies/1958-national-child-development-study/>).
The NCDS project is a continuing survey which follows the lives of over
17,000 people born in England, Scotland and Wales in a same week of the
year 1958.

## Usage

``` r
ncds_5
```

## Format

A data.frame with 365 participants (rows) and 6 variables

- ncdsid:

  the anonymised ncds identifier

- gender:

  the gender of the participant stored in a 2-levels factor: `1` for
  male, `2` for female

- RG91:

  the RG social class 91 scale coded as a 7-levels factor: `10` for
  professional educations, `20` for managerial and technical
  occupations, `31` for skilled non-manual occupations, `32` for skilled
  manual occupations, `40` for party-skilled occupations, `50` for
  unskilled occupations `50`, and `0` when the scale was not applicable
  to the participant. This variable is complete.

- health:

  the health status of the participant stored in a 4 ordered levels
  factor: `1` for excellent, `2` for good, `3` for fair, `4` for poor.
  This variable has 2 NAs.

- employ:

  the employment status at inclusion stored in a 7-levels factor: `1`
  for unemployed status, `2` for govt sheme, `3` for full-time
  education, `4` for housework or childcare, `5` for sick or
  handicapped, `6` for other, `7` if employed between 16 and 33. This
  variable has 58 NAs.

- study:

  a 2-level factor equals to `1` for participant with completed graduate
  studies or `2` otherwise

## Source

INSERM - This database is a sample of the National Child Development
Study

## Details

The ncds identifier have been voluntarily anonymized to allow their
availability for the package.

This sample has 365 participants included in the study during the 5th
waves of data collection.
