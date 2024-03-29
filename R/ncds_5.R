#' National Child Development Study: a sample of the fifth wave of data collection
#'
#' This database is a sample of the fifth wave of data collection of the National Child Development Study (NCDS)
#' started in 1958 (\url{https://cls.ucl.ac.uk/cls-studies/1958-national-child-development-study/}).
#' The NCDS project is a continuing survey which follows the lives of over 17,000 people born in England,
#' Scotland and Wales in a same week of the year 1958.
#'
#' The ncds identifier have been voluntarily anonymized to allow their availability for the package.
#'
#' This sample has 365 participants included in the study during the 5th waves of data collection.
#'
#'
#' @format A data.frame with 365 participants (rows) and 6 variables
#' \describe{
#'   \item{ncdsid}{the anonymised ncds identifier}
#'   \item{gender}{the gender of the participant stored in a 2-levels factor: \code{1} for male, \code{2} for female}
#'   \item{RG91}{the RG social class 91 scale coded as a 7-levels factor: \code{10} for professional educations,
#'   \code{20} for managerial and technical occupations, \code{31} for skilled non-manual occupations,
#'   \code{32} for skilled manual occupations, \code{40} for party-skilled occupations, \code{50} for unskilled occupations \code{50},
#'   and \code{0} when the scale was not applicable to the participant. This variable is complete.}
#'   \item{health}{the health status of the participant stored in a 4 ordered levels factor: \code{1} for excellent,
#'   \code{2} for good, \code{3} for fair, \code{4} for poor. This variable has 2 NAs.}
#'   \item{employ}{the employment status at inclusion stored in a 7-levels factor: \code{1} for unemployed status, \code{2} for govt sheme,
#'   \code{3} for full-time education, \code{4} for housework or childcare, \code{5} for sick or handicapped, \code{6} for other, \code{7}
#'   if employed between 16 and 33. This variable has 58 NAs.}
#'   \item{study}{a 2-level factor equals to \code{1} for participant with completed graduate studies or \code{2} otherwise}
#' }
#' @source INSERM - This database is a sample of the National Child Development Study
"ncds_5"
