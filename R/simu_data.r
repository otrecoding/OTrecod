#' A simulated dataset to test the functions of the OTrecod package.
#' The first 300 rows belong to the database A, while the next 400 rows belong to the database B.
#' The 5 covariates: \code{Gender},\code{Treatment},\code{Dosage},\code{Smoking} and \code{Age} are
#' so common to both databases (same encodings).
#' The variables \code{Yb1} and \code{Yb2} are the target variables
#' that summarize a same information stored in 2 distinct encodings, that is why, \code{Yb1} is
#' missing in the database B and \code{Yb2} is missing in the database A.
#'
#' The purpose of the functions contained in this package is to predict the missing information on \code{Yb1} and \code{Yb2}
#' in DB A and DB B using the Optimal Transportation Theory.
#'
#' Missing information has been simulated to some covariates following a simple MCAR process.
#'
#'
#' @format A data frame with 700 observations on the following 8 variables.
#' \describe{
#'   \item{DB}{database identification, a character with 2 possible classes: A or B}
#'   \item{Yb1}{the target variable encoded in 3 ordered classes in DB A, missing in database B, stored as factor}
#'   \item{Yb2}{the target variable encoded in 5 ordered classes in DB B, missing in database A, stored as integer}
#'   \item{Gender}{a complete 2 class covariate stored as factor}
#'   \item{Treatment}{a 3 class covariates stored as a character with 2\% of missing information}
#'   \item{Dosage}{a 4 class covariate, stored as a factor, with 5\% of missing information}
#'   \item{Smoking}{a 2 class covariate stored as a character: NO for non smoker, YES otherwise. 10\% of missing information}
#'   \item{Age}{a double corresponding to the age of participants in years. The last covariate of the database counts 5\% of missing information}
#' }
#' @source randomly generated
"simu_data"
