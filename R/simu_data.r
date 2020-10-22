#' A simulated dataset to test the functions of the OTrecod package
#'
#' The first 300 rows belong to the database A, while the next 400 rows belong to the database B.
#' Five covariates: \code{Gender}, \code{Treatment}, \code{Dosage}, \code{Smoking} and \code{Age} are
#' common to both databases (same encodings). \code{Gender} is the only complete covariate.
#' The variables \code{Yb1} and \code{Yb2} are the target variables of A and B respectively, summarizing a same information encoded in two different scales.
#' that summarize a same information saved in two distinct encodings, that is why, \code{Yb1} is
#' missing in the database B and \code{Yb2} is missing in the database A.
#'
#' The purpose of the functions contained in this package is to predict the missing information on \code{Yb1} and \code{Yb2}
#' in database A and database B using the Optimal Transportation Theory.
#'
#' Missing information has been simulated to some covariates following a simple MCAR process.
#'
#'
#' @format A data.frame made of 2 overlayed databases (A and B) with 700 observations on the following 8 variables.
#' \describe{
#'   \item{DB}{the database identifier, a character with 2 possible classes: \code{A} or \code{B}}
#'   \item{Yb1}{the target variable of the database A, stored as factor and encoded in 3 ordered levels: \code{[20-40]}, \code{[40-60[},\code{[60-80]} (the values related to the database B are missing)}
#'   \item{Yb2}{the target variable of the database B, stored as integer (an unknown scale from 1 to 5) in the database B (the values related to A are missing)}
#'   \item{Gender}{a factor with 2 levels (\code{Female} or \code{Male}) and no missing values}
#'   \item{Treatment}{a covariate of 3 classes stored as a character with 2\% of missing values: \code{Placebo}, \code{Trt A}, \code{Trt B}}
#'   \item{Dosage}{a factor with 4 levels and 5\% of missing values: from \code{Dos 1} to \code{dos 4}}
#'   \item{Smoking}{a covariate of 2 classes stored as a character and 10\% of missing values: \code{NO} for non smoker, \code{YES} otherwise}
#'   \item{Age}{a numeric corresponding to the age of participants in years. This variable counts 5\% of missing values}
#' }
#' @source randomly generated
"simu_data"
