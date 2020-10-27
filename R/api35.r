#' Student performance in California schools: the results of the county 35
#'
#'
#' This database is a sample of the API program \url{https://www.cde.ca.gov/re/pr/api.asp} that ended in 2018.
#' The sample is extracted from the data \code{\link[survey]{api}} of the package \pkg{survey}, related to the results of the
#' county 35 (San Benito).
#' The database contains information for the 362 schools of this county having at least 100 students.
#' Missing information has been randomly (and voluntary) added to the \code{awards} and \code{ell} variables (4\% and 7\% respectively).
#' Several variables have been voluntary categorized from their initial types.
#'
#'
#'
#' @format A data.frame with 362 schools (rows) and 12 variables
#' \describe{
#'   \item{cds}{the school identifier}
#'   \item{apicl_1999}{the API score in 1999 classed in 4 ordered levels: \code{G1},\code{G2},\code{G3}, \code{G4}}
#'   \item{stype}{the school type in a 3 ordered levels factor: \code{Elementary}, \code{Middle} or \code{High School}}
#'   \item{awards}{the school eligible for awards program ? Two possible answers: \code{No} or \code{Yes}. This variable counts 4\% of missing information.}
#'   \item{acs.core}{the number of core academic courses in the school}
#'   \item{api.stu}{the number of students tested in the school}
#'   \item{acs.k3.20}{the average class size years K-3 in the school. This variable is stored in a 3-levels factor: \code{Unknown}, \code{<=20}, \code{>20}.}
#'   \item{grad.sch}{the percentage of parents with postgraduate education stored in a 3 ordered levels factor of percents: \code{0}, \code{1-10}, \code{>10}}
#'   \item{ell}{the percentage of English language learners stored in a 4 ordered levels factor: \code{[0-10]},\code{(10-30]},\code{(30-50]},\code{(50-100]}. This variable counts 7\% of missing information.}
#'   \item{mobility}{the percentage of students for whom this is the first year at the school, stored in 2 levels: \code{1} and \code{2}}
#'   \item{meals}{the percentage of students eligible for subsidized meals stored in a 4 balanced levels factor (By quartiles): \code{[0-25]}, \code{(25-50]}, \code{(50-75]}, \code{(75-100]}}
#'   \item{full}{the percentage of fully qualified teachers stored in a 2-levels factor: \code{1}: For strictly less than 90\%, \code{2} otherwise}
#' }
#'
#' @source This database is a sample of the data \code{\link[survey]{api}} from the package \pkg{survey}.
"api35"
