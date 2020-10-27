#' A simulated dataset to test the library
#'
#' A dataset of 10000 rows containing 3 covariables and 2 outcomes.
#'
#' @format A data frame with 5000 rows and 6 variables:
#' \describe{
#'   \item{ident}{identifier, 1 or 2}
#'   \item{Y1}{outcome 1 with 2 levels, observed for ident=1 and unobserved for ident=2}
#'   \item{Y2}{outcome 2 with 4 levels, observed for ident=2 and unobserved for ident=1}
#'   \item{X1}{covariate 1, integer}
#'   \item{X2}{covariate 2, integer}
#'   \item{X3}{covariate 3, integer}
#' }
#' @source randomly generated
"tab_test"
