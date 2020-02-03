#' error_group()
#'
#' # This function gives an assesment of the proximity of two multinomial distributions
#' # with different number of classes.
#'
#' # By supposing that Y and Z are 2 categorical variables (ordered or not) where the
#' # number of levels of Y i bigger than the number of levels of Z.
#' # The error_group() function researches the optimal grouping of modalities of Y
#' # to approach at best the distribution of Z to give an assessment of the proximity
#' # between the two  distinct encodings.
#'
#' @param REF A factor with a reference number of levels.
#' @param Z A factor with a number of levels greater than the number of levels of the reference.
#' @param ord A boolean. If TRUE, only the subsequent levels of Z can be grouped together.
#'
#' @return A data.frame with 2 columns:
#'         \item{combi}{The first column enumerates all possible groups of modalities of Z to obtain the same number of levels as the reference.}
#'         \item{error_rate}{The second columns gives the corresponding rate error of reclassification.}
#'
#' @author Gregory Guernec
#' \email{gregory.guernec@@inserm.fr}
#'
#' @aliases error_group
#'
#' @export
#'
#' @examples
#'
#' ### using a sample of the tab_test object (3 complete covariates)
#' ### Y1 and Y2 are a same variable encoded in 2 different forms in DB 1 and 2:
#' ### (4 levels for Y1 and 3 levels for Y2)
#'
#' data(tab_test)
#' # Example with n1 = n2 = 70 and only X1 and X2 as covariates
#' tab_test2 = tab_test[c(1:75,5001:5075),1:5]
#'
#' ### An example of JOINT model (Manhattan distance)
#' # Suppose we want to impute the missing parts of Y1 in DB2 only ...
#' try1J = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
#'                  prep_choice = "M", norm = 1,which.DB = 2)
#'
#' # Error rates between Y2 and the predictions of Y1 in the DB 2
#' # by grouping the levels of Y1:
#' error_group(try1J$DATA2_OT$Z,try1J$DATA2_OT$OTpred)
#' table(try1J$DATA2_OT$Z,try1J$DATA2_OT$OTpred)
#'
#' \dontrun{
#' # Error rates between Y1 and the predictions of Y2 in the DB 1
#' error_group(try1J$DATA1_OT$OTpred,try1J$DATA1_OT$Y)
#' table(try1J$DATA1_OT$OTpred,try1J$DATA1_OT$Y)
#'
#' # Basic examples:
#' Z1 = as.factor(sample(1:3,50,replace = TRUE)); length(Z1)
#' Z3 = as.factor(sample(1:2,50,replace = TRUE)); length(Z3)
#' Z2 = as.factor(sample(c("A","B","C","D"),50, replace = TRUE)); length(Z2)
#' Z4 = as.factor(sample(c("A","B","C","D","E"),50, replace = TRUE)); length(Z4)
#'
#'  # By only grouping consecutive levels of Z1:
#' error_group(Z1,Z4)
#' # By only all possible levels of Z1, consecutive or not:
#' error_group(Z3,Z1,FALSE)
#' }
#'
error_group = function(REF,Z,ord = TRUE){

  if ((is.null(levels(REF)))|(is.null(levels(Z)))){

    stop("REF and Z must be factors")

  } else {}


  if (length(levels(REF))>length(levels(Z))){

    stop("The number of levels for Z must be greater than the number of levels of REF")

  } else {}

  cc = try_group(Z,REF,ordin = ord)

  error_g = vector(length = nrow(cc[[1]]))

  for (k in 1:nrow(cc[[1]])){

    Zbis = as.character(Z)

    for (j in 1:ncol(cc[[1]])){

      Zbis[Zbis %in% cc[[2]][[cc[[1]][k,j]]]] = levels(REF)[j]

    }

    Zbis = as.factor(Zbis)

    error_g[k] = 100 - round(sum(diag(table(REF,Zbis)))*100/sum(table(REF,Zbis)),1)
    error_combi    = data.frame(combi = row.names(cc[[1]]),error_rate= error_g)
    error_combi    = error_combi[sort.list(error_combi[,2]),]

  }

  return(error_combi)

}
