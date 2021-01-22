#' error_group()
#'
#' This function studies the association between two categorical distributions with different numbers of modalities.
#'
#' Assuming that \eqn{Y} and \eqn{Z} are categorical variables summarizing a same information, and that one of the two related encodings is unknown by user
#' because this latter is, for example, the result of predictions provided by a given model or algorithm, the function \code{error_group} searches for potential links between the modalities of \eqn{Y} to approach at best the distribution of \eqn{Z}.
#'
#' Assuming that \eqn{Y} and \eqn{Z} have \eqn{n_Y} and \eqn{n_Z} modalities respectively so that \eqn{n_Y > n_Z}, in a first step, the
#' function \code{error_group} combines modalities of \eqn{Y} to build all possible variables \eqn{Y'} verifying  \eqn{n_{Y'} = n_Z}.
#' In a second step, the association between \eqn{Z} and each new variable \eqn{Y'} generated is measured by studying the ratio of concordant pairs related to the confusion matrix but also using standard criterions:
#' the Cramer's V (1), the Cohen's kappa coefficient (2) and the Spearman's rank correlation coefficient.
#'
#' According to the type of \eqn{Y}, different combinations of modalities are tested:
#' \itemize{
#' \item If \eqn{Y} and \eqn{Z} are ordinal (\code{ord = TRUE}), only consecutive modalities of \eqn{Y} will be grouped to build the variables \eqn{Y'}.
#' \item If \eqn{Y} and \eqn{Z} are nominal (\code{ord = FALSE}), all combinations of modalities of \eqn{Y} (consecutive or not) will be grouped to build the variables \eqn{Y'}.
#' }
#'
#' All the associations tested are listed in output as a data.frame object.
#' The function \code{error_group} is directly integrated in the function \code{\link{verif_OT}} to evaluate the proximity of two multinomial distributions, when one of them is estimated from the predictions of an OT algorithm.
#'
#' Example:
#' Assuming that \eqn{Y = (1,1,2,2,3,3,4,4)} and \eqn{Z = (1,1,1,1,2,2,2,2)}, so \eqn{n_Y = 4} and \eqn{n_Z = 2} and the related coefficient of correlation \eqn{cor(Y,Z)} is 0.89.
#' Are there groupings of modalities of \eqn{Y} which contribute to improving the proximity between \eqn{Y} and \eqn{Z} ?
#' From \eqn{Y}, the function \code{error_group} gives an answer to this question by successively constructing the variables: \eqn{Y_1 = (1,1,1,1,2,2,2,2)}, \eqn{Y_2 = (1,1,2,2,1,1,2,2)}, \eqn{Y_3 = (1,1,2,2,2,2,1,1)}
#' and tests \eqn{\mbox{cor}(Z,Y_1) = 1}, \eqn{\mbox{cor}(Z,Y_2) = 0}, \eqn{\mbox{cor}(Z,Y_3) = 0}.
#' Here, the tests permit to conclude that the difference of encodings between \eqn{Y} and \eqn{Z} resulted in fact in a simple grouping of modalities.
#'
#'
#' @param REF a factor with a reference number of levels.
#' @param Z a factor with a number of levels greater than the number of levels of the reference.
#' @param ord a boolean. If TRUE, only neighboring levels of \eqn{Z} will be grouped and tested together.
#'
#' @return A data.frame with five columns:
#'         \item{combi}{the first column enumerates all possible groups of modalities of \eqn{Y} to obtain the same number of levels as the reference.}
#'         \item{error_rate}{the second column gives the corresponding rate error from the confusion matrix (ratio of non-diagonal elements)}
#'         \item{Kappa}{this column indicates the result of the Cohen's kappa coefficient related to each combination of \eqn{Y}}
#'         \item{Vcramer}{this column indicates the result of the Cramer's V criterion related to each combination of \eqn{Y}}
#'         \item{RankCor}{this column indicates the result of the Spearman's coefficient of correlation related to each combination of \eqn{Y}}
#'
#'
#' @author Gregory Guernec
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @aliases error_group
#'
#' @importFrom stats cor
#' @importFrom StatMatch pw.assoc
#' @importFrom vcd Kappa
#'
#' @references
#' \enumerate{
#' \item Cramér, Harald. (1946). Mathematical Methods of Statistics. Princeton: Princeton University Press.
#' \item McHugh, Mary L. (2012). Interrater reliability: The kappa statistic. Biochemia Medica. 22 (3): 276–282
#' }
#'
#' @export
#'
#' @examples
#'
#' # Basic examples:
#' Z1 = as.factor(sample(1:3,50,replace = TRUE)); length(Z1)
#' Z3 = as.factor(sample(1:2,50,replace = TRUE)); length(Z3)
#' Z2 = as.factor(sample(c("A","B","C","D"),50, replace = TRUE)); length(Z2)
#' Z4 = as.factor(sample(c("A","B","C","D","E"),50, replace = TRUE)); length(Z4)
#'
#' # By only grouping consecutive levels of Z1:
#' error_group(Z1,Z4)
#' # By only all possible levels of Z1, consecutive or not:
#' error_group(Z3,Z1,FALSE)
#'
#'
#' \donttest{
#'
#' ### using a sample of the tab_test object (3 complete covariates)
#' ### Y1 and Y2 are a same variable encoded in 2 different forms in DB 1 and 2:
#' ### (4 levels for Y1 and 3 levels for Y2)
#'
#' data(tab_test)
#' # Example with n1 = n2 = 70 and only X1 and X2 as covariates
#' tab_test2 = tab_test[c(1:70,5001:5070),1:5]
#'
#' ### An example of JOINT model (Manhattan distance)
#' # Suppose we want to impute the missing parts of Y1 in DB2 only ...
#' try1J = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
#'                  dist.choice = "M", which.DB = "B")
#'
#' # Error rates between Y2 and the predictions of Y1 in the DB 2
#' # by grouping the levels of Y1:
#' error_group(try1J$DATA2_OT$Z,try1J$DATA2_OT$OTpred)
#' table(try1J$DATA2_OT$Z,try1J$DATA2_OT$OTpred)
#'
#' }
#'
#'
error_group = function(REF,Z,ord = TRUE){

  if ((is.null(levels(REF)))|(is.null(levels(Z)))){

    stop("REF and Z must be factors")

  } else {}


  if (length(levels(REF))>length(levels(Z))){

    stop("The number of levels for Z must be greater than the number of levels for REF")

  } else {}

  cc      = try_group(Z,REF,ordin = ord)

  error_g = vcram = kap = rankor = vector(length = nrow(cc[[1]]))

  for (k in 1:nrow(cc[[1]])){

    Zbis = as.character(Z)

    for (j in 1:ncol(cc[[1]])){

      Zbis[Zbis %in% cc[[2]][[cc[[1]][k,j]]]] = levels(REF)[j]

    }

    Zbis = ordered(as.factor(Zbis),levels = levels(REF))

    stoc      = data.frame(REF,Zbis)
    vcram[k]  = round(suppressWarnings(StatMatch::pw.assoc(REF~Zbis,data = stoc)$V),2)
    kap[k]    = round(vcd::Kappa(table(REF,Zbis))[[1]][1],3)
    rankor[k] = round(stats::cor(rank(REF),rank(Zbis)),3)

    error_g[k]  = 100 - round(sum(diag(table(REF,Zbis)))*100/sum(table(REF,Zbis)),1)
    error_combi = data.frame(combi = row.names(cc[[1]]),error_rate= error_g,
                             Kappa = kap, Vcramer = vcram, RankCor = rankor)
    error_combi = error_combi[sort.list(error_combi[,2]),]

  }

  return(error_combi)

}
