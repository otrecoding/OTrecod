
#' transfo_target()
#'
#' This function prepares the encoding of the target variable before running an algorithm using optimal transportation theory.
#'
#' The function \code{transfo_target} is an intermediate function direcly implemented in the functions \code{\link{OT_outcome}} and \code{\link{OT_joint}},
#' two functions dedicated to data fusion (see (1) and (2) for details). Nevertheless, this function can also be used separately to assist user in the conversion
#' of a target variable (outcome) according to the following rules:
#' \itemize{
#' \item A character variable is converted in factor if the argument \code{levels_order} is set to NULL. In this case, the levels of the factor are assigned by order of appearance in the database.
#' \item A character variable is converted in ordered factor if the argument \code{levels_order} differs from NULL. In this case, the levels of the factor correspond to those assigned in the argument.
#' \item A factor stays unchanged if the argument \code{levels_order} is set to NULL. Otherwise the factor is converted in ordered factor and the levels are ordered according to the argument \code{levels_order}.
#' \item A numeric variable, discrete or continuous is converted in factor if the argument \code{levels_order} is set to NULL, and the related levels are the values assigned in ascending order.
#' \item A numeric variable, discrete or continuous is converted in ordered factor if the argument \code{levels_order} differed from NULL, and the related levels correspond to those assigned in the argument.
#' }
#'
#'
#' @param z            a factor variable (ordered or not). A variable of another type will be, by default, convert to a factor.
#' @param levels_order a vector corresponding to the values of the levels of z. When the target is ordinal, the levels can be sorted by ascending order.
#'                     By default, the initial order is remained.
#' @return The list returned is:
#' \item{NEW}{an object of class factor of the same length as z}
#' \item{LEVELS_NEW}{the levels (ordered or not) retained for z}
#'
#' @author Gregory Guernec
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @seealso \code{\link{compare_lists}}
#'
#' @aliases transfo_target
#'
#' @references
#' \enumerate{
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679. doi:10.1515/ijb-2018-0106
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association. \doi{10.1080/01621459.2020.1775615}
#' }
#'
#' @export
#'
#' @examples
#' y      = rnorm(100,30,10)
#' aa     = transfo_target(y)
#'
#' newlev  = unique(as.integer(y))
#' bb      = transfo_target(y,levels_order = newlev)
#' newlev2 = newlev[-1]
#' cc      = transfo_target(y,levels_order = newlev2)
#'
#' outco   = c(rep("A",25),rep("B",50),rep("C",25))
#' dd      = transfo_target(outco,levels_order = c("B","C","A"))
#' ee      = transfo_target(outco,levels_order = c("E","C","A","F"))
#' ff      = transfo_target(outco)
#'
#' outco2  = c(rep("A",25),NA,rep("B",50),rep("C",25),NA,NA)
#' gg      = transfo_target(outco2)
#' hh      = transfo_target(outco2,levels_order = c("B","C","A"))


transfo_target = function(z,levels_order=NULL){


  nlev = length(levels_order)


  if ((class(z)[1] %in% c("character","factor","ordered"))&(nlev == 0)){

    z = as.factor(z)

  } else if ((class(z)[1] %in% c("character","factor","ordered"))&(nlev != 0)){

    z = as.factor(z)

    # if (length(union(levels(z),levels_order)) == nlev){
      if (length(levels(z)) == nlev){

      z = ordered(z,levels = levels_order)

    } else {

      z = as.factor(z)
      message("Inappropriate number or declared labels of levels","\n")
      message("The default levels have been kept","\n")

    }

  } else {}



  if ((is.numeric(z))&(nlev == 0)){

    # cat("Your target was numeric ... By default, it has been converted in factor of integers","\n")
    message(paste("Your target",deparse(substitute(z)),"was numeric ... By default, it has been converted in factor of integers",sep=" "),"\n")
    z = as.factor(as.integer(z))
    message(paste(length(levels(z)),"remaining levels",sep=" "),"\n")


  } else if ((is.numeric(z))&(nlev!=0)){

    message(paste("Your target",deparse(substitute(z)),"was numeric ... By default, it has been converted in factor of integers",sep=" "),"\n")
    z = as.factor(as.integer(z))

    if (length(union(levels(z),levels_order)) == nlev){

      z = ordered(z,levels = levels_order)

    } else {

      z = as.factor(as.integer(z))
      message("Inappropriate number or declared labels of levels","\n")
      message("The default labels have been kept","\n")

    }

  } else {}


  return(list(NEW = z,LEVELS_NEW = levels(z)))

}







