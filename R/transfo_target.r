
#' transfo_target()
#'
#' This function prepares the encoding of the target variable before OT algorithm
#'
#' @param z            A factor variable (ordered or not). a variable of another type will be, by default, convert to a factor
#' @param levels_order A vector corresponding to the values of the levels of z. When the target is ordinal, the levels can be sorted by ascending order.
#'                     By default, the initial order is retained
#' @return The list returned is:
#' \item{NEW}{An object of class factor of the same length as z}
#' \item{LEVELS_NEW}{The levels (ordered or not) retained for z}
#'
#' @author Gregory Guernec
#' \email{otrecod.pkg@@gmail.com}
#'
#' @seealso \code{\link{compare_lists}}
#'
#' @aliases transfo_target
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

    if (length(union(levels(z),levels_order)) == nlev){

      z = ordered(z,levels = levels_order)

    } else {

      z = as.factor(z)
      cat("Inappropriate number or declared labels of levels","\n")
      cat("The default levels have been kept","\n")

    }

  } else {}



  if ((is.numeric(z))&(nlev == 0)){

    # cat("Your target was numeric ... By default, it has been converted in factor of integers","\n")
    cat(paste("Your target",deparse(substitute(z)),"was numeric ... By default, it has been converted in factor of integers",sep=" "),"\n")
    z = as.factor(as.integer(z))
    cat(paste(length(levels(z)),"remaining levels",sep=" "),"\n")


  } else if ((is.numeric(z))&(nlev!=0)){

    cat(paste("Your target",deparse(substitute(z)),"was numeric ... By default, it has been converted in factor of integers",sep=" "),"\n")
    z = as.factor(as.integer(z))

    if (length(union(levels(z),levels_order)) == nlev){

      z = ordered(z,levels = levels_order)

    } else {

      z = as.factor(as.integer(z))
      cat("Inappropriate number or declared labels of levels","\n")
      cat("The default labels have been kept","\n")

    }

  } else {}


  return(list(NEW = z,LEVELS_NEW = levels(z)))

}







