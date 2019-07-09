
#' transfo_target()
#'
#' A function which prepare the encoding of the target variable before OT algorithm
#'
#' @param z            A vector
#' @param levels_order A vector of levels. When the target is ordinal, the levels can be sorted by ascending order.
#'                     By default, the initial order is retained.
#' @return A list containing two elements:
#'         NEW An object of class factor of the same length as z
#'         LEVELS_NEW The levels (ordered or not) retained for z
#' @export
#'
#' @examples
#' y      = rnorm(100,30,10)
#' aa     = transfo_target(y)
#'
#' newlev = unique(as.integer(y))
#' bb     = transfo_target(y,levels_order = newlev)
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

    cat("Your target was numeric ... By default, it has been converted in factor of integers","\n")
    z = as.factor(as.integer(z))
    cat(paste(length(levels(z)),"remaining levels. See details in output",sep=" "),"\n")


  } else if ((is.numeric(z))&(nlev!=0)){

    cat("Your target was numeric ... By default, it has been converted in factor of integers","\n")
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







