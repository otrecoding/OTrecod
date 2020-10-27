#' transfo_quali()
#'
#' A function that transforms a factor of n(>1) levels in (n-1) binary variables.
#'
#' @aliases transfo_quali
#'
#' @param x a factor
#' @param labx a new label for the generated binary variables (By default the name of the factor is conserved)
#'
#' @return A matrix of (n-1) binary variables
#'
#' @author Gregory Guernec
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @export
#'
#' @aliases transfo_quali
#'
#' @examples
#'
#' treat     = as.factor(c(rep("A",10),rep("B",15),rep("C",12)))
#' treat_bin = transfo_quali(treat,"trt")
#'
#'
transfo_quali = function(x,labx = NULL){

  if (is.factor(x)==FALSE){
    stop("Your variable needs to be convert in factor !")
  }

  lev_x = levels(x)



  if (length(lev_x) == 1){

    if (lev_x == "0"){

      x_bin = as.matrix(as.numeric(x)-1)
      colnames(x_bin)[1] = labx

    } else if (lev_x != "0"){

      x_bin = as.matrix(rep(1,length(x)))
      colnames(x_bin)[1] = labx


    }

  } else {

    x_bin = matrix(nrow = length(x), ncol = length(lev_x) - 1)

    for (j in 1:(length(lev_x)-1)){

      x_bin[,j] = ifelse(x == lev_x[j+1],1,0)

    }


    # Names of the binary variables

    colnames(x_bin) = paste(labx,2:(length(lev_x)),sep="_")

  }

  return(x_bin)

}
