
#' simpond()
#'
#' A function that provides rounded values that respect the margins to the solution matrix from the cost() function
#'
#' @param x A matrix of double
#' @param marg A vector of row margins
#'
#' @return A matrix of integers
#' @export
#'
#' @examples
#' soluc  = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
#' prep   = transfo_dist(soluc[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")
#' xx     = cost(prep,"G")
#'
#' dat = soluc[[1]]
#'
#' n1    = nrow(dat[dat[,1]==1,])
#' n2    = nrow(dat[dat[,1]==2,])
#'
#' val1  = sort(unique(dat[dat[,1]==1,2]))
#' val2  = sort(unique(dat[dat[,1]==2,3]))
#'
#' mat = matrix(round(xx$solution*n1,4),nrow=length(val2),ncol=length(val1)); mat
#' simpond(mat,prep$Y2)
#'
#'
simpond = function(matx,marg){

  if (!(is.matrix(matx))){

    stop("Your object is not a matrix")

  } else {}

  y      = floor(matx)
  xx     = matx - y
  nbcand = table(marg)-rowSums(y)
  indic  = (1:length(nbcand))[nbcand!=0]
  ordx   = t(apply(xx,1,function(x){order(x,decreasing = TRUE)}))

  for (i in indic){
    y[i,ordx[i,1:nbcand[i]]] = y[i,ordx[i,1:nbcand[i]]] + 1
  }

  return(y)

}




