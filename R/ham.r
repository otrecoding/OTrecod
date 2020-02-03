
#' ham()
#'
#' This function computes a matrix distance using the Hamming metric as proximity measure.
#'
#' It returns the pairwise distance between rows (observations) of a single matrix if \code{mat1} equals \code{mat2}.
#' Otherwise it returns the matrix distance between rows of the two matrices \code{mat1} and \code{mat2} if this 2 matrices are different in entry.
#' Calculate the Hamming distance stay possible despite the presence of missing data by applying the following formula.
#'
#' Suppose that A and B are 2 matrices such as ncol(A) = ncol(B). The Hamming distance between the ith row of A and the kt row of B equals:
#'
#' ham(A_i,B_k) = (sum_j(A_ij != B_kj)/sum_j(1))*(sum_j(1)/sum_j(ifelse((!is.na(A_ij))&(!is.na(B_kj)),1,0))
#'
#' Where: i = 1,...,nrow(A) ; k = 1,...,nrow(B); And the term after the "*" corresponds to a weighting applies in presence of NAs in A_i and,or B_k.
#'
#' This option is not implemented in the \code{\link[rdist]{cdist}} function and the Hamming distance can not be computed using the \code{\link[proxy]{dist}} function either.
#'
#' The only 2 situations where the Hamming distance can not be calculated are:
#' 1. If a row of A or B has only missing values (ie on each of the columns of A or B respectively).
#' 2. The union of the indexes of the missing values in row i of A with the indexes of the missing values in row j of B concerns the indexes of all considered columns.
#'    Example: By supposing that ncol(A) = ncol(B) = 3, if A_i = (1,NA,0) and B_j = (NA,1,NA), for each columns, either the information in row i is missing in A,
#'    or the information is missing in B, which induces: ham(A_i,B_K) = NA.
#'
#' If mat_1 is a vector and mat_2 is a matrix (or data.frame) or vice versa, the length of mat_1 must be equal to the number of columns of mat_2.
#'
#' @param mat_1 A vector, a matrix or a data.frame of binary values that may contain missing data
#' @param mat_2 A vector, a matrix or a data.frame of binary values with the same number of columns as \code{mat1} that may contain missing data
#'
#' @return A matrix distance
#' @export
#'
#' @author Gregory Guernec
#' \email{gregory.guernec@@inserm.fr}
#'
#' @aliases ham
#'
#' @examples
#' set.seed(3010); aaa = sample(c(0,1),12,replace = TRUE)
#' set.seed(3007); bbb = sample(c(0,1),15,replace = TRUE)
#' A = matrix(aaa, ncol = 3)
#' B = matrix(bbb, ncol = 3)
#'
#' # These 2 matrices have no missing values
#'
#' # Matrix of pairwise distances with A:
#' ham(A,A)
#'
#' # Matrix of distances between the rows of A and the rows of B:
#' ham(A,B)
#'
#' # If mat_1 is a vector of binary values:
#' ham(c(0,1,0),B)
#'
#' # Now by considering A_NA and B_NA two matrices built from A and B respectively,
#' # where missing values have been manually added:
#' A_NA        = A
#' A_NA[3,1]   = NA
#' A_NA[2,2:3] = rep(NA,2)
#'
#' B_NA = B
#' B_NA[2,2] = NA
#'
#' ham(A_NA,B_NA)
#'
ham = function(mat_1,mat_2){

  if ((is.null(dim(mat_1)))&(!is.null(dim(mat_2)))){

    mat_1 = matrix(mat_1,nrow = 1)

  } else if ((!is.null(dim(mat_1)))&(is.null(dim(mat_2)))){

    mat_2 = matrix(mat_2,nrow = 1)

  } else if ((is.null(dim(mat_1)))&(is.null(dim(mat_2)))){

    mat_1 = matrix(mat_1,ncol = 1)
    mat_2 = matrix(mat_2,ncol = 1)

  } else {}

  d_fun = function(x_1, x_2) (sum(x_1 != x_2,na.rm=TRUE)/length(x_1))*length(x_1)/sum(is.na(x_1)+is.na(x_2)==0)
  matr  = apply(mat_2,1,function(x) sapply(1:nrow(mat_1),function(j) d_fun(x,mat_1[j,])))

  return(matr)

}
