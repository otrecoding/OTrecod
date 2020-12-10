
#' ham()
#'
#' This function computes a matrix distance using the Hamming distance as proximity measure.
#'
#' \code{ham} returns the pairwise distances between rows (observations) of a single matrix if \code{mat_1} equals \code{mat_2}.
#' Otherwise \code{ham} returns the matrix distance between rows of the two matrices \code{mat_1} and \code{mat_2} if this 2 matrices are different in input.
#' Computing the Hamming distance stays possible despite the presence of missing data by applying the following formula. Assuming that A and B are 2 matrices such as \code{ncol(A) = ncol(B)}.
#' The Hamming distance between the \eqn{i^{th}} row of A and the \eqn{k^{th}} row of B equals:
#'
#' \deqn{\mbox{ham}(A_i,B_k) = \frac{\sum_j 1_{\left\{A_{ij} \neq B_{kj}\right\}}}{\sum_j 1}\times\left(\frac{\sum_j 1}{\sum_j 1_{\left\{!\mbox{is.na}(A_{ij}) \& !\mbox{is.na}( B_{kj})\right\}}}\right)}
#'
#' where: \eqn{i = 1,\dots,\mbox{nrow}(A)} and  \eqn{k = 1,\dots,\mbox{nrow}(B)}; And the expression located to the right term of the multiplication corresponds to a specific weigh applied in presence of NAs in \eqn{A_i} and/or \eqn{B_k}.
#'
#' This specificity is not implemented in the \code{cdist} function and the Hamming distance can not be computed using the \code{\link[proxy]{dist}} function either.
#'
#' The Hamming distance can not be calculated in only two situations:
#' \enumerate{
#' \item If a row of A or B has only missing values (ie for each of the columns of A or B respectively).
#' \item The union of the indexes of the missing values in row i of A with the indexes of the missing values in row j of B concerns the indexes of all considered columns.
#' }
#' Example: Assuming that \eqn{\mbox{ncol}(A) = \mbox{ncol}(B) = 3}, if \eqn{A_i = (1,\mbox{NA},0)} and \eqn{B_j = (\mbox{NA},1,\mbox{NA})}, for each column, either the information in row i is missing in A,
#' or the information is missing in B, which induces: \eqn{\mbox{ham}(A_i,B_k) = \mbox{NA}}.
#'
#' If \code{mat_1} is a vector and \code{mat_2} is a matrix (or data.frame) or vice versa, the length of \code{mat_1} must be equal to the number of columns of \code{mat_2}.
#'
#' @param mat_1 a vector, a matrix or a data.frame of binary values that may contain missing data
#' @param mat_2 a vector, a matrix or a data.frame of binary values with the same number of columns as \code{mat_1} that may contain missing data
#'
#' @return A distance matrix
#' @export
#'
#' @author Gregory Guernec
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @references
#' Roth R (2006). Introduction to Coding Theory. Cambridge University Press.
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

  if (sum(is.na(matr))!=0){

    vec1 = as.vector(matr)
    vec1[which(is.na(vec1))] = NA
    matr = matrix(vec1, ncol = ncol(matr))

  } else{}

  return(matr)

}
