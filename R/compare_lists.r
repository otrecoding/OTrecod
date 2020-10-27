
#' compare_lists()
#'
#'    This function compares the elements of two lists of same length.
#'
#' @param listA a first list
#' @param listB a second list
#'
#' @return A boolean vector of same length as the two lists,
#' which ith element is \code{TRUE} if the ith element is different
#' between the 2 lists, or \code{FALSE} otherwise
#'
#' @author Gregory Guernec
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @export
#'
#' @aliases compare_lists
#'
#' @examples
#' data1 = data.frame(Gender = rep(c("m","f"),5),Age = rnorm(5,20,4))
#' data2 = data.frame(Gender = rep(c("m","f"),5),Age = rnorm(5,21,5))
#'
#' list1 = list(A = 1:4, B = as.factor(c("A","B","C")), C = matrix(1:6,ncol = 3))
#' list2 = list(A = 1:4, B = as.factor(c("A","B")), C = matrix(1:6,ncol = 3))
#' list3 = list(A = 1:4, B = as.factor(c("A","B","C")),C = matrix(c(1:5,7),ncol = 3))
#' list4 = list(A = 1:4, B = as.factor(c("A","B","C")), C = matrix(1:6,ncol = 2))
#' list5 = list(A = 1:4, B = as.factor(c("A","B")), C = matrix(1:6,ncol = 2))
#' list6 = list(A = 1:4, B = as.factor(c("A","B")), C = data1)
#' list7 = list(A = 1:4, B = as.factor(c("A","B")), C = data2)
#'
#' OTrecod::compare_lists(list1,list2)
#' OTrecod::compare_lists(list1,list3)
#' OTrecod::compare_lists(list1,list4)
#' OTrecod::compare_lists(list1,list5)
#' OTrecod::compare_lists(list6,list7)

compare_lists = function(listA,listB){

  stopifnot(is.list(listA), is.list(listB), length(listA) == length(listB))

  return(as.vector(!mapply(function(a, b) identical(as.matrix(a), as.matrix(b)), listA, listB)))

}




