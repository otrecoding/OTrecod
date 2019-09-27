
#' compare_lists()
#'
#'    This function compares the elements of two lists of same lengths
#'    The order of entry of the 2 lists have no importance
#'
#' @param listA A first list
#' @param listB A second list
#'
#' @return A boolean vector of same length as the two lists,
#' which ith element is \code{TRUE} if the ith element is different
#' between the 2 lists, or \code{FALSE} otherwise
#'
#' @author Gregory Guernec
#' \email{gregory.guernec@@inserm.fr}
#' @export
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
#' compare_lists(list1,list2)
#' compare_lists(list1,list3)
#' compare_lists(list1,list4)
#' compare_lists(list1,list5)
#' compare_lists(list6,list7)

compare_lists = function(listA,listB){

  if ((!is.list(listA))|(!is.list(listB))){

    stop("At least one of your two objects is not a list object!")} else {}


  nA  = length(listA); nB = length(listB)
  vec = vector(length = nA)

  if (nA == nB){

    for (i in 1:nA){

      if (class(listA[[i]])!= class(listB[[i]])){

        vec[i] = TRUE

      } else if (class(listA[[i]]) %in% c("matrix","data.frame")){

        if (any(dim(listA[[i]])!= dim(listB[[i]]))){

          vec[i] = TRUE

        } else if (all(dim(listA[[i]])== dim(listB[[i]]))){

          vec[i] = ifelse(any(listA[[i]]!= listB[[i]]),TRUE,FALSE)

        } else {}


      } else if (!(class(listA[[i]]) %in% c("matrix","data.frame"))){

        vec[i] = ifelse(length(listA[[i]])!= length(listB[[i]]),TRUE,
                        ifelse(any(listA[[i]]!= listB[[i]]),TRUE,FALSE))

      }  else {}

    }

  } else {stop("Your 2 lists have different lengths !!!")}

  return(vec)

}




