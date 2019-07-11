
#' family_part()
#'
#' Function that gives the power set of a set
#'
#' @param n integer. The cardinal of the set
#' @param ordinal boolean. If TRUE the power set is only composed of subsets of consecutive elements, FALSE (by default) otherwise
#'
#' @return A list of (2^n -1) subsets (The empty set is exclude)
#' @export
#'
#' @examples
#' fam  =  family_part(4)
#' fam2 =  family_part(4,ordinal = TRUE)

family_part = function(n, ordinal = FALSE){

  E        = 1:n

  BDD = rep(rep(1:n,each=n^(n-1)),n^0)
  for (k in 1:(n-1)){

    BDD = cbind(BDD,rep(rep(1:n,each=n^(n-(k+1))),n^k))

  }

  recup = list()

  for (i in 1:(n^n)){

    recup[[i]] = unique(BDD[i,])

  }

  recup2 = unique(lapply(recup,sort))


  if (ordinal == TRUE){

    recup3 = lapply(recup2,function(x) setdiff(E,x))

    recup2bis = lapply(recup2,diff)
    recup3bis = lapply(recup3,diff)

    indic2 = sapply(recup2bis,function(x) any(x>1))
    indic3 = sapply(recup2bis,function(x){length(x)==n-1})


    recup_new = recup2[(indic2==F) & (indic3==F)]
    recup_new[[length(recup_new)+1]] = 1:n



  } else {

    recup_new = recup2

  }

  return(recup_new)

}
