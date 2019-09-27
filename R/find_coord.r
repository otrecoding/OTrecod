

#' find_coord()
#'
#' An intermediate function that returns the index of a list from \code{family_part} output to compose the partition of a set
#'
#' @param all_ind A vector of integers from the count_pos() function
#' @param all_list A list of sublists of subsets of a power set (from family_part() by example) composed of 1,2,...,n elements
#' @param all_coord A matrix of integers which lines gives the coordinates of a list from family_part()
#'
#' @return A matrix of indexes of an output list from family_part(), which indexes of each row corresponds to subsets of a \code{family_part} output to give a partition of the starting set
#' @export
#' 
#' @author Gregory Guernec 
#' \email{gregory.guernec@@inserm.fr}
#' 
#' @seealso \code{\link{family_part}}, \code{\link{count_pos}}, \code{\link{try_group}}, \code{\link{error_group}}
#'
#' @examples
#' 
#' ## Suppose a set E of 4 elements. The set of parts of this set is given by:
#' fam      =  family_part(4)
#' 
#' ## We store in a vector the length of each elements of the previous object:
#' size_grp = sapply(fam,length)
#' 
#' ## We want to determine all the partitions of E made of 3 subsets (with no overlapping) 
#' a1       = count_pos(4,3)
#' 
#' The count_pos() function indicates that all the researched paritions will be made of 3 subsets: On of them must be of size 2, the other 2 will be of size 1
#' 
#' grp_list = list()   # The list of possible subsets of E of 1,2 or 3 elements
#' ind_list = list()   # The corresponding indexes of fam
#' ind      = 1:3
#' for (k in 1:3){
#'   grp_list[[k]] = fam[size_grp == ind[k]]
#'   ind_list[[k]] = (1:length(fam))[size_grp == ind[k]]
#' }
#' 
#' ## Each row of the list_coord object gives the indexes of the fam object to make a possible partition of E made of 3 subsets
#' list_coord = find_coord(unlist(a1[[1]]),grp_list,ind_list)
#' ## The 1st row corresponds to the partition of E made of 3 subsets given by: (1),(2),(3,4)
#' ## The 2nd row corresponds to the partition of E made of 3 subsets given by: (1),(2,4),(3) ... And so on ...
#'
find_coord = function(all_ind, all_list, all_coord) {
  
  
  if ((!(is.list(all_list)))|(!(is.list(all_coord)))){
    
    stop("The two options all_list and all_coord must be lists")
    
  } else {}
  
  
  
  list_glob = list_glob2 = list()
  
  eff1 = 1
  
  for (i in 1:length(all_ind)) {
    list_glob[[i]]  = all_list[[all_ind[i]]]
    list_glob2[[i]] = all_coord[[all_ind[i]]]
    eff1 = c(eff1, length(list_glob[[i]]))
    
  }
  
  nef     = length(eff1)
  eff2    = rev(eff1)[-length(eff1)]
  BDD     = rep(rep(list_glob[[1]], each = prod(eff2[-length(eff2)])), prod(eff1[1]))
  BDD_ind = rep(rep(list_glob2[[1]], each = prod(eff2[-length(eff2)])), prod(eff1[1]))
  
  for (k in 2:(length(eff1) - 1)) {
    BDD     = cbind(BDD, rep(rep(list_glob[[k]], each = prod(eff1[(2 + k):nef])), prod(eff1[1:k])))
    BDD_ind = cbind(BDD_ind, rep(rep(list_glob2[[k]], each = prod(eff1[(2 +
                                                                          k):nef])), prod(eff1[1:k])))
  }
  
  
  sol1 = sort(unlist(BDD[1, ]))
  sol2 = unlist(BDD_ind[1, ])
  sol3 = unlist(BDD[1, ])
  for (k in 2:nrow(BDD)) {
    sol1 = rbind(sol1, sort(unlist(BDD[k, ])))
    sol2 = rbind(sol2, unlist(BDD_ind[k, ]))
    sol3 = rbind(sol3, unlist(BDD[k, ]))
  }
  
  exemp = 1:(length(all_list) + 1)
  
  test1 = apply(sol1, 1, function(x) {
    setequal(x, exemp)
  })
  sol_coord = (1:nrow(sol1))[test1 == TRUE]
  
  sol4 = t(apply(sol2[test1 == TRUE, ], 1, sort))
  sol4 = sol4[!duplicated(sol4), ]
  
  sol5 = sol3[sol_coord, ]
  
  return(sol4)
}
