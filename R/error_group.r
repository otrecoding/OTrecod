
#' error_group()
#' 
#' By supposing that REF and Z are two factors (ordered of not) with a different number of levels (n_REF and n_Z respectiveley) such as n_Z > n_REF,
#' This function studies similarities between REF and Z when the levels of Z are grouped to obtain a number of levels equals to the number of levels of Z.
#' 
#' Applying this function to the results of the data integration
#'
#' @param REF A factor with a reference number of levels
#' @param Z A factor with a number of levels higher than the number of levels of the reference number
#' @param ord A boolean. If TRUE, only the subsequent levels of Z can be grouped together
#'
#' @return A data.frame ot 2 columns:
#' \item{grouping_levels}{Enumerates all possible groups of levels of Z to obtain the same number of levels as the reference}
#' \item{error_rate}{Gives the corresponding rate error after specific groupings of Z levels in n_REF levels}    
#' 
#' @export
#' 
#' @author Gregory Guernec
#' \email{gregory.guernec@@inserm.fr}
#' 
#' @seealso \code{\link{count_pos}}, \code{\link{find_coord}}, \code{\link{family_part}}, \code{\link{try_group}}
#'
#' @examples
#' Z1 = as.factor(sample(1:3,50,replace = TRUE)); length(Z1)
#' Z2 = as.factor(sample(c("A","B","C","D"),50, replace = TRUE)); length(Z2)
#' Z3 = as.factor(sample(1:2,50,replace = TRUE)); length(Z3)
#' Z4 = as.factor(sample(c("A","B","C","D","E"),50, replace = TRUE)); length(Z4)
#' error_group(Z1,Z4)
#' error_group(Z3,Z1,FALSE)
#' 
#' ## Using the "try1" object from the OT function example 
#' ## For database A, we have:
#' head(try1$DATA1_OT)
#' ##  ... where the Y variable summarizes an information in 7 levels and Z summarizes the same information in 4 unknown levels 
#' ## predictions using Optimal Transportation Theory are stored in the "OTpred" variable in a database A. So by doing:
#' 
#' Ypred = as.factor(try1$DATA1_OT$OTpred)
#' groups_similarities = error_group(Ypred,try1$DATA1_OT$Y, ord = TRUE)
#' 
#' 

error_group = function(REF,Z,ord = TRUE){
  
  # REF = as.factor(REF)
  # Z   = as.factor(Z)
  
  if ((!(is.factor(REF)))|(!(is.factor(Z)))){
    
    stop("REF and Z must be factors")
    
  } else {}
  
  if (!(is.logical(ord))){
    
    stop("The ord option is boolean: TRUE or FALSE expected")
    
    
  } else {}
  

  # cc = try_group(Z,REF,ordin = ord)
  cc = try_group(REF,Z,ordin = ord)

  error_g = vector(length = nrow(cc[[1]]))

  for (k in 1:nrow(cc[[1]])){

    Zbis = as.character(Z)

    for (j in 1:ncol(cc[[1]])){

      Zbis[Zbis %in% cc[[2]][[cc[[1]][k,j]]]] = levels(REF)[j]

    }

    Zbis = as.factor(Zbis)

    error_g[k] = 100 - round(sum(diag(table(REF,Zbis)))*100/sum(table(REF,Zbis)),1)
    error_combi    = data.frame(grouping_levels = row.names(cc[[1]]),error_rate= error_g)
    error_combi    = error_combi[sort.list(error_combi[,2]),]

  }

  return(error_combi)

}

