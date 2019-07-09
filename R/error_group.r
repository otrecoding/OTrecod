
#' error_group()
#'
#' @param REF A factor with a reference number of levels
#' @param Z A factor with a number of levels greater than the number of levels of the reference
#' @param ord A boolean. If TRUE, only the subsequent levels of Z can be grouped together
#'
#' @return A data.frame ot 2 columns:
#'         The first column enumerates all possible groups of modalities of Z to obtain the same number of levels as the reference
#'         The second columns gives the corresponding rate error of reclassification
#' @export
#'
#' @examples
#' Z1 = as.factor(sample(1:3,50,replace = TRUE)); length(Z1)
#' Z3 = as.factor(sample(1:2,50,replace = TRUE)); length(Z3)
#' Z2 = as.factor(sample(c("A","B","C","D"),50, replace = TRUE)); length(Z2)
#' Z4 = as.factor(sample(c("A","B","C","D","E"),50, replace = TRUE)); length(Z4)
#' error_group(Z1,Z4)
#' error_group(Z3,Z1,FALSE)
#'
error_group = function(REF,Z,ord = TRUE){

  cc = try_group(Z,REF,ordin = ord)

  error_g = vector(length = nrow(cc[[1]]))

  for (k in 1:nrow(cc[[1]])){

    Zbis = as.character(Z)

    for (j in 1:ncol(cc[[1]])){

      Zbis[Zbis %in% cc[[2]][[cc[[1]][k,j]]]] = levels(REF)[j]

    }

    Zbis = as.factor(Zbis)

    error_g[k] = 100 - round(sum(diag(table(REF,Zbis)))*100/sum(table(REF,Zbis)),1)
    error_combi    = data.frame(combi = row.names(cc[[1]]),error_rate= error_g)
    error_combi    = error_combi[sort.list(error_combi[,2]),]

  }

  return(error_combi)

}

