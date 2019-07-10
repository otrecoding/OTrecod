#' try_group()
#'
#' @param Y1 outcomes first database
#' @param Y2 outcomes second database
#' @param ordin family_part option
#'
#' @return list
#' @export
#'
#' @examples
#' Z1 = as.factor(sample(1:3,50,replace = TRUE)); length(Z1)
#' Z3 = as.factor(sample(1:2,50,replace = TRUE)); length(Z3)
#' Z2 = as.factor(sample(c("A","B","C","D"),50, replace = TRUE)); length(Z2)
#' Z4 = as.factor(sample(c("A","B","C","D","E"),50, replace = TRUE)); length(Z4)
#'
#' a1 = try_group(Z1,Z2,ordin = FALSE)
#' summary(a1); a1[[1]]
#' aa = try_group(Z3,Z2,ordin = FALSE)
#' bb = try_group(Z4,Z2,ordin = TRUE)
#' cc = try_group(Z4,Z1,ordin = TRUE)
#'
try_group = function(Y1,Y2,ordin = FALSE){

  ny1        = length(levels(Y1))
  ny2        = length(levels(Y2))

  n          = max(ny1,ny2)

  if (ny1 > ny2){  Y = Y1  } else {  Y = Y2  }

  fam_new    =  family_part(n,ordinal = ordin)


  # Decomposition possibilities

  nbpos     = lapply(count_pos(n,min(ny1,ny2)),sort,decreasing = TRUE)
  nb        = length(nbpos)
  nb2       = length(nbpos[[1]])

  # Constitution of the lists

  list_grp   = list()
  list_ind   = list()

  length_grp = sapply(fam_new,length)

  maxin      = (max(length_grp)-1)
  indic      = 1:maxin


  for (k in 1:length(indic)){

    list_grp[[k]] = fam_new[length_grp == indic[k]]
    list_ind[[k]] = (1:length(fam_new))[length_grp == indic[k]]

  }

  list_glob  = list()
  mat_coord  = matrix(,ncol = min(ny1,ny2))

  for (k in 1:nb){

    list_glob[[k]] = find_coord(unlist(nbpos[[k]]),list_grp,list_ind)
    mat_coord      = rbind(mat_coord,list_glob[[k]])

  }

  mat_coord  = mat_coord[-1,]
  ###
  fam2 = list()

  for (k in 1:length(fam_new)){
    fam2[[k]] = vector(length = length(fam_new[[k]]))

    for (i in 1:length(fam_new[[k]])){
      fam2[[k]][i] = levels(Y)[fam_new[[k]][i]]
    }
  }

  fam_new = fam2
  ###

  fam_new2   = lapply(fam_new,paste,collapse=" ")

  combin = vector(length = nrow(mat_coord))

  for (i in 1:nrow(mat_coord)){

    combin1 = fam_new2[[mat_coord[i,1]]]

    for (j in 2:nb2){

      combin1 = paste(combin1,fam_new2[[mat_coord[i,j]]],sep="/")
    }

    combin[i] = combin1

  }

  row.names(mat_coord) = combin

  return(list(COORD_COMBI = mat_coord,PART = fam_new))

}
