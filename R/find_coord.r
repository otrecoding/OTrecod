find_coord = function(all_ind,all_list,all_coord){


  list_glob = list_glob2 = list()

  eff1 = 1

  for (i in 1:length(all_ind)){

    list_glob[[i]]  = all_list[[all_ind[i]]]
    list_glob2[[i]] = all_coord[[all_ind[i]]]
    eff1 = c(eff1,length(list_glob[[i]]))

  }

  nef     = length(eff1)
  eff2    = rev(eff1)[-length(eff1)]
  BDD     = rep(rep(list_glob[[1]],each=prod(eff2[-length(eff2)])),prod(eff1[1]))
  BDD_ind = rep(rep(list_glob2[[1]],each=prod(eff2[-length(eff2)])),prod(eff1[1]))

  for (k in 2:(length(eff1)-1)){
    BDD     = cbind(BDD,rep(rep(list_glob[[k]],each=prod(eff1[(2+k):nef])),prod(eff1[1:k])))
    BDD_ind = cbind(BDD_ind,rep(rep(list_glob2[[k]],each=prod(eff1[(2+k):nef])),prod(eff1[1:k])))
  }


  sol1 = sort(unlist(BDD[1,]))
  sol2 = unlist(BDD_ind[1,])
  sol3 = unlist(BDD[1,])
  for (k in 2:nrow(BDD)){
    sol1 = rbind(sol1,sort(unlist(BDD[k,])))
    sol2 = rbind(sol2,unlist(BDD_ind[k,]))
    sol3 = rbind(sol3,unlist(BDD[k,]))
  }

  exemp = 1:(length(all_list)+1)

  test1 = apply(sol1,1,function(x){setequal(x,exemp)})
  sol_coord = (1:nrow(sol1))[test1==TRUE]

  sol4 = t(apply(sol2[test1 == TRUE,],1,sort))
  sol4 = sol4[!duplicated(sol4),]

  sol5 = sol3[sol_coord,]

  return(sol4)
}

