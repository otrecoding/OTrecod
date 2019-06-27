
#' OT
#' 
#' An original function for data integration using Optimal Transport theory
#'
#' @param DB_READY A data.frame on specific form obtained after transfo_dist() execution
#' @param which.dist A character (with quotes) corresponding to the chosen distance to evaluate similarities in OT. 
#'                   Three distances are actually available: The Gower distance adapted for mixed data ("G" by default),
#'                   The Manhattan distance ("M") specially useful with qualitative covariates, or the Euclidean Distance ("E") 
#'                   for continuous raw or transformed covariates
#' @param which.datab A character string. If "BOTH" (Default), the OT algorithm is used on the targets variable of the 2 stacked DBs. If "DB1" the imputation is only done on the target of the 1st DB.
#'                    If "DB2" the imputation is only done on the target of the 2nd DB
#' @param group_test A character string with quotes. "YES" (Default) if "BOTH" is selected in the which.datab option, compare the proximity between the repartition of the theorical target and those assessed with OT 
#' @param rate.resample A scalar or vector of proportions of training sets for estimating the accuracy of OT. Only if "BOTH" is selected in the which.datab option
#' @param repeat.resample An integer or a vector of integers (same length as rate.resample), corresponding to the numbers of repetitions for each resampling previously declared
#'
#' @return A list of 5 elements:
#'         COST_solution The estimated gamma corresponding to the optimal joint distribution of the 2 targets. Solution from the cost function.
#'         OT_AFFECT A list of 2 elements corresponding to the OT imputation in one of the 2 databases
#'         DIST A character corresponding to the distance function selected for OT
#'         ACCURACY_RESAMP A data.frame summaryzing the accuracies global and by database (mean accuracies if many repetitions by selected rates) of each resampling selected
#'         GRP_SCALE A data.frame summarizing the error rates of confusion matrix between the targets in same numbers of levels
#' @export
#'
#' @examples
#' soluc      = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
#' prep       = transfo_dist(soluc[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")
#' res_OT     = OT(prep,which.dist = "G",which.datab = "BOTH",group_test = "YES",rate.resample = c(0.40,0.30,0.10),repeat.resample = c(1,1,2))
#' summary(res_OT)

OT = function(DB_READY,which.dist = "G",which.datab = "BOTH",group_test = "NO",rate.resample = NULL,repeat.resample = NULL){


  stopifnot(is.data.frame(DB_READY))
  stopifnot(which.dist  %in% c("G","M","E"))
  stopifnot(which.datab %in% c("DB1","DB2","BOTH"))
  stopifnot(group_test  %in% c("NO","YES"))
  
  if (length(repeat.resample)!= length(rate.resample)){
    
    stop("repeat.resample and rate.resample options have different lengths")
    
  } else {} 
  
  
  if ((which.datab %in% c("DB1","DB2"))&(length(rate.resample)!=0)){
    
    stop("Accuracy is available on training set only if the 2 DBs are imputed, see which.datab option for more details.")
    
  } else {}
  
 # OT ALGORITHM
 #--------------
   
  cat("OT in progress. Please wait ...","\n")
  
  xx         = cost(DB_READY,which.dist)  
  yy         = affectation(xx$solution,DB_READY,dist_choice = which.dist,which.DB = which.datab) 
  

 # PROXIMITY OF THE ESTIMATED REPARTITIONS BY RESCALING
 #------------------------------------------------------
  
  if (group_test == "YES"){

    if (which.datab %in% c("DB1","BOTH")){
      
      zz     = yy[[1]]
      indcol = 2
      
    } else if (which.datab == "DB2"){
      
      zz     = yy[[2]]
      indcol = 3
      
    } else {}
      
      
    if(length(levels(as.factor(zz$OT)))<length(levels(as.factor(zz[,indcol])))){
          
          YA = zz$OT
          YB = zz[,indcol]
          
     } else {
          
          YA = zz[,indcol]
          YB = zz$OT
        
     }
        
     err_gp = error_group(YA,YB)
     
     cat("-------------","\n")
     cat("SCALING GROUPED: ok","\n")
     
    
  } else {err_gp = NULL}
 
  
  # ACCURACY OF THE METHOD BY RESAMPLING
  #--------------------------------------
  
  Nrate = length(rate.resample)
  
  
  if (Nrate!=0){
  
    cat("-------------","\n")
    cat("TRAINING ACCURACY in progress ...","\n")
    
    accur_res = rep(NA,11)
    
    for (k in  1:Nrate){
      
      print(paste(k,"/",Nrate))
      
      accur_test = rep(NA,8)
      
      for (j in 1:repeat.resample[k]){
        
          accur_test = rbind(accur_test,training_accuracy(yy,rate.resample[k])) 
      
      }
      
      accur_mean = apply(accur_test,2,mean,na.rm= TRUE)
      accur_sd   = apply(accur_test,2,sd,na.rm= TRUE)
      
      accur_res  = rbind(accur_res,c(accur_test[2,1:2],accur_mean[3],accur_sd[3],accur_test[2,4],accur_mean[5],accur_sd[5],
                         accur_test[2,6],accur_mean[7],accur_sd[7],repeat.resample[k]))
      

    }
    
    accur_res           = accur_res[-1,]
    accur_res           = round(accur_res,2)
    accur_res           = as.data.frame(accur_res)
    colnames(accur_res) = c("RATE_RESAMPLE","N","AC_GLOB_MEAN","AC_GLOB_SD","N1","AC_Y1_MEAN","AC_Y1_SD","N2","AC_Y2_MEAN","AC_Y2_SD","REP")
    accur_res = accur_res[sort.list(accur_res[,1]),]
    
  } else {accur_res = NULL}

  cat("TRAINING ACCURACY: ok","\n")
  
  
  return (list(COST_solution = xx$solution,OT_AFFECT = yy,DIST = which.dist,ACCURACY_RESAMP = accur_res,GRP_SCALE = err_gp))  

}  



  
















