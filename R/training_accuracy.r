
#' training_accuracy()
#' 
#' A function which gives an estimation of the accuracy of OT on training set of the OT algorithm
#'
#' @param listdata A list, output of the affectation() function.
#' @param rate A double between 0 and 1. It corresponds to the rate of resampling for the accuracy verification.
#' @param seed1 An integer used as argument by the set.seed() for offsetting the random number generator (Random integer by default)
#'
#' @return A vector of 8 doubles. 
#'         RATE_RESAMPLE The rate of database used for the accuracy estimation
#'         N Number of remaining rows in the database
#'         AC_GLOB The global accuracy in training set
#'         N1 Number of remaining rows to predict Y1
#'         AC_Y1 Rate of accuracy with OT to predict Y1
#'         N2 Number of remaining rows to predict Y2
#'         AC_Y2 Rate of accuracy with OT to predict Y2
#'         SEED The random number generator selected
#' @export
#'
#' @examples
#' soluc      = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
#'  prep       = transfo_dist(soluc[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")
#'  xx         = cost(prep,"G")
#'  yy         = affectation(xx$solution,prep,dist_choice = "G",which.DB = "BOTH")
#'  accur_test = training_accuracy(yy,0.30,2014)

training_accuracy = function(listdata,rate,seed1 = sample(1:999999,1)){
  
  if ((rate > 1)|(rate < 0)){
    
    stop("Inadequate rate in input. This value must be chosen between 0 and 1")
    
  } else {}
  
  
  if (!(is.list(listdata))){
    
    stop ("Your object must be a list")
    
  } else {}
  
  if (length(listdata)!=2){
    
    stop("Inappropriate dimension for the input list")
    
  } else {}
  
  if ((is.null(listdata[[1]]))|(is.null(listdata[[2]]))){
    
    stop("At least one of the 2 list is not declared")
    
  } else {}
  
  
  DB1_OT  = listdata[[1]]
  DB2_OT  = listdata[[2]]
  
  indicbd  = c(2,1)
  indicot  = c(3,2)
  
  OTDB     = list()
  Ytheo    = list()
  
  for (k in 1:2){
    
    DB_OT    = listdata[[k]]
    new_size = round(rate*nrow(DB_OT))
    
    set.seed(seed1); samp = sample(1:nrow(DB_OT),new_size,replace = FALSE)
    
    DBOT_new = DB_OT[samp,]
    
    Ytheo[[indicbd[k]]] = as.character(DBOT_new[,k+1])
    
    DBOT_new[,c(1:3,ncol(DBOT_new))] = apply(DBOT_new[,c(1:3,ncol(DBOT_new))],2,as.character)
    
    DBOT_new[,1] = rep(indicbd[k],nrow(DBOT_new))
    
    DBOT_new[,indicot[k]] = DBOT_new[,ncol(DBOT_new)]
    DBOT_new[,k+1] = rep(NA,nrow(DBOT_new))
    
    DBOT_new = DBOT_new[,-ncol(DBOT_new)]
    
    OTDB[[indicbd[k]]] = DBOT_new
    
    
  }
  
  OT_ready       = rbind(OTDB[[1]],OTDB[[2]])
  
  for (j in 2:3){
    
    DB_OT  = listdata[[k-1]]
    
    if (is.ordered(DB_OT[,j])){
      
      OT_ready[,j] = as.factor(OT_ready[,j])
      
    } else {
      
      OT_ready[,j] = as.character(OT_ready[,j])
      
    }
    
  }
  
  prep_new       = transfo_dist(OT_ready,quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")
  xx_new         = cost(prep_new,"G",PRINT1 = FALSE)
  listdata_new   = affectation(xx_new$solution,prep_new,dist_choice = "G",which.DB = "BOTH",PRINT2 = FALSE)
  
  yt1 = as.character(listdata_new[[1]][,ncol(listdata_new[[1]])])
  yt2 = as.character(listdata_new[[2]][,ncol(listdata_new[[2]])])
  
  glob_ok   = sum(diag(table(Ytheo[[1]],yt1))) + sum(diag(table(Ytheo[[2]],yt2)))
  glob_sum  = sum(table(Ytheo[[1]],yt1)) + sum(table(Ytheo[[2]],yt2))
  
  pred_glob = round(glob_ok*100 / glob_sum,2)
  
  predY2    = sum(diag(table(Ytheo[[1]],yt1)))*100 / sum(table(Ytheo[[1]],yt1))
  predY1    = sum(diag(table(Ytheo[[2]],yt2)))*100 / sum(table(Ytheo[[2]],yt2))
  
  return(c(RATE_RESAMPLE = rate,N = length(yt1) + length(yt2),AC_GLOB = pred_glob,N1 = length(yt2),AC_Y1 = predY1,N2 = length(yt1),AC_Y2 = predY2,SEED = seed1))
  
}

