
#' average_distance_to_closest()
#' 
#' This intermediate function of the OT algorithm computes the cost matrix where each cell symbolize the difficulty of moving from a given modality of the target variable in the 1st database (A) to another given modality of the target variable in the 2nd database (B)
#' These costs are assessed by calculating the average distances between covariations of given part of individuals (the percent closest neighbors) classified according to their belonging to modalities of the same target variable differently encoded in A and B
#'
#' @param inst An object corresponding to the output of the \code{Instance} function
#' @param percent_closest A value between 0 and 1 (by default) corresponding to the desired \code{percent closest} indivisuals taken in the computation of the distances
#'
#' @return A list of 3 matrix is returned:
#' \item{Davg}{The cost matrixwhich number of rows corresponds to the number of levels of the target variable in database A, and which number of columns corresponds to the number of levels of the target variable in database B.
#' This cost matrix can be interpreted as the cost of moving from one level of the target factor in the 1st database A to one level of the target factor in the 2nd database B.
#' Davg[P,Q] so refers to the distance between covariates vectors of individuals in database A having modality P and individuals in database B having modality Q}
#' \item{DindivA}{A matrix which number of rows corresponds to the number of rows of the 1st database A and the number of columns corresponds to the number of levels of the target variable in the 2nd database B.
#' DindivA[i,Q] refers to the average distance of the i^th individual of the 1st database and a given proportion of individuals (\code{percent_closest} set by the user) of the 2nd database having the modality Q as value of the target variable}
#' \item{DindivB}{A matrix which number of rows corresponds to the number of rows of the 2nd database B and the number of columns corresponds to the number of levels of the target variable in the 1st database A.
#' DindivB[k,P] refers to the average distance of the k^th individual of the 1st database and a given proportion of individuals (\code{percent_closest} set by the user) of the 2nd database having the modality P as value of the target variable}
#' 
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer 
#' \email{gregory.guernec@@inserm.fr}
#' 
#' @references
#' Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' 0, 20180106 (2019) | \url{https://doi.org/10.1515/ijb-2018-0106}
#' 
#' @seealso \code{\link{Instance}}
#'
#' @example
#' # See the example of the merge_dbs() function to obtain the soluc1 object
#' try1 = prep_dbs(soluc1$DB_READY,nominal = 6,ordinal = 1:5) 
#' res1 = Instance(try1)   # Using the Manhattan distance
#' 
#' res2 = average_distance_to_closest(res1,0.90) 

average_distance_to_closest=function(inst, percent_closest = 1){
  
  if (!is.list(inst)){
    
    stop("This object must be a list returned by the Instance function")
    
  } else {}
  
  if ((percent_closest>1)|(percent_closest <= 0)){
    
    stop("Incorrect value for the percent_closest option")
    
  } else {}
  
  
  # Redefine A and B for the model
  A = (1:inst$nA)
  B = (1:inst$nB)
  Y = inst$Y
  Z = inst$Z
  indY = inst$indY
  indZ = inst$indZ
  
  # Compute average distances
  Davg     = matrix(0,length(Y),length(Z));
  DindivA  = matrix(0,inst$nA,length(Z));
  DindivB  = matrix(0,inst$nB,length(Y));
  
  for (y in Y){
    for (i in indY[[y]]){
      for (z in Z){
        nbclose = max(round(percent_closest*length(indZ[[z]])),1);
        
        distance     = sort(inst$D[i,indZ[[z]]])
        DindivA[i,z] = sum(distance[1:nbclose])/nbclose;
        Davg[y,z]    = Davg[y,z] + sum(distance[1:nbclose])/nbclose/length(indY[[y]])/2.0;
      }
    }
  }
  for (z in Z){
    for (j in indZ[[z]]){
      for (y in Y){
        nbclose      = max(round(percent_closest*length(indY[[y]])),1);
        
        distance     = sort(inst$D[indY[[y]],j])
        DindivB[j,y] = sum(distance[1:nbclose])/nbclose;
        Davg[y,z]    = Davg[y,z] + sum(distance[1:nbclose])/nbclose/length(indZ[[z]])/2.0;
      }
    }
  }
  
  return(list(Davg=Davg, DindivA=DindivA, DindivB=DindivB))
}



