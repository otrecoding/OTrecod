
#' average_distance_to_closest(inst, percent_closest)
#'
#' @param inst todo list
#' @param percent_closest todo list
#'
#' @return todo list
#'
# @examples
average_distance_to_closest=function(inst, percent_closest){
  
  # Redefine A and B for the model
  A = (1:inst$nA)
  B = (1:inst$nB)
  Y = inst$Y
  Z = inst$Z
  indY = inst$indY
  indZ = inst$indZ
  
  # Compute average distances as described in the above
  Davg = matrix(0,length(Y),length(Z));
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
