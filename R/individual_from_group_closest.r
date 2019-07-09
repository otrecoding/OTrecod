
individual_from_group_closest=function(inst, jointprobaA, jointprobaB, percent_closest=1.0){
  
  # Redefine A and B for the model
  A = 1:inst$nA
  B = 1:inst$nB
  Y = inst$Y
  Z = inst$Z
  Yobserv = inst$Yobserv
  Zobserv = inst$Zobserv
  indY = inst$indY
  indZ = inst$indZ
  nbindY = numeric(0)
  for (y in Y){
    nbindY = c(nbindY,length(inst$indY[[y]]))}
  nbindZ= numeric(0)
  for (z in Z){
    nbindZ = c(nbindZ,length(inst$indZ[[z]]))}
  freqY= numeric(0)
  for (y in Y){
    freqY = c(freqY,nbindY[y] / length(A))}
  freqZ= numeric(0)
  for (z in Z){
    freqZ = c(freqZ,nbindZ[z] / length(B))}
  
  # In essence, assign to each individual the modality that is closest, 
  # where the distance from an individual to a modality is computed as the 
  # average distance to the individuals having this modality (in the other base)
  YAtrans =numeric(inst$nB);
  YBtrans =numeric(inst$nA);
  average_distance = average_distance_to_closest(inst, percent_closest);
  Davg=average_distance$Davg
  DindivA= average_distance$DindivA
  DindivB = average_distance$DindivB
  # DA = [(z,[sum(inst$D[i,j] for j in indZ[z])/nbindZ[z] for i in A]) for z in Z]
  # DB = [(y,[sum(inst$D[i,j] for i in indY[y])/nbindY[y] for j in B]) for y in Y]
  DA = DB =  list()
  
  for (z in Z){
    
    DA[[z]] = vector(length = length(A))
    
    for (i in A){
      
      DA[[z]][i] = DindivA[i,z]
      
    }
  } 
  
  for (y in Y){
    
    DB[[y]] = vector(length = length(B))
    
    for (i in B){
      
      DB[[y]][i] = DindivB[i,y]
      
    }
  } 
  
  for (y in Y){
    indtrans = indY[[y]]
    
    for (z in Z){
      
      nbtrans = min(round(jointprobaA[y,z]/freqY[y] * nbindY[y]),length(indtrans));
      
      distance = numeric(0)
      for (i in indtrans){
        
        distance = c(distance,DA[[z]][i])
        
      }
      
      names(distance) = indtrans
      distance        = sort(distance)
      ident_dist      = as.numeric(names(distance))
      
      for (k in 1:nbtrans){
        YBtrans[ident_dist[k]] = z;
        indtrans = setdiff(indtrans,ident_dist[k])
      }
    }
    
    
    # affect potential individuals that have not been transported due to
    # rounding
    for (i in indtrans){
      YBtrans[i] = Zobserv[which.min(inst$D[i,])]
    }
  }
  
  for (z in Z){
    indtrans = indZ[[z]]
    
    for (y in Y){
      
      nbtrans = min(round(jointprobaB[y,z]/freqZ[z] * nbindZ[z]), length(indtrans));
      distance = numeric(0)
      for (j in indtrans){
        
        distance = c(distance,DB[[y]][j])            
        
      }
      
      names(distance) = indtrans
      distance        = sort(distance)
      ident_dist      = as.numeric(names(distance))            
      
      for (k in 1:nbtrans){
        YAtrans[ident_dist[k]] = y;
        indtrans = indtrans = setdiff(indtrans,ident_dist[k])
      }
    }
    
    # affect potential individuals that have not been transported due to
    # rounding
    for (j in indtrans){
      YAtrans[j] = Yobserv[which.min(inst$D[,j])]
    }
  }
  
  return(list(YAtrans = YAtrans, YBtrans = YBtrans))
}    

