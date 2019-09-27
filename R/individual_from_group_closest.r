
#' individual_from_group_closest()
#' 
#' Assuming that Y and Z summarize a same information encoded in two distinct forms stored in two distinct databases A and B.
#' Using an estimation of the joint probability (Y,Z),this function sequentially assigns the missing modalities of the target variable in each encoding of the two databases. 
#' It corresponds to the execution of a nearest neighbor procedure as described in Algorithm 1 of the  reference article.
#' Please, consult the reference article for more details about this implemented method.
#'
#' @param inst An object corresponding to the output of the \code{Instance} function
#' @param jointprobaA A matrix which number of columns equals the number of modalities of the target variable Y in database A, and which number of rows equals the number of modalities of Z in database B. It gives an estimation of the joint probability (Y,Z) in DB A.
#' The sum of cells of this matrix must be equal to 1
#' @param jointprobaB A matrix which number of columns equals the number of modalities of the target variable Y in database A, and which number of rows equals the number of modalities of Z in database B. It gives an estimation of the joint probability (Y,Z) in DB B.
#' The sum of cells of this matrix must be equal to 1 
#' @param percent_closest A value between 0 and 1 (by default) corresponding to the desired \code{percent closest} indivisuals taken in the computation of the distances
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer 
#' \email{gregory.guernec@@inserm.fr}
#' 
#' @seealso \code{\link{Instance}},\code{\link{average_distance_to_closest}}
#'
#' @return A list of two vectors of numeric values:
#' \item{YAtrans}{A vector corresponding to the predicted values of Y in database B using the Optimal Transportation theory}
#' \item{ZBtrans}{A vector corresponding to the predicted values of Z in database A using the Optimal Transportation theory}
#' 
#' @references
#' Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' 0, 20180106 (2019) | \url{https://doi.org/10.1515/ijb-2018-0106}
#'
#' @example
#' # See the example of the merge_dbs() function to obtain the soluc1 object
#' try1 = prep_dbs(soluc1$DB_READY,nominal = 6,ordinal = 1:5) 
#' res1 = Instance(try1)                          # Using the Manhattan distance
#' 
#' 
#' ### Y and Z are the same variable encoded in 2 different forms (7 levels for Y and 4 levels for Z) stored in two distinct databases, A and B, respectively
#' ### The marginal distribution of Y in the database B is so unknown, as the marginal distribution of Z in the database A ...
#' 
#' # Empirical distribution of Y in database A
#' freqY = prop.table(table(try1$Y1)); freqY
#' 
#' # Empirical distribution of Z in database B
#' freqZ = prop.table(table(try1$Y2)); freqZ
#' 
#' # By supposing that the following matrix called transport symbolizes an estimation of the joint distribution L(Y,Z) ...
#' 
#' transport1 = matrix(c(0.19222222,0,0,0,0.20555556,0,0,0,0.01704523,0.1718437,0,0,0,0.1655457,0.01000983,0,0,0,0.10333333,0,0,0,0.08190020,0.005877581,0,0,0,0.046666667),ncol = 4,byrow = TRUE)
#' 
#' # ... So that the marginal distributions of this object corresponds to freqY and freqZ:
#' apply(transport1,1,sum)  # = freqY
#' apply(transport1,2,sum)  # = freqZ
#' 
#' # The affectation of the predicted values of Y in database B and Z in database A are stored in the following object:
#' 
#' res3      = individual_from_group_closest(res1,jointprobaA = transport1, jointprobaB = transport1, percent_closest= 0.90)
#' summary(res3)

individual_from_group_closest=function(inst, jointprobaA, jointprobaB, percent_closest = 1.0){
  
  if (!is.list(inst)){
    
    stop("This object must be a list returned by the Instance function")
    
  } else {}
  
  if ((!is.matrix(jointprobaA))|(!is.matrix(jointprobaA))){
    
    stop("The joint distributions must be store in matrix objects")
    
  } else {}  
  
  if ((ncol(jointprobaA) != ncol(jointprobaB))|(nrow(jointprobaA) != nrow(jointprobaB))){
    
    stop("The joint distributions must be store in matrix of same size")
    
  } else {} 
  
  
  if ((format(sum(jointprobaA)) != "1")| (format(sum(jointprobaB)) != "1")){
    
    stop("The sum of the jointprobaA matrix or the sum of the jointprobaB matrix differs from 1 !")
    
  } else {}   
  
  
  if ((percent_closest>1)|(percent_closest <= 0)){
    
    stop("Incorrect value for the percent_closest option")
    
  } else {}
  
  
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
  
  return(list(YAtrans = YAtrans, ZBtrans = YBtrans))
}    

