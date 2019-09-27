#' OT()
#' 
#' Implementation of the Optimal Transportation (OT) algorithm for data integration with or without relaxation on the constraints of marginal
#' distributions
#' 
#' Assuming that:
#' 
#' 1. Y and Z summarize a same information of interest encoded in two distinct forms stored in two independent databases A and B (no individuals or rows in common), 
#' 2. The two databases have a set of common covariates X (i.e the same variables stored in the same forms in the two bases),
#' 3. You would like to merge vertically A and B keeping this variable of interest in at least one of its two encodings.
#'  
#' So, this function gives you a possible answer to this recoding problem by using the Optimal Transportation Theory.
#' The models implemented in this function correspond to the OUTCOME and R-OUTCOME models described in the reference article.
#' 
#' The OUTCOME model performs data fusion of the two databases assuming that the transport of the target variable is equally distributed in the two databases of interest (ie L(Y) is equal in A and B, L(Z) is equal in A and B).
#' By default, the OUTCOME model is implemented by fixing the \code{maxrelax} option equals to 0.
#' 
#' When the \code{maxrelax} option > 0, a relaxation of the constraints of marginal is possible and the R-OUTCOME model is performed. 
#' This model allows to override the distributional assumption previously described (Here Y and Z can be differently distributed in A and B).
#' 
#' In this function, these models are implemented in two parts:
#' 1. The 1st part partially resolves the recoding problem using a linear programming model: It gives you an estimator of the joint distributions of Y and Z in database A (and B).
#' 2. The 2nd part uses the result of the 1st part and an independent procedure to gives you the individuals predictions of Y in B and Z in A.
#' 
#' 2 procedures are implemented in this function and can be chosen using the \code{indiv_method} option.
#' The 1st procedure (\code{indiv_method} = "sequential") corresponds to a nearest neighbor algorithm and can make some arbitrary decisions.
#' The 2nd procedure (\code{indiv_method} = "optimal") searchs for the individual predictions (of Z in A by example) that minimize the sum of the indivual distances in A with the modalities of Z in B, by solving an optimization problem (simplex algorithm).
#' Be a little patient please, the time execution of this last procedure may be quite long.
#' 
#'
#' @param datab A data.frame with a specific order of columns. The 1st comumn is a key column for individual identification,
#' the second and third columns corresponds to the target variable encoded with different number of levels in databases A(called Y) and B(called Z) (The levels must be beforehand convert in numeric, started from 1),
#' the following columns corresponds to common covariates (X matrix) between the 2 databases whatever their number, and whatever their order. Nevertheless, covariates must be factors beforehand convert in  numeric
#' @param nominal Vector of index of columns corresponding to nominal variables
#' @param ordinal Vector of index of columns corresponding to ordinal variables
#' @param logic Vector of index of columns corresponding to boolean variables
#' @param percent_c Percent of closest neighbors taken in the computation of the costs
#' @param maxrelax Maximum percentage of deviation from expected probability masses (0 for the OUTCOME model, a non-zero value otherwise). Please consult the reference article for more details
#' @param norm A value (with no quotes) corresponding to the distance function chosen to calculate the distances between covariates. Equals to 0 for the Hamming distance, 1 for the Manhattan distance (by default), or 2 for the euclidean distance
#' @param indiv_method A string of characters that specifies the method used to get individual transport from group joint probabilities, "sequential" or "optimal". Please consult the reference article for more details
#'
#' @return A list containing 7 elements:
#'     \item{TIME_EXE}{Running time of the function}
#'     \item{TRANSPORT_A}{Cost matrix corresponding to an estimation (gamma, see reference for more details) to the joint distribution of (YA,ZA)} 
#'     \item{TRANSPORT_B}{Cost matrix corresponding to an estimation to the joint distribution of (YB,ZB)} 
#'     \item{estimatorZA}{Estimates of the probability distribution of Z conditional to X and Y in database A from predictions stored in an array. The number of rows of each table corresponds to the total number of profiles of covariates.
#'     The number of columns of each table corresponds to the number of levels of Y. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of Z.}
#'     \item{estimatorYB}{Estimates of the probability distribution of Y conditional to X and Z in database B from predictions stored in an array. The number of rows of each table corresponds to the total number of profiles of covariates.
#'     The number of columns of each table corresponds to the number of levels of Z. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of Y.}
#'     \item{DATA1_OT}{database A with imputed individual prediction on Z using OT}
#'     \item{DATA2_OT}{database B with imputed individual prediction on Y using OT}
#'     
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer 
#' \email{gregory.guernec@@inserm.fr}
#' 
#' @references
#' Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' 0, 20180106 (2019),\url{https://doi.org/10.1515/ijb-2018-0106}
#' 
#' 
#' @seealso \code{\link{merge_dbs}},\code{\link{OT_joint}},\code{\link{Instance}}, \code{\link{average_distance_to_closest}}, \code{\link{individual_from_group_closest}}
#' 
#' @import rdist dplyr ROI ROI.plugin.glpk                           
#' @importFrom ompr sum_expr
#' @importFrom ompr.roi with_ROI
#' @export
#'
#'@examples
#'
#' # See the example of the merge_dbs() function to obtain the soluc1 object
#' 
#' 
#' ### An example of OUTCOME model with nearest neighbor procedure
#' try1 = OT(soluc1$DB_READY,nominal = 6,ordinal = 1:5, percent_c = 0.90, maxrelax = 0, norm = 1,indiv_method="sequential")
#' summary(try1)
#' head(try1$estimatorZA[,,1])
#' ### ... Corresponds to P[Z = 1|Y,age,hsize,sex_2]
#' ### The 1st line of this table corresponds to P[Z = 1|Y,age = (60,100], hsize = 1, sex_2 = 1)] because the 1st rox names indicates "311".
#' ### The "3" corresponds to the third class of "age", the follwing "1" corresponds to the class of "hsize", and the last "1" corresponds to the value of the binary sex_2 variable. 
#' ### The columns corresponds to the possible levels for Y so we can read that:
#' ### P[Z = 1|Y = (35,Inf], age = (60,100], hsize = 1, sex_2 = 1)]] = 1.
#' ### Thus, we can conclude that all individuals with this specific profile of covariates will be affected to the 1st level of Z in database A.
#' ### The reasoning will be the same for the estimatorYB object.
#' 
#' 
#' ### An example of OUTCOME model using an optimization procedure for the individual predictions
#' try2 = OT(soluc1$DB_READY,nominal = 6,ordinal = 1:5, percent_c = 0.90, maxrelax = 0, norm = 1,indiv_method="optimal")
#' summary(try2)
#' 
#' ### An example of R-OUTCOME model with nearest neighbor procedure
#' try3 = OT(soluc1$DB_READY,nominal = 6,ordinal = 1:5, percent_c = 0.90, maxrelax = 1, norm = 1,indiv_method="sequential")
#' summary(try3)
#' 
#' An example of R-OUTCOME model using an optimization procedure for the individual predictions
#' try4 = OT(soluc1$DB_READY,nominal = 6,ordinal = 1:5, percent_c = 0.90, maxrelax = 1, norm = 1,indiv_method="optimal")
#' # The following `from` values were not present in `x`: 3    
#' summary(try4)

OT = function(datab, nominal = NULL, ordinal = NULL,logic = NULL, percent_c = 1.0, maxrelax = 0, norm = 1, indiv_method){
  
  
  cat("---------------------------------------","\n")
  cat("OT PROCEDURE in progress ..."           ,"\n")
  cat("---------------------------------------","\n")
  cat("Type                  = ", ifelse(maxrelax == 0,"OUTCOME","R-OUTCOME"),"\n")
  cat("Distance              = ", ifelse(norm == 0,"Hamming",ifelse(norm == 1,"Manhattan","Euclidean")),"\n")
  cat("Percent closest       = ", 100.0*percent_c, "%","\n")
  cat("Relaxation term       = ", ifelse(maxrelax == 0,"NO","YES"),"\n")
  cat("Relaxation weight     = ", maxrelax,"\n")
  cat("Recoding method       = ", ifelse(indiv_method == "sequential","Sequential","Optimal"),"\n")
  cat("---------------------------------------","\n")
  
  
  
  tstart = Sys.time()
  
  if (!(indiv_method %in% c("sequential","optimal"))){
    
    stop("Incorrect syntax int the indiv_method option")
    
  } else {}
  
  dataB = prep_dbs(datab,nominal = nominal,ordinal = ordinal,logic = logic) 
  
  inst = Instance(dataB,norme = norm)
  
  
  # Redefine A and B for the model
  nA    = inst$nA
  nB    = inst$nB
  A     = 1:inst$nA
  B     = 1:inst$nB
  Y     = inst$Y
  Z     = inst$Z
  indY  = inst$indY
  indZ  = inst$indZ
  indXA = inst$indXA
  indXB = inst$indXB
  prof  = do.call(paste0,unique(inst$Xobserv))
  
  nbindY = numeric(0)
  
  for (y in Y){
    
    nbindY  = c(nbindY,length(inst$indY[[y]]))}
    nbindZ  = numeric(0)
  
  for (z in Z){
      
    nbindZ  = c(nbindZ,length(inst$indZ[[z]]))}
    freqY   = numeric(0)
  
  for (y in Y){
    
    freqY   = c(freqY,nbindY[y] / length(A))}
    freqZ   = numeric(0)
  
  for (z in Z){
      
    freqZ = c(freqZ,nbindZ[z] / length(B))
    
  }
  
  
  ###########################################################################
  # Compute data for aggregation of the individuals
  ###########################################################################
  
  nbX = length(indXA);
  
  # Computation of the cost matrix as average distances between the
  # individuals of two groups
  C = average_distance_to_closest(inst, percent_closest = percent_c)[[1]]
  
  
  if (maxrelax == 0){
    
    # Objective: minimize the distance between individuals of A and B
    
    Min = c(as.numeric(C),rep(0,2*length(Y)+2*length(Z)))
    
    
    result          = ompr::MIPModel() %>%
      
      ompr::add_variable(transport[y,z], y = Y, z = Z, type = "continuous") %>%
      ompr::add_variable(deviationA[y] , y = Y       , type = "continuous") %>%
      ompr::add_variable(absdevA[y]    , y = Y       , type = "continuous") %>%
      ompr::add_variable(deviationB[z]        , z = Z, type = "continuous") %>%
      ompr::add_variable(absdevB[z]           , z = Z, type = "continuous") %>%
      
      ompr::set_objective(sum_expr(C[y,z]*transport[y,z], y = Y, z = Z) , "min") %>%
      
      ompr::add_constraint(sum_expr(transport[y,z]       , z = Z) == freqY[y] + deviationA[y], y = Y) %>%
      ompr::add_constraint(sum_expr(transport[y,z], y = Y       ) == freqZ[z] + deviationB[z], z = Z) %>%
      ompr::add_constraint(sum_expr(deviationA[y] , y = Y) == 0 ) %>%
      ompr::add_constraint(sum_expr(deviationB[z]        ,z = Z )== 0) %>%
      ompr::add_constraint(deviationB[z]<= absdevB[z]    , z = Z) %>%
      ompr::add_constraint(deviationB[z]>= -absdevB[z]   , z = Z) %>%
      ompr::add_constraint(sum_expr(absdevB[z],z=Z)<= maxrelax/2.0) %>%
      ompr::add_constraint(deviationA[y] <= absdevA[y],y =Y) %>%
      ompr::add_constraint(deviationA[y] >= -absdevA[y],y =Y) %>%
      ompr::add_constraint(sum_expr(absdevA[y],y = Y)<= maxrelax/2.0) %>%
      
      ompr::solve_model(with_ROI(solver = "glpk"))
    
    # Solve the problem
    
    solution       = ompr::get_solution(result, transport[y,z]) 
    
    # Extract the values of the solution
    
    transportA_val = matrix(solution[,4], length(Y),length(Z))
    transportB_val = transportA_val
    
    
  } else {
    
    b2 = numeric(length(Z))
    for (z in Z){
      stc_sum2 = vector(length = nbX)
      for (i in 1:nbX){
        
        stc_sum2[i] = ifelse(length(indXB[[i]])==0,1/length(Z),length(indXB[[i]][inst$Zobserv[indXB[[i]]+nA] == z])/ length(indXB[[i]]))*length(indXA[[i]])/nA
        
      } 
      b2[z]= sum(stc_sum2)}
    
    result <-  MIPModel() %>%
      
      ompr::add_variable(transportA[y,z],y = Y, z= Z,type = "continuous",lb=0) %>%
      ompr::add_variable(deviationB[z], z = Z, type = "continuous") %>%
      ompr::add_variable(absdevB[z], z = Z, type = "continuous",lb=0) %>%
      
      ompr::set_objective(sum_expr(C[y,z]*transportA[y,z],y = Y,z=Z) , "min") %>%
      
      ompr::add_constraint(sum_expr(transportA[y,z], z = Z) == freqY[y], y =Y) %>%
      ompr::add_constraint(sum_expr(transportA[y,z],y = Y) - deviationB[z] == b2[z] , z = Z) %>%
      ompr::add_constraint(sum_expr(deviationB[z],z = Z)== 0) %>%
      ompr::add_constraint(deviationB[z]<= absdevB[z], z = Z) %>%
      ompr::add_constraint(deviationB[z]>= -absdevB[z], z = Z) %>%
      ompr::add_constraint(sum_expr(absdevB[z],z=Z)<= maxrelax/2.0) %>%
      
      ompr::solve_model(with_ROI(solver = "glpk"))
    
    solution       = ompr::get_solution(result, transportA[y,z]) 
    
    transportA_val = matrix(solution$value, length(Y),length(Z))
    
    
    b1 = numeric(length(Y))
    
    for (y in Y){
      
      stc_sum2 = vector(length = nbX)
      
      for (i in 1:nbX){
        
        stc_sum2[i] = ifelse(length(indXA[[i]])==0,1/length(Y),length(indXA[[i]][inst$Yobserv[indXA[[i]]] == y])/ length(indXA[[i]]))*length(indXB[[i]])/nB
        
      } 
      
      b1[y]= sum(stc_sum2)}
    
      result <-  MIPModel() %>%
        
        ompr::add_variable(transportB[y,z],y = Y, z= Z,type = "continuous",lb=0) %>%
        ompr::add_variable(deviationA[y], y = Y, type = "continuous") %>%
        ompr::add_variable(absdevA[y], y = Y, type = "continuous",lb=0) %>%
        
        ompr::set_objective(sum_expr(C[y,z]*transportB[y,z],y = Y,z=Z) , "min") %>%
        
        ompr::add_constraint(sum_expr(transportB[y,z], z = Z) - deviationA[y] == b1[y], y =Y) %>%
        ompr::add_constraint(sum_expr(transportB[y,z],y = Y) == freqZ[z], z = Z) %>%
        ompr::add_constraint(sum_expr(deviationA[y],y = Y)== 0) %>%
        ompr::add_constraint(deviationA[y] <= absdevA[y],y =Y) %>%
        ompr::add_constraint(deviationA[y] >= -absdevA[y],y =Y) %>%
        ompr::add_constraint(sum_expr(absdevA[y],y = Y)<= maxrelax/2.0) %>%
        
        ompr::solve_model(with_ROI(solver = "glpk"))
      
    
      solution       = ompr::get_solution(result, transportB[y,z]) 
      transportB_val = matrix(solution$value, length(Y),length(Z))
    
 
    }
  
  ####
  # Get the individual transport from the group transport
  
  if (indiv_method == "sequential"){
    
    YApred = individual_from_group_closest(inst, transportA_val, transportB_val, percent_closest = percent_c)$YAtrans
    YBpred = individual_from_group_closest(inst, transportA_val, transportB_val, percent_closest = percent_c)$ZBtrans
    
  } else if (indiv_method == "optimal"){
    
    YApred = individual_from_group_optimal(inst, transportA_val, transportB_val, percent_closest = percent_c)$YAtrans
    YBpred = individual_from_group_optimal(inst, transportA_val, transportB_val, percent_closest = percent_c)$ZBtrans
  }
  
  # Compute the estimated probability distributions from predictions
  indXA = inst$indXA; indXB = inst$indXB;
  nbX = length(indXA);
  
  estimatorZA = array(rep(0,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))
  # estimatorYB = estimatorZA
  estimatorYB = array(rep(0,nbX*length(Y)*length(Z)),dim = c(nbX,length(Z),length(Y)))
  
  
  for (x in 1:nbX){
    
    for (i in indXA[[x]]){
      
      estimatorZA[x,inst$Yobserv[i],YBpred[i]] = estimatorZA[x,inst$Yobserv[i],YBpred[i]] +  1/sum(inst$Yobserv[indXA[[x]]] == inst$Yobserv[i])
    }
    
    
    for (y in Y){
      
      if (sum(inst$Yobserv[indXA[[x]]] == y) == 0){
        
        estimatorZA[x,y,] = 1/length(Z)*rep(1,length(Z))
      }
    
    }
    
    for (i in indXB[[x]]){
      
      # estimatorYB[x,YApred[i],inst$Zobserv[i+nA]] = estimatorYB[x,YApred[i],inst$Zobserv[i+nA]] + 1/ sum(inst$Zobserv[indXB[[x]] + nA] == inst$Zobserv[i + nA])
      estimatorYB[x,inst$Zobserv[i+nA],YApred[i]] = estimatorYB[x,inst$Zobserv[i+nA],YApred[i]] + 1/ sum(inst$Zobserv[indXB[[x]] + nA] == inst$Zobserv[i + nA])
    }
    
    for (z in Z){
      
      if (sum(inst$Zobserv[indXB[[x]]+inst$nA] == z) == 0){
        
        # estimatorYB[x,,z] = 1/length(Y)*rep(1,length(Y))
        estimatorYB[x,z,] = 1/length(Y)*rep(1,length(Y))
      
      }
    }
  }
  
  row.names(estimatorZA) = row.names(estimatorYB) = prof
  colnames(estimatorZA)  = as.character(levels(datab[,2]))
  colnames(estimatorYB)  = as.character(levels(datab[,3]))

  DATA1_OT         = datab[datab[,1] == 1,]
  DATA1_OT$OTpred  = plyr::mapvalues(YBpred,from = 1:max(dataB[,3],na.rm = TRUE), to = levels(datab[,3]))
  
  DATA2_OT         = datab[datab[,1] == 2,]
  DATA2_OT$OTpred  = plyr::mapvalues(YApred,from = 1:max(dataB[,2],na.rm = TRUE), to = levels(datab[,2]))
  
  
  tend = Sys.time()
  
  return(list(TIME_EXE = difftime(tend,tstart),TRANSPORT_A = transportA_val,TRANSPORT_B =transportB_val,estimatorZA= estimatorZA,estimatorYB = estimatorYB,DATA1_OT = DATA1_OT,DATA2_OT  = DATA2_OT))
}


