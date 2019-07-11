
#' OT_joint(datab, 
#' maxrelax=0.0, 
#' lambda_reg = 0.0, 
#' percent_clos = 0.2, 
#' norm = 1, 
#' aggregate_tol=0.5, 
#' full_disp = FALSE, 
#' solver_disp = FALSE)
#'
#' @param datab todo list
#' @param maxrelax todo list
#' @param lambda_reg todo list
#' @param percent_clos todo list
#' @param norm todo list
#' @param aggregate_tol todo list
#' @param full_disp todo list
#' @param solver_disp todo list
#'
#' @return List with time and solution
#' 
#' @importFrom dplyr %>%
#' @importFrom ompr MIPModel get_solution
#' @importFrom ompr.roi with_ROI
#' 
#' @export
#'
# @examples
OT_joint = function(datab, maxrelax=0.0, lambda_reg = 0.0, percent_clos = 0.2, norm = 1, aggregate_tol=0.5, full_disp = FALSE, solver_disp = FALSE){
  
  cat("---------------------------------------","\n")
  cat("AGGREGATE INDIVIDUALS WRT COVARIATES","\n")
  cat("Reg. weight           =  ", lambda_reg,"\n")
  cat("Percent closest       = ", 100.0*percent_clos, "%","\n")
  cat("Aggregation tolerance = ", aggregate_tol,"\n")
  cat("---------------------------------------","\n")
  
  
  tstart = Sys.time()
  
  
  
  dataB = datab
  
  for (k in 1:ncol(datab)){
    
    dataB[,k] = as.numeric(dataB[,k])
    
  }
  
  inst = Instance(dataB,norme = norm)
  
  
  # Local redefinitions of parameters of  the instance
  nA      = inst$nA;
  nB      = inst$nB;
  A       = 1:nA;
  B       = 1:nB;
  Y       = inst$Y;
  Z       = inst$Z;
  indY    = inst$indY;
  indZ    = inst$indZ;
  Xobserv = inst$Xobserv;
  Yobserv = inst$Yobserv;
  Zobserv = inst$Zobserv;
  
  # Create a model for the optimal transport of individuals
  # modelA = Model(with_optimizer(Gurobi.Optimizer,LogToConsole=0,Method=2,Crossover=0));#ClpSolver(LogLevel=0)); 
  # modelB = Model(with_optimizer(Gurobi.Optimizer,LogToConsole=0,Method=2,Crossover=0));#Model(with_optimizer(Clp.Optimizer,LogLevel=0));
  
  
  ###########################################################################
  # Compute data for aggregation of the individuals
  ###########################################################################
  # println("... aggregating individuals")
  indXA = inst$indXA; indXB = inst$indXB;
  nbX = length(indXA);
  
  # compute the neighbors of the covariates for regularization
  Xvalues = unique(Xobserv)
  # norme = 1
  
  if (norm == 0){
    dist_X = rdist::cdist(Xvalues,Xvalues,"hamming")
    
  } else if (norm == 1){
    dist_X = rdist::cdist(Xvalues,Xvalues,"manhattan")
  }
  voisins_X = dist_X <= 1
  
  # println("... computing costs")
  C = average_distance_to_closest(inst, percent_closest = percent_clos)[1];
  
  ###########################################################################
  # Compute the estimators that appear in the model
  ###########################################################################
  
  estim_XA = estim_XB = estim_XA_YA =  estim_XB_ZB = list()
  
  for (x in 1:nbX){
    
    estim_XA[[x]] = length(indXA[[x]])/nA
    estim_XB[[x]] = length(indXB[[x]])/nB
    
  }
  
  for (x in 1:nbX){
    
    estim_XA_YA[[x]] = estim_XB_ZB[[x]] = numeric(0)
    
    
    for (y in Y){  
      estim_XA_YA[[x]][y] = length(indXA[[x]][Yobserv[indXA[[x]]] == y])/nA 
    }
    
    for (z in Z){
      estim_XB_ZB[[x]][z] = length(indXB[[x]][Zobserv[indXB[[x]] + nA] == z])/nB
    }
    
  }
  
  
  
  Cf <- function(y,z) {
    C$Davg[y,z]
  }
  
  estim_XBf <- function(x) {
    estim_XB[[x]]
  }
  
  voisin = function(x1){lambda_reg *(1/length(voisins_X[x1,]))}
  ind_voisins = list()
  
  for (x1 in 1:nrow(voisins_X)){
    ind_voisins[[x1]] = which(voisins_X[x1,])
  }
  ###########################################################################
  # Basic part of the model
  ###########################################################################
  
  # println("... creating model")
  # Variables
  # - gammaA[x,y,z]: joint probability of X=x, Y=y and Z=z in base A
  result <-  MIPModel() %>%
    ompr::add_variable(gammaA[x,y,z],x = 1:nbX, y=Y,z= Z,type = "continuous") %>%
    ompr::add_variable(errorA_XY[x,y],x=1:nbX, y = Y, type = "continuous") %>%
    ompr::add_variable(abserrorA_XY[x,y],x=1:nbX, y = Y, type = "continuous") %>%
    ompr::add_variable(errorA_XZ[x,z],x=1:nbX, z = Z, type = "continuous") %>%
    ompr::add_variable(abserrorA_XZ[x,z],x=1:nbX, z = Z, type = "continuous") %>%
    ompr::add_variable(reg_absA[x1, x2,y,z],x1=1:nbX, x2=1:nbX,y=Y,z= Z, type = "continuous") %>%
    #set_objective(sum_expr(C(y,z) * gammaA[x,y,z],y = Y,z=Z, x = 1:nbX) + lambda_reg * sum_expr(1/length(=voisin(x1)) *reg_absA[x1,x2,y,z]) , "min") %>%
    #set_objective(sum_expr(Cf(y,z) * gammaA[x,y,z],y = Y,z=Z, x = 1:nbX)  + lambda_reg * sum_expr(1/length(voisins_X[x1,]) *reg_absA[x1,x2,y,z],x1=1:nbX, x2= which(voisins_X[x1,]),y=Y,z= Z) , "min") %>%
    ompr::set_objective(sum_expr(Cf(y,z) * gammaA[x,y,z],y = Y,z=Z, x = 1:nbX) + sum_expr(voisin(x1)*reg_absA[x1,x2,y,z],x1=1:nbX, x2= ind_voisins[[x1]],y=Y,z= Z), "min") %>%
    ompr::add_constraint(sum_expr(gammaA[x,y,z], z = Z) - errorA_XY[x,y] == estim_XA_YA[[x]][y] , x = 1:nbX,y =Y) %>%
    ompr::add_constraint(estim_XBf(x)*sum_expr(gammaA[x,y,z],y = Y)  - estim_XBf(x)*errorA_XZ[x,z]== estim_XB_ZB[[x]][z] * estim_XA[[x]] , x = 1:nbX, z = Z) %>%
    ompr::add_constraint(errorA_XY[x,y] <= abserrorA_XY[x,y], x = 1:nbX,y =Y) %>%
    ompr::add_constraint(-errorA_XY[x,y] <= abserrorA_XY[x,y], x = 1:nbX,y =Y) %>%
    ompr::add_constraint(sum_expr(abserrorA_XY[x,y], x = 1:nbX,y = Y)<= maxrelax/2.0) %>%
    ompr::add_constraint(sum_expr(errorA_XY[x,y], x = 1:nbX,y = Y)== 0.0) %>%
    ompr::add_constraint(errorA_XZ[x,z] <= abserrorA_XZ[x,z], x = 1:nbX,z =Z) %>%
    ompr::add_constraint (-errorA_XZ[x,z] <= abserrorA_XZ[x,z], x = 1:nbX,z =Z) %>%
    ompr::add_constraint(sum_expr(abserrorA_XZ[x,z], x = 1:nbX,z = Z)<= maxrelax/2.0) %>%
    ompr::add_constraint(sum_expr(errorA_XZ[x,z], x = 1:nbX,z = Z)<= maxrelax/2.0) %>%
    ompr::solve_model(with_ROI(solver = "glpk"))
  
  solution  = ompr::get_solution(result, gammaA[x,y,z]) 
  gammaA_val= array(solution$value,dim = c(nbX,length(Y),length(Z)))
  
  # for (x in 1:nbX){
  # for (z in Z){
  #  for (y in Y){ 
  #    gammaA_sol[x,y,z] = solution$value[solution$y== y & solution$z==z & solution$x==x]
  #  }}}
  
  result <-  ompr::MIPModel() %>%
    ompr::add_variable(gammaB[x,y,z],x = 1:nbX, y=Y,z= Z,type = "continuous") %>%
    ompr::add_variable(errorB_XY[x,y],x=1:nbX, y = Y, type = "continuous") %>%
    ompr::add_variable(abserrorB_XY[x,y],x=1:nbX, y = Y, type = "continuous") %>%
    ompr::add_variable(errorB_XZ[x,z],x=1:nbX, z = Z, type = "continuous") %>%
    ompr::add_variable(abserrorB_XZ[x,z],x=1:nbX, z = Z, type = "continuous") %>%
    ompr::add_variable(reg_absB[x1, x2,y,z],x1=1:nbX, x2=1:nbX,y=Y,z= Z, type = "continuous") %>%
    ompr::set_objective(sum_expr(Cf(y,z) * gammaB[x,y,z],y = Y,z=Z, x = 1:nbX)  + sum_expr(voisin(x1)*reg_absB[x1,x2,y,z],x1=1:nbX, x2= ind_voisins[[x1]],y=Y,z= Z), "min")  %>%
    ompr::add_constraint(sum_expr(gammaB[x,y,z], y = Y) -errorB_XZ[x,z]  == estim_XB_ZB[[x]][z] , x = 1:nbX,z =Z) %>%
    ompr::add_constraint(estim_XA[[x]]*sum_expr(gammaB[x,y,z] ,z = Z) - estim_XA[[x]] * errorB_XY[x,y] == estim_XA_YA[[x]][y] * estim_XB[[x]] , x = 1:nbX, y = Y) %>%
    ompr::add_constraint(errorB_XY[x,y] <= abserrorB_XY[x,y], x = 1:nbX,y =Y) %>%
    ompr::add_constraint(-errorB_XY[x,y] <= abserrorB_XY[x,y], x = 1:nbX,y =Y) %>%
    ompr::add_constraint(sum_expr( abserrorB_XY[x,y], x = 1:nbX,y = Y)<= maxrelax/2.0) %>%
    ompr::add_constraint(sum_expr(errorB_XY[x,y], x = 1:nbX,y = Y)== 0.0) %>%
    ompr::add_constraint(errorB_XZ[x,z] <= abserrorB_XZ[x,z], x = 1:nbX,z =Z) %>%
    ompr::add_constraint( -errorB_XZ[x,z] <= abserrorB_XZ[x,z], x = 1:nbX,z =Z) %>%
    ompr::add_constraint(sum_expr(abserrorB_XZ[x,z], x = 1:nbX,z = Z)<= maxrelax/2.0) %>%
    ompr::add_constraint(sum_expr(errorB_XZ[x,z], x = 1:nbX,z = Z)<= maxrelax/2.0) %>%
    ompr::solve_model(with_ROI(solver = "glpk"))
  
  solution <- get_solution(result, gammaB[x,y,z]) 
  
  
  # Solve the problem
  # optimize!(modelA);
  # optimize!(modelB);
  
  # Extract the values of the solution
  # gammaA_val = [value(gammaA[x,y,z]) for x in 1:nbX, y in Y, z in Z];
  # gammaB_val = [value(gammaB[x,y,z]) for x in 1:nbX, y in Y, z in Z];
  gammaB_val= array(solution$value,dim = c(nbX,length(Y),length(Z)))
  
  # compute the resulting estimators for the distributions of Z conditional to X and Y in base A and of Y conditional to X and Z in base B
  estimatorZA = 1/length(Z) * array(rep(1,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))
  
  for (x in 1:nbX){
    
    for (y in Y){
      
      proba_c_mA = apply(gammaA_val,c(1,2),sum)[x,y]
      
      if (proba_c_mA > 1.0e-6){
        
        estimatorZA[x,y,] = 1/proba_c_mA * gammaA_val[x,y,];
        
      } else {}
    }
  }
  
  estimatorYB = 1/length(Y) * array(rep(1,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))
  
  for (x in 1:nbX){
    
    for (z in Z){
      
      proba_c_mB = apply(gammaB_val,c(1,3),sum)[x,z]
      
      if (proba_c_mB > 1.0e-6){
        
        estimatorYB[x,,z] = 1/proba_c_mB * gammaB_val[x,,z];
        
      } else {}
    }
  }
  
  #deduce the individual distributions of probability for each individual
  probaZindivA = matrix(0,nA,length(Z))
  probaYindivB = matrix(0,nB,length(Y))
  
  for (x in 1:nbX){
    for (i in indXA[[x]]){
      probaZindivA[i,] = estimatorZA[x,Yobserv[i],]
    }
    for (i in indXB[[x]]){
      probaYindivB[i,] = estimatorYB[x,,Zobserv[i+nA]]
    }
  }
  
  
  
  # Transport the Ylity that maximizes frequency
  
  predZA = numeric(0)
  for (i in A){
    
    predZA = c(predZA,which.max(probaZindivA[i,]))
    
  }
  
  predYB = numeric(0)
  for (j in B){
    
    predYB = c(predYB,which.max(probaYindivB[j,]))
    
  }
  
  # Display the solution
  # println("Solution of the joint probability transport");
  # println("Distance cost = ", sum(C[y,z] * (gammaA_val[x,y,z]+gammaB_val[x,y,z]) for y in Y, z in Z, x in 1:nbX));
  # println("Regularization cost = ", lambda_reg * value(regterm));
  
  
  DATA1_OT         = datab[datab[,1] == 1,]
  DATA1_OT$OTpred  = plyr::mapvalues(predZA,from = 1:max(dataB[,3],na.rm = TRUE), to = levels(datab[,3]))
  
  DATA2_OT         = datab[datab[,1] == 2,]
  DATA2_OT$OTpred  = plyr::mapvalues(predYB,from = 1:max(dataB[,2],na.rm = TRUE), to = levels(datab[,2]))
  
  
  tend = Sys.time()
  
  return(list(TIME_EXE = difftime(tend,tstart),
              GAMMA_A = apply(gammaA_val,c(2,3),sum),
              GAMMA_B = apply(gammaB_val,c(2,3),sum), 
              estimatorZA = estimatorZA, 
              estimatorYB = estimatorYB, 
              DATA1_OT = DATA1_OT, 
              DATA2_OT = DATA2_OT))
}

