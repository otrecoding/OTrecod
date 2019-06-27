
# OTjoint avec la partie optimisation à ajouter tout le reste de la fonction R ok

# require(rdist)


#############################################################################################################################################################
# Model where we directly compute the distribution of the outcomes for each
# individual or for sets of indviduals that similar values of covariates
# aggregate_tol: quantify how much individuals' covariates must be close for aggregation
# reg_norm: norm1, norm2 or entropy dep}ing on the type of regularization
# percent_closest: percent of closest neighbors taken into consideration in regularization
# lambda_reg: coefficient measuing the importance of the regularization term
# full_disp: if true, write the transported value of each individual; otherwise, juste write the number of missed transports
# solver_disp: if false, do not display the outputs of the solver
#############################################################################################################################################################

OT_joint=function(inst, maxrelax=0.0, lambda_reg:=0.0, percent_closest=0.2, norme=0, aggregate_tol=0.5, full_disp=false, solver_disp=false){

    # println("#################################################################")
    # println("AGGREGATE INDIVIDUALS WRT COVARIATES")
    # println("Reg. weight =  ", lambda_reg)
    # println("Percent closest = ", 100.0*percent_closest, "\%")
    # println("Aggregation tolerance = ", aggregate_tol,"\n")
    # println("#################################################################\n")


    tstart = time()

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
    
    if (norme == 0){
        dist_X = cdist(Xvalues),Xvalues,"hamming")}
    else if (norme == 1){
        dist_X = cdist(Xvalues,Xvalues,"manhattan")
    }
    voisins_X = dist_X <= 1

    # println("... computing costs")
    C = average_distance_to_closest(inst, percent_closest)[1];

    ###########################################################################
    # Compute the estimators that appear in the model
    ###########################################################################

    estim_XA = estim_XB = estim_XA_YA =  estim_XB_ZB = list()
    
    for (x in 1:nbX){
    
      estim_XA[[x]] = length(indXA[[x]])/nA
      estim_XB[[x]] = length(indXB[[x]])/nB
      
    }
    
    for (x in 1:nbX){
      
      for (y in Y){  
        estim_XA_YA[[x]][y] = length(indXA[[x]][Yobserv[indXA[[x]]] == y])/nA 
      }
      
      for (z in Z){
        estim_XB_ZB[[x]][y] = length(indXB[[x]][Zobserv[indXB[[x]] + nA] == z])/nB
      }
      
    }
      

    ###########################################################################
    # Basic part of the model
    ###########################################################################

    # println("... creating model")
    # Variables
    # - gammaA[x,y,z]: joint probability of X=x, Y=y and Z=z in base A
    @variable(modelA, gammaA[x in 1:nbX, y in Y, z in Z] >= 0, base_name="gammaA")

    # - gammaB[x,y,z]: joint probability of X=x, Y=y and Z=z in base B
    @variable(modelB, gammaB[x in 1:nbX, y in Y, z in Z] >= 0, base_name="gammaB")

    @variable(modelA, errorA_XY[x in 1:nbX, y in Y], base_name="errorA_XY");
    @variable(modelA, abserrorA_XY[x in 1:nbX, y in Y] >= 0, base_name="abserrorA_XY");
    @variable(modelA, errorA_XZ[x in 1:nbX, z in Z], base_name="errorA_XZ");
    @variable(modelA, abserrorA_XZ[x in 1:nbX, z in Z] >= 0, base_name="abserrorA_XZ");

    @variable(modelB, errorB_XY[x in 1:nbX, y in Y], base_name="errorB_XY");
    @variable(modelB, abserrorB_XY[x in 1:nbX, y in Y] >= 0, base_name="abserrorB_XY");
    @variable(modelB, errorB_XZ[x in 1:nbX, z in Z], base_name="errorB_XZ");
    @variable(modelB, abserrorB_XZ[x in 1:nbX, z in Z] >= 0, base_name="abserrorB_XZ");

    # Constraints
    # - assign sufficient probability to each class of covariates with the same outcome
    @constraint(modelA, ctYandXinA[x in 1:nbX, y in Y], sum(gammaA[x,y,z] for z in Z) == estim_XA_YA[x,y] + errorA_XY[x,y]);
    @constraint(modelB, ctZandXinB[x in 1:nbX, z in Z], sum(gammaB[x,y,z] for y in Y) == estim_XB_ZB[x,z] + errorB_XZ[x,z]);


    # - we impose that the probability of Y conditional to X is the same in the two databases
    # - the consequence is that the probability of Y and Z conditional to Y is also the same in the two bases
    @constraint(modelA, ctZandXinA[x in 1:nbX, z in Z], estim_XB[x]*sum(gammaA[x,y,z] for y in Y) == estim_XB_ZB[x,z] * estim_XA[x] + estim_XB[x]*errorA_XZ[x,z]);

    @constraint(modelB, ctYandXinB[x in 1:nbX, y in Y], estim_XA[x]*sum(gammaB[x,y,z] for z in Z) == estim_XA_YA[x,y] * estim_XB[x] + estim_XA[x] * errorB_XY[x,y]);

    # - recover the norm 1 of the error
    @constraint(modelA,[x in 1:nbX, y in Y], errorA_XY[x,y] <= abserrorA_XY[x,y]);
    @constraint(modelA,[x in 1:nbX, y in Y], -errorA_XY[x,y] <= abserrorA_XY[x,y]);
    @constraint(modelA, sum( abserrorA_XY[x,y] for x in 1:nbX, y in Y) <= maxrelax/2.0);
    @constraint(modelA, sum(errorA_XY[x,y] for x in 1:nbX, y in Y) == 0.0);
    @constraint(modelA,[x in 1:nbX, z in Z], errorA_XZ[x,z] <= abserrorA_XZ[x,z]);
    @constraint(modelA,[x in 1:nbX, z in Z], -errorA_XZ[x,z] <= abserrorA_XZ[x,z]);
    @constraint(modelA, sum( abserrorA_XZ[x,z] for x in 1:nbX, z in Z) <= maxrelax/2.0);
    @constraint(modelA, sum(errorA_XZ[x,z] for x in 1:nbX, z in Z) == 0.0);

    @constraint(modelB,[x in 1:nbX, y in Y], errorB_XY[x,y] <= abserrorB_XY[x,y]);
    @constraint(modelB,[x in 1:nbX, y in Y], -errorB_XY[x,y] <= abserrorB_XY[x,y]);
    @constraint(modelB, sum( abserrorB_XY[x,y] for x in 1:nbX, y in Y) <= maxrelax/2.0);
    @constraint(modelB, sum(errorB_XY[x,y] for x in 1:nbX, y in Y) == 0.0);
    @constraint(modelB,[x in 1:nbX, z in Z], errorB_XZ[x,z] <= abserrorB_XZ[x,z]);
    @constraint(modelB,[x in 1:nbX, z in Z], -errorB_XZ[x,z] <= abserrorB_XZ[x,z]);
    @constraint(modelB, sum( abserrorB_XZ[x,z] for x in 1:nbX, z in Z) <= maxrelax/2.0);
    @constraint(modelB, sum(errorB_XZ[x,z] for x in 1:nbX, z in Z) == 0.0);

    # - regularization
    @variable(modelA, reg_absA[x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z] >= 0);
    @constraint(modelA, [x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z], reg_absA[x1,x2,y,z] >= gammaA[x1,y,z]/(max(1,length(indXA[x1]))/nA)-gammaA[x2,y,z]/(max(1,length(indXA[x2]))/nA));
    @constraint(modelA, [x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z], reg_absA[x1,x2,y,z] >= gammaA[x2,y,z]/(max(1,length(indXA[x2]))/nA)-gammaA[x1,y,z]/(max(1,length(indXA[x1]))/nA));
    @expression(modelA, regterm, sum(1/length(voisins_X[x1,]) *reg_absA[x1,x2,y,z] for x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z));

    @variable(modelB, reg_absB[x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z] >= 0);
    @constraint(modelB, [x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z], reg_absB[x1,x2,y,z] >= gammaB[x1,y,z]/(max(1,length(indXB[x1]))/nB)-gammaB[x2,y,z]/(max(1,length(indXB[x2]))/nB));
    @constraint(modelB, [x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z], reg_absB[x1,x2,y,z] >= gammaB[x2,y,z]/(max(1,length(indXB[x2]))/nB)-gammaB[x1,y,z]/(max(1,length(indXB[x1]))/nB));
    @expression(modelB, regterm, sum(1/length(voisins_X[x1,]) *reg_absB[x1,x2,y,z] for x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z));

    # by default, the OT cost and regularization term are weighted to lie in the same interval
    @objective(modelA, Min, sum(C[y,z] * gammaA[x,y,z]  for y in Y, z in Z, x in 1:nbX) + lambda_reg *  sum(1/length(voisins_X[x1,]) *reg_absA[x1,x2,y,z] for x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z));

    @objective(modelB, Min, sum(C[y,z] * gammaB[x,y,z]  for y in Y, z in Z, x in 1:nbX) + lambda_reg *  sum(1/length(voisins_X[x1,]) *reg_absB[x1,x2,y,z] for x1 in 1:nbX, x2 in findall(voisins_X[x1,]), y in Y, z in Z));

   # Solve the problem
   # optimize!(modelA);
   # optimize!(modelB);

   # Extract the values of the solution
   # gammaA_val = [value(gammaA[x,y,z]) for x in 1:nbX, y in Y, z in Z];
   # gammaB_val = [value(gammaB[x,y,z]) for x in 1:nbX, y in Y, z in Z];

   
   # compute the resulting estimators for the distributions of Z conditional to X and Y in base A and of Y conditional to X and Z in base B
   estimatorZA = 1/length(Z) * array(rep(1,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))
   
   for (x = 1:nbX){
   
       for (y in Y){
           
          proba_c_mA = apply(gammaA_val,c(1,2),sum)[x,y]
          
          if (proba_c_mA > 1.0e-6){
          
                estimatorZA[x,y,] = 1/proba_c_mA * gammaA_val[x,y,];
          
          } else {}
       }
   }
   
   estimatorYB = 1/length(Y) * array(rep(1,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))
   
   for (x = 1:nbX){
   
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
   
     for (z in Z){
     
        predZA = c(predZA,which.max(probaZindivA[i,z])
        
     }
   }
   
   predYB = numeric(0)
   for (j in B){
   
     for (y in Y){
   
        predYB = c(predYB,which.max([probaYindivB[j,y])
        
     }
     
   }

   # Display the solution
   # println("Solution of the joint probability transport");
   # println("Distance cost = ", sum(C[y,z] * (gammaA_val[x,y,z]+gammaB_val[x,y,z]) for y in Y, z in Z, x in 1:nbX));
   # println("Regularization cost = ", lambda_reg * value(regterm));

   tend = Sys.time()

   return list(difftime(tend,tstart),apply(gammaA_val,c(2,3),sum),apply(gammaB_val,c(2,3),sum), estZA = estimatorZA, estYB = estimatorYB))
}
