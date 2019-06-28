OT = function(inst, percent_closest=0.2, maxrelax=0, norme=0, indiv_method, full_disp=false, solver_disp=false){
    
    tstart = Sys.time()
    
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
    
    
    
    ###########################################################################
    # Compute data for aggregation of the individuals
    ###########################################################################
    
    nbX = length(indXA);
    
    # Computation of the cost matrix as average distances between the
    # individuals of two groups
    C = average_distance_to_closest(inst, percent_closest)[[1]]
    
    # G: ----> Jusqu'ici OK !!
    #---------------------
    
    if (maxrelax == 0){
        
        # Objective: minimize the distance between individuals of A and B
        
        Min = c(as.numeric(C),rep(0,2*length(Y)+2*length(Z)))
        
        
        result <-  MIPModel() %>%
            add_variable(transport[y,z],y = Y, z= Z,type = "continuous") %>%
            add_variable(deviationA[y], y = Y, type = "continuous") %>%
            add_variable(absdevA[y], y = Y, type = "continuous") %>%
            add_variable(deviationB[z], z = Z, type = "continuous") %>%
            add_variable(absdevB[z], z = Z, type = "continuous") %>%
            set_objective(sum_expr(C[y,z]*transport[y,z],y = Y,z=Z) , "min") %>%
            add_constraint(sum_expr(transport[y,z], z = Z) == freqY[y] + deviationA[y], y =Y) %>%
            add_constraint(sum_expr(transport[y,z],y = Y) == freqZ[z] + deviationB[z], z = Z) %>%
            add_constraint(sum_expr(deviationA[y],y = Y)== 0) %>%
            add_constraint(sum_expr(deviationB[z],z = Z)== 0) %>%
            add_constraint(deviationB[z]<= absdevB[z], z = Z) %>%
            add_constraint(deviationB[z]>= -absdevB[z], z = Z) %>%
            add_constraint(sum_expr(absdevB[z],z=Z)<= maxrelax/2.0) %>%
            add_constraint(deviationA[y] <= absdevA[y],y =Y) %>%
            add_constraint(deviationA[y] >= -absdevA[y],y =Y) %>%
            add_constraint(sum_expr(absdevA[y],y = Y)<= maxrelax/2.0) %>%
            solve_model(with_ROI(solver = "glpk"))
        
        solution <- get_solution(result, transport[y,z])  
        
    
    else{
        
        b2 = numeric(length(Z))
        for (z in Z){
            stc_sum2 = vector(length = nbX)
            for (i in 1:nbX){
                
                stc_sum2[i] = ifelse(length(indXB[[i]])==0,1/length(Z),length(indXB[[i]][inst$Zobserv[indXB[[i]]+nA] == z])/ length(indXB[[i]]))*length(indXA[[i]])/nA
                
            } 
            b2[z]= sum(stc_sum2)}
        
        result <-  MIPModel() %>%
            add_variable(transportA[y,z],y = Y, z= Z,type = "continuous") %>%
            add_variable(deviationB[z], z = Z, type = "continuous") %>%
            add_variable(absdevB[z], z = Z, type = "continuous") %>%
            set_objective(sum_expr(C[y,z]*transportA[y,z],y = Y,z=Z) , "min") %>%
            add_constraint(sum_expr(transportA[y,z], z = Z) == freqY[y], y =Y) %>%
            add_constraint(sum_expr(transportA[y,z],y = Y) == b2[z] , z = Z) %>%
            add_constraint(sum_expr(deviationB[z],z = Z)== 0) %>%
            add_constraint(deviationB[z]<= absdevB[z], z = Z) %>%
            add_constraint(deviationB[z]>= -absdevB[z], z = Z) %>%
            add_constraint(sum_expr(absdevB[z],z=Z)<= maxrelax/2.0) %>%
            solve_model(with_ROI(solver = "glpk"))
        
        solution <- get_solution(result, transportA[y,z]) 
        transportA = matrix(solution, length(Y),length(Z))
        
        b1 = numeric(length(Z))
        for (y in Y){
            stc_sum2 = vector(length = nbX)
            for (i in 1:nbX){
                
                stc_sum2[i] = ifelse(length(indXA[[i]])==0,1/length(Y),length(indXA[[i]][inst$Zobserv[indXA[[i]]] == z])/ length(indXA[[i]]))*length(indXB[[i]])/nB
                
            } 
            b1[y]= sum(stc_sum2)}
        
        result <-  MIPModel() %>%
            add_variable(transportB[y,z],y = Y, z= Z,type = "continuous") %>%
            add_variable(deviationA[y], y = Y, type = "continuous") %>%
            add_variable(absdevA[y], y = Y, type = "continuous") %>%
            set_objective(sum_expr(C[y,z]*transportB[y,z],y = Y,z=Z) , "min") %>%
            add_constraint(sum_expr(transportB[y,z], z = Z) == b1[y], y =Y) %>%
            add_constraint(sum_expr(transportB[y,z],y = Y) == freqZ[z], z = Z) %>%
            add_constraint(sum_expr(deviationA[y],y = Y)== 0) %>%
            add_constraint(deviationA[y] <= absdevA[y],y =Y) %>%
            add_constraint(deviationA[y] >= -absdevA[y],y =Y) %>%
            add_constraint(sum_expr(absdevA[y],y = Y)<= maxrelax/2.0) %>%
            solve_model(with_ROI(solver = "glpk"))
        
        solution <- get_solution(result, transportB[y,z]) 
        transportB = matrix(solution, length(Y),length(Z))
        
    }
    
    
    
    
    # Get the individual transport from the group transport
    if (indiv_method == "sequential"){
        
        YApred = individual_from_group_closest(inst, transportA_val, transportB_val, percent_closest)$YAtrans
        YBpred = individual_from_group_closest(inst, transportA_val, transportB_val, percent_closest)$YBtrans
        
    } else if (indiv_method == "optimal"){
        
        YApred = individual_from_group_optimal(inst, transportA_val, transportB_val, percent_closest)$YAtrans
        YBpred = individual_from_group_optimal(inst, transportA_val, transportB_val, percent_closest)$YBtrans
    }
    
    # Compute the estimated probability distributions from predictions
    indXA = inst$indXA; indXB = inst$indXB;
    nbX = length(indXA);
    
    estimatorZA = array(rep(0,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))
    estimatorYB = estimatorZA
    
    for (x in 1:nbX){
        for (i in indXA[[x]]){
            estimatorZA[x,inst$Yobserv[i],YBpred[i]] = estimatorZA[x,inst$Yobserv[i],YBpred[i]] +  1/sum(inst$Yobserv[indXA[[x]]] == inst$Yobserv[i])
        }
        
        
        for (y in Y){
            if (sum(inst$Yobserv[indXA[x]] == y) == 0){
                estimatorZA[x,y,] = 1/length(Z)*rep(1,length(Z))
            }
        }
        for (i in indXB[[x]]){
            estimatorYB[x,YApred[i],inst$Zobserv[i+nA]] = estimatorYB[x,YApred[i],inst$Zobserv[i+nA]] + 1/ sum(inst$Zobserv[indXB[[x]] + nA] == inst$Zobserv[i + nA])
        }
        for (z in Z){
            if (sum(inst$Zobserv[indXB[x]] == z)) == 0){
                estimatorYB[x,,z] = 1/length(Y)*rep(1,length(Y))
            }
        }
    }
    
    tend = Sys.time()
    
    return(list(difftime(tend,tstart),transportA_val,transportB_val,estimatorZA,estimatorYB))
}