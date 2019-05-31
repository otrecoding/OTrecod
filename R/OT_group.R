#using JuMP
#using Gurobi
#using Cbc
#using Clp
# using Ipopt
#using Base
#include("utils.jl")

#@enum IndivFromGroup sequential optimal



###############################################################################
# Sequentially assign the modality of the individuals to that of the closest
# neighbor in the other base until the joint probability values are met
###############################################################################

individual_from_group_closest=function(inst, jointprobaA, jointprobaB, percent_closest=1.0){

    # Redefine A and B for the model
    A = 1:inst$nA
    B = 1:inst$nB
    Y = inst$Y
    Z = inst$Z
    indY = inst$indY
    indZ = inst$indZ
    nbindY = numeric(0)
    for (y in Y){
    nbindY = c(nbindY,length(indY[y]))}
    nbindZ= numeric(0)
    for (z in Z){
    nbindZ = c(nbindZ,length(indZ[z]))}
    freqY= numeric(0)
    for (y in Y){
    freqY = c(freqY,nbindY[y] / length(A))}
    freqZ= numeric(0)
    for (z in Z){
    freqZ = c(freqZ,nbindZ[z] / length(B))}

    # In essence, assign to each individual the modality that is closest, where the distance from an individual to a modality is computed as the average distance to the individuals having this modality (in the other base)
    YAtrans = c(undef,inst$nB);
    YBtrans = c(undef,inst$nA);
    average_distance = average_distance_to_closest(inst, percent_closest);
    Davg=average_distance$Davg
    DindivA= average_distance$DindivA
    DindivB = average_distance$DindivB
    # DA = [(z,[sum(inst$D[i,j] for j in indZ[z])/nbindZ[z] for i in A]) for z in Z]
    # DB = [(y,[sum(inst$D[i,j] for i in indY[y])/nbindY[y] for j in B]) for y in Y]
    DA = [(z,[DindivA[i,z] for i in A]) for z in Z];
    DB = [(y,[DindivB[j,y] for j in B]) for y in Y];
    for (y in Y){
        indtrans = indY[y]
        for(z in Z){
            nbtrans = min(round(Int,jointprobaA[y,z]/freqY[y] * nbindY[y]), length(indtrans));
            distance = [(i,DA[z][2][i]) for i in indtrans]
            sort!(distance, by = x->x[2])
            for (k = 1:nbtrans){
                YBtrans[distance[k][1]] = z;
                indtrans = subset(indtrans,indtrans != distance[k][1])
            }
        }

        # affect potential individuals that have not been transported due to
        # rounding
        for (i in indtrans){
            YBtrans[i] = inst$Zobserv[which.min([inst$D[i,j] for j in B])[2]]
        }
    }

    for (z in Z){
        indtrans = indZ[z];
        for (y in Y){
            nbtrans = min(round(Int,jointprobaB[y,z]/freqZ[z] * nbindZ[z]), length(indtrans));
            distance = [(j,DB[y][2][j]) for j in indtrans];
            sort!(distance, by = x->x[2]);
            for (k in 1:nbtrans){
                YAtrans[distance[k][1]] = y;
                indtrans = subset(indtrans,indtrans != distance[k][1])
            }
        }

        # affect potential individuals that have not been transported due to
        # rounding
        for (j in indtrans){
            YAtrans[j] = inst$Yobserv[which.min([inst$D[i,j] for i in A])[2]];
        }
    }

    return(list(YAtrans, YBtrans))
}

###############################################################################
# Solve an optimization problem to get the individual transport that minimizes
# total distance while satisfying the joint probability computed by the model by
# group
###############################################################################

individual_from_group_optimal=function (inst, jointprobaA, jointprobaB, percent_closest=1.0){


    # Redefine A and B for the model
    A = 1:inst$nA
    B = 1:inst$nB
    Y = inst$Y
    Z = inst$Z
    indY = inst$indY
    indZ = inst$indZ
    nbindY = numeric(0)
    for (y in Y){
    nbindY = c(nbindY,length(indY[y]))}
    nbindZ = numeric(0)
    for (z in Z){
    nbindZ = c(nbindZ,length(indZ[z]))

    # Create a model for the optimal transport of individuals
    # indiv = Model(solver=IpoptSolver(print_level=4))
    indiv = Model(with_optimizer(Clp.Optimizer,LogLevel=0))
    # indiv = Model(solver=GurobiSolver(Method=2,LogToConsole=0)); #Presolve=0,Method=2,Crossover=0))


    # Variables
    # - assignA[i][z] : fraction of individual i assigned to modality z
    @variable(indiv, assignA[i in A, z in Z] >= 0, base_name="assignA")
    # - assignB[j][y] : fraction of individual j assigned to modality y
    @variable(indiv, assignB[j in B, y in Y] >= 0, base_name="assignB")

    # compute the average distance between the individuals and the modalities of
    # the other base
    CA = zeros(inst$nA,length(Z))
    CB = zeros(inst$nB, length(Y))
    for (i in A){
        for (z in Z){
            nbclose = round(Int,percent_closest*nbindZ[z])
            distance = numeric(0)
            for (j in indZ[z]){
            distance = c(distance,sort(inst$D[i,j]))}
            sort(distance)
            CA[i,z] = sum(distance[1:nbclose])/nbclose
        }
    }
    for (j in B){
        for (y in Y){
            nbclose = round(Int,percent_closest*nbindY[y])
            distance = numeric(0)
            for (i in indY[y]){
            distance = c(distance,sort(inst$D[i,j]))}
            sort(distance)
            CB[j,y] = CB[j,y]+ sum(distance[1:nbclose])/nbclose
        }
    }

    # Objective: minimize the distance between individuals of A and B
    @objective(indiv, Min, sum(CA[i,z]*assignA[i,z] for i in A, z in Z)
                            + sum(CB[j,y]*assignB[j,y] for j in B, y in Y))

    # Constraints
    # - assign the individuals so as to satisfy the joint probability computed
    #   with the model by group
    @constraint(indiv, ctjointprobaA[y in Y, z in Z],
        sum(assignA[i,z] for i in indY[y]) == jointprobaA[y,z])
    @constraint(indiv, ctjointprobaB[y in Y, z in Z],
        sum(assignB[j,y] for j in indZ[z]) == jointprobaB[y,z])

    # - assign sufficient modality to each individual
    @constraint(indiv, ctassignA[i in A], sum(assignA[i,z] for z in Z) == 1/(length(A)))
    @constraint(indiv, ctassignB[j in B], sum(assignB[j,y] for y in Y) == 1/(length(B)))


    # Solve the problem
    optimize!(indiv)

    # Extract the values of the solution
    assignA_val = [value(assignA[i,z]) for i in A, z in Z]
    assignB_val = [value(assignB[j,y]) for j in B, y in Y]

    # Transport the modality that maximizes frequency
    YBtrans = [findmax([assignA_val[i,z]  for z in Z])[2] for i in A]
    YAtrans = [findmax([assignB_val[j,y]  for y in Y])[2] for j in B]

    return(list(YAtrans, YBtrans))
}


###############################################################################
# Model of group transport
# percent_closest: percent of closest neighbors taken in the computation of the
#   costs
# maxrelax: maximum percentage of deviation from expected probability masses
# indiv_method: specifies the method used to get individual transport from
#   group joint probabilities
# full_disp: if true, write the transported value of each individual;
#       otherwise, juste write the number of missed transports
# solver_disp: if false, do not display the outputs of the solver
###############################################################################

OT_group=function(inst, percent_closest=0.2, maxrelax=0.0, norme=0, indiv_method, full_disp=false, solver_disp=false){

    tstart = time()

    # Redefine A and B for the model
    nA = inst$nA
    nB = inst$nB
    A = 1:inst$nA
    B = 1:inst$nB
    Y = inst$Y
    Z = inst$Z
    indY = inst$indY
    indZ = inst$indZ
    nbindY = numeric(0)
    for(y in Y){
    nbindY = c(nbindY,length(indY[y]))
    nbindZ = numeric(0)
    for (z in Z){
    nbindZ =  c(nbindZ,length(indZ[z]))}
    freqY= numeric(0)
    for (y in Y){
    freqY =  c(freqY,nbindY[y] / length(A)) }
    freqZ= numeric(0)
    for (z in Z){
    freqZ =  c(freqZ,nbindZ[z] / length(B)) }

    ###########################################################################
    # Compute data for aggregation of the individuals
    ###########################################################################
    indXA = inst$indXA; indXB = inst$indXB;
    nbX = length(indXA);

    # Computation of the cost matrix as average distances between the
    # individuals of two groups
    C = average_distance_to_closest(inst, percent_closest)[1];

    if (maxrelax == 0.0){
        # Create a model for the optimal transport of individuals
        # group = Model(solver=GurobiSolver(Method=2, LogToConsole=0));
        group = Model(with_optimizer(Clp.Optimizer));

        # Variables
        # - transport[y,z] : joint probability of modalities y and z
        @variable(group, transport[y in Y, z in Z] >= 0, base_name="transport");
        # - deviationA[y]: deviation of the probability mass of Y from
        #   that observed in base A
        @variable(group, deviationA[y in Y], base_name="deviationA")
        # - absdevA[y]: absolute value of the deviation of the probability mass of
        #   Y from that observed in base A
        @variable(group, absdevA[y in Y] >= 0, base_name="absdevA")
        # - deviationB[z]: deviation of the probability mass of Z from
        #   that observed in base B
        @variable(group, deviationB[z in Z], base_name="deviationB")
        # - absdevB[z]: absolute value of the deviation of the probability mass of
        #  Z from that observed in base B
        @variable(group, absdevB[z in Z] >= 0, base_name="absdevB")

        # Objective: minimize the distance between individuals of A and B
        @objective(group, Min, sum(C[y,z]*transport[y,z] for y in Y, z in Z));

        # Constraints
        # - satisfy probabilities of modality A
        @constraint(group, cttransportA[y in Y], sum(transport[y,z] for z in Z) == freqY[y] + deviationA[y]);

        # - satisfy probabilities of modality B
        @constraint(group, cttransportB[z in Z], sum(transport[y,z] for y in Y) == freqZ[z] + deviationB[z]);

        # - the deviations must sum to zero to conserve a well-defined probability measure
        @constraint(group, ctsumdeviationA, sum(deviationA[y] for y in Y) == 0)
        @constraint(group, ctsumdeviationB, sum(deviationB[z] for z in Z) == 0)

        # - bound the norm 1 of deviations
        @constraint(group, ctabsdevBplus[z in Z], deviationB[z] <= absdevB[z]);
        @constraint(group, ctabsdevBmoins[z in Z], deviationB[z] >= -absdevB[z]);
        @constraint(group, ctbounddevB, sum([absdevB[z] for z in Z]) <= maxrelax/2.0);
        @constraint(group, ctabsdevAplus[y in Y], deviationA[y] <= absdevA[y]);
        @constraint(group, ctabsdevAmoins[y in Y], deviationA[y] >= -absdevA[y]);
        @constraint(group, ctbounddevA, sum([absdevA[y] for y in Y]) <= maxrelax/2.0);

        # Solve the problem
        optimize!(group);

        # Extract the values of the solution
        transportA_val = [value(transport[y,z]) for y in Y, z in Z];
        transportB_val = transportA_val}


    # Otherwise, the empirical distribution constraints are relaxed and the
    # model can be decomposed by database
    else{
       # Create a model for the optimal transport of individuals
       groupA = Model(with_optimizer(Clp.Optimizer))
       # groupA = Model(solver=GurobiSolver(Method=2,LogToConsole=0))
       # Create a model for the optimal transport of individuals
       groupB = Model(with_optimizer(Clp.Optimizer))
       # groupB = Model(solver=GurobiSolver(Method=2,LogToConsole=0))


       # - transportA[y,z] : joint probability of modalities y and z if in base A
       @variable(groupA, transportA[y in Y, z in Z] >= 0, base_name="transportA")
       # - transportB[y,z] : joint probability of modalities y and z if in base B
       @variable(groupB, transportB[y in Y, z in Z] >= 0, base_name="transportB")

       # - deviationA[y]: deviation of the probability mass of y in base B from
       #   that observed in base A
       @variable(groupB, deviationA[y in Y], base_name="deviationA")
       # - absdevA[y]: absolute value of the deviation of the probability mass of
       #   y in base B from that observed in base A
       @variable(groupB, absdevA[y in Y] >= 0, base_name="absdevA")
       # - deviationB[z]: deviation of the probability mass of z in base A from
       #   that observed in base B
       @variable(groupA, deviationB[z in Z], base_name="deviationB")

       # - absdevB[z]: absolute value of the deviation of the probability mass of
       #   z in base A from that observed in base B
       @variable(groupA, absdevB[z in Z] >= 0, base_name="absdevB")

       # Objective: minimize the distance between individuals of A and B
       @objective(groupA, Min, sum(C[y,z]*transportA[y,z] for y in Y, z in Z))
       @objective(groupB, Min, sum(C[y,z]*transportB[y,z] for y in Y, z in Z))


       # Constraints
       # - satisfy probabilities of modality A for individuals in base A
       @constraint(groupA, cttransportAinA[y in Y], sum(transportA[y,z] for z in Z) == freqY[y])
       # - satisfy probabilities of modality B with possible deviation for individuals in base A
       @constraint(groupA, cttransportBinA[z in Z],
           sum(transportA[y,z] for y in Y) == sum((length(indXB[c])==0 ? 1.0 / length(Z) :  length(indXB[c][findall(inst$Zobserv[indXB[c]nA] == z)])/length(indXB[c]) ) * length(indXA[c])/nA for c in 1:nbX)+ deviationB[z]);

       # - satisfy probabilities of modality A with possible deviation for individuals in base B
       @constraint(groupB, cttransportAinB[y in Y],
           sum(transportB[y,z] for z in Z) == sum( (length(indXA[c])==0 ? 1.0 / length(Y) :  length(indXA[c][findall(inst$Yobserv[indXA[c]] == y)])/length(indXA[c]) ) * length(indXB[c])/nB for c in 1:nbX)+ deviationA[y]);
       # - satisfy probabilities of modality B for individuals in base B
       @constraint(groupB, cttransportBinB[z in Z],
           sum(transportB[y,z] for y in Y) == freqZ[z])

       # - the deviations must sum to zero to conserve a well-defined probability measure
       @constraint(groupB, ctsumdeviationA, sum(deviationA[y] for y in Y) == 0)
       @constraint(groupA, ctsumdeviationB, sum(deviationB[z] for z in Z) == 0)

       # - bound the norm 1 of deviations
       @constraint(groupA, ctabsdevBplus[z in Z], deviationB[z] <= absdevB[z])
       @constraint(groupA, ctabsdevBmoins[z in Z], deviationB[z] >= -absdevB[z])
       @constraint(groupA, ctbounddevB, sum([absdevB[z] for z in Z]) <= maxrelax/2.0)
       @constraint(groupB, ctabsdevAplus[y in Y], deviationA[y] <= absdevA[y])
       @constraint(groupB, ctabsdevAmoins[y in Y], deviationA[y] >= -absdevA[y])
       @constraint(groupB, ctbounddevA, sum([absdevA[y] for y in Y]) <= maxrelax/2.0)

       # Solve the problem
       optimize!(groupA)
       optimize!(groupB)

       # Extract the values of the solution
       transportA_val = [value(transportA[y,z]) for y in Y, z in Z]
       transportB_val = [value(transportB[y,z]) for y in Y, z in Z]
    }



    # Get the individual transport from the group transport
    if (indiv_method == sequential){
        YApred, YBpred = individual_from_group_closest(inst, transportA_val, transportB_val, percent_closest)}
    else if (indiv_method == optimal){
        YApred, YBpred = individual_from_group_optimal(inst, transportA_val, transportB_val, percent_closest);
    }

    # Compute the estimated probability distributions from predictions
    indXA = inst$indXA; indXB = inst$indXB;
    nbX = length(indXA);
    estimatorZA = zeros(nbX,length(Y),length(Z));
    estimatorYB = zeros(nbX,length(Y),length(Z));

    for (x in 1:nbX){
        for (i in indXA[x]){
            estimatorZA[x,inst$Yobserv[i],YBpred[i]] += 1/ length(findall(inst$Yobserv[indXA[x]] == inst$Yobserv[i]));
        }
        for (y in Y){
            if (length(findall(inst$Yobserv[indXA[x]] == y)) == 0){
                estimatorZA[x,y,] = 1/length(Z)*ones(length(Z),1);
            }
        }
        for (i in indXB[x]){
            estimatorYB[x,YApred[i],inst$Zobserv[i+nA]] += 1/ length(findall(inst$Zobserv[indXB[x]nA] == inst$Zobserv[i+nA]));
        }
        for (z in Z){
            if (length(findall(inst$Zobserv[indXB[x]] == z)) == 0){
                estimatorYB[x,:,z] = 1/length(Y)*ones(length(Y),1);
            }
        }
    }

    return(list(Solution(time()-tstart,transportA_val, transportB_val,estimatorZA,estimatorYB))
}
