#using JuMP
#using Distances
#using Printf
#using DelimitedFiles
#using DataFrames

#enum DataBase baseA baseB

# Definition and initialization of an Instance structure
#struct Instance
#    name::AbstractString
#    nA::Int64
#    nB::Int64
#    Xobserv::Array{Float64,2}
#    Yobserv::Array{Int64,1}
#    Zobserv::Array{Int64,1}
#    D::Array{Float64,2}
#    Y::Array{Int64,1}
#    Z::Array{Int64,1}
#    indY::Dict{Int64,Array{Int64,1}}
#    indZ::Dict{Int64,Array{Int64,1}}
#    indXA:: Dict{Int64,Array{Int64}}; # indexes of subjects of A with given X value
#    indXB:: Dict{Int64,Array{Int64}}; # indexes of subjects of B with given X value
#    DA::Array{Float64,2}
#    DB::Array{Float64,2}
setwd("C:\\Users\\vagares\\Google Drive\\3_Recherche\\R1_Recherche\\2_INSA\\DataMerging_4\\Programmes\\Data\\data\\SR-0.8")
data_file = read.table("tab1.txt",sep =" ",header=TRUE)
Instance  = function(data_file,norme){
      data = data_file
      # number of covariables
      nbcvar = dim(data)[2] - 3

      # recover the sets of individuals in base 1 and 2
      base = data[1:dim(data)[1],1]
      indA = which(base == 1)
      indB = which(base == 2)
      nA = length(indA)
      nB = length(indB)

      # recover the input data
      Xobserv = data[1:dim(data)[1], 4:dim(data)[2]]
      Yobserv = data[1:dim(data)[1], 2]
      Zobserv = data[1:dim(data)[1], 3]

      # modify order so that base A comes first and then base B
      Xobserv = cbind(Xobserv[indA,],Xobserv[indB,])
      Yobserv = c(Yobserv[indA],Yobserv[indB])
      Zobserv = c(Zobserv[indA],Zobserv[indB])
      indA = 1:nA;
      indB = nA+1:nA+nB;

      # Modify Y and Z so that they go from 1 to the number of modalities
      Y = sort(unique(Yobserv[Yobserv != -1]));
      Z = sort(unique(Zobserv[Zobserv != -1]));
      for (i in 1:length(Y))
          {Yobserv[Yobserv == Y[i]]= i}

      #Y = [i for i in 1:length(Y)];
      Y=1:length(Y)
      for (i in 1:length(Z))
          {Zobserv[Zobserv == Z[i]] = i}

      #Z = [i for i in 1:length(Z)];
      Z=1:length(Z)
      # list the distinct modalities in A and B
      indY = numeric(0);indZ = numeric(0)
      for (m in Y){indY = c(indY,which(Yobserv[1:nA] == m))}
      for (m in Z){indY = c(indZ,which(Zobserv[(nA+1):(nA+nB)] == m))}

      # compute the distance between pairs of individuals in different bases
      # devectorize all the computations to go about twice faster
      # only compute norm 1 here
      a = t(Xobserv[indA,])
      b = t(Xobserv[indB,])
      if (norme == 1)
         {D = pairwise(Cityblock(), a, b, dims=2)
          DA = pairwise(Cityblock(), a, a, dims=2)
          DB = pairwise(Cityblock(), b, b, dims=2)}
      else if (norme == 2)
          {D = pairwise(Euclidean(), a, b, dims=2)
          DA = pairwise(Euclidean(), a, a, dims=2)
          DB = pairwise(Euclidean(), b, b, dims=2)}
      else if (norme == 0)
          {D = pairwise(Hamming(), a, b, dims=2)
          DA = pairwise(Hamming(), a, a, dims=2)
          DB = pairwise(Hamming(), b, b, dims=2) }
      }

      # Compute the indexes of individuals with same covariates
      A = 1:nA;
      B = 1:nB;
      nbX = 0;
      indXA = numeric(dim(Xval[1]);
      indXB = numeric(dim(Xval[1]);

      X1val = sort(unique(Xobserv[,1]));
      X2val = sort(unique(Xobserv[,2]));
      X3val = sort(unique(Xobserv[,3]));

      Xval = as.matrix(Xobserv)

      # aggregate both bases
      for (i in  (1:dim(Xval[1])))
          {nbX = nbX + 1;
          x = matrix(0,dim(Xval[2]),1);
          x[,1] = Xval[i,(1:dim(Xval)[2])];
          if (norme == 1)
              {distA = pairwise(Cityblock(), x, t(Xobserv[A,]), dims=2);
              distB = pairwise(Cityblock(), x, t(Xobserv[B + nA,]), dims=2)}
          else {if (norme == 2)
              distA = pairwise(Euclidean(), x, t(Xobserv[A,]), dims=2);
              distB = pairwise(Euclidean(), x, t(Xobserv[B + nA,]), dims=2)}
          else if (norme == 0)
              {distA = pairwise(Hamming(), x, t(Xobserv[A,]), dims=2);
              distB = pairwise(Hamming(), x, t(Xobserv[B + nA,]), dims=2);}
          indXA[nbX] = (distA[1,] < 0.1);
          indXB[nbX] = (distB[1,] < 0.1)}


      file_name = base_name(data_file)
      return(list(file_name=file_name,nA=nA, nB=nB, Xobserv=Xobserv, Yobserv=Yobserv, Zobserv=Zobserv, D=D, Y=Y, Z=Z, 
                  indY=indY, indZ=indZ, indXA=indXA, indXB=indXB, DA=DA, DB=DB))
}
#revoir acolades

aggregate_per_covar_mixed = function(Inst, norme=1, aggregate_tol)
{
    A = 1:(inst$nA);
    B = 1:(inst$nB);

    # initialization of the structures
    nbX = 0;
    #indXA = Dict{Int64,Array{Int64}}();
    #indXB = Dict{Int64,Array{Int64}}();
    notaggA = A;
    notaggB = B;

    # aggregate until every individual in base A is aggregated
    while (length(notaggA)>0){
        nbX = nbX+1;
        ind = notaggA[1];
        isinset =  which(inst$DA[ind,notaggA] < aggregate_tol);
        indXA[nbX] = notaggA[isinset];
        notaggA = subset(notaggA, !sinset);
        isinset =  which(inst$D[ind,notaggB] < aggregate_tol);
        indXB[nbX] = notaggB[isinset];
        notaggB = subset(notaggB, !isinset);
    }

    # complete the aggregation with the individuals of base B that are not aggregated yet
    while (length(notaggB)>0){
        nbX = nbX+1;
        ind = notaggB[1];
        isinset = which(inst$DB[ind,notaggB] < aggregate_tol);
        indXB[nbX] = notaggB[isinset];
        indXA[nbX] = numeric(0);
        notaggB = subset(notaggB, !isinset);
    }

    return(list(indXA=indXA,indXB=indXB))
}

###############################################################################
# Compute a bound on the average prediction error in each base
# The bound is computed as the expected prediction error assuming that the
# distribution of Z in base A (and that of Y in base B) is known, and the
# prediction done with the value that maximizes the probability
###############################################################################
bound_prediction_error = function(Instance, norme=1, aggregate_tol=0.5)
{
    # Local redefinitions of parameters of  the instance
    nA = inst$nA;
    nB = inst$nB;
    Y = inst$Y;
    Z = inst$Z;

    # compute the bound in base A
    boundpredZA = 0.0;
    for (x  in (1:length(inst$indXA))){
        for (mA in Y){
            indwithmA = inst$indXA[x][inst$Yobserv[inst$indXA[x]] == mA];
            nindiv = length(indwithmA);}
            if (nindiv == 0)
                {
            for (mB in Z){
              indwithmB = indwithmA[inst$Zobserv[indwithmA] == mB];
                boundpredZA = length(indwithmB)/nindiv * (1-length(indwithmB)/nindiv) * nindiv/ nA;
            }
            }
    }
    # print("Bound on average prediction error in A : %.1f %%\n", 100.*boundpredZA);

    # compute the bound in base B
    boundpredYB = 0.0;
    for (x  in 1:length(inst$indXB){
        for (mB in Z)
            {indwithmB = inst$indXB[x][inst$Zobserv[inst$indXB[x][,nA] == mB]+ nA;
            nindiv = length(indwithmB)}
            if (nindiv == 0)
            {for (mA in Y){
                indwithmA = indwithmB[inst$Yobserv[indwithmB] == mA];
                boundpredYB = boundpredYB+ length(indwithmA)/nindiv * (1-length(indwithmA)/nindiv) * nindiv/ nB;
            }
        }
    }
    # print("Bound on average prediction error in B : %.1f %%\n", 100.*boundpredYB);

    return(list(boundpredZA=boundpredZA, boundpredYB= boundpredYB))
}

###############################################################################
# Return the empirical cardinality of the joint occurrences of (C=x,Y=mA,Z=mB)
# in both bases
###############################################################################
empirical_distribution=function(inst, norme=0, aggregate_tol=0.5){

    # Local redefinitions of parameters of  the instance
    nA = inst$nA;
    nB = inst$nB;
    A = (1:inst$nA);
    B = (1:inst$nB);
    Y = inst$Y;
    Z = inst$Z;

    # aggregate the individuals per covariate
    nbX = length(inst$indXA);

    # count the cardinality of occurrence of each triplet (x,mA,mB) in both bases
    cardA_c_mA_mB = zeros(nbX, length(Y), length(Z));
    for (x = 1:nbX){
        for (i in inst$indXA[x])
            {cardA_c_mA_mB[x,inst$Yobserv[i],inst$Zobserv[i]] += 1;
          }
    }
    cardB_c_mA_mB = zeros(nbX, length(Y), length(Z));
    for (x = 1:nbX){
        for (in inst$indXB[x]){
            cardB_c_mA_mB[x,inst$Yobserv[i+nA],inst$Zobserv[i+nA]] += 1;
        }
    }

    return(list(cardA_c_mA_mB=cardA_c_mA_mB, cardB_c_mA_mB=cardB_c_mA_mB))
}


###############################################################################
# Display information about the distance between the modalities
###############################################################################

disp_inst_info = function(inst)

    # local definitions
    nA = inst$nA
    nB = inst$nB
    A = (1:inst$nA)
    B = (1:inst$nB)
    Y = inst$Y
    Z = inst$Z
    indY = inst$indY
    indZ = inst$indZ


    println("\n#################################################################")
    println("INFORMATION ABOUT THE INSTANCE")
    println("#################################################################\n")

    # return indicators about the original density of the modalities
    print("Average distance between objects of base 1: ")
    print("%.2f\n", 10*sum([DA[i,j] for i in A, j in A])/(nA^2))
    print("Average distance between objects of base 2: ")
    print("%.2f\n", 10*sum([DB[i,j] for i in B, j in B])/(nB^2))
    print("Crossed average distance between objects of base 1 and 2: ")
    print("%.2f\n", 10*sum([inst$D[i,j] for i in A, j in B])/(nA*nB))

    # restrict the average distance to the 10% closest individuals
    println("\nAverage distance between objects per modality")
    println("Modalities of Y, individuals of A:")
    percent_closest = 0.1
    for (y1 in Y){
        for (y2 in Y){
            if (y1 > y2){
               
            avg = avg_distance_closest(inst,baseA,baseA,baseA,y1,y2,1.0)
            print("\tModalities %d and %d : %.2f\n", y1, y2, 10*avg)
            avg = avg_distance_closest(inst,baseA,baseA,baseA,y1,y2,percent_closest)
            print("\t\trestricted to the %.1f %% closest: %.2f\n", 100*percent_closest, 10*avg)
        }}}
    
    println("\nModalities of Z, individuals of B:")
    for (z1 in Z){
        for (z2 in Z){
            if (z1 > z2){
                
            avg = avg_distance_closest(inst,baseB,baseB,baseB,z1,z2,1.0)
            print(paste("Modalities ",z1," and ",z2," : ", 10*avg))
            avg = avg_distance_closest(inst,baseB,baseB,baseB,z1,z2,percent_closest)
            print(paste("trestricted to the",100*percent_closest," closest: ",10*avg)
        }}}
    
    println("\nModalities of Y, crossed bases:")
    for (y1 in Y){
        for (y2 in Y){
            avg = avg_distance_closest(inst,baseA,baseB,baseA,y1,y2,1.0)
            print(paste("Modalities ",y1 ," and ",y2," : ", 10*avg))
            avg = avg_distance_closest(inst,baseA,baseB,baseA,y1,y2,percent_closest)
            print(paste("restricted to the ",100*percent_closest," closest: ",10*avg))
        }
    }
    println("\nModalities of Z, crossed bases:")
    for (z1 in Z){
        for (z2 in Z){
            avg = avg_distance_closest(inst,baseA,baseB,baseB,z1,z2,1.0)
            print(paste("Modalities ", z1," and ", z2," : ", 10*avg)
            avg = avg_distance_closest(inst,baseA,baseB,baseB,z1,z2,percent_closest)
            print("restricted to the ",100.0*percent_closest," closest: ", 10*avg)
        }
    }
}

###############################################################################
# Compute the average distance between individuals of base1 with modality m1
# for outcome and individuals of base2 with modality m2 for outcome
# Consider only the percent_closest individuals in the computation of the
# distance
###############################################################################

avg_distance_closest = function(inst, base1, base2, outcome, m1, m2, percent_closes){
    #inst=instance
    #base1=DataBase
     #base2=DataBase
   #outcome=DataBase
    # indices of individuals of base A with given outcomes in base B (and reciprocally)
  indZinA=numeric(0)
  for (z in inst$Z){
    indZinA = c(indZinA,which(inst$Zobserv[1:inst$nA] == z))}
  indYinB=numeric(0)
  for (y in inst$Y){
    indYinB =c(indYinB,which(inst$Yobserv[inst$nA+1:n1+nB] == y))}
    ind1 = base1 == baseA ? (outcome == baseA ? inst$indY[m1] : indZinA[m1]) : (outcome == baseA ? indYinB[m1] : inst$indZ[m1])
    ind2 = base2 == baseA ? (outcome == baseA ? inst$indY[m2] : indZinA[m2]) : (outcome == baseA ? indYinB[m2] : inst$indZ[m2])

    # select the distance matrix deping on the base
    # println(base1)
    # println(base2)
    # println(outcome)
    # println(typeof(Dscaled))
    D = base1 == baseA ? ( base2==baseA ? inst$DA : inst$D) : ( base2==baseA ? inst$D : inst$DB)

    # swap the two sets of indices if base1=baseB and base2=baseA
    if (base1==baseB && base2==baseA){
        ind = ind1
        ind1 = ind2
        ind2 = ind1
    }

    # compute the average distance between the individuals in ind1 and the
    # percent_closest in ind2 and reciprocally
    avg = 0.0
    for (i in ind1){
        nbclose = round(Int,percent_closest*length(ind2))
        distance=numeric(0)
        for (j in ind2){
        distance = c(distance,sort([D[i,j])}
        avg = avg  + sum(distance[1:nbclose])/nbclose/length(ind1)/2
    }
    for (j in ind2){
        nbclose = round(Int,percent_closest*length(ind1))
        distance =numeric(0)
        for (i in ind1){
          distance = c(distance,sort([D[i,j])}
        avg += sum(distance[1:nbclose])/nbclose/length(ind2)/2
    }

    return(avg)
}



#mutable struct Solution
#tsolve::Float64  # solution time
#    jointYZA::Array{Float64,2} # joint distribution of Y and Z in A
#    jointYZB::Array{Float64,2} # joint distribution of Y and Z in B
#    estimatorZA::Array{Float64,3} # estimator of probability of Z for individuals in base A
#    estimatorYB::Array{Float64,3} # estimator of probability of Y for individuals in base B
#    errorpredZA::Float64
#    errorpredYB::Float64
#    errorpredavg::Float64
#    errordistribZA::Float64
#    errordistribYB::Float64
#    errordistribavg::Float64

    Solution(t,jointYZA,jointYZB) = new(t,jointYZA,jointYZB);
    Solution(t,jointYZA,jointYZB,estimatorZA,estimatorYB) = new(t,jointYZA,jointYZB,estimatorZA,estimatorYB);



###############################################################################
# Compute prediction errors in a solution
###############################################################################
compute_pred_error= function(inst, sol, proba_disp=true, mis_disp=false, full_disp=false){

    A = 1:inst$nA
    B = 1:inst$nB

    # display the transported and real modalities
    if (full_disp){
        println("Modalities of base 1 individuals:")
        for (i in A){
            print(paste("Index: " , i, ", real value: ", inst$Zobserv[i], ", transported value: ", sol.predZA[i]))
        }
        # display the transported and real modalities
        print("Modalities of base 2 individuals:")
        for (j in B){
            print(paste("Index: " , j, ", real value: ", inst$Yobserv[inst$nA+j], ", transported value: ", sol.predYB[j]))
        }
    }

    # Count the number of mistakes in the transport

    # Base 1
    nbmisA = 0
    misA = numeric(0)
    for (i in A){
        if (sol.predZA[i] != inst$Zobserv[i]){
          nbmisA = nbmisA+ 1
          push!(misA, i )
        }
    }

    # Base 2
    nbmisB = 0
    misB = numeric(0)
    for (j in B){
        if (sol.predYB[j] != inst$Yobserv[inst$nA+j]){
          nbmisB = nbmisB+ 1
          push!(misB, j)
        }
    }


    if (proba_disp){
        if (nbmisA == 0){
            print("No mistake in the transport of base A")}
        else
            {print(paste("Probability of error in base A: ", 100.0 * nbmisA / inst$nA));
            if (mis_disp){
                print("Indices with mistakes in base A:", misA)
            }
        }

        if (nbmisB == 0){
            print("No mistake in the transport of base B")}
        else
            {print(paste("Probability of error in base B: ", 100.0 * nbmisB/inst$nB));
            if (mis_disp)
                {print(paste("Indices with mistakes in base 2:", misB))}
            }
        }
    

    sol.errorpredZA, sol.errorpredYB = nbmisA/inst$nA, nbmisB/inst$nB;
    sol.errorpredavg = (inst$nA*sol.errorpredZA + inst$nB*sol.errorpredYB)/(inst$nA+inst$nB);

    return(sol)
}

###############################################################################
# Compute errors in the conditional distributions of a solution
###############################################################################
compute_distrib_error = function(inst, sol, empiricalZA, empiricalYB){

    nA = inst$nA ;  nB = inst$nB ; Y = inst$Y; Z = inst$Z;
    nbX = length(inst$indXA);

    for (x in 1:nbX){
      for (y in Y){
    sol.errordistribZA = sum(length(inst$indXA[x][inst$Yobserv[inst$indXA[x]] == y)])/nA * sum(max(sol.estimatorZA[x,y,] - empiricalZA[x,y,],0))}} 
    for (x in 1:nbX){
      for (z in Z){
    sol.errordistribYB = sum(length(inst$indXB[x][inst$Zobserv[inst$indXB[x]+nA] == z)])/nB * sum(max.(sol.estimatorYB[x,,z] - empiricalYB[x,,z],0))}}
    sol.errordistribavg = (nA * sol.errordistribZA + nB * sol.errordistribYB)/(nA+nB);

    return(sol)
}

###############################################################################
# Compute the cost between pairs of outcomes as the average distance between
# covariations of individuals with these outcomes, but considering only the
# percent closest neighbors
###############################################################################

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
        for (i in indY[y]){
            for (z in Z){
                nbclose = max(round(Int,percent_closest*length(indZ[z])),1);
                distance=numeric(0)
                for (j in indZ[z]{}
                distance = c(distance,sort(inst$D[i,j]))
                DindivA[i,z] =  sum(distance[1:nbclose])/nbclose;
                Davg[y,z] += sum(distance[1:nbclose])/nbclose/length(indY[y])/2.0;
            }
        }
    }
    for (z in Z){
        for (j in indZ[z]){
            for (y in Y){
                nbclose = max(round(Int,percent_closest*length(indY[y])),1);
                distance = numeric(0)
                for (i in indY[y]]){
                distance = c(distance,sort(inst$D[i,j])) }
                DindivB[j,y] = sum(distance[1:nbclose])/nbclose;
                Davg[y,z] += sum(distance[1:nbclose])/nbclose/length(indZ[z])/2.0;
            }
        }
    }

    return(list(Davg=Davg, DindivA=DindivA, DindivB=DindivB)
}
