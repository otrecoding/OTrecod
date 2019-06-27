
#-------------------------------------------------------------------------------
# PB opimisation dans la fonction "individual_from_group_optimal"
# (la derni?re de ce fichier)
#-------------------------------------------------------------------------------


setwd("F:\\inserm1027\\cloeDM\\OT_new")
tab1      = read.table("tab1.txt",sep =" ",header=TRUE)

# Pour que ?a tourne plus vite on ne conserve que 400 sujets (200 base A/200 base B)
tab1      = tab1[c(1:200,5001:5200),]


jointprobaA = jointprobaB = matrix(c(0.0834,0.0834,0.0832,0.0884,0.0826,0.0790,
                                     0.0908,0.0786,0.0806,0.0872,0.0816,0.0812),
                                     ncol = 3,byrow = T)
                                     
# Packages  -----------------------------
library(rdist)
#--> Nouveaux packages pour l'optimisation (remplace linprog qui posait pb)
library(dplyr)
library(ROI)
library(ROI.plugin.symphony)
library(ompr)
library(ompr.roi)
#----------------------------------------


Instance  = function(data_file,norme){

      dat    = data_file
      # number of covariables
      nbcvar = dim(dat)[2] - 3

      # recover the sets of individuals in base 1 and 2
      Base = dat[1:dim(dat)[1],1]
      indA = which(Base == 1)
      indB = which(Base == 2)
      nA = length(indA)
      nB = length(indB)

      # recover the input data
      Xobserv = dat[1:nrow(dat),4:ncol(dat)]
      Yobserv = dat[1:nrow(dat),2]
      Zobserv = dat[1:nrow(dat),3]

      # modify order so that base A comes first and then base B
      # Xobserv = cbind(Xobserv[indA,],Xobserv[indB,])
      Xobserv = rbind(Xobserv[indA,],Xobserv[indB,])
      Yobserv = c(Yobserv[indA],Yobserv[indB])
      Zobserv = c(Zobserv[indA],Zobserv[indB])
      indA = 1:nA;
      indB = (nA+1):(nA+nB)

      # Modify Y and Z so that they go from 1 to the number of modalities
      Y = sort(unique(Yobserv[Yobserv != -1]));
      Z = sort(unique(Zobserv[Zobserv != -1]));
      for (i in 1:length(Y)){

        Yobserv[Yobserv == Y[i]]= i

      }

      #Y = [i for i in 1:length(Y)];
      Y=1:length(Y)
      for (i in 1:length(Z)){

        Zobserv[Zobserv == Z[i]] = i

      }

      #Z = [i for i in 1:length(Z)];
      Z=1:length(Z)


      # list the distinct modalities in A and B
      indY = indZ = list()
      for (m in Y){indY[[m]] = which(Yobserv[1:nA] == m)}
      for (m in Z){indZ[[m]] = which(Zobserv[(nA+1):(nA+nB)] == m)}

      # compute the distance between pairs of individuals in different bases
      # devectorize all the computations to go about twice faster
      # only compute norm 1 here
      a = Xobserv[indA,]
      b = Xobserv[indB,]

      stopifnot(norme %in% c(0,1,2))


      if (norme == 1){

        #   D = pairwise(Cityblock(), a, b, dims=2)
        #  DA = pairwise(Cityblock(), a, a, dims=2)
        #  DB = pairwise(Cityblock(), b, b, dims=2)

        D  = cdist(a,b,"manhattan")
        DA = cdist(a,a,"manhattan")
        DB = cdist(b,b,"manhattan")

      } else if (norme == 2){

        D  = cdist(a,b,"euclidean")
        DA = cdist(a,a,"euclidean")
        DB = cdist(b,b,"euclidean")

      } else if (norme == 0){

        D  = cdist(a,b,"hamming")
        DA = cdist(a,a,"hamming")
        DB = cdist(b,b,"hamming")

      }

      # Compute the indexes of individuals with same covariates
      A     = 1:nA;
      B     = 1:nB;
      nbX   = 0;
      # indXA = numeric(dim(Xval)[1]);
      # indXB = numeric(dim(Xval)[1]));

      Xval  = unique(Xobserv)
      indXA = indXB = rep(0,nrow(Xval))


      X1val = sort(unique(Xobserv[,1]));
      X2val = sort(unique(Xobserv[,2]));
      X3val = sort(unique(Xobserv[,3]));



      # aggregate both bases

      indXA = indXB =  list()

      for (i in  (1:nrow(Xval))){

          # if (i %in% seq(0,nrow(Xval),50)){print(i)} else {}
	    print(i)

          nbX = nbX + 1;
          # x = matrix(0,dim(Xval[2]),1);
          # plut?t: x = rep(0,ncol(Xval)) mais inutile

          # x[,1] = Xval[i,(1:dim(Xval)[2])];
          #x[:,1] = [Xval[i,j] for j in 1:size(Xval,2)];

	    x = Xval[i,]


          if (norme == 1){

            # distA = pairwise(Cityblock(), x, t(Xobserv[A,]), dims=2)
            # distB = pairwise(Cityblock(), x, t(Xobserv[B + nA,]), dims=2)

            distA  = cdist(x,Xobserv[A,],"manhattan")
            distB  = cdist(x,Xobserv[B + nA,],"manhattan")



          } else if (norme == 2){

            distA  = cdist(x,Xobserv[A,],"euclidean")
            distB  = cdist(x,Xobserv[B + nA,],"euclidean")

          } else if (norme == 0){

            distA  = cdist(x,Xobserv[A,],"hamming")
            distB  = cdist(x,Xobserv[B + nA,],"hamming")

          }

          # indXA[nbX] = (distA[1,] < 0.1)
          # indXB[nbX] = (distB[1,] < 0.1)

          indXA[[nbX]] = which(distA < 0.1)
          indXB[[nbX]] = which(distB < 0.1)

      }

      # file_name = base_name(data_file)
      file_name = deparse(substitute(data_file))

      return(list(FILE_NAME =file_name,nA=nA, nB=nB, Xobserv=Xobserv,
                  Yobserv=Yobserv, Zobserv=Zobserv, D=D, Y=Y, Z=Z,
                  indY=indY, indZ=indZ, indXA=indXA, indXB=indXB, DA=DA, DB=DB))
}

#--------------------------------------
# Exemple OK
stock_res = Instance(tab1,norme = 1)
#--------------------------------------



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

#--------------------------------------
# Exemple OK
try1 = average_distance_to_closest(stock_res,1)
#--------------------------------------




inst = stock_res   # Pour pouvoir travailler directement dans la fonction
                   # ci-dessous sans avoir ? l'ex?cuter


# --> PB pour utiliser le code relatif ? l'optimisation via les nouveaux
#     packages

###############################################################################
# Solve an optimization problem to get the individual transport that minimizes
# total distance while satisfying the joint probability computed by the model by
# group
###############################################################################

individual_from_group_optimal=function(inst, jointprobaA, jointprobaB, percent_closest=1.0){


    # Redefine A and B for the model
    A = 1:inst$nA
    B = 1:inst$nB
    Y = inst$Y
    Z = inst$Z
    indY = inst$indY
    indZ = inst$indZ
    nbindY = numeric(0)

    for (y in Y){
    
      nbindY = c(nbindY,length(inst$indY[[y]]))

    }

    nbindZ= numeric(0)
    for (z in Z){
    
      nbindZ = c(nbindZ,length(inst$indZ[[z]]))

    }

    # Create a model for the optimal transport of individuals
    # indiv = Model(solver=IpoptSolver(print_level=4))
    # indiv = Model(with_optimizer(Clp.Optimizer,LogLevel=0))
    # indiv = Model(solver=GurobiSolver(Method=2,LogToConsole=0)); #Presolve=0,Method=2,Crossover=0))


    # Variables
    # - assignA[i][z] : fraction of individual i assigned to modality z
    # @variable(indiv, assignA[i in A, z in Z] >= 0, base_name="assignA")
    # - assignB[j][y] : fraction of individual j assigned to modality y
    # @variable(indiv, assignB[j in B, y in Y] >= 0, base_name="assignB")

    # compute the average distance between the individuals and the modalities of
    # the other base
    CA = matrix(rep(0,inst$nA * length(Z)), nrow = inst$nA,ncol = length(Z))
    CB = matrix(rep(0,inst$nB * length(Y)), nrow = inst$nB,ncol = length(Y))

    for (i in A){
        for (z in Z){
            nbclose      = round(percent_closest*nbindZ[z])
            distance     = sort(inst$D[i,indZ[[z]]])
            CA[i,z]      = sum(distance[1:nbclose])/nbclose
        }
    }
    for (j in B){
        for (y in Y){
            nbclose      = round(percent_closest*nbindY[y])
            distance     = sort(inst$D[indY[[y]],j])
            CB[j,y]      = CB[j,y] + sum(distance[1:nbclose])/nbclose
        }
    }

 # OPTIMISATION JULIA------------------------------------------------------------------
    # Objective: minimize the distance between individuals of A and B
    # @objective(indiv, Min, sum(CA[i,z]*assignA[i,z] for i in A, z in Z)
    #                           + sum(CB[j,y]*assignB[j,y] for j in B, y in Y))

    # Constraints
    # - assign the individuals so as to satisfy the joint probability computed
    #   with the model by group
    #@constraint(indiv, ctjointprobaA[y in Y, z in Z],
    #    sum(assignA[i,z] for i in indY[y]) == jointprobaA[y,z])
    #@constraint(indiv, ctjointprobaB[y in Y, z in Z],
    #    sum(assignB[j,y] for j in indZ[z]) == jointprobaB[y,z])

    # - assign sufficient modality to each individual
    #@constraint(indiv, ctassignA[i in A], sum(assignA[i,z] for z in Z) == 1/(length(A)))
    #@constraint(indiv, ctassignB[j in B], sum(assignB[j,y] for y in Y) == 1/(length(B)))
  #-------------------------------------------------------------------------------------
  
  # --> Equivalent optimisation R (qui ne fonctionne pas pour le moment)
  
  result =  MIPModel() %>%
  
    add_variable(assignA[i,z],  i = A, z = Z,type = "continuous",lb=0) %>%
    add_variable(assignB[j,y],  j = B, y = Y,type = "continuous",lb=0) %>%
    set_objective(sum_expr(CA[i,z]*assignA[i,z], i = A,z=Z) + sum_expr(CB[j,y]*assignB[j,y], j = B,y=Y), "min") %>%
    #set_objective(sum_expr(CA[i,z]*assignA[i,z], i = A,z=Z), "min") %>%
    #add_constraint(sum_expr(assignA[i,z], i=indY[y])   == jointprobaA[y,z], y = Y, z = Z) %>%
    #add_constraint(sum_expr(assignB[j,y], j = indZ[z]) == jointprobaB[y,z], y = Y, z = Z) %>%
    add_constraint(sum_expr(assignA[i,z], z = Z) == 1/(length(A)),i = A) %>%
    add_constraint(sum_expr(assignB[j,y], y = Y) == 1/(length(B)),j = B) %>%
    solve_model(with_ROI(solver = "glpk"))

    solution = get_solution(result, assignA[i,z])
    solution = get_solution(result, assignB[j,y])
  
    # Extract the values of the solution
    # assignA_val = [value(assignA[i,z]) for i in A, z in Z]
    # assignB_val = [value(assignB[j,y]) for j in B, y in Y]
    assignA_val = matrix(solution[1:(length(A)*length(Z))], length(A),length(Z))
    assignB_val = matrix(solution[(length(A)*length(Z)+1):(length(A)*length(Z)+length(B)*length(Y))], length(B),length(Z))

    # Transport the modality that maximizes frequency
    # YBtrans = [findmax([assignA_val[i,z]  for z in Z])[2] for i in A]
    # YAtrans = [findmax([assignB_val[j,y]  for y in Y])[2] for j in B]
    YBtrans = apply(assignA_val,1,which.max)
    YAtrans = apply(assignB_val,1,which.max)

    return(list(YAtrans, YBtrans))
}

