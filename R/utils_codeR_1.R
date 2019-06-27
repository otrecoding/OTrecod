
# Espace de travail ? charger depuis la dropbox qui contient tous les objets
setwd("C:\\Users\\vagares\\Documents\\Dropbox\\Conversion Julia")
load("otgroup_julia.RData")
ls()
#tab1 = tab1[c(1:200,5001:5200),]
# Packages utiles
library(rdist); library(linprog); library(mvtnorm) ; library(lpSolveAPI) ;



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
# Example OK
# tab1bis   = tab1[c(1:500,5001:5500),]
# stock_res = Instance(tab1bis,norme = 0)
setwd("F:\\inserm1027\\cloeDM\\OT_new")
tab1      = read.table("tab1.txt",sep =" ",header=TRUE)
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
average_distance = try1
#--------------------------------------



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
    
    return(list(YAtrans, YBtrans))
}    
    

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
    nbindY = c(nbindY,length(inst$indY[[y]]))}
    nbindZ= numeric(0)
    for (z in Z){
      nbindZ = c(nbindZ,length(inst$indZ[[z]]))}
    
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

    # Objective: minimize the distance between individuals of A and B
    #@objective(indiv, Min, sum(CA[i,z]*assignA[i,z] for i in A, z in Z)
    #                        + sum(CB[j,y]*assignB[j,y] for j in B, y in Y))
    
    Min = c(as.numeric(t(CA)),as.numeric(t(CB)))
    

    # Constraints
    # - assign the individuals so as to satisfy the joint probability computed
    #   with the model by group
    # @constraint(indiv, ctjointprobaA[y in Y, z in Z],
    #    sum(assignA[i,z] for i in indY[y]) == jointprobaA[y,z])
        
    
    c1 = matrix(0,length(Z)*length(Y),length(A)*length(Z)+length(B)*length(Y))
    b1 = numeric(length(Z)*length(Y))
    for (y in Y){
        for (z in Z){
           for (i in indY[y]){
                c1[length(Z)*(y-1)+ z,length(Z)*(i-1)+ z] = 1
                b1[length(Z)*(y-1)+ z] = jointprobaA[y,z]
           }
        }
    }

    s1 = rep("=",length(Z)*length(Y))
        
        
    # @constraint(indiv, ctjointprobaB[y in Y, z in Z],
    #    sum(assignB[j,y] for j in indZ[z]) == jointprobaB[y,z])

    c2 = matrix(0,length(Z)*length(Y),length(A)*length(Z)+length(B)*length(Y))
    b2 = numeric(length(Z)*length(Y))
    for (y in Y){
         for (z in Z){
            for (j in indZ[z]){
                 c1[length(Z)*(y-1)+ z,length(A)*length(Z) + length(Y)*(j-1)+ y] = 1
                 b1[length(Z)*(y-1)+ z] = jointprobaB[y,z]
            }
         }
    }

    s2 = rep("=",length(Z)*length(Y))



    # - assign sufficient modality to each individual
    # @constraint(indiv, ctassignA[i in A], sum(assignA[i,z] for z in Z) == 1/(length(A)))
    
    
    c3 = matrix(0,length(A),length(A)*length(Z)+length(B)*length(Y))
    b3 = numeric(length(A))
    for (i in A){
      for (z in Z){
          c3[i,length(Z)*(i-1)+ z] = 1
          b3[i] = 1/(length(A))
        }
      }
    
    
    s3 = rep("=",length(A))
    

    
    # @constraint(indiv, ctassignB[j in B], sum(assignB[j,y] for y in Y) == 1/(length(B)))

    
    c4 = matrix(0,length(B),length(A)*length(Z)+length(B)*length(Y))
    b4 = numeric(length(B))
    for (j in B){
      for (y in Y){
        c4[j,length(Y)*(j-1)+ y] = 1
        b4[i] = 1/(length(B))
      }
    }
    
    
    s4 = rep("=",length(B))
    nc = dim(c1)[1] + dim(c2)[1] + dim(c3)[1]+ dim(c4)[1]
    nv = length(A)*length(Z)+length(B)*length(Y)
    lps.model <- make.lp(nc, nv)
    for (i in (1:dim(c1)[1])){
    xt = c1[i,]
    st = s1[i]
    bt = b1[i]
    add.constraint(lps.model, xt, st, bt)}
   
    for (i in (1:dim(c2)[1])){
      xt = c2[i,]
      st = s2[i]
      bt = b2[i]
      add.constraint(lps.model, xt, st, bt)}
    
    for (i in (1:dim(c3)[1])){
      xt = c3[i,]
      st = s3[i]
      bt = b3[i]
      add.constraint(lps.model, xt, st, bt)}
    
    for (i in (1:dim(c4)[1])){
      xt = c4[i,]
      st = s4[i]
      bt = b4[i]
      add.constraint(lps.model, xt, st, bt)}
    set.objfn(lps.model,Min)
    solution = solve(lps.model)
    solution = get.primal.solution(lps.model)
    # Solve the problem
    # optimize!(indiv)
    
    c     = rbind(c1,c2,c3,c4)
    b     = c(b1,b2,b3,b4)
    s     = c(s1,s2,s3,s4)
    indiv = solveLP(cvec=Min, bvec=b, Amat=c, const.dir = s,lpSolve=TRUE)
 

    # Extract the values of the solution
    # assignA_val = [value(assignA[i,z]) for i in A, z in Z]
    # assignB_val = [value(assignB[j,y]) for j in B, y in Y]
    assignA = matrix(solution[1:length(A)*length(Z)], length(A),length(Z))
    assignB=  matrix(solution[(length(A)*length(Z)+1):(length(A)*length(Z)+length(B)*length(Y))], length(B),length(Y))

    # Transport the modality that maximizes frequency
    # YBtrans = [findmax([assignA_val[i,z]  for z in Z])[2] for i in A]             
    # YAtrans = [findmax([assignB_val[j,y]  for y in Y])[2] for j in B]
    YBtrans = apply(assignA,1,which.max)
    YAtrans = apply(assignB,1,which.max)

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
    C = average_distance_to_closest(inst, percent_closest)[[1]];
    
    # G: ----> Jusqu'ici OK !!
    #---------------------
    
    if (maxrelax == 0.0){
      
      # G: Peut-on supprimer tous les commentaires ci-dessous ???
      #----------------------------------------------------------
    
      # min = function(transport){sum(as.vector(C)*transport)}
      # 
      # maxi = c(length(val1),length(val2))
      # 
      # 
      # S = matrix(ncol = maxi[1]*maxi[2],nrow = maxi[1]+maxi[2])
      # 
      # 
      # seq1 = rep(1,maxi[2])
      # seq0 = rep(0,maxi[2])
      # seq2 = c(1,rep(0,maxi[2]-1))
      # 
      # 
      # S[1,]         = L1 = c(seq1,rep(seq0,maxi[1]-1))
      # S[maxi[1]+1,] = L2 = rep(seq2,maxi[1])
      # 
      # 
      # for (i in setdiff(2:nrow(S),maxi[1]+1)){
      #   
      #   if (i <= maxi[1]){
      #     
      #     L1    = c(seq0,L1)
      #     S[i,] = L1[1:ncol(S)]
      #   }
      #   
      #   else {
      #     
      #     L2    = c(0,L2)
      #     S[i,] = L2[1:ncol(S)]
      #     
      #   }
      # }
      # 
      # 
      # CVal = as.numeric(c(table(dat[dat[,1]==1,2])/nrow(dat[dat[,1]==1,]),
      #                     table(dat[dat[,1]==2,3])/nrow(dat[dat[,1]==2,])))
      # 
      # res = solveLP(cvec=Risk, bvec=CVal, Amat=S, const.dir = rep( "==",
      #                                                              length(CVal)),lpSolve=TRUE)
      # transport = matrix(res$solution, lenth(Y),lenght(Z))
      # Extract the values of the solution
     
        # Create a model for the optimal transport of individuals
        # group = Model(solver=GurobiSolver(Method=2, LogToConsole=0));
       



        # Objective: minimize the distance between individuals of A and B
        
        # G:---> Ajout de algo1 ici
        
        Min = c(as.numeric(C),rep(0,2*length(Y)+2*length(Z)))
        #VERIFIER QUE C SOIT DNS LE MEME ORDRE que la solution x transport
        # --> G: je pense que c'est bon
        
        
        # Constraints
        # - satisfy probabilities of modality A
        #constraint(group, cttransportA[y in Y], sum(transport[y,z] for z in Z) == freqY[y] + deviationA[y]);

        S = matrix(0,length(Y),length(Y)*length(Z))
        b1 = numeric(length(Y))
           for (y in Y){
             for (z in Z){
            S[y,length(Z)*(y-1)+ z] = 1
            b1[y] = freqY[y]
          }
        }
        
        c1 = cbind(S,diag(length(Y))*(-1),matrix(0,length(Y),length(Y)),matrix(0,length(Y),length(Z)),matrix(0,length(Y),length(Z)))
        s1 = rep("==",length(Y))
        
        # - satisfy probabilities of modality B
        #@constraint(group, cttransportB[z in Z], sum(transport[y,z] for y in Y) == freqZ[z] + deviationB[z]);
        
        S = matrix(0,length(Z),length(Y)*length(Z))
        b2 = numeric(length(Z))
        for (z in Z){
          for (y in Y){
            S[z,length(Z)*(y-1)+ z] = 1
            b2[y] = freqZ[z]
          }
        }
        
       
        c2 = cbind(S,matrix(0,length(Z),length(Y)),matrix(0,length(Z),length(Y)),diag(length(Z))*(-1),matrix(0,length(Z),length(Z)))
        b2 = freqZ
        s2 = rep("==",length(Z))

        # - the deviations must sum to zero to conserve a well-defined probability measure
        #@constraint(group, ctsumdeviationA, sum(deviationA[y] for y in Y) == 0)
        # PB de dimension ici !!! c3 = cbind(matrix(0,1,length(Y)*length(Z)),rep(1,length(Y)),matrix(0,1,length(Y)),matrix(0,1,length(Z)),matrix(0,1,length(Z)))
        # --> G: Je corrige par ligne ci-dessou: (est-ce ok ?)
        c3 = c(matrix(0,1,length(Y)*length(Z)),rep(1,length(Y)),matrix(0,1,length(Y)),matrix(0,1,length(Z)),matrix(0,1,length(Z)))
        b3 = rep(0,1)
        s3 = rep("==",1)
        
        #@constraint(group, ctsumdeviationB, sum(deviationB[z] for z in Z) == 0)
        # IDEM ICI c4 = cbind(matrix(0,1,length(Y)*length(Z)),matrix(0,1,length(Y)),matrix(0,1,length(Y)),rep(1,length(Z)),matrix(0,1,length(Z)))
        c4 = c(matrix(0,1,length(Y)*length(Z)),matrix(0,1,length(Y)),matrix(0,1,length(Y)),rep(1,length(Z)),matrix(0,1,length(Z)))
        b4 = rep(0,1)
        s4 = rep("==",1)
        # - bound the norm 1 of deviations
        #@constraint(group, ctabsdevBplus[z in Z], deviationB[z] <= absdevB[z]);
        c5 = cbind(matrix(0,length(Z),length(Y)*length(Z)),matrix(0,length(Z),length(Y)),matrix(0,length(Z),length(Y)),diag(length(Z)),diag(length(Z))*(-1))
        b5 = rep(0,length(Z))
        s5 = rep("<=",length(Z))
        #@constraint(group, ctabsdevBmoins[z in Z], deviationB[z] >= -absdevB[z]);
        c6 = cbind(matrix(0,length(Z),length(Y)*length(Z)),matrix(0,length(Z),length(Y)),matrix(0,length(Z),length(Y)),diag(length(Z)),diag(length(Z)))
        b6 = rep(0,length(Z))
        s6 = rep(">=",length(Z))
        #@constraint(group, ctbounddevB, sum([absdevB[z] for z in Z]) <= maxrelax/2.0);
        c7 = c(matrix(0,1,length(Y)*length(Z)),matrix(0,1,length(Y)),matrix(0,1,length(Y)),matrix(0,1,length(Z)),rep(1,length(Z)))
        b7 = rep(maxrelax/2.0,1)
        s7 = rep("<=",1)
        #@constraint(group, ctabsdevAplus[y in Y], deviationA[y] <= absdevA[y]);
        c8 = cbind(matrix(0,length(Y),length(Y)*length(Z)),diag(length(Y)),diag(length(Y))*(-1),matrix(0,length(Y),length(Z)),matrix(0,length(Y),length(Z)))
        b8 = rep(0,length(Y))
        s8 = rep("<=",length(Y))
        #@constraint(group, ctabsdevAmoins[y in Y], deviationA[y] >= -absdevA[y]);
        c9 = cbind(matrix(0,length(Y),length(Y)*length(Z)),diag(length(Y)),diag(length(Y)),matrix(0,length(Y),length(Z)),matrix(0,length(Y),length(Z)))
        b9 = rep(0,length(Y))
        s9 = rep(">=",length(Y))
        #@constraint(group, ctbounddevA, sum([absdevA[y] for y in Y]) <= maxrelax/2.0);
        c10 = c(matrix(0,1,length(Y)*length(Z)),matrix(0,1,length(Y)),rep(1,length(Y)),matrix(0,1,length(Z)),matrix(0,1,length(Z)))
        b10 = rep(maxrelax/2.0,1)
        s10 = rep("<=",1)
        # Solve the problem
        c= rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)
        b= c(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
        s= c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
        #optimize!(group);
        group = solveLP(cvec=Min, bvec=b, Amat=c, const.dir = s,lpSolve=TRUE)
                                                              
        # Variables
        # - transport[y,z] : joint probability of modalities y and z
        #@variable(group, transport[y in Y, z in Z] >= 0, base_name="transport");
        transport = matrix(group$solution[1:length(Y)*lenght(Z)], length(Y),lenght(Z))
        # - deviationA[y]: deviation of the probability mass of Y from
        #   that observed in base A
        #@variable(group, deviationA[y in Y], base_name="deviationA")
        deviationA = group$solution[(length(Y)*lenght(Z)+1):(length(Y)*lenght(Z)+length(Y))]
        # - absdevA[y]: absolute value of the deviation of the probability mass of
        #   Y from that observed in base A
        #@variable(group, absdevA[y in Y] >= 0, base_name="absdevA")
        absdevA = group$solution[(length(Y)*lenght(Z)+length(Y)+1):(length(Y)*lenght(Z)+2*length(Y))]
        # - deviationB[z]: deviation of the probability mass of Z from
        #   that observed in base B
        #@variable(group, deviationB[z in Z], base_name="deviationB")
        deviationB = group$solution[(length(Y)*lenght(Z)+2*length(Y)+1):(length(Y)*lenght(Z)+2*length(Y)+length(Z))]
        
        # - absdevB[z]: absolute value of the deviation of the probability mass of
        #  Z from that observed in base B
        #@variable(group, absdevB[z in Z] >= 0, base_name="absdevB")
        absdevB = group$solution[(length(Y)*lenght(Z)+2*length(Y)+length(Z)+1):(length(Y)*lenght(Z)+2*length(Y)+2*length(Z))]

        
        
        
        
        
        
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

        }


    # Otherwise, the empirical distribution constraints are relaxed and the
    # model can be decomposed by database
    else{
      #groupA = Model(with_optimizer(Clp.Optimizer))
      # groupA = Model(solver=GurobiSolver(Method=2,LogToConsole=0))
      # Create a model for the optimal transport of individuals
      
      #@variable(groupA, transportA[y in Y, z in Z] >= 0, base_name="transportA")
      
      
      
      #@variable(groupA, deviationB[z in Z], base_name="deviationB")
      
      # - absdevB[z]: absolute value of the deviation of the probability mass of
      #   z in base A from that observed in base B
      #@variable(groupA, absdevB[z in Z] >= 0, base_name="absdevB")
      
      # Objective: minimize the distance between individuals of A and B
      #@objective(groupA, Min, sum(C[y,z]*transportA[y,z] for y in Y, z in Z))
      Min = c(as.numeric(C),rep(0,2*length(Y)+2*length(Z)))
      
      #@constraint(groupA, cttransportAinA[y in Y], sum(transportA[y,z] for z in Z) == freqY[y])
      
      S = matrix(0,length(Y),length(Y)*length(Z))
      b1 = numeric(length(Y))
      for (y in Y){
        for (z in Z){
          S[y,length(Z)*(y-1)+ z] = 1
          b1[i] = freqY[y]
        }
      }
      
      
      c1 = cbind(S,diag(length(Y))*(-1),matrix(0,length(Y),length(Y)),matrix(0,length(Y),length(Z)),matrix(0,length(Y),length(Z)))
      s1 = rep("== ",length(Y))
      # - satisfy probabilities of modality B with possible deviation for individuals in base A
      #@constraint(groupA, cttransportBinA[z in Z],
      #            sum(transportA[y,z] for y in Y) == sum((length(indXB[c])==0 ? 1.0 / length(Z) :  length(indXB[c][findall(inst$Zobserv[indXB[c]nA] == z)])/length(indXB[c]) ) * length(indXA[c])/nA for c in 1:nbX)+ deviationB[z]);
      
      
      S = matrix(0,length(Z),length(Y)*length(Z))
      b2 = numeric(length(Z))
      for (z in Z){
        for (y in Y){
          S[z,length(Z)*(y-1)+ z] = 1
          b2[i] = freqZ[z]
        }
      }
      
      
      c2 = cbind(S,matrix(0,length(Z),length(Y)),matrix(0,length(Z),length(Y)),diag(length(Z))*(-1),matrix(0,length(Z),length(Z)))
      b2 = numeric(length(Z))
      for (z in Z){
        stc_sum2 = vector(length = nbX)
        for (i in 1:nbX){
          
          stc_sum2[i] = ifelse(length(indXB[[i]])==0,1/length(Z),length(indXB[[i]][inst$Zobserv[indXB[[i]]+nA] == z])/ length(indXB[[i]]))*length(indXA[[i]])/nA
          
        } 
        b2[z]= sum(stc_sum2)}
      s2 = rep("== ",length(Z))
      #@constraint(groupA, ctsumdeviationB, sum(deviationB[z] for z in Z) == 0)
      c3 = c(matrix(0,1,length(Y)*length(Z)),matrix(0,1,length(Y)),matrix(0,1,length(Y)),rep(1,length(Z)),matrix(0,1,length(Z)))
      b3 = rep(0,1)
      s3 = rep("==",1)
      # - bound the norm 1 of deviations
      #@constraint(groupA, ctabsdevBplus[z in Z], deviationB[z] <= absdevB[z])
      c4 = cbind(matrix(0,length(Z),length(Y)*length(Z)),matrix(0,length(Z),length(Y)),matrix(0,length(Z),length(Y)),diag(length(Z)),diag(length(Z))*(-1))
      b4 = rep(0,length(Z))
      s4 = rep("<=",length(Z))
      #@constraint(groupA, ctabsdevBmoins[z in Z], deviationB[z] >= -absdevB[z])
      c5 = cbind(matrix(0,length(Z),length(Y)*length(Z)),matrix(0,length(Z),length(Y)),matrix(0,length(Z),length(Y)),diag(length(Z)),diag(length(Z)))
      b5 = rep(0,length(Z))
      s5 = rep(">= ",length(Z))
      #@constraint(groupA, ctbounddevB, sum([absdevB[z] for z in Z]) <= maxrelax/2.0)
      c6 = c(matrix(0,1,length(Y)*length(Z)),matrix(0,1,length(Y)),matrix(0,1,length(Y)),matrix(0,1,length(Z)),rep(1,length(Z)))
      b6 = rep(maxrelax/2.0,1)
      s6 = rep("<=",1)
      c= rbind(c1,c2,c3,c4,c5,c6)
      b= c(b1,b2,b3,b4,b5,b6)
      s= c(s1,s2,s3,s4,s5,s6)
      groupA = solveLP(cvec=Min, bvec=b, Amat=c, const.dir = s,lpSolve=TRUE)
      transportA = matrix(groupA$solution[1:length(Y)*lenght(Z)], length(Y),lenght(Z))
      
      
      
      # groupB = Model(solver=GurobiSolver(Method=2,LogToConsole=0))
      
      
      # - transportB[y,z] : joint probability of modalities y and z if in base B
      #@variable(groupB, transportB[y in Y, z in Z] >= 0, base_name="transportB")
      
      # - deviationA[y]: deviation of the probability mass of y in base B from
      #   that observed in base A
      #@variable(groupB, deviationA[y in Y], base_name="deviationA")
      # - absdevA[y]: absolute value of the deviation of the probability mass of
      #   y in base B from that observed in base A
      #@variable(groupB, absdevA[y in Y] >= 0, base_name="absdevA")
      # - deviationB[z]: deviation of the probability mass of z in base A from
      
      
      
      #@objective(groupB, Min, sum(C[y,z]*transportB[y,z] for y in Y, z in Z))
      Min = c(as.numeric(C),rep(0,2*length(Y)+2*length(Z)))
      
      
      
      
      # - satisfy probabilities of modality A with possible deviation for individuals in base B
      #@constraint(groupB, cttransportAinB[y in Y],
      #            sum(transportB[y,z] for z in Z) == sum( (length(indXA[c])==0 ? 1.0 / length(Y) :  length(indXA[c][findall(inst$Yobserv[indXA[c]] == y)])/length(indXA[c]) ) * length(indXB[c])/nB for c in 1:nbX)+ deviationA[y]);
  
      S = matrix(0,length(Y),length(Y)*length(Z))
      for (y in Y){
        for (z in Z){
          S[y,length(Z)*(y-1)+ z] = 1
        }
      }
      c1 = cbind(S,matrix(0,length(Y),length(Y)),diag(length(Y))*(-1),matrix(0,length(Y),length(Z)),matrix(0,length(Y),length(Z)))
      b1 = numeric(length(Z))
      for (y in Y){
        stc_sum2 = vector(length = nbX)
        for (i in 1:nbX){
          
          stc_sum2[i] = ifelse(length(indXA[[i]])==0,1/length(Y),length(indXA[[i]][inst$Zobserv[indXA[[i]]] == z])/ length(indXA[[i]]))*length(indXB[[i]])/nB
          
        } 
        b1[y]= sum(stc_sum2)}
      s1 = rep("==",length(Y))
      
      # - satisfy probabilities of modality B for individuals in base B
      #@constraint(groupB, cttransportBinB[z in Z],
      #            sum(transportB[y,z] for y in Y) == freqZ[z])

      S = matrix(0,length(Z),length(Y)*length(Z))
       for (z in Z){
         for (y in Y){
            S[z,length(Z)*(y-1)+ z] = 1
        }
      }
      c2 = cbind(S,matrix(0,length(Z),length(Y)),matrix(0,length(Z),length(Y)),matrix(0,length(Z),length(Z)),matrix(0,length(Z),length(Z)))
      b2 = numeric(length(Z))
      for (z in Z){
        b2[z]= freqZ[z]}
      s2 = rep("==",length(Z))
      # - the deviations must sum to zero to conserve a well-defined probability measure
      #@constraint(groupB, ctsumdeviationA, sum(deviationA[y] for y in Y) == 0)
      c3 = c(matrix(0,1,length(Y)*length(Z)),rep(1,length(Y)),matrix(0,1,length(Y)),matrix(0,1,length(Z)),matrix(0,1,length(Z)))
      b3 = rep(0,1)
      s3 = rep("==",1)
      #@constraint(groupB, ctabsdevAplus[y in Y], deviationA[y] <= absdevA[y])
      c4 = cbind(matrix(0,length(Y),length(Y)*length(Z)),diag(length(Y)),diag(length(Y))*(-1),matrix(0,length(Y),length(Z)),matrix(0,length(Y),length(Z)))
      b4 = rep(0,length(Y))
      s4 = rep("<=",length(Y))
      #@constraint(groupB, ctabsdevAmoins[y in Y], deviationA[y] >= -absdevA[y])
      c5 = cbind(matrix(0,length(Y),length(Y)*length(Z)),diag(length(Y)),diag(length(Y)),matrix(0,length(Y),length(Z)),matrix(0,length(Y),length(Z)))
      b5 = rep(0,length(Y))
      s5 = rep(">=",length(Y))
      #@constraint(groupB, ctbounddevA, sum([absdevA[y] for y in Y]) <= maxrelax/2.0)
      c6 = c(matrix(0,1,length(Y)*length(Z)),matrix(0,1,length(Y)),rep(1,length(Y)),matrix(0,1,length(Z)),matrix(0,1,length(Z)))
      b6 = rep(maxrelax/2.0,1)
      s6 = rep("<=",1)
      # Solve the problem
      optimize!(groupA)
      optimize!(groupB)
      c= rbind(c1,c2,c3,c4,c5,c6)
      b= c(b1,b2,b3,b4,b5,b6)
      s= c(s1,s2,s3,s4,s5,s6)
      #optimize!(group);
      
      groupB = solveLP(cvec=Min, bvec=b, Amat=c, const.dir = s,lpSolve=TRUE)
      transportB = matrix(groupB$solution[1:length(Y)*lenght(Z)], length(Y),lenght(Z))
      
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
            estimatorZA[x,inst$Yobserv[i],YBpred[i]] = estimatorZA[x,inst$Yobserv[i],YBpred[i]] +  1/ length(findall(inst$Yobserv[indXA[x]] == inst$Yobserv[i]));
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