<<<<<<< HEAD
#' Instance()
#' 
#' This function computes the distances between pairs of individuals in 2 distinct databases using the information contained in the common variables
#' 
#' This function allows covariates in ordered or unordered factors only
#' 
#' @param data_file A data.frame with a specific order of columns. The 1st comumn is a key column for individual identification,
#' the second and third columns corresponds to the target variable encoded with different number of levels in databases 1 and 2 (The levels must be beforehand convert in numeric, started from 1),
#' the following columns corresponds to common covariates between the 2 databases whatever their number, and whatever their order. Nevertheless, covariates must be factors beforehand convert in  numeric
#' @param norme An integer that can only take the values 0,1, or 2, corresponding to the choice of the distance function. This latest depends on the type of the common covariates.The Hamming distance
#' for binary covariates (norme = 0) useful for disjunctive forms by example, or the Manhattan distance (1, by default) or the euclidean distance (2)
#' 
#' @return A list of 15 elements is returned:
#' \item{FILE_NAME}{A simple reminder of the name of the raw database}
#' \item{nA}{The length of the 1st database}
#' \item{nB}{The length of the 2nd database}
#' \item{Xobserv}{The subset of the two merged databases composed of the common variables}
#' \item{Yobserv}{The values of the target variable in the 1st database}
#' \item{Zobserv}{The values of the target variable in the 2nd database}
#' \item{D}{A distance matrix corresponding to the distances computed between individuals of the 2 databases}
#' \item{Y}{The levels of the target variable in numeric form, in the 1st database}
#' \item{Z}{The levels of the target variable in numeric form, in the 2nd database}
#' \item{indY}{A list of groups of individual index related to each level of the target variable in the 1st database}
#' \item{indZ}{A list of groups of individual index related to each level of the target variable in the 2nd database}
#' \item{indXA}{A list of groups of individual index related to all possible profiles of covariates in the 1st database}
#' \item{indXB}{A list of groups of individual index related to all possible profiles of covariates in the 2nd database}
#' \item{DA}{A distance matrix corresponding to the pairwise distances between individuals of the 1st database}
#' \item{DB}{A distance matrix corresponding to the pairwise distances between individuals of the 2nd database}
#'  
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer 
#' \email{gregory.guernec@@inserm.fr}
#'
#' @importFrom  rdist cdist
#' @export
#'
#' @examples 
#'### Ex 1
#'data(tab)
#'stock_res = Instance(tab)
#'summary(stock_res)
#'
#'
#'### Ex 2
#'# Please, see the example of the merge_dbs function to obtain data3 and data4
#'soluc1  = merge_dbs(data3,data4,NAME_Y = "c.neti",NAME_Z = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9),impute = "MICE",R_MICE = 2, seed_func = 4036)
#'summary(soluc1$DB_READY)
#'
#'soluc1$DB_READY$hsize5 = as.numeric(soluc1$DB_READY$hsize5)
#'soluc1$DB_READY$sex = as.numeric(soluc1$DB_READY$sex)
#'
#'res1 = Instance(soluc1$DB_READY,norme = 1)
#'
#'data_test       = soluc1$DB_READY
#'data_test$poids = rnorm(nrow(data_test),mean = 65,sd = 10)
#'res2            = Instance(data_test,norme = 1)
#'
#'# Preliminary step: Discretize the continuous covariates if necessary, like "age" in this example
#'soluc1$DB_READY$age = cut(soluc1$DB_READY$age,c(0,30,60,100)); summary(soluc1$DB_READY)
#' 
#'# Step2: Use the prep_dbs() function to prepare the database for calculating the distances
#'try1 = prep_dbs(soluc1$DB_READY,nominal = 6,ordinal = 1:5) 
#'
#'res1 = Instance(try1,norme = 0)   # Using the Hamming distance
#'
#'try2 = try1[,-6]
#'res2 = Instance(try2,norme = 0) 
#'
Instance  = function(data_file,norme = 1){
  
  
  if (!is.data.frame(data_file)){
    
    stop("This object must be a data.frame !")
    
  } else {
  }
  
  if (!(norme %in% c(0,1,2))){
    
    stop("Incorrect value for the choice of the distance: Only 0,1 or 2 is allowed with no quotes")  
    
  } else {}
  
  
  dat    = data_file
  
=======
#' Instance(data_file,norme)
#'
#' @param data_file todo list
#' @param norme todo list
#'
#' @return todo list
#' @export
#'
# @examples
Instance  = function(data_file,norme){
  
  dat    = data_file
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
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
<<<<<<< HEAD
=======
  # Xobserv = cbind(Xobserv[indA,],Xobserv[indB,])
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
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
  
<<<<<<< HEAD
  # Compute the distance between pairs of individuals in different bases
  a = Xobserv[indA,]
  b = Xobserv[indB,]
  
  # stopifnot(norme %in% c(0,1,2))
=======
  # compute the distance between pairs of individuals in different bases
  # devectorize all the computations to go about twice faster
  # only compute norm 1 here
  a = Xobserv[indA,]
  b = Xobserv[indB,]
  
  stopifnot(norme %in% c(0,1,2))
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
  
  
  if (norme == 1){
    
<<<<<<< HEAD
    D  = dist(a,b, method = "manhattan")
    DA = dist(a,a, method = "manhattan")
    DB = dist(b,b, method = "manhattan")
    
  } else if (norme == 2){
    
    D  = dist(a,b, method = "euclidean")
    DA = dist(a,a, method = "euclidean")
    DB = dist(b,b, method = "euclidean")
    
  } else if (norme == 0){
    
    if ((sum(is.na(a))!=0)|(sum(is.na(b))!=0)){
      
      stop("Covariates in at leat one of the two DB have missing values: Please impute them or delete corresponding individuals before using this distance")
      
    }
    
=======
    #   D = pairwise(Cityblock(), a, b, dims=2)
    #  DA = pairwise(Cityblock(), a, a, dims=2)
    #  DB = pairwise(Cityblock(), b, b, dims=2)
    
    D  = rdist::cdist(a,b,"manhattan")
    DA = rdist::cdist(a,a,"manhattan")
    DB = rdist::cdist(b,b,"manhattan")
    
  } else if (norme == 2){
    
    D  = rdist::cdist(a,b,"euclidean")
    DA = rdist::cdist(a,a,"euclidean")
    DB = rdist::cdist(b,b,"euclidean")
    
  } else if (norme == 0){
    
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
    D  = rdist::cdist(a,b,"hamming")
    DA = rdist::cdist(a,a,"hamming")
    DB = rdist::cdist(b,b,"hamming")
    
  }
  
  # Compute the indexes of individuals with same covariates
  A     = 1:nA;
  B     = 1:nB;
  nbX   = 0;
<<<<<<< HEAD
  
=======
  # indXA = numeric(dim(Xval)[1]);
  # indXB = numeric(dim(Xval)[1]));
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
  
  Xval  = unique(Xobserv)
  indXA = indXB = rep(0,nrow(Xval))
  
  
<<<<<<< HEAD
  # X1val = sort(unique(Xobserv[,1]));
  # X2val = sort(unique(Xobserv[,2]));
  # X3val = sort(unique(Xobserv[,3]));
=======
  X1val = sort(unique(Xobserv[,1]));
  X2val = sort(unique(Xobserv[,2]));
  X3val = sort(unique(Xobserv[,3]));
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
  
  
  
  # aggregate both bases
  
  indXA = indXB =  list()
  
  for (i in  (1:nrow(Xval))){
    
<<<<<<< HEAD
    nbX = nbX + 1;
    x   = Xval[i,]
=======
    # if (i %in% seq(0,nrow(Xval),50)){print(i)} else {}
    # print(i)
    
    nbX = nbX + 1;
    # x = matrix(0,dim(Xval[2]),1);
    # plut?t: x = rep(0,ncol(Xval)) mais inutile
    
    # x[,1] = Xval[i,(1:dim(Xval)[2])];
    #x[:,1] = [Xval[i,j] for j in 1:size(Xval,2)];
    
    x = Xval[i,]
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
    
    
    if (norme == 1){
      
<<<<<<< HEAD
      distA  = dist(x,Xobserv[A,]     , method =  "manhattan")
      distB  = dist(x,Xobserv[B + nA,], method =  "manhattan")
      
    } else if (norme == 2){
      
      distA  = dist(x,Xobserv[A,]     , method = "euclidean")
      distB  = dist(x,Xobserv[B + nA,], method = "euclidean")
=======
      # distA = pairwise(Cityblock(), x, t(Xobserv[A,]), dims=2)
      # distB = pairwise(Cityblock(), x, t(Xobserv[B + nA,]), dims=2)
      
      distA  = rdist::cdist(x,Xobserv[A,],"manhattan")
      distB  = rdist::cdist(x,Xobserv[B + nA,],"manhattan")
      
      
      
    } else if (norme == 2){
      
      distA  = rdist::cdist(x,Xobserv[A,],"euclidean")
      distB  = rdist::cdist(x,Xobserv[B + nA,],"euclidean")
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
      
    } else if (norme == 0){
      
      distA  = rdist::cdist(x,Xobserv[A,],"hamming")
      distB  = rdist::cdist(x,Xobserv[B + nA,],"hamming")
      
    }
    
<<<<<<< HEAD
=======
    # indXA[nbX] = (distA[1,] < 0.1)
    # indXB[nbX] = (distB[1,] < 0.1)
    
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
    indXA[[nbX]] = which(distA < 0.1)
    indXB[[nbX]] = which(distB < 0.1)
    
  }
  
  # file_name = base_name(data_file)
  file_name = deparse(substitute(data_file))
  
  return(list(FILE_NAME =file_name,nA=nA, nB=nB, Xobserv=Xobserv,
              Yobserv=Yobserv, Zobserv=Zobserv, D=D, Y=Y, Z=Z,
              indY=indY, indZ=indZ, indXA=indXA, indXB=indXB, DA=DA, DB=DB))
}
<<<<<<< HEAD


=======
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
