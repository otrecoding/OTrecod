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
    
    D  = rdist::cdist(a,b,"manhattan")
    DA = rdist::cdist(a,a,"manhattan")
    DB = rdist::cdist(b,b,"manhattan")
    
  } else if (norme == 2){
    
    D  = rdist::cdist(a,b,"euclidean")
    DA = rdist::cdist(a,a,"euclidean")
    DB = rdist::cdist(b,b,"euclidean")
    
  } else if (norme == 0){
    
    D  = rdist::cdist(a,b,"hamming")
    DA = rdist::cdist(a,a,"hamming")
    DB = rdist::cdist(b,b,"hamming")
    
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
    # print(i)
    
    nbX = nbX + 1;
    # x = matrix(0,dim(Xval[2]),1);
    # plut?t: x = rep(0,ncol(Xval)) mais inutile
    
    # x[,1] = Xval[i,(1:dim(Xval)[2])];
    #x[:,1] = [Xval[i,j] for j in 1:size(Xval,2)];
    
    x = Xval[i,]
    
    
    if (norme == 1){
      
      # distA = pairwise(Cityblock(), x, t(Xobserv[A,]), dims=2)
      # distB = pairwise(Cityblock(), x, t(Xobserv[B + nA,]), dims=2)
      
      distA  = rdist::cdist(x,Xobserv[A,],"manhattan")
      distB  = rdist::cdist(x,Xobserv[B + nA,],"manhattan")
      
      
      
    } else if (norme == 2){
      
      distA  = rdist::cdist(x,Xobserv[A,],"euclidean")
      distB  = rdist::cdist(x,Xobserv[B + nA,],"euclidean")
      
    } else if (norme == 0){
      
      distA  = rdist::cdist(x,Xobserv[A,],"hamming")
      distB  = rdist::cdist(x,Xobserv[B + nA,],"hamming")
      
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
