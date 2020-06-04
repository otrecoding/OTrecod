#' proxim_dist()
#'
#' \code{proxim_dist} computes pairwise distance matrice of a database and cross-distance matrix between two databases according to various distance functions in the context of data fusion.
#'
#'
#' This function is the first step of a family of algorithms that solve recoding problems of data fusion using Optimal Transportation theory (\code{OUTCOME},\code{R_OUTCOME},\code{JOINT} and \code{R_JOINT}, in the two following references of Gares and Al.)
#' \code{proxim_dist} is directly implemented in the \code{OT_outcome} and \code{OT_JOINT} functions but can also be used separately as long as the input database has as suitable structure.
#' The function also provides in output the list of all the profiles of covariates encountered.
#'
#' A. REQUIRED STRUCTURE OF THE DATABASE
#'
#' The input database is a data.frame that must be saved in a specific form by users: Two overlayed databases containing a common column of databases' identifiers (A and B, 1 or 2, by examples, encoded in numeric or factor form),
#' a column corresponding to the target variable with its specific encoding in A (By example a factor Y encoded in nY levels, ordered or not, with NAs in the corresponding rows of B), a column corresponding to the same variable with its specific endoded in B (By example a factor Z in nZ levels,
#' with NAs in rows of A), and a set of shared covariates (at least one).
#'
#' The order of these variables in the database have no importance but the indexes of the columns related to the 3rd columns previously described (ie ID, Y and Z) must be specified in the \code{indx_DB_Y_Z} option.
#' Please do not hesitate to refer to the related examples for a better understanding of the structure expected for your database in argument.
#'
#' Missing values are allowed on covariates only (and in presence of more than one covariate), and are excluded from all computations involving the rows within which they occur.
#' If some columns are excluded when computing a Euclidean, Manhattan, or Hamming distance, the sum is scaled up proportionally to the number of columns used as proposed by the standard (\code{\link[stats]{dist}}) function.
#' If all pairs are excluded when computing a particular distance, instead of putting NA in the corresponding cell of the distance matrix, the process stops and an object listing the problematic rows is proposed in output.
#' It suggests users to remove these rows before running the process again or impute NAs related to these rows (See ex 6 for more details).
#'
#' In the particular case where only one covariate with NAs is used, we recommend working with imputed or complete case only to avoid the presence of NA in the distance matrix.
#'
#' If the database counts many covariates and some of them have missing data, user can keep them or apply beforehand the \code{\link{imput_cov}} function on data.frame to deal with this problem.
#'
#'
#' B. DISTANCE FUNCTIONS AND TYPES OF COVARIATES
#'
#' The choice of the distance function must be made according to the type of the covariates.
#' It requires the covariates to have been previously stored in some specific forms in the database, and the use of the \code{\link{transfo_dist}} function is there to help the user in this task.
#' The Euclidean or Manhattan distance (\code{norm} = "E" or "M") only requires numeric covariates with finite or infinite number of values (ie continuous variables with or without decimals values). Binary ciovariates are allowed.
#' The Hamming distance (\code{norm} = "H") only requires binary variables. In this case, continuous variables could have been converted in factor of k levels (k>2) beforehand, and then transformed in disjunctive table (containing the (k-1) corresponding binary variables) before use.
#' The \code{\link{transfo_quali}} function can help users in this task.
#' Using the Hamming distance could be quite long, in presence NAs on covariates, so we suggerate users to be a little patient in this case.
#' The Gower distance (\code{norm} = G") uses the (\code{\link[StatMatch]{gower.dist}}) function and so allows logical, factor, character and numeric variables as described in its documentation.
#'
#'
#' C. PROFILES OF COVARIATES AND OUTPUT DETAILS
#'
#' By example, assuming that a data.frame (composed of 2 overlayed data.frame A and B) have 3 shared binary covariates (ie identically encoded in A and B) so \code{011} and \code{101} will be considered as two distinct profiles of covariates.
#' Now, by supposing that this data.frame has 3 shared covariates,and that each covariate is a factor of n1, n2 and n3 levels respectively, so it exists at most n1*n2*n3 possible profiles.
#' This estimation is considered as a maximum here because, in fact, only the profiles of covariates present in at least one of the two databases will be kept for the other step of an OT algorithm.
#'
#' Whatever the type (mixed or not) and the number of covariates in the data.frame of interest, \code{proxim_dist} firstly detect all the possible profiles (or combinations) of covariates from the 2 databases.
#' \code{proxim_dist} classifies individuals from the 2 databases according to their proximities to each profile of covariates and saves the corresponding indexes of rows from A and B in 2 lists \code{indXA} and \code{indXB} respectively.
#' \code{indXA} and \code{indXB} thus contained as many objects as covariates profiles and the proximity between a given profile and a given individual is defined as follows.
#'
#' As a decision rule, for a given profile of covariates \code{P_j}, an individual i will be considered as a neighbor of \code{P_j} if \eqn{dist(i,P_j) < \code{prox} * max(dist(i,P_j))} where \code{prox} will be fixed by user.
#'
#'
#' @param data_file A data.frame corresponding ideally to the output of \code{\link{transfo_dist}}. Ohterwise this data.frame must have an ID variable for the identification of the 2 superimposed databases (1 and 2 or "A" and "B" by example), a target variable (called Y by example) only encoded in the 1st database, a target variable (Z) only stored in the 2nd database, such that Y and Z summarize a same information differently encoded in the 2 databases and set of common covariates (at least one) of any type.
#' The order of the variables in the data.frame have no importance.The type of the covariates must be in accordance with the distance measurement chosen in the \code{norm} option.
#' @param indx_DB_Y_Z A vector of 3 column numbers corresponding to the colum indexes of the ID variable, the target variable in the 1st database and the target variable in the 2nd database. The indexes must be declared in this specific order.
#' @param norm A character (with quotes) indicating the choice of the distance function. This latest depends on the type of the common covariates.The Hamming distance
#' for binary covariates only (\code{norm} = "H"), the Manhattan distance ("M", by default) and the euclidean distance ("E") for continuous covariates only, or the Gower distance for mixed covariates ("G").
#' @param prox A percentage (betwen 0 and 1) uses to calculate the distance threshold below which an individual (a row) is considered as a neighbor of a given profile of covariates.
#'
#' @return A list of 15 elements (the first 15 detailed below) is returned containing various distance matrices and lists useful for the OT algorithm. A list of 2 objects (The last 2 of the following list) is returned if distance matrices contained NAs.
#' \item{FILE_NAME}{A simple reminder of the name of the raw database}
#' \item{nA}{The row numbers of the 1st database}
#' \item{nB}{The row numbers of the 2nd database}
#' \item{Xobserv}{The subset of the two merged databases composed of the common variables only}
#' \item{Yobserv}{The values of the target variable in the 1st database}
#' \item{Zobserv}{The values of the target variable in the 2nd database}
#' \item{D}{A distance matrix corresponding to the distances computed between individuals of the 2 databases}
#' \item{Y}{The nY levels of the target variable in numeric form, in the 1st database}
#' \item{Z}{The nZ levels of the target variable in numeric form, in the 2nd database}
#' \item{indY}{A list of nY groups of individual (or row) numbers where each group corresponds to the individuals indexes related to a given level of Y in the 1st database}
#' \item{indZ}{A list of nZ groups of individual (or row) numbers where each group corresponds to the individuals indexes related to a given level of Z in the 2nd database}
#' \item{indXA}{A list of individual (row) indexes from the 1st database, sorted by profiles of covariates according to their proximities. See the \code{Details} part for more information}
#' \item{indXB}{A list of individual (row) indexes from the 2nd database, sorted by profiles of covariates according to their proximities. See the \code{Details} part for more information}
#' \item{DA}{A distance matrix corresponding to the pairwise distances between individuals of the 1st database}
#' \item{DB}{A distance matrix corresponding to the pairwise distances between individuals of the 2nd database}
#' \item{ROWS_TABLE}{Combinations of row numbers of the 2 databases that generate NAs in D}
#' \item{ROWS_TO_RM}{Number of times a row of the 1st or 2nd database is involved in the NA process of D}
#'
#'
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#' \email{gregory.guernec@@inserm.fr}
#'
#' @importFrom  rdist cdist
#' @importFrom  proxy dist
#' @importFrom  stats na.omit
#' @importFrom  StatMatch gower.dist
#'
#' @seealso \code{\link{transfo_dist}}, \code{\link{imput_cov}}
#'
#' @references
#' ### About OT algorithms for data integration:
#' Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679 | \url{https://doi.org/10.1515/ijb-2018-0106}
#'
#' Gares V, Omer J. Regularized optimal transport of covariates and outcomes in datarecoding(2019).hal-02123109 \url{https://hal.archives-ouvertes.fr/hal-02123109/document}
#'
#'
#' ### About the Gower distance:
#' Gower, J. C. (1971), “A general coefficient of similarity and some of its properties”. Biometrics, 27, 623--637.
#'
#' ### About distance measures:
#' Anderberg, M.R. (1973), Cluster analysis for applications, 359 pp., Academic Press, New York, NY, USA.
#'
#' Borg, I. and Groenen, P. (1997) Modern Multidimensional Scaling. Theory and Applications. Springer.
#'
#'
#' @aliases proxim_dist
#'
#'
#' @export
#'
#' @examples
#'
#' data(simu_data)
#' ### The covariates of the data are prepared according to the distance chosen
#' ### using the transfo_dist function
#'
#' ### Ex 1: The Manhattan distance
#'
#' try1 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),
#'                     ordinal = c(2,6), logic = NULL, prep_choice = "M")
#' res1 = proxim_dist(try1,norm = "M")  # try1 compatible with norm = "E" for Euclidean
#'
#'
#' ### Ex 2: The Euclidean and Manhattan distance applied on coordinates from FAMD
#'
#' try2 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),
#'                     ordinal = c(2,6), logic = NULL, prep_choice = "FAMD",info = 0.80)
#' res2_E  = proxim_dist(try2,norm = "E")
#' \dontrun{
#' res2_M  = proxim_dist(try2,norm = "M")
#' }
#'
#' ### Ex 3: The Gower distance with mixed covariates
#'
#' try3 = transfo_dist(simu_data[c(1:100,301:400),],quanti = c(3,8), nominal = c(1,4:5,7),
#'                     ordinal = c(2,6), logic = NULL, prep_choice = "G")
#' res3 = proxim_dist(try3,norm = "G")
#'
#' \dontrun{
#' ### Ex 4a: The Hamming distance with binary (but incomplete) covariates only
#'
#' # categorization of the continuous covariates age by tertiles
#' try4 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),ordinal = c(2,6),
#'                     convert_num = 8, convert_clss = 3, prep_choice = "H")
#' res4 = proxim_dist(try4,norm = "H")
#' # Be patient ... It could take few minutes
#'
#' ### Ex 4b: The Hamming distance with complete cases on nominal and ordinal covariates only
#' simu_data_CC = simu_data[(!is.na(simu_data[,5]))&(!is.na(simu_data[,6]))
#'                           &(!is.na(simu_data[,7])),1:7]
#' try4b = transfo_dist(simu_data_CC,quanti = 3, nominal = c(1,4:5,7),ordinal = c(2,6),
#'                      prep_choice = "H")
#' res4b = proxim_dist(try4b,norm = "H")
#' }
#'
#' ### Ex 5: PARTICULAR CASE, If only one covariate with no NAs
#'
#' try5 = try1[,c(1:3,7)]           # Only Smoking variable
#' try6 = try5[!is.na(try5[,4]),]   # Keep complete case
#' res6_M = proxim_dist(try6,norm = "M",prox = 0.01)
#' \dontrun{
#' res7_H = proxim_dist(try6,norm = "H") # Hamming
#' }
#' \dontrun{
#' ### Ex 6: PARTICULAR CASE, many covariates but NAs in distance matrix
#'
#' # We generated NAs in the try1 object so that:
#' # dist(A4,B102) and dist(A122,B102) returns NA whatever the norm chosen:
#' try7 = try1
#' try7[4,7:9]   = NA
#' try7[122,6:9] = NA
#' try7[300+102,4:6] = NA
#' res7 = proxim_dist(try7,norm = "M")
#' # The process stopped indicates 2 NAs and the corresponding row numbers
#' # The 2nd output of res7 indicates that removing first the 102th row of the database
#' # B is enough to solve the pb:
#' try8 = try7[-402,]
#' res8 = proxim_dist(try8,norm = "M")
#' }
#'
proxim_dist  = function(data_file,indx_DB_Y_Z = 1:3,norm = "E", prox = 0.80){


  if (!is.data.frame(data_file)){

    stop("This object must be a data.frame !")

  } else {
  }

  if (!(norm %in% c("H","M","E","G"))){

    stop("Improper character in the norm argument")

  } else {}

  if (length(indx_DB_Y_Z)!=3){

    stop("Invalid length for indx_DB_Y_Z: This option must contain the index of the DB identification variable, and the indexes of columns related to the targets in A and B")

  } else {}

  if (max(indx_DB_Y_Z)>ncol(data_file)){

    stop("Invalid index of column")

  } else {}


  dat     = data_file[,c(indx_DB_Y_Z,setdiff(1:ncol(data_file),indx_DB_Y_Z))]
  dat[,1] = as.factor(dat[,1])
  dat[dat[,1] == levels(dat[,1])[1],3] = NA
  dat[dat[,1] == levels(dat[,1])[2],2] = NA

  if (unique(dat[,1])[1] != levels(dat[,1])[1]){

    stop("Please change the name of your databases in the ID column so that the names of the 2 databases will be alphanumerically ranked in ascending order. You can use A for data1 and B for data2 or simply 1 and 2 by example.")

  } else {}


  dat[,1] = as.numeric((dat[,1]))

  if (sum(is.na(dat[,1]))!=0){

    stop("The ID variable has at least one NA, please correct it and return")

  } else {}

  if (nrow(dat[dat[,1] == 1,1:2])!= nrow(na.omit(dat[dat[,1] == 1,1:2]))){

    stop("Y has at least one NA in the 1st database, please correct it and return")

  } else if (nrow(dat[dat[,1] == 2,c(1,3)])!= nrow(na.omit(dat[dat[,1] == 2,c(1,3)]))){

    stop("Z has at least one NA in the 2nd database, please correct it and return")

  } else {}


  # number of covariates
  nbcvar  = ncol(dat) - 3

  if (nbcvar < 0){

    stop("Invalid number of variables in your database: At least 4")

  } else {}

  if (nbcvar == 0){

    stop("No covariate in your database. At least 1 covariate is required")

  } else {}

  if ((nbcvar == 1)&(sum(is.na(dat[,4]))!=0)){

    stop("Only one covariate with NAs: Please impute NAs or work with complete case only")

  } else {}



  # recover the sets of individuals in base 1 and 2
  Base    = dat[1:nrow(dat),1]
  indA    = which(Base == 1)
  indB    = which(Base == 2)
  nA      = length(indA)
  nB      = length(indB)

  # recover the input data
  Xobserv = dat[1:nrow(dat),4:ncol(dat)]
  Yobserv = dat[1:nrow(dat),2]
  Zobserv = dat[1:nrow(dat),3]

  # modify order so that base A comes first and then base B

  if (nbcvar == 1){

    Xobserv = c(Xobserv[indA],Xobserv[indB])

  } else {

    Xobserv = rbind(Xobserv[indA,],Xobserv[indB,])

  }

  Yobserv = c(Yobserv[indA],Yobserv[indB])
  Zobserv = c(Zobserv[indA],Zobserv[indB])
  indA = 1:nA
  indB = (nA+1):(nA+nB)

  # Modify Y and Z so that they go from 1 to the number of modalities
  Y = sort(unique(Yobserv[Yobserv != -1]));
  Z = sort(unique(Zobserv[Zobserv != -1]));
  for (i in 1:length(Y)){

    Yobserv[Yobserv == Y[i]]= i

  }

  #Y = [i for i in 1:length(Y)];
  Y = 1:length(Y)
  for (i in 1:length(Z)){

    Zobserv[Zobserv == Z[i]] = i

  }

  #Z = [i for i in 1:length(Z)];
  Z = 1:length(Z)


  # list the distinct modalities in A and B
  indY = indZ = list()
  for (m in Y){indY[[m]] = which(Yobserv[1:nA] == m)}
  for (m in Z){indZ[[m]] = which(Zobserv[(nA+1):(nA+nB)] == m)}


  # Compute the distance between pairs of individuals in different bases

  if (nbcvar == 1){

    a = Xobserv[indA]
    b = Xobserv[indB]

  } else {

    a = Xobserv[indA,]
    b = Xobserv[indB,]

  }

  if (norm == "M"){

    D  = proxy::dist(a,b, method = "manhattan")
    DA = proxy::dist(a,a, method = "manhattan")
    DB = proxy::dist(b,b, method = "manhattan")

  } else if (norm == "E"){

    D  = proxy::dist(a,b, method = "euclidean")
    DA = proxy::dist(a,a, method = "euclidean")
    DB = proxy::dist(b,b, method = "euclidean")

  } else if (norm == "H"){

    if (nrow(as.data.frame(dat[,4:ncol(dat)]))== nrow(na.omit(as.data.frame(dat[,4:ncol(dat)])))){

      D  = rdist::cdist(a,b, metric = "hamming")
      DA = rdist::cdist(a,a, metric = "hamming")
      DB = rdist::cdist(b,b, metric = "hamming")

    } else {

      D  = ham(a,b)
      print("Hamming 1/3")
      DA = ham(a,a)
      print("Hamming 2/3")
      DB = ham(b,b)
      print("Hamming 3/3")

    }

  } else if (norm == "G"){

    D  = StatMatch::gower.dist(a,b)
    DA = StatMatch::gower.dist(a,a)
    DB = StatMatch::gower.dist(b,b)

  }

  # D = as.data.frame(D)

  if (nrow(D) != nrow(na.omit(D))){

    rowD    = (1:nrow(D))[apply(is.na(D),1,sum)!=0]
    colD    = (1:ncol(D))[apply(is.na(D),2,sum)!=0]
    coord   = cbind(rep(rowD,each= length(colD)),rep(colD,length(rowD)))
    tab_pb  = coord[is.na(apply(coord,1,function(x){D[x[1],x[2]]})),]

    if (length(nrow(tab_pb)) == 0){

      names(tab_pb) = c("A","B")
      freq_pb       = sort(table(c(paste("A",tab_pb[1],sep=""),paste("B",tab_pb[2],sep=""))),decreasing = TRUE)

    } else {

      colnames(tab_pb) = c("A","B")
      freq_pb = sort(table(c(paste("A",tab_pb[,1],sep=""),paste("B",tab_pb[,2],sep=""))),decreasing = TRUE)

    }

    warning("THE PROCESS STOPPED")

    cat("!!! Because of the presence of NAs in distance matrix, the process stopped",
        "Combinations of rows of A and B with NAs cause pbs and have to be removed or imputed",
        "To help you, the indexes of rows are listed in the returned object",sep="\n")

    return(list(ROWS_TABLE = tab_pb, ROWS_TO_RM = freq_pb))

  } else {}


  # Compute the indexes of individuals with same covariates
  A     = 1:nA;
  B     = 1:nB;
  nbX   = 0;


  Xval  = unique(Xobserv)
  # indXA = indXB = rep(0,nrow(Xval))


  # X1val = sort(unique(Xobserv[,1]));
  # X2val = sort(unique(Xobserv[,2]));
  # X3val = sort(unique(Xobserv[,3]));



  # aggregate both bases

  indXA = indXB =  list()

  n_Xval = ifelse(nbcvar != 1,nrow(Xval),length(Xval))

  for (i in  (1:n_Xval)){

    nbX = nbX + 1;

    if (nbcvar != 1){

      x      = Xval[i,]
      Xobs_A = Xobserv[A,]
      Xobs_B = Xobserv[B + nA,]

    } else {

      x      = Xval[i]
      Xobs_A = Xobserv[A]
      Xobs_B = Xobserv[B + nA]

    }

    if (norm == "M"){

      distA  = proxy::dist(x, Xobs_A, method =  "manhattan")
      distB  = proxy::dist(x, Xobs_B, method =  "manhattan")

    } else if (norm == "E"){

      distA  = proxy::dist(x, Xobs_A, method = "euclidean")
      distB  = proxy::dist(x, Xobs_B, method = "euclidean")

    } else if (norm == "H"){

      if (nrow(as.data.frame(dat[,4:ncol(dat)]))== nrow(na.omit(as.data.frame(dat[,4:ncol(dat)])))){

        distA  = rdist::cdist(x, Xobs_A, metric = "hamming")
        distB  = rdist::cdist(x, Xobs_B, metric = "hamming")

      } else {

        distA  = ham(x, Xobs_A)
        distB  = ham(x, Xobs_B)

      }

    } else if (norm == "G"){

      distA  = StatMatch::gower.dist(x,Xobs_A)
      distB  = StatMatch::gower.dist(x,Xobs_B)

    } else {}

    indXA[[nbX]] = which(distA < prox*max(distA,na.rm=TRUE)); names(indXA[[nbX]]) = NULL
    indXB[[nbX]] = which(distB < prox*max(distA,na.rm=TRUE)); names(indXB[[nbX]]) = NULL

  }

  # file_name = base_name(data_file)
  file_name = deparse(substitute(data_file))

  return(list(FILE_NAME = file_name, nA = nA, nB = nB, Xobserv = Xobserv,
              Yobserv   = Yobserv  , Zobserv = Zobserv, D = D, Y = Y, Z = Z,
              indY=indY, indZ=indZ, indXA=indXA, indXB=indXB, DA=DA, DB=DB))
}


