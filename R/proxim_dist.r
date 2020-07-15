#' proxim_dist()
#'
#' \code{proxim_dist} computes the pairwise distance matrix of a database and cross-distance matrix between two databases according to various distance functions in the context of data fusion.
#'
#'
#' This function is the first step of a family of algorithms that solve recoding problems of data fusion using Optimal Transportation theory (See the details of these corresponding models \code{OUTCOME}, \code{R_OUTCOME}, \code{JOINT} and \code{R_JOINT} in (1) and (2))
#' The function \code{proxim_dist} is directly implemented in the functions \code{\link{OT_outcome}} and \code{\link{OT_joint}} but can also be used separately as long as the input database has as suitable structure. Nevertheless, its preparation will have to be rigorously made in two steps detailled in the following sections.
#'
#' A. REQUIRED STRUCTURE FOR THE DATABASE
#'
#' Firsly, the initial database required is a data.frame that must be prepared in a specific form by users. From two separate databases, the function \code{\link{merge_dbs}} available in this package can assist users in this initial merging, nevertheless notice that this preliminary transformation can also be made directly by following the imposed structure described below:
#' two overlayed databases containing a common column of database identifiers (A and B for examples, encoded in numeric or factor form),
#' a column corresponding to the target variable with its specific encoding in A (for example a factor \eqn{Y} encoded in \eqn{n_Y} levels, ordered or not, with NAs in the corresponding rows of B), a column corresponding to the same variable with its specific endoded in B (for example a factor \eqn{Z} in \eqn{n_Z} levels,
#' with NAs in database A), and a set of shared covariates (at least one) between the two databases.
#'
#' The order of these variables in the database have no importance but the column indexes related to database identifier, \eqn{Y} and \eqn{Z}, must be specified in the \code{indx_DB_Y_Z} option.
#' Users can refer to the structure of the table \code{\link{simu_data}} available in the package to adapt their databases to the inital format required.
#'
#' Missing values are allowed on covariates only (and in presence of more than one covariate), and are excluded from all computations involving the rows within which they occur.
#' In the particular case where only one covariate with NAs is used, we recommend working with imputed or complete case only to avoid the presence of NA in the distance matrix that will be computed a posteriori.
#'
#' If the database counts many covariates and some of them have missing data, user can keep them or apply beforehand the \code{\link{imput_cov}} function on data.frame to deal with this problem.
#'
#'
#' B. DISTANCE FUNCTIONS AND TYPES OF COVARIATES
#'
#' In a second step, before including the data.frame in the function, the shared variables of the merged database will have to be encoded according to the choice of the distance function fixed by user, knowing that it is also frequent that it is the type of the variables which fixes the distance function to choose.
#' The function \code{\link{transfo_dist}} is available in the package to assist users in this task but a user can also decide to make this preparation by themselves.
#' Thus, with the Euclidean or Manhattan distance ((3), \code{norm} = "E" or "M"), if all types of variables will be allowed, logical variables will have to be previously transformed in binary variables, and categorical variables (factors ordered or not) will have to be transformed by their related disjunctive tables (the function \code{\link{transfo_quali}} can make these specific transformations).
#' The Hamming distance (\code{norm} = "H") only requires binary variables (all other forms are not allowed). In this context, continuous variables could have been converted in factor of k levels (\eqn{k>2}) beforehand. The categorical covariates are then transformed in disjunctive tables (containing the (\eqn{k-1}) corresponding binary variables) before use. With this distance, categorical variables are also transformed in disjunctive tables.
#' Notice that, using the Hamming distance could be quite long in presence NAs on covariates, so please be a little patient and let the function runs few minutes if necessary.
#' Finally, the Gower distance ((4), \code{norm} = "G") uses the (\code{\link[StatMatch]{gower.dist}}) function (5) and so allows logical, categorical and numeric variables with no preliminary transformations.
#'
#' In conclusion, the structure of the data.frame required in input of the function \code{proxim_dist} corresponds to two overlayed databases with two target outcomes and a set of shared covariates whose encodings depend on the distance function choosen by user.
#'
#' If some columns are excluded when computing an Euclidean, Manhattan, or Hamming distance between two rows, the sum is scaled up proportionally to the number of columns used in the computation as proposed by the standard (\code{\link[stats]{dist}}) function.
#' If all pairs are excluded when computing a particular distance, instead of putting NA in the corresponding cell of the distance matrix, the process stops and an object listing the problematic rows is proposed in output.
#' It suggests users to remove these rows before running the process again or impute NAs related to these rows (see (6) for more details).
#'
#' C. PROFILES OF COVARIATES AND OUTPUT DETAILS
#'
#' Whatever the type (mixed or not) and the number of covariates in the data.frame of interest, the function \code{proxim_dist} firstly detects all the possible profiles (or combinations) of covariates from the two databases, and saves them in the output \code{profile}.
#' For example, assuming that a data.frame in input (composed of two overlayed data.frames A and B) have three shared binary covariates (ie identically encoded in A and B) so the sequences \code{011} and \code{101} will be considered as two distinct profiles of covariates.
#' If each covariate is a factor of \eqn{n_1}, \eqn{n_2} and \eqn{n_3} levels respectively, so it exists at most \eqn{n_1 \times n_2 \times n_3} possible profiles of covariates.
#' This number is considered as a maximum here because only the profiles of covariates met in at least one of the two databases will be kept for the study.
#'
#' \code{proxim_dist} classifies individuals from the two databases according to their proximities to each profile of covariates and saves the corresponding indexes of rows from A and B in two lists \code{indXA} and \code{indXB} respectively.
#' \code{indXA} and \code{indXB} thus contain as many objects as covariates profiles and the proximity between a given profile and a given individual is defined as follows.
#' The function also provides in output the list of all the encountered profiles of covariates.
#' As a decision rule, for a given profile of covariates \eqn{P_j}, an individual i will be considered as a neighbor of \eqn{P_j} if \eqn{dist(i,P_j) < prox \times max(dist(i,P_j))} where \code{prox} will be fixed by user.
#'
#'
#' @param data_file A data.frame corresponding ideally to an output of the function \code{\link{transfo_dist}}. Otherwise this data.frame is the result of two overlayed databases with a column of database identifier ("A" and "B", 1 and 2, for example), a target variable (called \eqn{Y} by example) only known in the first database, a target variable (\eqn{Z}) only stored in the second database, such that \eqn{Y} and \eqn{Z} summarize a same information differently encoded in the two databases and set of common covariates (at least one) of any type.
#' The order of the variables in the data.frame have no importance. The type of the covariates must be in accordance with the chosen distance function in the \code{norm} option.
#' @param indx_DB_Y_Z A vector of three column numbers corresponding to the colum indexes of the identifier, the target variable in the first database and the target variable in the second database. The indexes must be declared in this specific order.
#' @param norm A character (with quotes) indicating the choice of the distance function. This latest depends on the type of the common covariates: the Hamming distance
#' for binary covariates only (\code{norm} = "H"), the Manhattan distance ("M", by default) and the euclidean distance ("E") for continuous covariates only, or the Gower distance for mixed covariates ("G").
#' @param prox A percentage (betwen 0 and 1) used to calculate the distance threshold below which an individual (a row) is considered as a neighbor of a given profile of covariates.
#'
#' @return A list of 16 elements (the first 16 detailed below) is returned containing various distance matrices and lists useful for the algorithms that used Optimal Transportation theory. Two more objects (the last two of the following list) will be returned if distance matrices contain NAs.
#' \item{FILE_NAME}{A simple reminder of the name of the raw database}
#' \item{nA}{The row numbers of the first database (A)}
#' \item{nB}{The row numbers of the second database (B)}
#' \item{Xobserv}{The subset of the two merged databases composed of the common variables only}
#' \item{profile}{The different encountered profiles of covariates according to the data.frame}
#' \item{Yobserv}{The values of the target variable in the first database}
#' \item{Zobserv}{The values of the target variable in the second database}
#' \item{D}{A distance matrix corresponding to the computed distances between individuals of the two databases}
#' \item{Y}{The \eqn{n_Y} levels of the target variable in numeric form, in the first database}
#' \item{Z}{The \eqn{n_Z} levels of the target variable in numeric form, in the second database}
#' \item{indY}{A list of \eqn{n_Y} groups of individual (or row) numbers where each group corresponds to the individuals indexes related to a given level of \eqn{Y} in the first database}
#' \item{indZ}{A list of \eqn{n_Z} groups of individual (or row) numbers where each group corresponds to the individuals indexes related to a given level of \eqn{Z} in the second database}
#' \item{indXA}{A list of individual (row) indexes from the first database, sorted by profiles of covariates according to their proximities. See the \code{Details} part for more information}
#' \item{indXB}{A list of individual (row) indexes from the second database, sorted by profiles of covariates according to their proximities. See the \code{Details} part for more information}
#' \item{DA}{A distance matrix corresponding to the pairwise distances between individuals of the first database}
#' \item{DB}{A distance matrix corresponding to the pairwise distances between individuals of the second database}
#' \item{ROWS_TABLE}{Combinations of row numbers of the two databases that generate NAs in D}
#' \item{ROWS_TO_RM}{Number of times a row of the first or second database is involved in the NA process of D}
#'
#'
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @importFrom  rdist cdist
#' @importFrom  proxy dist
#' @importFrom  stats na.omit
#' @importFrom  StatMatch gower.dist
#'
#' @seealso \code{\link{transfo_dist}}, \code{\link{imput_cov}}, \code{\link{merge_dbs}}, \code{\link{simu_data}}
#'
#' @references
#' \enumerate{
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679 | \url{https://doi.org/10.1515/ijb-2018-0106}
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association, DOI: 10.1080/01621459.2020.1775615
#' \item Anderberg, M.R. (1973), Cluster analysis for applications, 359 pp., Academic Press, New York, NY, USA.
#' \item Gower, J. C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 27, 623--637.
#' \item D'Orazio M. (2015). Integration and imputation of survey data in R: the StatMatch package. Romanian Statistical Review, vol. 63(2)
#' \item Borg, I. and Groenen, P. (1997) Modern Multidimensional Scaling. Theory and Applications. Springer.
#' }
#'
#' @aliases proxim_dist
#'
#'
#' @export
#'
#' @examples
#'
#' data(simu_data)
#' ### The covariates of the data are prepared according to the chosen distance
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
proxim_dist  = function(data_file, indx_DB_Y_Z = 1:3, norm = "E", prox = 0.80){


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

  return(list(FILE_NAME = file_name, nA = nA, nB = nB, Xobserv = Xobserv, profile = unique(Xobserv),
              Yobserv   = Yobserv  , Zobserv = Zobserv, D = D, Y = Y, Z = Z,
              indY=indY, indZ=indZ, indXA=indXA, indXB=indXB, DA=DA, DB=DB))
}


