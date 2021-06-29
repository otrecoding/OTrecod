#' OT_joint()
#'
#' The function \code{OT_joint} integrates two algorithms called (\code{JOINT}) and (\code{R-JOINT}) dedicated to the solving of recoding problems in data fusion
#' using optimal transportation of the joint distribution of outcomes and covariates.
#'
#'
#' A. THE RECODING PROBLEM IN DATA FUSION
#'
#' Assuming that \eqn{Y} and \eqn{Z} are two target variables which refered to the same target population in two separate databases A and B respectively (no overlapping rows),
#' so that \eqn{Y} and \eqn{Z} are never jointly observed. Assuming also that A and B share a subset of common covariates \eqn{X} of any types (same encodings in A and B)
#' completed or not. Merging these two databases often requires to solve a recoding problem by creating an unique database where
#' the missing information of \eqn{Y} and \eqn{Z} is fully completed.
#'
#'
#' B. INFORMATIONS ABOUT THE ALGORITHM
#'
#' As with the function \code{\link{OT_outcome}}, the function \code{OT_joint} provides a solution to the recoding problem by proposing an
#' application of optimal transportation which aims is to search for a bijective mapping between the joint distributions of \eqn{(Y,X)} and \eqn{(Z,X)} in A and B (see (2) for more details).
#' The principle of the algorithm is also based on the resolution of an optimization problem, which provides a solution \eqn{\gamma} (as called in (1) and (2)), estimate
#' of the joint distribution of \eqn{(X,Y,Z)} according to the database to complete (see the argument \code{which.DB} for the choice of the database). While the algorithms \code{OUTCOME} and \code{R_OUTCOME} integrated in
#' the function \code{\link{OT_outcome}} require post-treatment steps to provide individual predictions, the algorithm \code{JOINT} directly uses estimations of the conditional distributions \eqn{(Y|Z,X)} in B and
#' \eqn{(Z|Y,X)} in A to predict the corresponding incomplete individuals informations of \eqn{Y} and/or \eqn{Z} respectively.
#' This algorithm supposes that the conditional distribution \eqn{(Y|X)} must be identical in A and B. Respectively, \eqn{(Z|X)} is supposed identical in A and B.
#' Estimations a posteriori of conditional probabilities \eqn{P[Y|X,Z]} and \eqn{P[Z|X,Y]} are available for each profiles of covariates in output (See the objects \code{estimatorYB} and \code{estimatorZA}).
#' Estimations of \eqn{\gamma} are also available according to the chosen transport distributions (See the arguments \code{gamma_A} and \code{gamma_B}).
#'
#' The algorithm \code{R-JOINT} gathers enrichments of the algorithm \code{JOINT} and is also available via the function \code{OT_joint}. It allows users to add a relaxation term in the algorithm to relax distributional assumptions (\code{maxrelax}>0),
#' and (or) add also a positive regularization term (\code{lamdba.reg}>0) expressing that the transportation map does not vary to quickly with respect of covariates \eqn{X}.
#' Is is suggested to users to calibrate these two parameters a posteriori by studying the stability of the individual predictions in output.
#'
#'
#' C. EXPECTED STRUCTURE FOR THE INPUT DATABASE
#'
#' The input database is a data.frame that must satisfy a specific form:
#' \itemize{
#' \item Two overlayed databases containing a common column of databases identifiers (A and B, 1 or 2, by examples, encoded in numeric or factor form)
#' \item A column corresponding to the target variable with its specific encoding in A (For example a factor \eqn{Y} encoded in \eqn{n_Y} levels, ordered or not, with NAs in the corresponding rows of B)
#' \item A column corresponding to another target outcome summarizing the same latent information with its specific encoding in B (By example a factor \eqn{Z} with \eqn{n_Z} levels, with NAs in rows of A)
#' \item The order of the variables in the database have no importance but the column indexes related to the three columns previously described (ie ID, \eqn{Y} and \eqn{Z}) must be rigorously specified
#' in the argument \code{index_DB_Y_Z}.
#' \item A set of shared common categorical covariates (at least one but more is recommended) with or without missing values (provided that the number of covariates exceeds 1) is required. On the contrary to the
#' function \code{OT_outcome}, please notice, that the function \code{OT_joint} does not accept continuous covariates therefore these latters will have to be categorized beforehand or using the provided input process (see \code{convert.num}).
#' }
#' The function \code{\link{merge_dbs}} is available in this package to assist user in the preparation of their databases.
#'
#' Remarks about the target variables:
#' \itemize{
#' \item A target variable can be of categorical type, but also discrete, stored in factor, ordered or not. Nevertheless, notice that, if the variable is stored in numeric it will be automatically converted in ordered factors.
#' \item If a target variable is incomplete, the corresponding rows will be automatically dropped during the execution of the function.
#' }
#' The type of each variables (including \eqn{ID}, \eqn{Y} and \eqn{Z}) of the database must be rigorously specified, in one of the four arguments \code{quanti}, \code{nominal}, \code{ordinal} and \code{logic}.
#'
#'
#' D. TRANSFORMATIONS OF CONTINUOUS COVARIATES
#'
#' Continuous shared variables (predictors) with infinite numbers of values have to be categorized before being introduced in the function.
#' To assist users in this task, the function \code{OT_joint} integrates in its syntax a process dedicated to the categorization of continuous covariates. For this, it is necessary to rigorously fill in
#' the arguments \code{quanti} and \code{convert.clss}.
#' The first one informs about the column indexes of the continuous variables to be transformed in ordered factor while the second one specifies the corresponding number of desired balanced levels (for unbalanced levels, users must do transformations by themselves).
#' Therefore \code{convert.num} and \code{convert.clss} must be vectors of same length, but if the length of \code{quanti} exceeds 1, while the length of \code{convert.clss} is 1, then, by default, all the covariates to convert will have the same number of classes (transformation by quantiles),
#' that corresponds to the value specified in the argument \code{convert.clss}.
#' Notice that only covariates can be transformed (not target variables) and that any incomplete information must have been taken into account beforehand (via the dedicated functions \code{\link{merge_dbs}} or \code{\link{imput_cov}} for examples).
#' Moreover, all the indexes informed in the argument \code{convert.num} must also be informed in the argument \code{quanti}.
#' Finally, it is recommended to declare all discrete covariates as ordinal factors using the argument \code{ordinal}.
#'
#'
#' E. INFORMATIONS ABOUT DISTANCE FUNCTIONS AND RELATED PARAMETERS
#'
#' Each individual (or row) of a given database is here characterized by a vector of covariates, so the distance between two individuals or groups of individuals depends on similarities between covariates
#' according to the distance function chosen by user (via the argument \code{dist.choice}). Actually four distance functions are implemented in \code{OT_joint} to take into account the most frequently encountered situation (see (3)):
#' \itemize{
#' \item the Manhattan distance ("M")
#' \item the Euclidean distance ("E")
#' \item the Gower distance for mixed data (see (4): "G")
#' \item the Hamming distance for binary data ("H")
#' }
#'
#' Finally, two profiles of covariates \eqn{P_1} (\eqn{n_1} individuals) and \eqn{P_2} (\eqn{n_2} individuals) will be considered as neighbors if \eqn{dist(P_1,P_2) < prox.X \times max(dist(P_i,P_j))} where \eqn{prox.X} must be fixed by user (\eqn{i = 1,\dots,n_1} and \eqn{j = 1,\dots,n_2}). This choice is used in the computation of the \code{JOINT} and \code{R_JOINT} algorithms.
#' The \code{prox.X} argument influences a lot the running time of the algorithm. The greater, the more the value will be close to 1, the more the convergence of the algorithm will be difficult or even impossible.
#'
#' Each individual \eqn{i} from A or B is here considered as a neighbor of only one profile of covariates \eqn{P_j}.
#'
#'
#' F. INFORMATIONS ABOUT THE SOLVER
#'
#' The argument \code{solvR} permits user to choose the solver of the optimization algorithm. The default solver is "glpk" that corresponds to the GNU Linear Programming Kit (see (5) for more details).
#' Moreover, the function actually uses the \code{R} optimization infrastructure of the package \pkg{ROI} which offers a wide choice of solver to users by easily loading the associated plugins of \pkg{ROI} (see (6)).
#'
#' For more details about the algorithms integrated in \code{OT_joint}, please consult (2).
#'
#'
#' @param datab a data.frame made up of two overlayed databases with at least four columns sorted in a random order. One column must be a column dedicated to the identification of the two databases ranked in ascending order
#' (For example: 1 for the top database and 2 for the database from below, or more logically here A and B  ...But not B and A!). One column (\eqn{Y} here but other names are allowed)
#' must correspond to the target variable related to the information of interest to merge with its specific encoding in the database A (corresponding encoding should be missing in the database B). In the same way,
#' one column (\eqn{Z} here) corresponds to the second target variable with its specific encoding in the database B (corresponding encoding should be missing in the database A).
#' Finally, the input database must have at least one shared covariate with same encoding in A and B. Please notice that, if your data.frame has only one shared covariate (four columns) with missing values (because no imputation is desired)
#' then a warning will appear and the algorithm will only run with complete cases.
#' @param index_DB_Y_Z a vector of three indexes of variables. The first index must correspond to the index of the databases identifier column. The second index corresponds
#' to the index of the target variable in the first database (A) while the third index corresponds to the column index related to the target variable in the second database (B).
#' @param nominal a vector of column indexes of all the nominal (not ordered) variables (database identifier and target variables included if it is the case for them).
#' @param ordinal a vector of column indexes of all the ordinal variables (database identifier and target variables included if it is the case for them).
#' @param logic a vector of column indexes of all the boolean variables of the data.frame.
#' @param convert.num indexes of the continuous (quantitative) variables. They will be automatically converted in ordered factors. By default, no continuous variables is assumed in the database.
#' @param convert.clss a vector indicating for each continuous variable to convert, the corresponding desired number of levels. If the length of the argument \code{convert_num} exceeds 1 while the length of \code{convert_clss} equals 1 (only one integer),
#' each discretization will count the same number of levels (quantiles).
#' @param dist.choice a character string (with quotes) corresponding to the distance function chosen between: the euclidean distance ("E", by default), the Manhattan distance ("M"),
#' the Gower distance ("G"), and the Hamming distance ("H") for binary covariates only.
#' @param percent.knn the ratio of closest neighbors involved in the computations of the cost matrices. 1 is the default value that includes all rows in the computation.
#' @param maxrelax the maximum percentage of deviation from expected probability masses. It must be equal to 0 (default value) for the \code{JOINT} algorithm, and equal to a strictly positive value for the R-JOINT algorithm.
#' @param lambda.reg a coefficient measuring the importance of the regularization term. It corresponds to the \code{R-JOINT} algorithm for a value other than 0 (default value).
#' @param prox.X a probability (betwen 0 and 1) used to calculate the distance threshold below which two covariates' profiles are supposed as neighbors.
#' If \code{prox.X = 1}, all profiles are considered as neighbors.
#' @param solvR a character string that specifies the type of method selected to solve the optimization algorithms. The default solver is "glpk".
#' @param which.DB a character string indicating the database to complete ("BOTH" by default, for the prediction of \eqn{Y} and \eqn{Z} in the two databases), "A" only for the imputation of \eqn{Z} in A, "B" only for the imputation of \eqn{Y} in B.
#'
#'
#' @return A "otres" class object of 9 elements:
#'     \item{time_exe}{running time of the function}
#'     \item{gamma_A}{estimate of \eqn{\gamma} for the completion of A. A matrix that corresponds to the joint distribution of \eqn{(Y,Z,X)} in A}
#'     \item{gamma_B}{estimate of \eqn{\gamma} for the completion of B. A matrix that corresponds to the joint distribution of \eqn{(Y,Z,X)} in B}
#'     \item{profile}{a data.frame that gives all details about the remaining P profiles of covariates. These informations can be linked to the \code{estimatorZA} and the \code{estimatorYB} objects for a better interpretation of the results.}
#'     \item{res_prox}{a \code{proxim_dist} object}
#'     \item{estimatorZA}{an array that corresponds to estimates of the probability distribution of \eqn{Z} conditional to \eqn{X} and \eqn{Y} in database A. The number of rows of each table corresponds to the total number of profiles of covariates.
#'     The first dimension of each table (rownames) correspond to the profiles of covariates sorted by order of appearance in the merged database. The second dimension of the array (columns of the tables) corresponds to the levels of \eqn{Y} while the third element corresponds to the levels of \eqn{Z}.}
#'     \item{estimatorYB}{an array that corresponds to estimates of the probability distribution of \eqn{Y} conditional to \eqn{X} and \eqn{Z} in database B. The number of rows of each table corresponds to the total number of profiles of covariates.
#'     The first dimension of each table (rownames) correspond to the profiles of covariates sorted by order of appearance in the merged database. The second dimension of the array (columns of the tables) corresponds to the levels of \eqn{Z} while the third element corresponds to the levels of \eqn{Y}.}
#'     \item{DATA1_OT}{the database A with the individual predictions of \eqn{Z} using an optimal transportation algorithm (\code{JOINT}) or \code{R-JOINT}}
#'     \item{DATA2_OT}{the database B with the individual predictions of \eqn{Y} using an optimal transportation algorithm (\code{JOINT}) or \code{R-JOINT}}
#'
#' @import ROI ROI.plugin.glpk
#'
#' @importFrom dplyr %>%
#' @importFrom ompr MIPModel get_solution
#' @importFrom ompr.roi with_ROI
#' @importFrom rdist cdist
#' @importFrom proxy dist
#' @importFrom plyr mapvalues
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @aliases OT_joint
#'
#' @references
#' \enumerate{
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679. doi:10.1515/ijb-2018-0106
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association. \doi{10.1080/01621459.2020.1775615}
#' \item Anderberg, M.R. (1973), Cluster analysis for applications, 359 pp., Academic Press, New York, NY, USA.
#' \item Gower J.C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 27, 623–637
#' \item Makhorin A (2011). GNU Linear Programming Kit Reference Manual Version 4.47.\url{http://www.gnu.org/software/glpk/}
#' \item Theussl S, Schwendinger F, Hornik K (2020). ROI: An Extensible R Optimization Infrastructure.Journal of Statistical Software,94(15), 1-64. \doi{10.18637/jss.v094.i15}
#' }
#'
#' @seealso \code{\link{merge_dbs}}, \code{\link{OT_outcome}}, \code{\link{proxim_dist}}, \code{\link{avg_dist_closest}}
#'
#' @export
#'
#' @examples
#'
#' ### An example of JOINT algorithm with:
#' #-----
#' # - A sample of the database tab_test
#' # - Y1 and Y2 are a 2 outcomes encoded in 2 different forms in DB 1 and 2:
#' #   4 levels for Y1 and 3 levels for Y2
#' # - n1 = n2 = 40
#' # - 2 discrete covariates X1 and X2 defined as ordinal
#' # - Distances estimated using the Gower function
#' # Predictions are assessed for Y1 in B only
#' #-----
#'
#' data(tab_test)
#' tab_test2 = tab_test[c(1:40,5001:5040),1:5]
#'
#'
#' try1J = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
#'                  dist.choice = "G", which.DB = "B")
#'
#' \donttest{
#'
#' ### An example of R-JOINT algorithm using the previous database,
#' ### and keeping the same options excepted for:
#' #-----
#' # - The distances are estimated using the Gower function
#' # - Inclusion of an error term in the constraints on
#' #   the marginals (relaxation term)
#' # Predictions are assessed for Y1 AND Y2 in A and B respectively
#' #-----
#'
#' try1RJ = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
#'                   dist.choice = "G", maxrelax = 0.4,
#'                   which.DB = "BOTH")
#'
#' ### The previous example of R-JOINT algorithm with:
#' # - Adding a regularization term
#' # Predictions are assessed for Y1 and Y2 in A and B respectively
#' #-----
#'
#' try2RJ = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
#'                   dist.choice = "G", maxrelax = 0.4, lambda.reg = 0.9,
#'                   which.DB = "BOTH")
#'
#'
#'
#' ### Another example of JOINT algorithm with:
#' #-----
#' # - A sample of the database simu_data
#' # - Y1 and Y2 are a 2 outcomes encoded in 2 different forms in DB A and B:
#' #   (3 levels for Y and 5 levels for Z)
#' # - n1 = n2 = 100
#' # - 3 covariates: Gender, Smoking and Age in a qualitative form
#' # - Complete Case study
#' # - The Hamming distance
#' # Predictions are assessed for Y1 and Y2 in A and B respectively
#' #-----
#'
#' data(simu_data)
#' simu_data2 = simu_data[c(1:100,401:500),c(1:4,7:8)]
#' simu_data3 = simu_data2[!is.na(simu_data2$Age),]
#'
#' try2J = OT_joint(simu_data3, convert.num = 6, convert.clss = 3,
#'                  nominal = c(1,4:5), ordinal = 2:3,
#'                  dist.choice = "H", which.DB = "BOTH")
#'
#' }

OT_joint = function(datab, index_DB_Y_Z = 1:3,
                    nominal = NULL, ordinal = NULL,logic = NULL,
                    convert.num = NULL, convert.clss = NULL, dist.choice = "E", percent.knn = 1,
                    maxrelax = 0, lambda.reg = 0.0, prox.X = 0.10, solvR = "glpk", which.DB = "BOTH"){

  if (dist.choice %in% c("M","Manhattan","manhattan")){

    dist.choice = "M"

  } else if (dist.choice %in% c("E","Euclidean","euclidean")){

    dist.choice = "E"

  } else if (dist.choice %in% c("G","Gower","gower")){

    dist.choice = "G"

  } else if (dist.choice %in% c("H","Hamming","hamming")){

    dist.choice = "H"

  } else {

    stop("Invalid dist.choice argument")

  }


  if (!(which.DB %in% c("A","B","BOTH"))){

    stop("Invalid which.DB argument")

  } else {}


  if (percent.knn > 1){

    stop("Improper value for percent.knn")

  } else {}

  if (prox.X > 1){

    stop("Improper value for prox.X")

  } else {}

  if (length(convert.num) < length(convert.clss)){

    stop("The arguments quanti and convert.num must be equal because the algorithm does not handle continuous covariates")

  } else {}

  message("---------------------------------------","\n")
  message("OT JOINT PROCEDURE in progress ..."     ,"\n")
  message("---------------------------------------","\n")
  message("Type                  = ", ifelse((maxrelax == 0)&(lambda.reg == 0),"JOINT","R-JOINT"),"\n")
  message("Distance              = ", ifelse(dist.choice == "H","Hamming",ifelse(dist.choice == "M","Manhattan",ifelse(dist.choice == "E","Euclidean","Gower"))),"\n")
  message("Percent closest       = ", 100.0*percent.knn, "%","\n")
  message("Relaxation term       = ", maxrelax, "\n")
  message("Regularization term   = ", lambda.reg  ,"\n")
  message("Aggregation tol cov   = ", prox.X,"\n")
  message("DB imputed            = ", which.DB,"\n")
  message("---------------------------------------","\n")


  tstart = Sys.time()

  dataB = transfo_dist(datab,index_DB_Y_Z = index_DB_Y_Z,
                       quanti = convert.num, nominal = nominal, ordinal = ordinal, logic = logic,
                       convert_num = convert.num, convert_clss = convert.clss,
                       prep_choice = dist.choice)

  if (dist.choice == "H"){

    datac  = dataB[,-index_DB_Y_Z]

    test_H = apply(as.data.frame(datac),2,function(x){length(table(x))})

    if (max(test_H)>2){stop("With Hamming distance, all your covariates must be binaries !")} else {}

  } else {}


  inst = proxim_dist(dataB, norm = dist.choice, prox = 0)


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
  # prof    = do.call(paste0,unique(inst$Xobserv))
  prof    = as.data.frame(unique(inst$Xobserv))
  ID_prof = paste(rep("P",nrow(prof)),1:nrow(prof),sep="_")


  # Create a model for the optimal transport of individuals
  # modelA = Model(with_optimizer(Gurobi.Optimizer,LogToConsole=0,Method=2,Crossover=0));#ClpSolver(LogLevel=0));
  # modelB = Model(with_optimizer(Gurobi.Optimizer,LogToConsole=0,Method=2,Crossover=0));#Model(with_optimizer(Clp.Optimizer,LogLevel=0));


  ###########################################################################
  # Compute data for aggregation of the individuals
  ###########################################################################
  # println("... aggregating individuals")

  indXA = inst$indXA
  indXB = inst$indXB
  nbX   = length(indXA)

  # compute the neighbors of the covariates for regularization
  Xvalues = unique(Xobserv)

  if (dist.choice == "M"){

    dist_X  = proxy::dist(Xvalues,Xvalues, method = "manhattan")

  } else if (dist.choice == "E"){

    dist_X  = proxy::dist(Xvalues,Xvalues, method = "euclidean")

  } else if (dist.choice == "H"){

    if (nrow(dataB[,4:ncol(dataB)])== nrow(na.omit(dataB[,4:ncol(dataB)]))){

      dist_X  = rdist::cdist(Xvalues,Xvalues, metric = "hamming")

    } else {

      dist_X  = ham(Xvalues,Xvalues)

    }

  } else if (dist.choice == "G"){

    dist_X  = StatMatch::gower.dist(Xvalues,Xvalues)

  }

  tol.X        = prox.X * max(dist_X)
  voisins_X    = dist_X <= tol.X
  C            = avg_dist_closest(inst, percent_closest = percent.knn)[[1]]

  ###########################################################################
  # Compute the estimators that appear in the model
  ###########################################################################

  estim_XA = estim_XB = estim_XA_YA =  estim_XB_ZB = list()

  for (x in 1:nbX){

    estim_XA[[x]] = length(indXA[[x]])/nA
    estim_XB[[x]] = length(indXB[[x]])/nB

    # estim_XA[[x]] = length(indXA[[x]])/length(unlist(indXA))
    # estim_XB[[x]] = length(indXB[[x]])/length(unlist(indXB))

  }

  for (x in 1:nbX){

    estim_XA_YA[[x]] = estim_XB_ZB[[x]] = numeric(0)


    for (y in Y){
      estim_XA_YA[[x]][y] = length(indXA[[x]][Yobserv[indXA[[x]]] == y])/nA
      # estim_XA_YA[[x]][y] = length(indXA[[x]][Yobserv[indXA[[x]]] == y])/length(unlist(indXA))
    }

    for (z in Z){
      estim_XB_ZB[[x]][z] = length(indXB[[x]][Zobserv[indXB[[x]] + nA] == z])/nB
      # estim_XB_ZB[[x]][z] = length(indXB[[x]][Zobserv[indXB[[x]] + nA] == z])/length(unlist(indXB))
    }

  }


  Cf <- function(y,z) {
    return(C[y,z])
  }

  estim_XBf <- function(x) {
    return(estim_XB[[x]])
  }

  voisin = function(x1){lambda.reg *(1/length(voisins_X[x1,]))}

  ind_voisins = list()
  for (x1 in 1:nrow(voisins_X)){

    ind_voisins[[x1]] = which(voisins_X[x1,])

  }


  ###########################################################################
  # Basic part of the model
  ###########################################################################

  if (which.DB %in% c("A","BOTH")){

    # COMPLETE Z IN DATABASE A

    result <-  MIPModel() %>%
      # DEFINE VARIABLES ----------------------------------------------------------------------
    # gammaA[x,y,z]: joint probability of X=x, Y=y and Z=z in base A
    ompr::add_variable(gammaA[x,y,z]       , x = 1:nbX, y = Y, z= Z, type = "continuous", lb = 0, ub =1) %>%
      ompr::add_variable(errorA_XY[x,y]      , x = 1:nbX, y = Y,       type = "continuous") %>%
      ompr::add_variable(abserrorA_XY[x,y]   , x = 1:nbX, y = Y,       type = "continuous", lb = 0, ub =1) %>%
      ompr::add_variable(errorA_XZ[x,z]      , x = 1:nbX, z = Z,       type = "continuous") %>%
      ompr::add_variable(abserrorA_XZ[x,z]   , x = 1:nbX, z = Z,       type = "continuous", lb = 0, ub =1) %>%
      # REGULARIZATION ---------------------------------------------------------------------------------
    # ompr::add_variable(reg_absA[x1,x2,y,z] , x1 = 1:nbX, x2 = 1:nbX, y= Y, z= Z, type = "continuous") %>%
    ompr::add_variable(reg_absA[x1,x2,y,z]   , x1 = 1:nbX, x2 = ind_voisins[[x1]], y= Y, z= Z, type = "continuous", lb = 0) %>%
      # OBJECTIVE ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ompr::set_objective(sum_expr(Cf(y,z) * gammaA[x,y,z], y = Y, z = Z, x = 1:nbX) +
                          sum_expr(voisin(x1)*reg_absA[x1,x2,y,z], x1 = 1:nbX, x2 = ind_voisins[[x1]],y=Y,z= Z), "min") %>%
      # CONSTRAINTS ----------------------------------------------------------------------------------------------------
    ompr::add_constraint(sum_expr(gammaA[x,y,z], z = Z) - errorA_XY[x,y] == estim_XA_YA[[x]][y] , x = 1:nbX,y =Y) %>%
      ompr::add_constraint(estim_XBf(x)*sum_expr(gammaA[x,y,z],y = Y)  - estim_XBf(x)*errorA_XZ[x,z]== estim_XB_ZB[[x]][z] * estim_XA[[x]] ,
                           x = 1:nbX, z = Z) %>%

      ompr::add_constraint( errorA_XY[x,y] <= abserrorA_XY[x,y], x = 1:nbX,y =Y) %>%
      ompr::add_constraint(-errorA_XY[x,y] <= abserrorA_XY[x,y], x = 1:nbX,y =Y) %>%

      ompr::add_constraint(sum_expr(abserrorA_XY[x,y], x = 1:nbX, y = Y)<= maxrelax/2.0) %>%
      ompr::add_constraint(sum_expr(errorA_XY[x,y]   , x = 1:nbX, y = Y) == 0.0) %>%

      ompr::add_constraint( errorA_XZ[x,z] <= abserrorA_XZ[x,z], x = 1:nbX, z =Z) %>%
      ompr::add_constraint(-errorA_XZ[x,z] <= abserrorA_XZ[x,z], x = 1:nbX, z =Z) %>%

      ompr::add_constraint(sum_expr(abserrorA_XZ[x,z], x = 1:nbX, z = Z)<= maxrelax/2.0) %>%
      ompr::add_constraint(sum_expr(errorA_XZ[x,z]   , x = 1:nbX, z = Z) == 0) %>%
      # ompr::add_constraint(sum_expr(errorA_XZ[x,z]   , x = 1:nbX, z = Z)<= maxrelax/2.0) %>%

      # Constraints regularization - Ajout GG
      ompr::add_constraint(reg_absA[x1,x2,y,z] + gammaA[x2,y,z]/(max(1,length(indXA[[x2]]))/nA) >= gammaA[x1,y,z]/(max(1,length(indXA[[x1]]))/nA),
                           x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z) %>%
      ompr::add_constraint(reg_absA[x1,x2,y,z] + gammaA[x1,y,z]/(max(1,length(indXA[[x1]]))/nA) >= gammaA[x2,y,z]/(max(1,length(indXA[[x2]]))/nA),
                           x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z) %>%

    # SOLUTION -------------------------------------------------------
    ompr::solve_model(with_ROI(solver = solvR))
    solution  = ompr::get_solution(result, gammaA[x,y,z])
    gammaA_val= array(solution$value,dim = c(nbX,length(Y),length(Z)))
    #------------ END OPTIMIZATION STEP ------------------------------


    ### Compute the resulting estimators for the distributions of Z conditional to X and Y in base A

    estimatorZA = 1/length(Z) * array(rep(1,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))

    for (x in 1:nbX){

      for (y in Y){

        proba_c_mA = apply(gammaA_val,c(1,2),sum)[x,y]

        if (proba_c_mA > 1.0e-6){

          estimatorZA[x,y,] = 1/proba_c_mA * gammaA_val[x,y,];

        } else {}
      }
    }

    row.names(estimatorZA) =  ID_prof
    colnames(estimatorZA)  = as.character(levels(dataB[,2]))


    ### Deduce the individual distributions of probability for each individual

    probaZindivA = matrix(0,nA,length(Z))

    for (x in 1:nbX){

      for (i in indXA[[x]]){

        probaZindivA[i,] = estimatorZA[x,Yobserv[i],]

      }

    }


    # Transport the probability that maximizes frequency

    predZA = apply(probaZindivA,1,function(x)which.max(x))

    DATA1_OT         = dataB[dataB[,1] == unique(dataB[,1])[1],]


    if (index_DB_Y_Z[3] %in% nominal){

      DATA1_OT$OTpred  = factor(plyr::mapvalues(predZA,
                                                from    = sort(unique(predZA)),
                                                to      = levels(dataB[,3])[sort(unique(predZA))]),
                                levels  = levels(dataB[,3])[sort(unique(predZA))])
    } else {

      DATA1_OT$OTpred  = ordered(plyr::mapvalues(predZA,
                                                 from   = sort(unique(predZA)),
                                                 to     = levels(dataB[,3])[sort(unique(predZA))]),
                                 levels = levels(dataB[,3])[sort(unique(predZA))])

    }

  } else {}

  #--->  END FOR DATABASE A


  if (which.DB %in% c("B","BOTH")){

    # COMPLETE Y IN DATABASE B

    result <-  ompr::MIPModel() %>%
      # DEFINE VARIABLES -----------------------------------------------------------------------------------
    # gammaA[x,y,z]: joint probability of X=x, Y=y and Z=z in base B

    ompr::add_variable(gammaB[x,y,z]    , x = 1:nbX, y = Y, z = Z,type = "continuous", lb = 0, ub = 1) %>%

      ompr::add_variable(errorB_XY[x,y]   , x = 1:nbX, y = Y,       type = "continuous") %>%
      ompr::add_variable(abserrorB_XY[x,y], x = 1:nbX, y = Y,       type = "continuous", lb = 0, ub = 1) %>%
      ompr::add_variable(errorB_XZ[x,z]   , x = 1:nbX,        z = Z,type = "continuous") %>%
      ompr::add_variable(abserrorB_XZ[x,z], x = 1:nbX,        z = Z,type = "continuous", lb = 0, ub = 1) %>%

      # REGULARIZATION ------------------------------------------------------------------------------------------

    # ompr::add_variable(reg_absB[x1, x2,y,z], x1 = 1:nbX, x2 = 1:nbX, y = Y, z = Z, type = "continuous") %>%
    ompr::add_variable(reg_absB[x1,x2,y,z], x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z, type = "continuous", lb = 0) %>%

      # OBJECTIVE ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ompr::set_objective(sum_expr(Cf(y,z) * gammaB[x,y,z],y = Y, z = Z, x = 1:nbX)  + sum_expr(voisin(x1)*reg_absB[x1,x2,y,z],x1=1:nbX, x2= ind_voisins[[x1]],y=Y,z= Z), "min")  %>%

      # CONSTRAINTS ----------------------------------------------------------------------------------------------------

    ompr::add_constraint(sum_expr(gammaB[x,y,z], y = Y) -errorB_XZ[x,z]  == estim_XB_ZB[[x]][z] , x = 1:nbX, z = Z) %>%

      ompr::add_constraint(estim_XA[[x]]*sum_expr(gammaB[x,y,z] ,z = Z) - estim_XA[[x]] * errorB_XY[x,y] == estim_XA_YA[[x]][y] * estim_XB[[x]] , x = 1:nbX, y = Y) %>%

      ompr::add_constraint( errorB_XY[x,y] <= abserrorB_XY[x,y], x = 1:nbX,y =Y) %>%
      ompr::add_constraint(-errorB_XY[x,y] <= abserrorB_XY[x,y], x = 1:nbX,y =Y) %>%

      ompr::add_constraint(sum_expr(abserrorB_XY[x,y], x = 1:nbX, y = Y)<= maxrelax/2.0) %>%
      ompr::add_constraint(sum_expr(errorB_XY[x,y]   , x = 1:nbX, y = Y) == 0.0) %>%

      ompr::add_constraint( errorB_XZ[x,z] <= abserrorB_XZ[x,z], x = 1:nbX,z =Z) %>%
      ompr::add_constraint(-errorB_XZ[x,z] <= abserrorB_XZ[x,z], x = 1:nbX,z =Z) %>%

      ompr::add_constraint(sum_expr(abserrorB_XZ[x,z], x = 1:nbX, z = Z)<= maxrelax/2.0) %>%
      ompr::add_constraint(sum_expr(errorB_XZ[x,z]   , x = 1:nbX, z = Z) == 0) %>%
      # ompr::add_constraint(sum_expr(errorB_XZ[x,z]   , x = 1:nbX, z = Z)<= maxrelax/2.0) %>%

      # Constraints regularization - Ajout GG
      ompr::add_constraint(reg_absB[x1,x2,y,z] + gammaB[x2,y,z]/(max(1,length(indXB[[x2]]))/nB) >= gammaB[x1,y,z]/(max(1,length(indXB[[x1]]))/nB), x1 = 1:nbX, x2 = ind_voisins[[x1]], y= Y, z = Z) %>%
      ompr::add_constraint(reg_absB[x1,x2,y,z] + gammaB[x1,y,z]/(max(1,length(indXB[[x1]]))/nB) >= gammaB[x2,y,z]/(max(1,length(indXB[[x2]]))/nB), x1 = 1:nbX, x2 = ind_voisins[[x1]], y= Y, z = Z) %>%

      # SOLUTION ------------------------------------------------------
    ompr::solve_model(with_ROI(solver = solvR))

    solution   = get_solution(result, gammaB[x,y,z])
    gammaB_val = array(solution$value,dim = c(nbX,length(Y),length(Z)))
    #------------ END OPTIMIZATION STEP -------------------------------


    ### compute the resulting estimators for the distributions of Y conditional to X and Z in base B

    estimatorYB = 1/length(Y) * array(rep(1,nbX*length(Y)*length(Z)),dim = c(nbX,length(Z),length(Y)))

    for (x in 1:nbX){

      for (z in Z){

        proba_c_mB = apply(gammaB_val,c(1,3),sum)[x,z]

        if (proba_c_mB > 1.0e-6){

          estimatorYB[x,z,] = 1/proba_c_mB * gammaB_val[x,,z];

        } else {}
      }
    }

    row.names(estimatorYB) = ID_prof
    colnames(estimatorYB)  = as.character(levels(dataB[,3]))


    ### Deduce the individual distributions of probability for each individual

    probaZindivA = matrix(0,nA,length(Z))
    probaYindivB = matrix(0,nB,length(Y))

    for (x in 1:nbX){

      for (i in indXB[[x]]){

        probaYindivB[i,] = estimatorYB[x,Zobserv[i+nA],]

      }
    }


    ### Transport the Ylity that maximizes frequency

    predYB = apply(probaYindivB,1,function(x)which.max(x))


    DATA2_OT         = dataB[dataB[,1] == unique(dataB[,1])[2],]

    if (index_DB_Y_Z[2] %in% nominal){

      DATA2_OT$OTpred  = factor(plyr::mapvalues(predYB,from = sort(unique(predYB)),
                                                to = levels(dataB[,2])[sort(unique(predYB))]),
                                levels = levels(dataB[,2])[sort(unique(predYB))])
    } else {

      DATA2_OT$OTpred  = ordered(plyr::mapvalues(predYB,from = sort(unique(predYB)),
                                                 to = levels(dataB[,2])[sort(unique(predYB))]),
                                 levels = levels(dataB[,2])[sort(unique(predYB))])
    }


  } else {}

  # ---> END DATABASE B


  # OUTPUT OBJECTS -------------------------------------------------

  if (which.DB == "A"){

    GAMMA_B            = NULL
    estimatorYB        = NULL
    DATA2_OT           = dataB[dataB[,1] == unique(dataB[,1])[2],]
    GAMMA_A            = apply(gammaA_val,c(2,3),sum)
    colnames(GAMMA_A)  = levels(dataB[,3])
    row.names(GAMMA_A) = levels(dataB[,2])

  } else if (which.DB == "B"){

    GAMMA_A            = NULL
    estimatorZA        = NULL
    DATA1_OT           = dataB[dataB[,1] == unique(dataB[,1])[1],]
    GAMMA_B            = apply(gammaB_val,c(2,3),sum)
    colnames(GAMMA_B)  = levels(dataB[,3])
    row.names(GAMMA_B) = levels(dataB[,2])

  } else {

    GAMMA_A            = apply(gammaA_val,c(2,3),sum)
    GAMMA_B            = apply(gammaB_val,c(2,3),sum)
    colnames(GAMMA_A)  = colnames(GAMMA_B)  = levels(dataB[,3])
    row.names(GAMMA_A) = row.names(GAMMA_B) = levels(dataB[,2])

  }

  tend = Sys.time()

  res_OT = list(time_exe    = difftime(tend,tstart),
                gamma_A     = GAMMA_A,
                gamma_B     = GAMMA_B,
                profile     = data.frame(ID = ID_prof,prof),
                res_prox    = inst,
                estimatorZA = estimatorZA,
                estimatorYB = estimatorYB,
                DATA1_OT    = DATA1_OT,
                DATA2_OT    = DATA2_OT)

  # otres class object
  class(res_OT) = "otres"

  return(res_OT)

}

