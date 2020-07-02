#' OT_joint()
#'
#' The function \code{OT_joint} integrates two models called (\code{JOINT}) and (\code{R-JOINT}) dedicated to the solving of recoding problems in data integration
#' using optimal transportation of the joint distribution of outcomes and covariates.
#'
#'
#' A. THE RECODING PROBLEM IN DATA FUSION
#'
#' Assuming that \eqn{Y} and \eqn{Z} are two variables that summarize a same latent information in two separate (no overlapping rows) databases A and B respectively,
#' so that \eqn{Y} and \eqn{Z} are never jointly observed in A and B. Assuming also that A and B share a subset of common covariates X of any types (with same encodings in A and B)
#' completed or not. Integrating these two databases often requires to solve the recoding problem observed between \eqn{Y} and \eqn{Z} by creating an unique database where
#' the missing information of \eqn{Y} and/or \eqn{Z} is fully completed.
#'
#'
#' B. INFORMATIONS ABOUT THE ALGORITHM
#'
#' As with the function \code{\link{OT_outcome}}, the function \code{OT_joint} provides a solution to the recoding problem by proposing an
#' application of optimal transportation which aims is to search for a bijective mapping between the joint distributions of \eqn{(Y,X)} and \eqn{(Z,X)} in A and B (see (2) for more details).
#' The principle of the algorithm is also based on the resolution of an optimization problem, which provides a \eqn{\gamma} solution (as called in (1) and (2)) which is an estimation
#' of the joint distribution of \eqn{(X,Y,Z)} according to the database to complete (see the argument \code{which.DB} for the choice of the database). While the models \code{OUTCOME} and \code{R_OUTCOME} integrated in
#' the function \code{\link{OT_outcome}} require post-treatment steps to provide individual predictions, the algorithm \code{JOINT} directly uses estimations of the conditional distributions \eqn{(Y|Z,X)} in B and
#' \eqn{(Z|Y,X)} in A to predict the corresponding incomplete individuals informations of \eqn{Y} and/or \eqn{Z} respectively.
#' This algorithm supposes that the conditional distribution \eqn{(Y|X)} must be identical in A and B. Respectively, \eqn{(Z|X)} is supposed identical in A and B.
#' Estimations a posteriori of conditional probabilities \eqn{P[Y|X,Z]} and \eqn{P[Z|X,Y]} are available by profiles of covariates in output (See the objects \code{estimatorYB} and \code{estimatorZA}).
#' Estimations of \eqn{\gamma} are also available according to the chosen transport distributions (See the arguments \code{gamma_A} and \code{gamma_B}).
#'
#' The model \code{R-JOINT} gathers enrichments of the model \code{JOINT} and is also available via the function \code{OT_joint} It allows users to add a relaxation term in the algorithm to relax distributional assumptions (\code{maxrelax}>0),
#' and (or) add also a positive regularization term (\code{lamdba.reg}>0) expressing that the transportation map does not vary to quickly with respect of covariates \eqn{X}.
#' Is is suggested to users to calibrate these two parameters a posteriori by studying the stability of the individual predictions in output.
#'
#'
#' C. TYPE OF REQUIRED INPUT DATABASE
#'
#' The input database is a data.frame that must be saved in a specific form by users:
#' \itemize{
#' \item Two overlayed databases containing a common column of databases' identifiers (A and B, 1 or 2, by examples, encoded in numeric or factor form)
#' \item A column corresponding to the target variable with its specific encoding in A (For example a factor \eqn{Y} encoded in \eqn{n_Y} levels, ordered or not, with NAs in the corresponding rows of B)
#' \item A column corresponding to another target outcome summarizing the same latent information with its specific endoded in B (By example a factor \eqn{Z} in \eqn{n_Z} levels, with NAs in rows of A)
#' \item The order of the variables in the database have no importance but the indexes of the columns related to the 3rd columns previously described (ie ID, \eqn{Y} and \eqn{Z}) must be rigorously specified
#' in the argument \code{index_DB_Y_Z}.
#' \item A set of shared common categorical covariates (at least one but more is recommended) complete or not (provided that the number of covariates exceeds 1) is required. On the contrary to the
#' function \code{OT_outcome}, please notice, that the function \code{OT_joint} does not accept continuous covariates therefore these latters will have to be categorized beforehand or using the provided input process (see \code{quanti}).
#' }
#' The function \code{\link{merge_dbs}} is available in this package to assist user in the preparation of their databases, so please, do not hesitate to use it beforehand if necessary.
#'
#' Remarks about the outcomes:
#' \itemize{
#' \item A target outcome can be categorical, in factor, ordered or not, discrete (with a finite number of values ONLY) but, notice that, if they are stored in numeric they will be automatically converted in ordered factors.
#' \item If a target outcome is incomplete, the corresponding rows will be automatically dropped during the execution of the function.
#' }
#' The type of each variables (including \eqn{ID}, \eqn{Y} and \eqn{Z}) of the database must be rigorously specified, in one of the four arguments \code{quanti}, \code{nominal}, \code{ordinal} and \code{logic}.
#'
#'
#' D. TRANSFORMATIONS OF CONTINUOUS COVARIATES
#'
#' Continuous covariates with infinite numbers of values will have to be categorized beforehand their inclusions in the function.
#' To assist users in this task, the function \code{OT_joint} integrates in its syntax a process dedicated to the categorization of continuous covariates. For this, it is necessary to rigorously fill in
#' the arguments \code{quanti} and \code{convert.clss}.
#' The first one informs about the indexes in database of the continuous variables to transform in ordered factor while the second one specifies the corresponding number of desired balanced levels (for unbalanced levels, users must do transformations by themselves).
#' Therefore \code{convert.num} and \code{convert.clss} must be vectors of same length, but if the length of \code{quanti} exceeds 1, while the length of \code{convert.clss} is 1, then, by default, all the covariates to convert will have the same number of classes (transformation by quantiles),
#' that corresponds to the value specified in the argument \code{convert.clss}.
#' Please notice that only covariates can be transformed (not outcomes) and missing informations are not taken into account for the transformations.
#' Moreover, all the indexes informed in the argument \code{convert.num} must also be informed in the argument \code{quanti}.
#' Finally, it is suggested to declare all discrete covariates as ordinal factors using the argument \code{ordinal}.
#'
#'
#' E. INFORMATIONS ABOUT DISTANCE FUNCTIONS AND RELATED PARAMETERS
#'
#' Each individual (or row) of a given database is here characterized by his covariates, so the distance between 2 individuals or groups of individuals depends on similarities between covariates
#' according to the distance function chosen by user (via the argument \code{dist.choice}). Actually four distance functions are implemented in \code{OT_joint} to take into account the most frequently encountered situation (see (3)):
#' \itemize{
#' \item the Manhattan distance ("M")
#' \item the Euclidean distance ("E")
#' \item the Gower distance for mixed data (see (4): "G")
#' \item the Hamming distance for binary data ("H")
#' }
#'
#' Let two profiles of covariates \eqn{P_1} (\eqn{n_1} individuals) and \eqn{P_2} (\eqn{n_2} individuals) will be considered as neighbors if \eqn{dist(P_1,P_2) < prox.X \times max(dist(P_i,P_j))} where \eqn{prox.X} must be fixed by user (\eqn{i = 1,\dots,n_1} and \eqn{j = 1,\dots,n_2}). This choice is used in the calculation of the \code{JOINT} and \code{R_JOINT} algorithms.
#' In the same way, for a given profile of covariates \eqn{Pj}, an individual i will be considered as a neighbor of \eqn{Pj} if \eqn{dist(i,P_j) < prox.dist \times max(dist(i,P_j))} where \eqn{prox.dist} will be fixed by user.
#'
#' For more details about the related algorithms integrated in \code{OT_joint}, please consult (2).
#'
#'
#' @param datab A data.frame with at least four columns sorted in a non-specific order. One column must be a column of identifier for the two databases, where the names of the two databases
#' must be ranked in ascending order (By example: 1 for the top database and 2 for the database from below, or more logically here A and B  ...But NOT B and A!). One column (\eqn{Y} here but other names are allowed)
#' must correspond to the target variable related to the information of interest to merge with its specific encoding in the database A (corresponding encoding should be so missing in the database B). In the same way,
#' one column (\eqn{Z} here) corresponds to the second target variable with its specific encoding in the database B (corresponding encoding should be so missing in the database A).
#' Finally, the input database must have at least one shared covariate with same encoding in A and B. Please notice that, if your data.frame has only four columns, that is to say with only one covariate, if this latter has NA, and unless user
#' has previously imputed the corresponding missing information, the OT algorithm will only run complete cases.
#' @param index_DB_Y_Z A vector of three column indexes. The first index must correspond to the index of the databases identifier column (DB identifier). The second index must correspond
#' to the index of the target variable in the first database (A) while the third index corresponds to the index of column related to the target variable in the second database (B).
#' @param nominal A vector of indexes that corresponds to the column indexes of all the nominal (not ordered) variables (DB identification and target variables included if it is the case for them).
#' @param ordinal A vector of indexes that corresponds to the column indexes of all the ordinal variables (DB identification and target variables included if it is the case for them).
#' @param logic A vector of indexes that corresponds to the column indexes of all the boolean variables of the data.frame.
#' @param convert.num Indexes of the continuous (quantitative) variables. They will be automatically converted in ordered factors. By default, no continuous variables is assumed in the database.
#' @param convert.clss A vector indicating for each continuous variable to convert, the corresponding desired number of levels. If the length of the argument \code{convert_num} exceeds 1 while the length of \code{convert_clss} equals 1 (only one integer),
#' each discretization will count the same number of levels.
#' @param dist.choice A character (with quotes) corresponding to the distance function chosen between: the euclidean distance ("E", by default), the Manhattan distance ("M"),
#' the Gower distance ("G"), and the Hamming distance ("H") for binary covariates only.
#' @param percent.knn Percent of closest neighbors taken in the computation of the cost matrix.
#' @param maxrelax The maximum percentage of deviation from expected probability masses. It must be equal to 0 (default value) for the \code{JOINT} model, and equal to a strictly positive value for the R-JOINT model
#' @param lambda.reg A coefficient measuring the importance of the regularization term. In the related reference (2), it corresponds to the \code{R-JOINT} model for a value other than 0 (default value))
#' @param prox.dist A percentage (betwen 0 and 1) used to calculate the distance threshold below which an individual (a row) is considered as a neighbor of a given profile of covariates.
#' @param prox.X A percentage (betwen 0 and 1) used to calculate the distance threshold below which two covariates' profiles are supposed as neighbors.
#' If \code{prox.X = 1}, all profiles are considered as neighbors.
#' @param which.DB A character indicating the database to complete ("BOTH" by default, for the prediction of \eqn{Y} and \eqn{Z} in the two databases), "A" only for the imputation of \eqn{Z} in A, "B" only for the imputation of \eqn{Y} in B.
#'
#'
#' @return A list of 7 elements containing:
#'     \item{time_exe}{Running time of the function}
#'     \item{gamma_A}{Estimation of \eqn{\gamma} for the completion of A. A cost matrix that corresponds to the joint distribution of \eqn{(Y,Z,X)} in A}
#'     \item{gamma_B}{Estimation of \eqn{\gamma} for the completion of B. A cost matrix that corresponds to the joint distribution of \eqn{(Y,Z,X)} in B}
#'     \item{profile}{A data.frame that gives all details about the remaining P profiles of covariates. These informations can be linked to the \code{estimatorZA} and the \code{estimatorYB} objects for a better interpretation of the results}
#'     \item{res_prox}{The outputs of the function \code{proxim_dist}}
#'     \item{estimatorZA}{An array that corresponds to estimates of the probability distribution of \eqn{Z} conditional to \eqn{X} and \eqn{Y} in database A. The number of rows of each table corresponds to the total number of profiles of covariates.
#'     The number of columns of each table corresponds to the number of levels of \eqn{Y}. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of \eqn{Z}.}
#'     \item{estimatorYB}{An array that corresponds to estimates of the probability distribution of \eqn{Y} conditional to \eqn{X} and \eqn{Z} in database B. The number of rows of each table corresponds to the total number of profiles of covariates.
#'     The number of columns of each table corresponds to the number of levels of \eqn{Z}. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of \eqn{Y}.}
#'     \item{DATA1_OT}{database A with imputed individual prediction on \eqn{Z} using OT}
#'     \item{DATA2_OT}{database B with imputed individual prediction on \eqn{Y} using OT}
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
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N. (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679 | \url{https://doi.org/10.1515/ijb-2018-0106}
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association, DOI: 10.1080/01621459.2020.1775615
#' \item Anderberg, M.R. (1973), Cluster analysis for applications, 359 pp., Academic Press, New York, NY, USA.
#' \item Gower J.C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 27, 623–637.
#' }
#'
#' @seealso \code{\link{merge_dbs}}, \code{\link{OT_outcome}}, \code{\link{proxim_dist}}, \code{\link{avg_dist_closest}}
#'
#' @export
#'
#' @examples
#'
#' ### An example of JOINT model with:
#' #-----
#' # - A sample of the database tab_test
#' # - Y1 and Y2 are a 2 outcomes encoded in 2 different forms in DB 1 and 2:
#' #   4 levels for Y1 and 3 levels for Y2
#' # - n1 = n2 = 75
#' # - 2 discrete covariates X1 and X2 defined as ordinal
#' # - Distances estimated using the Manhattan function
#' # Predictions are assessed for Y1 in B only
#' #-----
#'
#' data(tab_test)
#' tab_test2 = tab_test[c(1:75,5001:5075),1:5]
#'
#' try1J = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
#'                  dist.choice = "M", which.DB = "B")
#'
#' \dontrun{
#'
#' ### An example of R-JOINT model using the previous database,
#' ### and keeping the same options excepted for:
#' #-----
#' # - The distances are estimated using the Gower function
#' # - Inclusion of an error term in the constraints on
#'     the marginals (relaxation term)
#' # Predictions are assessed for Y1 AND Y2 in A and B respectively
#' #-----
#'
#' try1RJ = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
#'                   dist.choice = "G", maxrelax = 0.4,
#'                   which.DB = "BOTH")
#'
#' ### The same example of R-JOINT model as previously:
#' # - Adding a regularization term
#' # Predictions are assessed for Y1 AND Y2 in A and B respectively
#' #-----
#'
#' try2RJ = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
#'                   dist.choice = "G", maxrelax = 0.4, lambda.reg = 0.9,
#'                   which.DB = "BOTH")
#'
#'
#' ### Another example of JOINT model with:
#' #-----
#' # - A sample of the database simu_data
#' # - Y1 and Y2 are a 2 outcomes encoded in 2 different forms in DB A and B:
#' #   (3 levels for Y and 5 levels for Z)
#' # - n1 = n2 = 100
#' # - 2 categorical covariates: Gender and Dosage
#' # - 1 continuous covariate (Age) converted in ordered factor of 3 balanced
#' #   classes
#' # - Distances estimated using the Euclidean function
#' # Predictions are assessed for Y1 in B only
#' #-----
#'
#' simu_data2 = simu_data[c(1:100,401:500),c(1:4,6,8)]
#'
#' try2J = OT_joint(simu_data2, convert.num = 6, convert.clss = 3,
#'                 nominal = c(1,4), ordinal = c(2:3,5),
#'                 dist.choice = "E", percent.knn = 0.90,
#'                 which.DB = "A")
#'
#'
#' ### A last example of JOINT model with:
#' #-----
#' # - Another sample of the database simu_data using
#' # - n1 = n2 = 100
#' # - 3 covariates: Gender, Smoking and Age in a qualitative form
#' # - Complete Case study
#' # - The Hamming distance
#' # Predictions are assessed for Y1 AND Y2 in A and B respectively
#' #-----
#'
#' simu_data2 = simu_data[c(1:100,401:500),c(1:4,7:8)]
#' simu_data3 = simu_data2[!is.na(simu_data2$Age),]
#'
#' try3J = OT_joint(simu_data3, convert.num = 6, convert.clss = 3,
#'                  nominal = c(1,4:5), ordinal = 2:3,
#'                  dist.choice = "H", which.DB = "BOTH")
#'
#' }

OT_joint = function(datab, index_DB_Y_Z = 1:3,
                    nominal = NULL, ordinal = NULL,logic = NULL,
                    convert.num = NULL, convert.clss = NULL, dist.choice = "E", percent.knn = 1,
                    maxrelax = 0, lambda.reg = 0.0, prox.dist = 0.80, prox.X = 0.80, which.DB = "BOTH"){

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

  if (length(convert.num) < length(convert.clss)){

    stop("The arguments quanti and convert.num must be equal because the algorithm does not handle continuous covariates")

  } else {}

  cat("---------------------------------------","\n")
  cat("OT JOINT PROCEDURE in progress ..."     ,"\n")
  cat("---------------------------------------","\n")
  cat("Type                  = ", ifelse((maxrelax == 0)&(lambda.reg == 0),"JOINT","R-JOINT"),"\n")
  cat("Distance              = ", ifelse(dist.choice == "H","Hamming",ifelse(dist.choice == "M","Manhattan",ifelse(dist.choice == "E","Euclidean","Gower"))),"\n")
  cat("Percent closest       = ", 100.0*percent.knn, "%","\n")
  cat("Relaxation term       = ", ifelse(maxrelax == 0,"NO","YES"),"\n")
  cat("Regularization term   = ", lambda.reg  ,"\n")
  cat("Aggregation tol cov   = ", prox.X,"\n")
  cat("DB imputed            = ", which.DB,"\n")
  cat("---------------------------------------","\n")


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


  inst = proxim_dist(dataB, norm = dist.choice, prox = prox.dist)


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

  indXA = inst$indXA; indXB = inst$indXB;
  nbX = length(indXA);

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

   tol.X = prox.X * max(dist_X)
   voisins_X = dist_X <= tol.X
   C = avg_dist_closest(inst, percent_closest = percent.knn)[1];

  ###########################################################################
  # Compute the estimators that appear in the model
  ###########################################################################

      estim_XA = estim_XB = estim_XA_YA =  estim_XB_ZB = list()

      for (x in 1:nbX){

        estim_XA[[x]] = length(indXA[[x]])/nA
        estim_XB[[x]] = length(indXB[[x]])/nB

      }

      for (x in 1:nbX){

        estim_XA_YA[[x]] = estim_XB_ZB[[x]] = numeric(0)


        for (y in Y){
          estim_XA_YA[[x]][y] = length(indXA[[x]][Yobserv[indXA[[x]]] == y])/nA
        }

        for (z in Z){
          estim_XB_ZB[[x]][z] = length(indXB[[x]][Zobserv[indXB[[x]] + nA] == z])/nB
        }

      }


      Cf <- function(y,z) {
        return(C$Davg[y,z])
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
          ompr::add_variable(gammaA[x,y,z]       , x = 1:nbX, y = Y, z= Z, type = "continuous") %>%
          ompr::add_variable(errorA_XY[x,y]      , x = 1:nbX, y = Y,       type = "continuous") %>%
          ompr::add_variable(abserrorA_XY[x,y]   , x = 1:nbX, y = Y,       type = "continuous") %>%
          ompr::add_variable(errorA_XZ[x,z]      , x = 1:nbX, z = Z,       type = "continuous") %>%
          ompr::add_variable(abserrorA_XZ[x,z]   , x = 1:nbX, z = Z,       type = "continuous") %>%
          # REGULARIZATION ---------------------------------------------------------------------------------
          # ompr::add_variable(reg_absA[x1,x2,y,z] , x1 = 1:nbX, x2 = 1:nbX, y= Y, z= Z, type = "continuous") %>%
          ompr::add_variable(reg_absA[x1,x2,y,z]   , x1 = 1:nbX, x2 = ind_voisins[[x1]], y= Y, z= Z, type = "continuous") %>%
          # OBJECTIVE ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
          ompr::set_objective(sum_expr(Cf(y,z) * gammaA[x,y,z], y = Y, z = Z, x = 1:nbX) + sum_expr(voisin(x1)*reg_absA[x1,x2,y,z], x1 = 1:nbX, x2 = ind_voisins[[x1]],y=Y,z= Z), "min") %>%
          # CONSTRAINTS ----------------------------------------------------------------------------------------------------
          ompr::add_constraint(sum_expr(gammaA[x,y,z], z = Z) - errorA_XY[x,y] == estim_XA_YA[[x]][y] , x = 1:nbX,y =Y) %>%
          ompr::add_constraint(estim_XBf(x)*sum_expr(gammaA[x,y,z],y = Y)  - estim_XBf(x)*errorA_XZ[x,z]== estim_XB_ZB[[x]][z] * estim_XA[[x]] , x = 1:nbX, z = Z) %>%

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
          ompr::add_constraint(reg_absA[x1,x2,y,z] + gammaA[x2,y,z]/(max(1,length(indXA[[x2]]))/nA) >= gammaA[x1,y,z]/(max(1,length(indXA[[x1]]))/nA), x1 = 1:nbX, x2 = which(voisins_X[x1,]), y= Y, z = Z) %>%
          ompr::add_constraint(reg_absA[x1,x2,y,z] + gammaA[x1,y,z]/(max(1,length(indXA[[x1]]))/nA) >= gammaA[x2,y,z]/(max(1,length(indXA[[x2]]))/nA), x1 = 1:nbX, x2 = which(voisins_X[x1,]), y= Y, z = Z) %>%

          # SOLUTION -------------------------------------------------------
        ompr::solve_model(with_ROI(solver = "glpk"))
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

        predZA = numeric(0)
        for (i in A){

          predZA = c(predZA,which.max(probaZindivA[i,]))

        }

        DATA1_OT         = dataB[dataB[,1] == unique(dataB[,1])[1],]
        # DATA1_OT$OTpred  = as.factor(plyr::mapvalues(predZA,from = sort(unique(predZA)), to = levels(dataB[,3])[sort(unique(predZA))]))


        if (index_DB_Y_Z[3] %in% nominal){

          DATA1_OT$OTpred  = factor(plyr::mapvalues(predZA,from = sort(unique(predZA)),
                                                    to = levels(dataB[,3])[sort(unique(predZA))]),
                                                    levels = levels(dataB[,3])[sort(unique(predZA))])
        } else {

          DATA1_OT$OTpred  = ordered(plyr::mapvalues(predZA,from = sort(unique(predZA)),
                                                     to = levels(dataB[,3])[sort(unique(predZA))]),
                                                     levels = levels(dataB[,3])[sort(unique(predZA))])

        }



        #if (is.ordered(dataB[,3])){

        #  DATA1_OT$OTpred = as.ordered(DATA1_OT$OTpred)

        #} else {}

      } else {}


      if (which.DB %in% c("B","BOTH")){

        # COMPLETE Y IN DATABASE B

        result <-  ompr::MIPModel() %>%
          # DEFINE VARIABLES ----------------------------------------------------------------------
        # gammaA[x,y,z]: joint probability of X=x, Y=y and Z=z in base B

          ompr::add_variable(gammaB[x,y,z]    , x = 1:nbX, y = Y, z = Z,type = "continuous") %>%

          ompr::add_variable(errorB_XY[x,y]   , x = 1:nbX, y = Y,       type = "continuous") %>%
          ompr::add_variable(abserrorB_XY[x,y], x = 1:nbX, y = Y,       type = "continuous") %>%
          ompr::add_variable(errorB_XZ[x,z]   , x = 1:nbX,        z = Z,type = "continuous") %>%
          ompr::add_variable(abserrorB_XZ[x,z], x = 1:nbX,        z = Z,type = "continuous") %>%

        # REGULARIZATION ----------------------------------------------------------------------------

          # ompr::add_variable(reg_absB[x1, x2,y,z], x1 = 1:nbX, x2 = 1:nbX, y = Y, z = Z, type = "continuous") %>%
          ompr::add_variable(reg_absB[x1,x2,y,z], x1 = 1:nbX, x2 = ind_voisins[[x1]], y = Y, z = Z, type = "continuous") %>%

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

          # SOLUTION -------------------------------------------------------
        ompr::solve_model(with_ROI(solver = "glpk"))

        solution  = get_solution(result, gammaB[x,y,z])
        gammaB_val= array(solution$value,dim = c(nbX,length(Y),length(Z)))
        #------------ END OPTIMIZATION STEP ------------------------------


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

        predYB = numeric(0)

        for (j in B){

          predYB = c(predYB,which.max(probaYindivB[j,]))

        }

        # Display the solution
        # println("Solution of the joint probability transport");
        # println("Distance cost = ", sum(C[y,z] * (gammaA_val[x,y,z]+gammaB_val[x,y,z]) for y in Y, z in Z, x in 1:nbX));
        # println("Regularization cost = ", lambda_reg * value(regterm));


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



        # DATA2_OT$OTpred  = as.factor(plyr::mapvalues(predYB,from = sort(unique(predYB)), to = levels(dataB[,2])[sort(unique(predYB))]))

        # if (is.ordered(dataB[,2])){

        #  DATA2_OT$OTpred = as.ordered(DATA2_OT$OTpred)

        #} else {}

      } else {}

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

        GAMMA_A = apply(gammaA_val,c(2,3),sum)
        GAMMA_B = apply(gammaB_val,c(2,3),sum)
        colnames(GAMMA_A)  = colnames(GAMMA_B)  = levels(dataB[,3])
        row.names(GAMMA_A) = row.names(GAMMA_B) = levels(dataB[,2])

      }

      tend = Sys.time()

      return(list(time_exe    = difftime(tend,tstart),
                  gamma_A     = GAMMA_A,
                  gamma_B     = GAMMA_B,
                  profile     = data.frame(ID = ID_prof,prof),
                  res_prox    = inst,
                  estimatorZA = estimatorZA,
                  estimatorYB = estimatorYB,
                  DATA1_OT    = DATA1_OT,
                  DATA2_OT    = DATA2_OT))
}
