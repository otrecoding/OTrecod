#' OT_joint()
#'
#' The function \code{OT_joint} integrates an original algorithm (\code{JOINT}) dedicated to the solving of recoding problems in data integration
#' using optimal transportation of the joint distribution of outcomes and covariates.
#'
#'
#' A. THE RECODING PROBLEM IN DATA FUSION
#'
#' Assuming that Y and Z are two variables that summarize a same latent information in two separate (no overlapping rows) databases A and B respectively,
#' so that Y and Z are never jointly observed in A and B. Assuming also that A and B share a subset of common covariates X of any types (but same encodings in A and B)
#' completed or not. Integrating these two databases often requires to solve the recoding problem observed between Y and Z by creating an unique database where
#' the missing information of Y and Z is completed.
#'
#'
#' B. INFORMATIONS ABOUT THE ALGORITHM
#'
#' As with the function \code{\link{OT_outcome}}, the function \code{OT_joint} provides a solution to the recoding problem previously described by considering the recoding as an
#' application of optimal transportation which aims at searching for a bijective mapping between the conditional distributions of (Y|X) and (Z|X) in A and B (See the 2nd referenced article for more details).
#' The principle of the algorithm is also based on the resolution of an optimization problem, provides a gamma solution (as called in the related articles) which is now an estimation
#' of the distribution of (X,Y,Z) according to the database to complete (see the argument \code{which.DB}). While the algorithms \code{OUTCOME} and \code{R_OUTCOME} integrated in
#' the function \code{\link{OT_outcome}} require post-treatment steps to provide individual predictions, the algorithm \code{JOINT} direcly use estimations of the conditional distributions (Y|Z,X) in B and
#' (Z|Y,X) in A to predict the corresponding incomplete individuals informations of Y,Z or both of them.
#' Estimations a posteriori of conditional probabilities \code{P[Y|X,Z]} and \code{P[Z|X,Y]} are available by profiles of covariates in output (See the objects \code{estimatorYB} and \code{estimatorZA}).
#' Estimations of gamma are also available according to the direction of the distributions transport chosen (See \code{TRANSPORT_A} and \code{TRANSPORT_B}).
#'
#' Enrichments of the algorithm \code{JOINT} is available with the function \code{OT_joint} by allowing users to add a relaxation term in the algorithm to relax distributional assumptions (\code{mawrelax}>0),
#' and (or) add also a positive regularization term (\code{lamdba.reg}>0) expressing that the transportation map does not vary to quickly with respect of X.
#' These enrichments corresponds to an extension of the raw algorithm called \code{R_joint} in the 2nd article referenced.
#' Is is suggested to users to calibrate these 2 parameters a posteriori by studying the stability of the individual predictions in output.
#'
#'
#' C. TYPE OF INPUT DATABASE REQUIRED
#'
#' The input database is a data.frame that must be saved in a specific form by users:
#' \itemize{
#' \item Two overlayed databases containing a common column of databases' identifiers (A and B, 1 or 2, by examples, encoded in numeric or factor form)
#' \item A column corresponding to the target variable with its specific encoding in A (By example a factor Y encoded in nY levels, ordered or not, with NAs in the corresponding rows of B)
#' \item A column corresponding to another target outcome summarizing the same latent information with its specific endoded in B (By example a factor Z in nZ levels, with NAs in rows of A)
#' \item The order of the variables in the database have no importance but the indexes of the columns related to the 3rd columns previously described (ie ID, Y and Z) must be rigorously specified
#' in the argument \code{index_DB_Y_Z}.
#' \item A set of shared common categorical covariates (at least one but more is recommended) complete or not (provided that the number of covariates exceeds 1) is required. On the contrary to the
#' function \code{OT_outcome}, please notice, that the function{OT_joint} does not accept continuous covariates therefore these latters will have to be categorized beforehand or using the input process provided (see \code{quanti}).
#' }
#' The function \code{\link{merge_dbs}} is available in this package to assist user in the preparation of their databases, so please, do not hesitate to use it beforehand if necessary.
#'
#' Remarks about the outcomes:
#' \itemize{
#' \item A target outcome can be categorical, in factor, ordered or not, discrete (with a finite number of values ONLY) but, notice that, if they are stored in numeric they will be automatically converted in ordered factors.
#' \item If a target outcome is incomplete, the corresponding rows will be automatically dropped during the execution of the function.
#' }
#' The type of each variables (including ID, Y and Z) of the database must be rigorously specified once, in one of the 4 arguments \code{quanti},\code{nominal}, \code{ordinal} and \code{logic}.
#'
#'
#' D. TRANSFORMATIONS OF CONTINUOUS COVARIATES
#'
#' Continuous covariates with infinite numbers of values will have to be categorized beforehand their inclusions in the function.
#' To assist users in this task, the function \code{OT_joint} integrates in is syntax a process dedicated to the categorization of continuous covariates. For this, it is necessary to rigorously filled in
#' the arguments \code{quanti} and \code{convert.clss}.
#' The first one informs about the indexes in database of the continuous variables to transform in ordered factor while the second one specifies the corresponding number of balanced levels desired (for unbalanced levels, users must do transformations by themselves).
#' Therefore \code{quanti} and \code{convert.clss} must be vectors of same length, but if the length of \code{quanti} exceeds 1, while the length of \code{convert.clss} is 1, then, by default, all the covariates to convert will have the same number of classes,
#' that corresponds to the value specified in the argument \code{convert.clss}.
#' Please notice that only covariates can be transformed (not outcomes) and missing informations are not taken into account for the transformations.
#' Moreover, all the indexes informed in the argument \code{convert.num} must also be informed in the argument \code{quanti}.
#' Finally, it is suggested to declare all discrete covariates as ordinal factors using the argument \code{ordinal}.
#'
#'
#' E. INFORMATIONS ABOUT DISTANCE FUNCTIONS AND RELATED PARAMETERS
#'
#' Each individual (or row) of a given database is here characterized by his covariates, so the distance between 2 individuals or groups of individuals depends on similarities between covariates
#' according to the distance function chosen by user. Users can specify the distance measure chosen via the argument \code{dist.choice}. Actually 4 distance functions are implemented in {OT_outcome} to take into account the most frequently encountered situation (see Anderberg 1973):
#' \itemize{
#' \item The Manhattan distance ("M")
#' \item The Euclidean distance ("E")
#' \item The Gower distance for mixed data (See Gower(1971):"G")
#' \item The Hamming distance for binary data ("H")
#' }
#'
#' Assuming that P1 and P2 are 2 profiles of covariates, they will be considered as neighbors if \eqn{dist(P1,P2) < \code{prox.X} * max(dist(Pi,Pj))} where \code{prox.X} must be also fixed by user. This choice comes into the calculation of the \code{JOINT} and \code{R_joint} algorithms.
#' In the same way, for a given profile of covariates \code{P_j}, an individual i will be considered as a neighbor of \code{P_j} if \eqn{dist(i,P_j) < \code{prox.dist} * max(dist(i,P_j))} where \code{prox.dist} will be fixed by user.
#'
#' For more details about the related algorithms integrated in \code{OT_joint}, please consult, the 2nd article referenced.
#'
#'
#' @param datab A data.frame that must have at least 4 columns sorted in a non-specific order. One column must be a key column of 2 classes for the databases identification, where the names of the two databases
#' must be alphanumerically ranked in ascending order (By examples: 1 for the top database and 2 for the database from below, or more logically here A and B  ...But NOT B and A!).One column (Y here but other names are allowed)
#' must correspond to the target variable related to the information of interest to merge with its specific encoding in the database A (corresponding encoding should be so missing in the database B). In the same way,
#' one column (Z here) corresponds to the target variable that summarizes the same information as Y but with its specific encoding in the database B (corresponding encoding should be so missing in the database A).
#' Finally, your database must have at least one shared covariate with same encoding in A and B. Please not that, if your data.frame has only 4 columns, that is to say, only one covariate, if this latter has NA, and unless user
#' has previously imputed the corresponding missing information, the OT algorithm will only run with complete cases.
#' @param index_DB_Y_Z A vector of 3 indexes of columns. The 1st index must correspond to the index of the databases identification column (DB identifier). The 2nd index corresponds
#' to the index of the target variable in the 1st database (A) while the 3rd index corresponds to the index of column related to the target variable in the 2nd database (B).
#' @param nominal A vector of indexes that corresponds to the indexes of columns of all the nominal (not ordered) variables (DB identification and target variables included if it is the case for them).
#' @param ordinal A vector of indexes that corresponds to the indexes of columns of all the ordinal variables (DB identification and target variables included if it is the case for them).
#' @param logic A vector of indexes that corresponds to the indexes of columns of all the boolean variables of the data.frame.
#' @param convert.num Indexes of the continuous (quantitative) variables. They will be automatically converted in ordered factors. By default, no continuous variables is assumed in the database.
#' @param convert.clss A vector indicating for each continuous variable to convert, the corresponding number of levels desired.If the length of the argument \code{convert_num} exceeds 1 while the length of \code{convert_clss} equals 1 (only one integer),
#' each discretization will count the same number of levels.
#' @param dist.choice A character (with quotes) corresponding to the distance function chosen between: The euclidean distance ("E", by default), The Manhattan distance ("M"),
#' the Gower distance ("G"), and the Hamming distance ("H") for binaries covariates only.
#' @param percent.knn Percent of closest neighbors taken in the computation of the cost matrix.
#' @param maxrelax Maximum percentage of deviation from expected probability masses. It must be equal to 0 (default value) for the JOINT model, and equal to a strictly positive value for the R-JOINT model
#' @param lambda.reg A coefficient measuring the importance of the regularization term. In the related reference, it corresponds to the R-JOINT model for a vue other than 0 (Default value))
#' @param prox.dist A percentage (betwen 0 and 1) uses to calculate the distance threshold below which an individual (a row) is considered as a neighbor of a given profile of covariates.
#' @param prox.X A percentage (betwen 0 and 1) uses to calculate the distance threshold below which 2 covariates' profiles are supposed as neighbors.
#' If \code{prox.X = 1}, all profiles are considered as neighbors.
#' @param which.DB A character indicating which database completed ("BOTH" by default), "A" for the imputation of Z in A, "B" for the imputation of Y in B.
#'
#'
#' @return A list of 7 elements containing:
#'     \item{time_exe}{Running time of the function}
#'     \item{gamma_A}{Estimation of gamma for the completion of A. A cost matrix that corresponds to the joint distribution of (YA,ZA)}
#'     \item{gamma_B}{Estimation of gamma for the completion of B. A cost matrix that corresponds to the joint distribution of (YB,ZB)}
#'     \item{profile}{A data.frame that gives all details about the remaining P profiles of covariates. These informations can be linked to the \code{estimatorZA} and the \code{estimatorYB} objects for a better interpretation of the results}
#'     \item{res_prox}{The outputs of the function \code{proxim_dist}}
#'     \item{estimatorZA}{Estimates of the probability distribution of Z conditional to X and Y in database A from predictions stored in an array. The number of rows of each table corresponds to the total number of profiles of covariates.
#'     The number of columns of each table corresponds to the number of levels of Y. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of Z.}
#'     \item{estimatorYB}{Estimates of the probability distribution of Y conditional to X and Z in database B from predictions stored in an array. The number of rows of each table corresponds to the total number of profiles of covariates.
#'     The number of columns of each table corresponds to the number of levels of Z. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of Y.}
#'     \item{DATA1_OT}{database A with imputed individual prediction on Z using OT}
#'     \item{DATA2_OT}{database B with imputed individual prediction on Y using OT}
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
#' \email{otrecod.pkg@@gmail.com}
#'
#' @aliases OT_joint
#'
#' @references
#' \enumerate{
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679 | \url{https://doi.org/10.1515/ijb-2018-0106}
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association, DOI: 10.1080/01621459.2020.1775615
#' }
#' # For the Gower distance:
#' Gower J. C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 27, 623–637.
#'
#' # About the other distance functions:
#' Anderberg, M.R. (1973), Cluster analysis for applications, 359 pp., Academic Press, New York, NY, USA.
#'
#'
#' @seealso \code{\link{OT_outcome}},\code{\link{proxim_dist}}, \code{\link{avg_dist_closest}}
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
#' # - The inclusion if of a term of error in the constraints on
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
          ompr::add_variable(reg_absA[x1,x2,y,z] , x1= 1:nbX, x2= 1:nbX, y= Y, z= Z, type = "continuous") %>%
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
          ompr::add_constraint(sum_expr(errorA_XZ[x,z]   , x = 1:nbX, z = Z)<= maxrelax/2.0) %>%

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
        DATA1_OT$OTpred  = as.factor(plyr::mapvalues(predZA,from = sort(unique(predZA)), to = levels(dataB[,3])[sort(unique(predZA))]))

        if (is.ordered(dataB[,3])){

          DATA1_OT$OTpred = as.ordered(DATA1_OT$OTpred)

        } else {}

      } else {}


      if (which.DB %in% c("B","BOTH")){

        # COMPLETE Y IN DATABASE B

        result <-  ompr::MIPModel() %>%
          # DEFINE VARIABLES ----------------------------------------------------------------------
        # gammaA[x,y,z]: joint probability of X=x, Y=y and Z=z in base A
        ompr::add_variable(gammaB[x,y,z]    , x = 1:nbX, y = Y, z= Z   ,type = "continuous") %>%
          ompr::add_variable(errorB_XY[x,y]   , x = 1:nbX, y = Y,       type = "continuous") %>%
          ompr::add_variable(abserrorB_XY[x,y], x = 1:nbX, y = Y,       type = "continuous") %>%
          ompr::add_variable(errorB_XZ[x,z]   , x = 1:nbX,        z = Z,type = "continuous") %>%
          ompr::add_variable(abserrorB_XZ[x,z], x = 1:nbX,        z = Z,type = "continuous") %>%
          # REGULARIZATION ----------------------------------------------------------------------------
        ompr::add_variable(reg_absB[x1, x2,y,z],x1=1:nbX, x2=1:nbX,y=Y,z= Z, type = "continuous") %>%
          # OBJECTIVE ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        ompr::set_objective(sum_expr(Cf(y,z) * gammaB[x,y,z],y = Y,z=Z, x = 1:nbX)  + sum_expr(voisin(x1)*reg_absB[x1,x2,y,z],x1=1:nbX, x2= ind_voisins[[x1]],y=Y,z= Z), "min")  %>%
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
          ompr::add_constraint(sum_expr(errorB_XZ[x,z]   , x = 1:nbX, z = Z)<= maxrelax/2.0) %>%

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
        DATA2_OT$OTpred  = as.factor(plyr::mapvalues(predYB,from = sort(unique(predYB)), to = levels(dataB[,2])[sort(unique(predYB))]))

        if (is.ordered(dataB[,2])){

          DATA2_OT$OTpred = as.ordered(DATA2_OT$OTpred)

        } else {}

      } else {}

      if (which.DB == "A"){

        GAMMA_B            = NULL
        estimatorYB        = NULL
        DATA2_OT           = NULL
        GAMMA_A            = apply(gammaA_val,c(2,3),sum)
        colnames(GAMMA_A)  = levels(dataB[,3])
        row.names(GAMMA_A) = levels(dataB[,2])

      } else if (which.DB == "B"){

        GAMMA_A            = NULL
        estimatorZA        = NULL
        DATA1_OT           = NULL
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
