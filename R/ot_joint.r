#' OT_joint()
#'
#' Implementation of the Optimal Transportation (OT) algorithm for data integration by directly estimating the joint distributions (Y,Z,X) in the two databases
#'
#' The context of use of this function is the same as for the \code{\link{OT}} function.Nevertheless, the algorithm related to the \code{OT_joint} function does not allowed the use of continuous quantitative covariates.
#' Bbe sure that continuous quantitatives covariates have been deleted from your analyses or discretize beforehand using this function.
#'
#' Assuming that:
#'
#' \enumerate{
#' \item Y and Z summarize a same information of interest encoded in two distinct forms stored in two independent databases A and B (no individuals or rows in common),
#' \item The two databases have a set of common covariates X (i.e the same variables stored in the same forms in the two bases),
#' \item You would like to merge vertically A and B keeping this variable of interest in at least one of its two encodings.
#' }
#'
#' So, this function gives you a possible answer to this recoding problem by using the Optimal Transportation Theory.
#' The models implemented in this function correspond to the JOINT and R-JOINT models described in the referenced article 2.
#'
#' The algorithm linked to this function transports the joint distributions of target variables AND covariates (X,Y,Z) while the \code{\link{OT}} function only derives the joint distribution of Y and Z.
#' By direclty solving the recoding problem (given the individual predictions of Z in A and Y in B), an independent post-treatment step like the nearest neighbor procedure implemented in the \code{\link{OT}} function is no longer required,
#' and the individual predictions of Z in A and Y in B are given at once.
#'
#' As for the \code{\link{OT}} function, this function allows to relax the constraints on marginal distributions (using the \code{maxrelax} option) and gives the possibility to add an L^1 regularization term (using the \code{lambda_reg} option) expressing that the transportation map should not vary too quickly
#' with respect to X.
#'
#' For more details, please consult the reference article 2.
#'
#' @param datab A data.frame that must have at least 4 columns sorted in a non-specific order. One column must be a key column of 2 classes for the databases identification, where the names of the two databases
#' must be alphanumerically ranked in ascending order (By examples: 1 for the top database and 2 for the database from below, or more logically here A and B  ...But NOT B and A!).One column (Y here but other names are allowed)
#' must correspond to the target variable related to the information of interest to merge with its specific encoding in the database A (corresponding encoding should be so missing in the database B). In the same way,
#' one column (Z here) corresponds to the target variable that summarizes the same information as Y but with its specific encoding in the database B (corresponding encoding should be so missing in the database A).
#' Finally, your database must have at least one covariate with same encoding in A and B. Please not that, if your data.frame has only 4 columns, that is to say, only one covariate, if this latter has NA, and unless user
#' has previously imputed the corresponding missing information,the OT algorithm will only run with complete cases.
#' @param index_DB_Y_Z A vector of exactly 3 integers. The 1st integer must correspond to the index of the databases identification column (DB identifier). The 2nd integer corresponds
#' to the index of the target variable in the 1st database (A) while the 3rd integer corresponds to the index of column related to the target variable in the 2nd database (B).
#' @param nominal A vector of integers that corresponds to the indexes of columns of all the nominal (not ordered) variables (DB identification and target variables included if it is the case for them).
#' @param ordinal A vector of integers that corresponds to the indexes of columns of all the ordinal variables (DB identification and target variables included if it is the case for them).
#' @param logic A vector of integers that corresponds to the indexes of columns of all the boolean variables of the data.frame.
#' @param prep_choice A character (with quotes) corresponding to the distance function chosen between: The euclidean distance ("E", by default), The Manhattan distance ("M"),
#' the Gower distance ("G"), and the Hamming distance ("H") for binaries covariates only.
#' @param infoFAMD A percent value (between 0 and 100) that corresponds to the part of variability taken into account by the principal components of the FAMD when this option is required.
#' @param percent_c Percent of closest neighbors taken in the computation of the cost matrix.
#' @param maxrelax Maximum percentage of deviation from expected probability masses. It must be equal to 0 (default value) for the JOINT model, and equal to a strictly positive value for the R-JOINT model
#' @param lambda_reg A coefficient measuring the importance of the regularization term. In the related reference, it corresponds to the R-JOINT model for a vue other than 0 (Default value))
#' @param norm A decimal number to choose between 0 and 1. This value specifies the distance function chosen to compute the neighbors of the covariates for regularization. Equals to 0 for the Hamming distance, to 1 for the Manhattan distance and to 2 for the Euclidean distance respectively
#' @param prox_dist A value between 0 and 1 that corresponds to a threshold (0.1 by default) below which an individual is considered significantly close (neighbor) to a given profile of covariates.
#' @param aggregate_tol A value between 0 and 1 that quantify how much individuals' covariates must be closed for aggregation
#' @param which.DB A value equals to 1 if the user only wants to impute Z in the 1st database (A), equals to 2 if the user only wants to impute Y in the 2nd database (B) or equals 3 for imputation of Y in Z in the both databases (A and B)
#'
#'
#' @return A list of 7 elements containing:
#'     \item{TIME_EXE}{Running time of the function}
#'     \item{GAMMA_A}{Cost matrix corresponding to an estimation (gamma, see reference for more details) to the joint distribution of (YA,ZA)}
#'     \item{GAMMA_B}{Cost matrix corresponding to an estimation to the joint distribution of (YB,ZB)}
#'     \item{profile}{A data.frame that gives all details about the remaining P profiles of covariates. These informations can be linked to the \code{estimatorZA} and the \code{estimatorYB} objects for a better interpretation of the results}
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
#' \email{gregory.guernec@@inserm.fr}
#'
#' @aliases OT_joint OT_JOINT ot_joint
#'
#' @references
#' # Article 1:
#' Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' 0, 20180106 (2019) | \url{https://doi.org/10.1515/ijb-2018-0106}
#'
#' # Article 2:
#' Gares V, Omer J. Regularized optimal transport of covariates and outcomes in datarecoding(2019).hal-02123109 \url{https://hal.archives-ouvertes.fr/hal-02123109/document}
#'
#' @seealso \code{\link{OT}},\code{\link{proxim_dist}}, \code{\link{avg_dist_closest}}, \code{\link{indiv_grp_closest}}, \code{\link{indiv_grp_optimal}}
#'
#' @export
#'
#'@examples
#'
#' ### Using the \code{simu_data} object
#' ### Y and Z are a same variable encoded in 2 different forms in DB A and B:
#' ### (3 levels for Y and 5 levels for Z)
#' #--------
#' data(simu_data)
#'
#' ### using a sample of the tab_test object (3 complete covariates)
#' ### Y1 and Y2 are a same variable encoded in 2 different forms in DB 1 and 2:
#' ### (4 levels for Y1 and 3 levels for Y2)
#'
#' data(tab_test)
#' # Example with n1 = n2 = 75 and only X1 and X2 as covariates
#' tab_test2 = tab_test[c(1:75,5001:5075),1:5]
#'
#' ### An example of JOINT model (Manhattan distance)
#' # Suppose we want to impute the missing parts of Y1 in DB2 only ...
#' try1J = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
#'                  prep_choice = "M", norm = 1, which.DB = 2)
#'
#' \dontrun{
#' ### An example of R-JOINT model (Gower distance) by including a term of error
#' ### in the constraints on the marginals (relaxation)
#' # Suppose we want to impute the missing parts of Y2 in DB1 only ...
#' try2J = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3,6),
#'                  prep_choice = "G", maxrelax = 0.4, norm = 3, which.DB = 1)
#'
#' ### An example of R-JOINT model (Euclidean distance) with relaxation and regularization terms
#' try3J = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3,6),
#'                  prep_choice = "E", maxrelax = 0.4, lambda_reg = 0.9, norm = 2, which.DB = 3)
#'
#'
#' ### Using the \code{simu_data} object
#' ### Y and Z are a same variable encoded in 2 different forms in DB A and B:
#' ### (3 levels for Y and 5 levels for Z)
#' #--------
#' data(simu_data)
#'
#' # By only keeping 3 covariates: Gender, Smoking and Age in a qualitative form:
#' simu_data2 = simu_data[c(1:100,401:500),c(1:4,7:8)]
#'
#' # The variable Age has to be discretized previously:
#' simu_data2$Age_class = cut(simu_data2$Age,breaks = c(34,47,53,66))
#' simu_data2           = simu_data2[,-(ncol(simu_data2)-1)]
#'
#' try4J      = OT_joint(simu_data2, nominal = c(1,4:5), ordinal = c(2:3,6),
#'                       prep_choice = "E", norm = 2, which.DB = 1)
#'
#' # Complete Case study and only one remaining covariate:
#' # (Manhattan distance, DB 2 only)
#' simu_data3 = na.omit(simu_data2[,c(1:3,6)])
#' try5J      = OT_joint(simu_data3, nominal = 1, ordinal = 2:4,
#'                       prep_choice = "M", norm = 2, which.DB = 2)
#'
#' # (Hamming distance, DB2 only)
#' try6J      = OT_joint(simu_data3, nominal = 1, ordinal = 2:4,
#'                       prep_choice = "H", norm = 0, which.DB = 2)
#'
#' }

OT_joint = function(datab, index_DB_Y_Z = 1:3, nominal = NULL, ordinal = NULL,logic = NULL, prep_choice = "E",infoFAMD = 80, percent_c = 1.0, maxrelax = 0.0,
                    lambda_reg = 0.0, norm = 1, prox_dist = 1, aggregate_tol = 1, which.DB = 3){

  if (prep_choice =="FAMD"){

    stop("The FAMD option can not be used here, because continuous quantitative covariates are not allowed here")

  } else {}


  cat("---------------------------------------","\n")
  cat("OT JOINT PROCEDURE in progress ..."     ,"\n")
  cat("---------------------------------------","\n")
  cat("Type                  = ", ifelse((maxrelax == 0)&(lambda_reg == 0),"JOINT","R-JOINT"),"\n")
  cat("Distance              = ", ifelse(norm == 0,"Hamming",ifelse(norm == 1,"Manhattan",ifelse(norm == 2,"Euclidean","Gower"))),"\n")
  cat("Percent closest       = ", 100.0*percent_c, "%","\n")
  cat("Relaxation term       = ", ifelse(maxrelax == 0,"NO","YES"),"\n")
  cat("Regularization weight = ", lambda_reg  ,"\n")
  cat("Aggregation tolerance = ", aggregate_tol,"\n")
  cat("IMPUTE ?              = ", ifelse(which.DB == 1,"1st DB",ifelse(which.DB == 2,"2nd DB","Both DBs")),"\n")
  cat("---------------------------------------","\n")

  if (!(which.DB %in% 1:3)){

    stop("Invalid value for the which.DB option")

  } else {}

  if (infoFAMD > 100){

    stop("Invalid value for the infoFAMD option")

  } else {}

  if (percent_c > 100){

    stop("Invalid value for the percent related to the definition of a significant neighborood")

  } else {}


  tstart = Sys.time()


  dataB = transfo_dist(datab,index_DB_Y_Z = index_DB_Y_Z,quanti = NULL, nominal = nominal, ordinal = ordinal,
                       logic = logic, prep_choice = prep_choice, info = infoFAMD)


  inst = proxim_dist(dataB, norme = norm, prox = prox_dist)


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
  # norme = 1

  if (norm == 1){

    dist_X  = proxy::dist(Xvalues,Xvalues, method = "manhattan")

  } else if (norm == 2){

    dist_X  = proxy::dist(Xvalues,Xvalues, method = "euclidean")

  } else if (norm == 0){

    if (nrow(dataB[,4:ncol(dataB)])== nrow(na.omit(dataB[,4:ncol(dataB)]))){

      dist_X  = rdist::cdist(Xvalues,Xvalues, metric = "hamming")


    } else {

      dist_X  = ham(Xvalues,Xvalues)

    }

  } else if (norm == 3){

    dist_X  = StatMatch::gower.dist(Xvalues,Xvalues)

  }


  voisins_X = dist_X <= aggregate_tol

  # println("... computing costs")
  C = avg_dist_closest(inst, percent_closest = percent_c)[1];

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
    C$Davg[y,z]
  }

  estim_XBf <- function(x) {
    estim_XB[[x]]
  }

  voisin = function(x1){lambda_reg *(1/length(voisins_X[x1,]))}
  ind_voisins = list()

  for (x1 in 1:nrow(voisins_X)){
    ind_voisins[[x1]] = which(voisins_X[x1,])
  }
  ###########################################################################
  # Basic part of the model
  ###########################################################################

  if (which.DB %in% c(1,3)){

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
    ompr::add_variable(reg_absA[x1,x2,y,z] , x1= 1:nbX, x2= 1:nbX, y=Y, z= Z, type = "continuous") %>%
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


  if (which.DB %in% c(2,3)){

  # COMPLETE Y IN DATABASE B

  result <-  ompr::MIPModel() %>%
    # DEFINE VARIABLES ----------------------------------------------------------------------
    # gammaA[x,y,z]: joint probability of X=x, Y=y and Z=z in base A
    ompr::add_variable(gammaB[x,y,z]    , x = 1:nbX, y = Y, z= Z, type = "continuous") %>%
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

  if (which.DB == 1){

    GAMMA_B            = NULL
    estimatorYB        = NULL
    DATA2_OT           = NULL
    GAMMA_A            = apply(gammaA_val,c(2,3),sum)
    colnames(GAMMA_A)  = levels(dataB[,3])
    row.names(GAMMA_A) = levels(dataB[,2])

  } else if (which.DB == 2){

    GAMMA_A            = NULL
    estimatorZA        = NULL
    DATA1_OT           = NULL
    GAMMA_B            = apply(gammaB_val,c(2,3),sum)
    colnames(GAMMA_B)  = levels(dataB[,3])
    row.names(GAMMA_B) = levels(dataB[,2])

  } else {

    GAMMA_A = apply(gammaA_val,c(2,3),sum)
    GAMMA_B = apply(gammaB_val,c(2,3),sum)
    colnames(GAMMA_A) = colnames(GAMMA_B)   = levels(dataB[,3])
    row.names(GAMMA_A) = row.names(GAMMA_B) = levels(dataB[,2])

  }

  tend = Sys.time()

  return(list(TIME_EXE = difftime(tend,tstart),
              GAMMA_A = GAMMA_A,
              GAMMA_B = GAMMA_B,
              profile = data.frame(ID = ID_prof,prof),
              estimatorZA = estimatorZA,
              estimatorYB = estimatorYB,
              DATA1_OT = DATA1_OT,
              DATA2_OT = DATA2_OT))
}
