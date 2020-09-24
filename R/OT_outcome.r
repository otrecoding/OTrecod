#' OT_outcome()
#'
#' The function \code{OT_outcome} integrates two models called (\code{OUTCOME}) and \code{R-OUTCOME} dedicated to the solving of recoding problems in data integration
#' using optimal transportation of outcomes.
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
#' The algorithm integrated in the function \code{OT_outcome} provides a solution to the recoding problem previously described by proposing an
#' application of optimal transportation which aims is to search for a bijective mapping between the distributions of of \eqn{Y} in A and \eqn{Z} in B.
#' Mathematically, the principle of the algorithm is based on the resolution of an optimization problem which provides a \eqn{\gamma} solution (as called in the related articles)
#' that pushes the distribution of \eqn{Y} forward to the distribution of \eqn{Z} and can be so interpreted as an estimator of the joint distribution
#' \eqn{(Y,Z)} in A (or B respetively). According to this result, a second step of the algorithm provides individual predictions of \eqn{Y} in B (of \eqn{Z} in A, or both, depending on the choice
#' specified by user in the argument \code{which.DB}) two possible approaches are available depending on the argument \code{indiv.method}:
#' \itemize{
#' \item When \code{indiv.method = "sequential"}, a nearest neighbor procedure is applied. This corresponds to the use of the function \code{\link{indiv_grp_closest}} directly
#' implemented in the function \code{OT_outcome}.
#' \item When \code{indiv.method = "optimal"}, a linear optimization problem is solved to determine the individual predictions that minimize the sum of the individual distances
#' in A (B) with the modalities of \eqn{Z} in B (\eqn{Y} in A). This uses directly the function \code{\link{indiv_grp_optimal}} in the function \code{OT_outcome}.
#' }
#' This algorithm supposes the respect of the two following assumptions:
#' \enumerate{
#' \item \eqn{Y} must follow the same distribution in A and B. In the same way \eqn{Z} follows the same distribution in the two databases.
#' \item The conditional distribution \eqn{(Y|X)} must be identical in A and B. Respectively, \eqn{(Z|X)} is supposed identical in A and B.
#' }
#' Because the first assumption can appear too strong in some situations, a relaxation of the constraints of marginal distribution is possible using the argument \code{maxrelax}.
#' When \code{indiv.method = "sequential"} and \code{maxrelax = 0}, the first version of the algorithm (called \code{OUTCOME} in (2) and used in all the article (1))
#' is applied. In all other situations, the algorithm applied corresponds to a model called \code{R_OUTCOME} (only described in (2)).
#' Estimations a posteriori of conditional probabilities \eqn{P[Y|X,Z]} and \eqn{P[Z|X,Y]} are available by profiles of covariates in output (See the objects \code{estimatorYB} and \code{estimatorZA}).
#' Estimations of \eqn{\gamma} are also available according to the direction of the distributions transport chosen (See \eqn{\gamma_A} and \eqn{\gamma_B}).
#'
#'
#' C. TYPE OF REQUIRED INPUT DATABASE
#'
#' The input database is a data.frame that must be saved in a specific form by users:
#' \itemize{
#' \item Two overlayed databases containing a common column of database identifiers (A and B, 1 or 2, by examples, encoded in numeric or factor form)
#' \item A column corresponding to the target variable with its specific encoding in A (For example a factor \eqn{Y} encoded in \eqn{n_Y} levels, ordered or not, with NAs in the corresponding rows of B)
#' \item A column corresponding to the second target outcome with its specific endoded in B (For example a factor \eqn{Z} in \eqn{n_Z} levels, with NAs in rows of A)
#' \item The order of the variables in the database have no importance but the column indexes related to the three columns previously described (ie ID, \eqn{Y} and \eqn{Z}) must be rigorously specified
#' in the argument \code{index_DB_Y_Z}.
#' \item A set of shared common covariates (at least one but more is recommended) of any type, complete or not (provided that the number of covariates exceeds 1) is required.
#' }
#' The function \code{\link{merge_dbs}} is available in this package to assist user in the preparation of their databases, so please, do not hesitate to use it beforehand if necessary.
#'
#' Remarks about the outcomes:
#' \itemize{
#' \item A target outcome can be categorical, in factor, ordered or not, discrete (with a finite number of values ONLY) but, notice that, if they are stored in numeric they will be automatically converted in ordered factors.
#' \item If a target outcome is incomplete, the corresponding rows will be automatically dropped during the execution of the function.
#' }
#' The type of each variables (including \eqn{ID}, \eqn{Y} and \eqn{Z}) of the database must be rigorously specified once, in one of the four arguments \code{quanti},\code{nominal}, \code{ordinal} and \code{logic}.
#'
#'
#' D. TRANSFORMATIONS OF CONTINUOUS COVARIATES
#'
#' The function \code{OT_outcome} integrates in its syntax a process dedicated to the categorization of continuous covariates. For this, it is necessary to rigorously fill in the arguments \code{convert.num} and \code{convert.clss}.
#' The first one informs about the indexes in database of the continuous variables to transform in ordered factor while the second one specifies the corresponding number of desired balanced levels (for unbalanced levels, users must do transformations by themselves).
#' Therefore \code{convert.num} and \code{convert.clss} must be vectors of same length, but if the length of \code{convert.num} exceeds 1, while the length of \code{convert.clss} is 1, then, by default, all the covariates to convert will have the same number of classes,
#' that corresponds to the value specified in the argument \code{convert.clss}.
#' Please notice that only covariates can be transformed (not outcomes) and missing informations are not taken into account for the transformations.
#' Moreover, all the indexes informed in the argument \code{convert.num} must also be informed in the argument \code{quanti}.
#'
#'
#' E. INFORMATIONS ABOUT DISTANCE FUNCTIONS
#'
#' Each individual (or row) of a given database is here characterized by their covariates, so the distance between two individuals or groups of individuals depends on similarities between covariates
#' according to the distance function chosen by user (via the argument \code{dist.choice}). Actually four distance functions are implemented in \code{OT_outcome} to take into account the most frequently encountered situation (see (3)):
#' \itemize{
#' \item the Manhattan distance ("M")
#' \item the Euclidean distance ("E")
#' \item the Gower distance for mixed data (see (4): "G")
#' \item the Hamming distance for binary data ("H")
#' }
#' Moreover, it is also possible to directly apply the first three distances mentioned on coordinates extracted from a multivariate analysis (Factor Analysis for Mixed Data, see (5)) applied on raw covariates using the arguments \code{FAMD.coord} and \code{FAMD.perc}.
#' This method is used in the article (1).
#'
#' As a decision rule, for a given profile of covariates \eqn{P_j}, an individual i will be considered as a neighbor of \eqn{P_j} if \eqn{dist(i,P_j) < \mbox{prox.dist} \times max(dist(i,P_j))} where \eqn{prox.dist} must be fixed by user.
#'
#'
#' F. INFORMATIONS ABOUT THE SOLVER
#'
#' The argument \code{solvR} permits user to choose the solver of the optimization algorithm. The default solver is "glpk" that corresponds to the GNU Linear Programming Kit (see (6) for more details). The solver "clp" (see (7)) for Coin-or Linear Programming, convenient in linear and quadratic situations, is also directly integrated in the function.
#' Moreover, the function actually uses the \code{R} optimization infrastructure of the package \pkg{ROI} which offers a wide choice of solver to users by easily loading the associated plugins of \pkg{ROI} (see (8)).
#'
#'
#' @aliases OT_outcome ot_outcome OT
#'
#' @param datab A data.frame with at least four columns sorted in a non-specific order. One column must be a column of identifier for the two databases, where the names of the two databases
#' must be alphanumerically ranked in ascending order (for example: 1 for the top database and 2 for the database from below, or more logically here A and B  ...But NOT B and A!). One column (\eqn{Y} here but other names are allowed)
#' must correspond to the target variable related to the information of interest to merge with its specific encoding in the database A (corresponding encoding should be so missing in the database B). In the same way,
#' one column (\eqn{Z} here) corresponds to the second target variable with its specific encoding in the database B (corresponding encoding should be so missing in the database A).
#' Finally, the input database must have at least one shared covariate with same encoding in A and B. Please notice that, if your data.frame has only four columns, that is to say, only one covariate, if this latter has NA, and unless user
#' has previously imputed the corresponding missing information, the OT algorithm will only run complete cases.
#' @param index_DB_Y_Z A vector of three column indexes. The first index must correspond to the index of the databases identifier column (DB identifier). The second index must correspond
#' to the index of the target variable in the first database (A) while the third index corresponds to the index of column related to the target variable in the second database (B).
#' @param quanti A vector of indexes that corresponds to the column indexes of all the quantitative variables (DB identification and target variables included if it is the case for them).
#' @param nominal A vector of indexes that corresponds to the column indexes of all the nominal (not ordered) variables (DB identification and target variables included if it is the case for them).
#' @param ordinal A vector of indexes that corresponds to the column indexes of all the ordinal variables (DB identification and target variables included if it is the case for them).
#' @param logic A vector of indexes that corresponds to the column indexes of all the boolean variables.
#' @param convert.num Indexes of the continuous (quantitative) variables to convert in ordered factors. All declared indexes in this argument must have been declared in the argument \code{quanti} (no conversion by default).
#' @param convert.clss A vector indicating for each continuous variable to convert, the desired corresponding number of levels. If the length of the argument \code{convert_num} exceeds 1 while the length of \code{convert_clss} equals 1 (only one integer),
#' each discretization will count the same number of levels (quantiles).
#' @param dist.choice A character (with quotes) corresponding to the distance function chosen between: The euclidean distance ("E", by default), The Manhattan distance ("M"),
#' the Gower distance ("G"), the Hamming (also called binary) distance and the Euclidean or Manhattan distance, calculated from principal components of a factor analysis of mixed data ("FAMD").
#' @param FAMD.coord A logical that must be equal to TRUE when user decides to work with principal components of a factor analysis for mixed data (FAMD) instead of the set of raw covariates (FALSE is the default value).
#' @param FAMD.perc A percent value (between 0 and 100) linked to the \code{FAMD.coord} argument. When this latter equals TRUE, this argument corresponds to the minimum part of variability that must be taken into account by the principal components of the FAMD method.
#' This option fixes the remaining number of principal components for the rest of the study.
#' @param percent.knn Percent of closest neighbors taken in the computation of the cost matrices.
#' @param maxrelax The maximum percentage of deviation from expected probability masses (0 for the \code{OUTCOME} model, a non-zero value otherwise). Please consult (2) for more details.
#' @param prox.dist A percentage (betwen 0 and 1) used to calculate the distance threshold below which an individual (a row) is considered as a neighbor of a given profile of covariates. This choice does not influence the estimate of \eqn{\gamma},
#' but the estimates of the conditional posterior probabilities stored, in output, in the objects \code{estimatorZA} and \code{estimatorYB}.
#' @param indiv.method A character string indicating the chosen method to get individual predictions from the joint probabilities assessed, "sequential" by default, or "optimal". Please consult the section \code{details} and the corresponding article 2 for more details.
#' @param solvR A character string that specifies the type of method selected to solve the optimization algorithms. The default solver is "glpk".
#' @param which.DB A character string (with quotes) that indicates which individual predictions compute: Only the individual predictions of \eqn{Y} in B ("B"), only those of \eqn{Z} in A ("A") or the both ("BOTH" by default).
#'
#' @return A list containing 9 elements:
#' \item{time_exe}{The running time of the function}
#' \item{gamma_A}{A matrix corresponding to an estimation of the joint distribution of \eqn{(Y,Z)} in A}
#' \item{gamma_B}{A matrix corresponding to an estimation of the joint distribution of \eqn{(Y,Z)} in B}
#' \item{profile}{A data.frame that gives all details about the remaining \eqn{P} profiles of covariates. These informations can be linked to the \code{estimatorZA} and the \code{estimatorYB} objects for a better interpretation of the results.}
#' \item{res_prox}{The outputs of the function \code{proxim_dist}}
#' \item{estimatorZA}{An array that corresponds to estimates of the probability distribution of \eqn{Z} conditional to \eqn{X} and \eqn{Y} in database A. The number of rows of each table corresponds to the total number of the P profiles of covariates.
#' The number of columns of each table corresponds to the number of levels of \eqn{Y}. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of \eqn{Z}}
#' \item{estimatorYB}{An array that corresponds to estimates of the probability distribution of \eqn{Y} conditional to \eqn{X} and \eqn{Z} in database B. The number of rows of each table corresponds to the total number of profiles of covariates.
#' The number of columns of each table corresponds to the number of levels of \eqn{Z}. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of \eqn{Y}}
#' \item{DATA1_OT}{The database A with imputed individual prediction on \eqn{Z} using OT}
#' \item{DATA2_OT}{The database B with imputed individual prediction on \eqn{Y} using OT}
#'
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @references
#' \enumerate{
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N. (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679 | \url{https://doi.org/10.1515/ijb-2018-0106}
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association, DOI: 10.1080/01621459.2020.1775615
#' \item Anderberg, M.R. (1973), Cluster analysis for applications, 359 pp., Academic Press, New York, NY, USA.
#' \item Gower J.C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 27, 623–637.
#' \item Pages J. (2004). Analyse factorielle de donnees mixtes. Revue Statistique Appliquee. LII (4). pp. 93-111.
#' \item Makhorin A (2011). GNU Linear Programming Kit Reference Manual Version 4.47.\url{http://www.gnu.org/software/glpk}
#' \item Forrest J, de la Nuez D, Lougee-Heimer R (2004). Clp User Guide. \url{http://www.coin-or.org/Clp/userguide/index.html}
#' \item Theussl S, Schwendinger F, Hornik K (2020). ROI: An Extensible R Optimization Infrastructure.Journal of Statistical Software,94(15), 1-64. doi: 10.18637/jss.v094.i15 \url{https://doi.org/10.18637/jss.v094.i15)}
#' }
#'
#'
#' @seealso \code{\link{transfo_dist}},\code{\link{proxim_dist}}, \code{\link{avg_dist_closest}}, \code{\link{indiv_grp_closest}}, \code{\link{indiv_grp_optimal}}
#'
#' @import ROI ROI.plugin.glpk ROI.plugin.clp
#' @importFrom ompr MIPModel sum_expr
#' @importFrom ompr.roi with_ROI
#' @importFrom plyr mapvalues
#' @importFrom dplyr %>%
#'
#' @export
#'
#'@examples
#'
#' ### Using a sample of simu_data dataset
#' ### Y and Z are a same variable encoded in 2 different forms:
#' ### (3 levels for Y and 5 levels for Z)
#' #--------
#' data(simu_data)
#' simu_dat  = simu_data[c(1:200,301:500),]
#'
#' ### An example of OUTCOME model that uses:
#' #-----
#' # - A nearest neighbor procedure for the estimation of individual predictions
#' # - The Manhattan distance function
#' # - 90% of individuals from each modalities to calculate average distances
#' #   between individuals and modalities
#' # Predictions are assessed for Y in B and Z in A
#' #-----
#'
#' try1 = OT_outcome(simu_dat, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                   dist.choice = "M",percent.knn = 0.90, maxrelax = 0,
#'                   indiv.method = "sequential")
#' head(try1$DATA1_OT)  # Part of the completed database A
#' head(try1$DATA2_OT)  # Part of the completed database B
#'
#' head(try1$estimatorZA[,,1])
#' # ... Corresponds to P[Z = 1|Y,P1] when P1 corresponds to the 1st profile of covariates (P_1)
#' # detailed in the 1st row of the profile object:
#' try1$profile[1,]   # Details of P_1
#'
#' # So estimatorZA[1,1,1]= 0.2 corresponds to an estimation of:
#' # P[Z = 1|Y=[20-40],Gender_2=0,Treatment_2=1,Treatment_3=0,Smoking_2=1,Dosage=3,Age=65.44]
#' # Thus, we can conclude that all individuals with the P_1 profile of covariates have
#' # 20% of chance to be affected to the 1st level of Z in database A.
#' # ... And so on, the reasoning is the same for the estimatorYB object.
#'
#'
#' \dontrun{
#'
#' ### An example of OUTCOME model with same conditions as the previous example, excepted that;
#' # - Only the individual predictions of Y in B are required
#' # - The continuous covariates "age" (related index = 8) will be converted in an ordinal factors
#' #   of 3 balanced classes (tertiles)
#' # - The Gower distance is now used
#' ###-----
#'
#' try2 = OT_outcome(simu_dat, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                   dist.choice = "G",percent.knn = 0.90, maxrelax = 0,
#'                   convert.num = 8, convert.clss = 3,
#'                   indiv.method = "sequential", which.DB = "B")
#'
#'
#' ### An example of OUTCOME model with same conditions as the first example, excepted that;
#' # - Only the individual predictions of Z in A are required
#' # - The continuous covariates "age" (related index = 8) will be converted in an ordinal factors
#' #   of 3 balanced classes (tertiles)
#' # - Here, the Hamming distance can be applied because, after conversion, all covariates are factors.
#' #   Disjunctive tables of each covariates will be automatically used to work with a set of binary
#' #   variables.
#' ###-----
#'
#' try3 = OT_outcome(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                   dist.choice = "H",percent.knn = 0.90, maxrelax = 0,
#'                   convert.num = 8, convert.clss = 3,
#'                   indiv.method = "sequential",which.DB = "B")
#'
#'
#' ### An example of R-OUTCOME model using:
#' # - An optimization procedure for individual predictions on the 2 databases
#' # - The Manhattan distance
#' # - Raw covariates
#' ###-----
#' try4 = OT_outcome(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                   dist.choice = "M", percent.knn = 0.90, maxrelax = 0,
#'                   indiv.method = "optimal")
#'
#'
#' ### An example of R-OUTCOME model with:
#' # - An optimization procedure for individual predictions on the 2 databases
#' # - The use of Euclidean distance on coordinates from FAMD
#' # - Raw covariates
#' ###-----
#'
#' try5 = OT_outcome(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                   dist.choice = "E", percent.knn = 0.90,
#'                   FAMD.coord = "YES", FAMD.perc = 0.8,
#'                   indiv.method = "optimal")
#'
#'
#' ### An example of R-OUTCOME model with relaxation on marginal distributions and:
#' # - An optimization procedure for individual predictions on the 2 databases
#' # - The use of the euclidean distance
#' # - An arbitrary coefficient of relaxation
#' # - Raw covariates
#' #-----
#'
#' try6 = OT_outcome(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                   dist.choice = "E", percent.knn = 0.90,
#'                   maxrelax = 0.4, indiv.method = "optimal")
#'
#'
#' }
#'
OT_outcome = function(datab, index_DB_Y_Z = 1:3,
                      quanti = NULL, nominal = NULL, ordinal = NULL,logic = NULL,
                      convert.num = NULL, convert.clss = NULL, FAMD.coord = "NO", FAMD.perc = 0.8,
                      dist.choice = "E", percent.knn = 1,
                      maxrelax = 0, indiv.method = "sequential", prox.dist = 0.80,
                      solvR = "glpk", which.DB = "BOTH"){



  tstart = Sys.time()

  cat("---------------------------------------","\n")
  cat("OT PROCEDURE in progress ..."           ,"\n")
  cat("---------------------------------------","\n")
  cat("Type                     = ", ifelse((maxrelax == 0) & (indiv.method == "sequential"),"OUTCOME","R-OUTCOME"),"\n")
  cat("Distance                 = ", ifelse(dist.choice == "H","Hamming",ifelse(dist.choice == "M","Manhattan",ifelse(dist.choice == "E","Euclidean","Gower"))),"\n")
  cat("Percent closest knn      = ", 100.0*percent.knn,"%","\n")
  cat("Relaxation parameter     = ", ifelse(maxrelax == 0,"NO","YES"),"\n")
  cat("Relaxation value         = ", maxrelax,"\n")
  cat("Individual pred process  = ", ifelse(indiv.method == "sequential","Sequential","Optimal"),"\n")
  cat("DB imputed               = ", which.DB,"\n")
  cat("---------------------------------------","\n")


  if (index_DB_Y_Z[2] %in% quanti){

    datab[,index_DB_Y_Z[2]] = ordered(datab[,index_DB_Y_Z[2]])
    quanti                  = setdiff(quanti,index_DB_Y_Z[2])
    ordinal                 = sort(c(ordinal,index_DB_Y_Z[2]))

  } else {}

  if (index_DB_Y_Z[3] %in% quanti){

    datab[,index_DB_Y_Z[3]] = ordered(datab[,index_DB_Y_Z[3]])
    quanti                  = setdiff(quanti,index_DB_Y_Z[3])
    ordinal                 = sort(c(ordinal,index_DB_Y_Z[3]))

  } else {}



  if (!(which.DB %in% c("A","B","BOTH"))){

    stop("Invalid which.DB argument")

  } else {}

  if (!(dist.choice %in% c("E","M","H","G"))){

    stop("Invalid dist.choice argument")

  } else {}

  if (!(indiv.method %in% c("sequential","optimal"))){

    stop("Invalid indiv.method argument")

  } else {}

  if (!(FAMD.coord %in% c("YES","NO","Y","N"))){

    stop("Invalid FAMD.coord argument")

  } else {}

  prep.choice = ifelse(FAMD.coord %in% c("YES","Y"), "FAMD", dist.choice)

  dataB = transfo_dist(datab,index_DB_Y_Z = index_DB_Y_Z,
                       quanti = quanti, nominal = nominal, ordinal = ordinal, logic = logic,
                       convert_num = convert.num, convert_clss = convert.clss,
                       prep_choice = dist.choice, info = FAMD.perc)


  if (prep.choice == "H"){

    datac  = dataB[,-index_DB_Y_Z]

    test_H = apply(as.data.frame(datac),2,function(x){length(table(x))})

    if (max(test_H)>2){stop("With Hamming distance, all your covariates must be binaries !")} else {}

  } else {}


  inst = proxim_dist(dataB, norm = dist.choice, prox = prox.dist)


  # Redefine A and B for the model
  nA      = inst$nA
  nB      = inst$nB
  A       = 1:inst$nA
  B       = 1:inst$nB
  Y       = inst$Y
  Z       = inst$Z
  indY    = inst$indY
  indZ    = inst$indZ
  indXA   = inst$indXA
  indXB   = inst$indXB
  nbX     = length(indXA)

  #prof  = do.call(paste0,unique(inst$Xobserv))
  prof    = as.data.frame(unique(inst$Xobserv))
  ID_prof = paste(rep("P",nrow(prof)),1:nrow(prof),sep="_")

  nbindY = numeric(0)

  for (y in Y){

    nbindY  = c(nbindY,length(inst$indY[[y]]))}
  nbindZ  = numeric(0)

  for (z in Z){

    nbindZ  = c(nbindZ,length(inst$indZ[[z]]))}
  freqY   = numeric(0)

  for (y in Y){

    freqY   = c(freqY,nbindY[y] / length(A))}
  freqZ   = numeric(0)

  for (z in Z){

    freqZ = c(freqZ,nbindZ[z] / length(B))

  }


  ###########################################################################
  # Compute data for aggregation of the individuals
  ###########################################################################

  nbX = length(indXA);

  # Computation of the cost matrix as average distances between the
  # individuals of two groups
  C = avg_dist_closest(inst, percent_closest = percent.knn)[[1]]

  # absdevA = absdevB = deviationA = deviationB = transport = transportA = transportB = NULL

  if (maxrelax == 0){

    # Objective: minimize the distance between individuals of A and B

    Min = c(as.numeric(C),rep(0,2*length(Y)+2*length(Z)))


    result          = ompr::MIPModel() %>%

      ompr::add_variable(transport[y,z], y = Y, z = Z, type = "continuous") %>%
      ompr::add_variable(deviationA[y] , y = Y       , type = "continuous") %>%
      ompr::add_variable(absdevA[y]    , y = Y       , type = "continuous") %>%
      ompr::add_variable(deviationB[z]        , z = Z, type = "continuous") %>%
      ompr::add_variable(absdevB[z]           , z = Z, type = "continuous") %>%

      ompr::set_objective(sum_expr(C[y,z]*transport[y,z], y = Y, z = Z) , "min") %>%

      ompr::add_constraint(sum_expr(transport[y,z]       , z = Z) == freqY[y] + deviationA[y], y = Y) %>%
      ompr::add_constraint(sum_expr(transport[y,z], y = Y       ) == freqZ[z] + deviationB[z], z = Z) %>%
      ompr::add_constraint(sum_expr(deviationA[y] , y = Y) == 0 ) %>%
      ompr::add_constraint(sum_expr(deviationB[z]        ,z = Z )== 0) %>%
      ompr::add_constraint(deviationB[z]<= absdevB[z]    , z = Z) %>%
      ompr::add_constraint(deviationB[z]>= -absdevB[z]   , z = Z) %>%
      ompr::add_constraint(sum_expr(absdevB[z],z=Z)<= maxrelax/2.0) %>%
      ompr::add_constraint(deviationA[y] <= absdevA[y],y =Y) %>%
      ompr::add_constraint(deviationA[y] >= -absdevA[y],y =Y) %>%
      ompr::add_constraint(sum_expr(absdevA[y],y = Y)<= maxrelax/2.0) %>%

      ompr::solve_model(with_ROI(solver = solvR))

    # Solve the problem

    solution       = ompr::get_solution(result, transport[y,z])

    # Extract the values of the solution

    transportA_val = matrix(solution[,4], length(Y),length(Z))
    transportB_val = transportA_val


  } else {

    b2 = numeric(length(Z))

    for (z in Z){

      stc_sum2 = vector(length = nbX)

      for (i in 1:nbX){

        stc_sum2[i] = ifelse(length(indXB[[i]])==0,1/length(Z),length(indXB[[i]][inst$Zobserv[indXB[[i]]+nA] == z])/ length(indXB[[i]]))*length(indXA[[i]])/nA

      }
      b2[z]= sum(stc_sum2)
    }
    b2   = b2/sum(b2)


    result <-  MIPModel() %>%

      ompr::add_variable(transportA[y,z],y = Y, z = Z, type = "continuous",lb=0) %>%
      ompr::add_variable(deviationB[z]  ,       z = Z, type = "continuous") %>%
      ompr::add_variable(absdevB[z]     ,       z = Z, type = "continuous",lb=0) %>%

      ompr::set_objective(sum_expr(C[y,z]*transportA[y,z],y = Y,z=Z) , "min") %>%

      ompr::add_constraint(sum_expr(transportA[y,z], z = Z) == freqY[y], y =Y) %>%
      ompr::add_constraint(sum_expr(transportA[y,z],y = Y) - deviationB[z] == b2[z], z = Z) %>%
      ompr::add_constraint(sum_expr(deviationB[z],z = Z)== 0) %>%
      ompr::add_constraint(deviationB[z]<= absdevB[z], z = Z) %>%
      ompr::add_constraint(deviationB[z]>= -absdevB[z], z = Z) %>%
      ompr::add_constraint(sum_expr(absdevB[z],z=Z)<= maxrelax/2.0) %>%

      ompr::solve_model(with_ROI(solver = solvR))

    solution       = ompr::get_solution(result, transportA[y,z])

    transportA_val = matrix(solution$value, length(Y),length(Z))


    b1 = numeric(length(Y))

    for (y in Y){

      stc_sum2 = vector(length = nbX)

      for (i in 1:nbX){

        stc_sum2[i] = ifelse(length(indXA[[i]])==0,1/length(Y),length(indXA[[i]][inst$Yobserv[indXA[[i]]] == y])/ length(indXA[[i]]))*length(indXB[[i]])/nB

      }

      b1[y]= sum(stc_sum2)}
    b1   = b1/sum(b1)

    result <-  MIPModel() %>%

      ompr::add_variable(transportB[y,z],y = Y, z= Z, type = "continuous",lb=0) %>%
      ompr::add_variable(deviationA[y]  , y = Y     , type = "continuous") %>%
      ompr::add_variable(absdevA[y]     , y = Y     , type = "continuous",lb=0) %>%

      ompr::set_objective(sum_expr(C[y,z]*transportB[y,z],y = Y,z=Z) , "min") %>%

      ompr::add_constraint(sum_expr(transportB[y,z], z = Z) - deviationA[y] == b1[y], y =Y) %>%
      ompr::add_constraint(sum_expr(transportB[y,z],y = Y) == freqZ[z], z = Z) %>%
      ompr::add_constraint(sum_expr(deviationA[y],y = Y)== 0) %>%
      ompr::add_constraint(deviationA[y] <= absdevA[y],y =Y) %>%
      ompr::add_constraint(deviationA[y] >= -absdevA[y],y =Y) %>%
      ompr::add_constraint(sum_expr(absdevA[y],y = Y)<= maxrelax/2.0) %>%

      ompr::solve_model(with_ROI(solver = solvR))


    solution       = ompr::get_solution(result, transportB[y,z])
    transportB_val = matrix(solution$value, length(Y),length(Z))


  }

  ####
  # Get the individual transport from the group transport

  if (indiv.method == "sequential"){

    indpred = indiv_grp_closest(inst, transportA_val, transportB_val, percent_closest = percent.knn, which.DB = which.DB)

  } else if (indiv.method == "optimal"){

    indpred = indiv_grp_optimal(inst, transportA_val, transportB_val, percent_closest = percent.knn, solvr = solvR, which.DB = which.DB)

  } else {}

  if (which.DB == "A"){

    YApred      = NULL
    YBpred      = indpred$ZBtrans
    estimatorZA = array(rep(0,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))
    estimatorYB = NULL

    DATA1_OT         = dataB[dataB[,1] == unique(dataB[,1])[1],]
    DATA2_OT         = dataB[dataB[,1] == unique(dataB[,1])[2],]

    if (index_DB_Y_Z[3] %in% nominal){

      DATA1_OT$OTpred  = factor(plyr::mapvalues(YBpred,from = sort(unique(YBpred)),
                                                to = levels(dataB[,3])[sort(unique(YBpred))]),
                                                levels = levels(dataB[,3])[sort(unique(YBpred))])
    } else {

      DATA1_OT$OTpred  = ordered(plyr::mapvalues(YBpred,from = sort(unique(YBpred)),
                                                 to = levels(dataB[,3])[sort(unique(YBpred))]),
                                                 levels = levels(dataB[,3])[sort(unique(YBpred))])

    }


    # if (is.ordered(dataB[,3])){

    #    DATA1_OT$OTpred = as.ordered(DATA1_OT$OTpred)

    # } else {}


  } else if (which.DB == "B"){

    YApred = indpred$YAtrans
    YBpred = NULL
    # Compute the estimated probability distributions from predictions
    estimatorZA = NULL
    estimatorYB = array(rep(0,nbX*length(Y)*length(Z)),dim = c(nbX,length(Z),length(Y)))

    DATA1_OT         = dataB[dataB[,1] == unique(dataB[,1])[1],]
    DATA2_OT         = dataB[dataB[,1] == unique(dataB[,1])[2],]

    if (index_DB_Y_Z[2] %in% nominal){

        DATA2_OT$OTpred  = factor(plyr::mapvalues(YApred,from = sort(unique(YApred)),
                                                  to = levels(dataB[,2])[sort(unique(YApred))]),
                                                  levels = levels(dataB[,2])[sort(unique(YApred))])
    } else {

        DATA2_OT$OTpred  = ordered(plyr::mapvalues(YApred,from = sort(unique(YApred)),
                                                   to = levels(dataB[,2])[sort(unique(YApred))]),
                                                   levels = levels(dataB[,2])[sort(unique(YApred))])
    }

    #if (is.ordered(dataB[,2])){

    #  DATA2_OT$OTpred = as.ordered(DATA2_OT$OTpred)

    # } else {}

  } else {

    YApred = indpred$YAtrans
    YBpred = indpred$ZBtrans
    # Compute the estimated probability distributions from predictions
    estimatorZA = array(rep(0,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))
    estimatorYB = array(rep(0,nbX*length(Y)*length(Z)),dim = c(nbX,length(Z),length(Y)))

    DATA1_OT         = dataB[dataB[,1] == unique(dataB[,1])[1],]
    DATA2_OT         = dataB[dataB[,1] == unique(dataB[,1])[2],]

    if (index_DB_Y_Z[3] %in% nominal){

      DATA1_OT$OTpred  = factor(plyr::mapvalues(YBpred,from = sort(unique(YBpred)),
                                                to = levels(dataB[,3])[sort(unique(YBpred))]),
                                                levels = levels(dataB[,3])[sort(unique(YBpred))])
    } else {

      DATA1_OT$OTpred  = ordered(plyr::mapvalues(YBpred,from = sort(unique(YBpred)),
                                                 to = levels(dataB[,3])[sort(unique(YBpred))]),
                                                 levels = levels(dataB[,3])[sort(unique(YBpred))])

    }

    if (index_DB_Y_Z[2] %in% nominal){

      DATA2_OT$OTpred  = factor(plyr::mapvalues(YApred,from = sort(unique(YApred)),
                                                to = levels(dataB[,2])[sort(unique(YApred))]),
                                levels = levels(dataB[,2])[sort(unique(YApred))])
    } else {

      DATA2_OT$OTpred  = ordered(plyr::mapvalues(YApred,from = sort(unique(YApred)),
                                                 to = levels(dataB[,2])[sort(unique(YApred))]),
                                 levels = levels(dataB[,2])[sort(unique(YApred))])
    }

    # DATA1_OT$OTpred  = as.factor(plyr::mapvalues(YBpred,from = sort(unique(YBpred)), to = levels(dataB[,3])[sort(unique(YBpred))]))
    # DATA2_OT$OTpred  = as.factor(plyr::mapvalues(YApred,from = sort(unique(YApred)), to = levels(dataB[,2])[sort(unique(YApred))]))

    # if (is.ordered(dataB[,3])){

      # DATA1_OT$OTpred = as.ordered(DATA1_OT$OTpred)

    #} else {}

    #if (is.ordered(dataB[,2])){

    #  DATA2_OT$OTpred = as.ordered(DATA2_OT$OTpred)

    # } else {}

  }

  for (x in 1:nbX){

    if (which.DB %in% c("BOTH","A")){

      for (i in indXA[[x]]){

        estimatorZA[x,inst$Yobserv[i],YBpred[i]] = estimatorZA[x,inst$Yobserv[i],YBpred[i]] +  1/sum(inst$Yobserv[indXA[[x]]] == inst$Yobserv[i])
      }


      for (y in Y){

        if (sum(inst$Yobserv[indXA[[x]]] == y) == 0){

          estimatorZA[x,y,] = 1/length(Z)*rep(1,length(Z))
        }

      }

      row.names(estimatorZA) = ID_prof
      colnames(estimatorZA)  = as.character(levels(dataB[,2]))

    } else {}

    if (which.DB %in% c("BOTH","B")){

      for (i in indXB[[x]]){

        # estimatorYB[x,YApred[i],inst$Zobserv[i+nA]] = estimatorYB[x,YApred[i],inst$Zobserv[i+nA]] + 1/ sum(inst$Zobserv[indXB[[x]] + nA] == inst$Zobserv[i + nA])
        estimatorYB[x,inst$Zobserv[i+nA],YApred[i]] = estimatorYB[x,inst$Zobserv[i+nA],YApred[i]] + 1/ sum(inst$Zobserv[indXB[[x]] + nA] == inst$Zobserv[i + nA])
      }

      for (z in Z){

        if (sum(inst$Zobserv[indXB[[x]]+inst$nA] == z) == 0){

          # estimatorYB[x,,z] = 1/length(Y)*rep(1,length(Y))
          estimatorYB[x,z,] = 1/length(Y)*rep(1,length(Y))

        }
      }

      row.names(estimatorYB) = ID_prof
      colnames(estimatorYB)  = as.character(levels(dataB[,3]))

    } else {}

  }

  tend = Sys.time()

  return(list(time_exe = difftime(tend,tstart), gamma_A = transportA_val,gamma_B = transportB_val,
              profile = data.frame(ID = ID_prof,prof),res_prox = inst, estimatorZA= estimatorZA,estimatorYB = estimatorYB,DATA1_OT = DATA1_OT,DATA2_OT  = DATA2_OT))
}



