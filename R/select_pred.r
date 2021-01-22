#' select_pred()
#'
#' Selection of a subset of non collinear predictors having relevant relationships with a given target outcome using a random forest procedure.
#'
#' The \code{select_pred} function provides several tools to identify, on the one hand, the relationships between predictors, by detecting especially potential problems of collinearity, and, on the other hand, proposes a parcimonious subset of relevant predictors (of the outcome) using appropriate random forest procedures.
#' The function which can be used as a preliminary step of prediction in regression areas is particularly adapted to the context of data fusion by providing relevant subsets of predictors (the matching variables) to algorithms dedicated to the solving of recoding problems.
#'
#'
#' A. REQUIRED STRUCTURE FOR THE DATABASE
#'
#' The expected input database is a data.frame that especially requires a specific column of row identifier and a target variable (or outcome) having a finite number of values or classes (ordinal, nominal or discrete type). Notice that if the chosen outcome is in numeric form, it will be automatically converted in ordinal type.
#' The number of predictors is not a constraint for \code{select_pred} (even if, with less than three variables a process of variables selection has no real sense...), and can exceed the number of rows (no problem of high dimensionality here).
#' The predictors can be continuous (quantitative), boolean, nominal or ordinal with or without missing values.
#' In presence of numeric variables, users can decide to discretize them or a part of them by themselves beforehand. They can also choose to use the internal process directly integrated in the function. Indeed, to assist users in this task, two arguments called \code{convert_num} and \code{convert_clss} dedicated to these transformations are available in input of the function.
#' These options make the function \code{select_pred} particularly adapted to the function \code{\link{OT_joint}} which only allows data.frame with categorical covariates.
#' With the argument \code{convert_num}, users choose the continuous variables to convert and the related argument \code{convert_clss} specifies the corresponding number of classes chosen for each discretization.
#' It is the reason why these two arguments must be two vectors of indexes of same length. Nevertheless, an unique exception exists when \code{convert_clss} is equalled to a scalar \eqn{S}. In this case, all the continuous predictors selected for conversion will be discretized with a same number of classes S.
#' By example, if \code{convert_clss = 4}, all the continuous variables specified in the \code{convert_num} argument will be discretized by quartiles. Moreover, notice that missing values from incomplete predictors to convert are not taken into account during the conversion, and that each predictor specified in the argument \code{convert_num}
#' must be also specified in the argument \code{quanti}.
#' In this situation, the label of the outcome must be entered in the argument \code{Y}, and the arguments \code{Z} and \code{OUT} must keep their default values.
#' Finally, the order of the column indexes related to the identifier and the outcome have no importance.
#'
#' For a better flexibility, the input database can also be the result of two overlayed databases.
#' In this case, the structure of the database must be similar to those observed in the datasets \code{\link{simu_data}} and \code{\link{tab_test}} available in the package with a column of database identifier, one target outcome by database (2 columns), and a subset of shared predictors.
#' Notice that, overlaying two separate databases can also be done easily using the function \code{\link{merge_dbs}} beforehand.
#' The labels of the two outcomes will have to be specified in the arguments \code{Y} for the top database, and in \code{Z} for the bottom one.
#' Notice also that the function \code{select_pred} deals with only one outcome at a time that will have to be specified in the argument \code{OUT} which must be equalled to "Y" for the study of the top database or "Z" for the study of the bottom one.
#'
#' Finally, whatever the structure of the database declared in input, each column index related to the database variable must be entered once (and only once) in one of the following four arguments: \code{quanti}, \code{nominal}, \code{ordinal}, \code{logic}.
#'
#'
#' B. PAIRWISE ASSOCIATIONS BETWEEN PREDICTORS
#'
#' In a first step of process, \code{select_pred} calculates standard pairwise associations between predictors according to their types.
#' \enumerate{
#' \item{Between categorical predictors (ordinal, nominal and logical):
#' Cramer's V (and Bias-corrected Cramer's V, see (1) for more details) are calculated between categorical predictors and the argument \code{thres_cat} fixed the associated threshold beyond which two predictors can be considered as redundant.
#' A similar process is done between the target variable and the subset of categorical variables which provides in output a first table ranking the top scoring predictors. This table summarizes the ability of each variable to predict the target outcome.}
#' \item{Between continuous predictors:
#' If the \code{ordinal} and \code{logic} arguments differ from NULL, all the corresponding predictors are beforehand converted in rank values.
#' For numeric (quantitative), logical and ordinal predictors, pairwise correlations between ranks (Spearman) are calculated and the argument \code{thresh_num} fixed the related threshold beyond which two predictors can be considered as redundant.
#' A similar process is done between the outcome and the subset of discrete variables which provides in output, a table ranking the top scoring predictor variates which summarizes their abilities to predict the target.
#' In addition, the result of a Farrar and Glauber test is provided. This test is based on the determinant of the correlation matrix of covariates and the related null hypothesis of the test corresponds to an absence of collinearity between them (see (2) for more details about the method).
#' In presence of a large number of numeric covariates and/or ordered factors, the approximate Farrar-Glauber test, based on the normal approximation of the null distribution is more adapted and its result is also provided in output.
#' These two tests are highly sensitive and, by consequence, it suggested to consider these results as simple indicators of collinearity between predictors rather than an essential condition of acceptability.}
#' }
#' If the initial number of predictors is not too important, these informations can be sufficient to the user for the visualization of potential problems of collinearity and for the selection of a subset of predictors (\code{RF = FALSE}).
#' It is nevertheless often necessary to complete this visualization by an automatical process of selection like the Random Forest approach (see Breiman 2001, for a better understanding of the method) linked to the function \code{select_pred} (\code{RF = TRUE}).
#'
#'
#' C. RANDOM FOREST PROCEDURE
#'
#' As a final step of the process, a random forest approach (RF(3)) is here prefered (to regression models) for two main reasons: RF methods allow notably the number of variables to exceed the number of rows and remain applicable whatever the types of covariates considered.
#' The function \code{select_pred} integrates in its algorithm the functions \code{\link[party]{cforest}} and \code{\link[party]{varimp}} of the package \pkg{party} (Hothorn, 2006) and so gives access to their main arguments.
#'
#' A RF approach generally provides two types of measures for estimating the mean variable importance of each covariate in the prediction of an outcome: the Gini importance and the permutation importance. These measurements must be used with caution, by taking into account the following constraints:
#' \enumerate{
#' \item{The Gini importance criterion can produce bias in favor of continuous variables and variables with many categories. To avoid this problem, only the permutation criterion is available in the function.}
#' \item{The permutation importance criterion can overestimate the importance of highly correlated predictors.}
#' }
#' The function \code{select_pred} proposes three different scenarios according to the types of predictors:
#' \enumerate{
#' \item{The first one consists in boiling down to a set of categorical variables (ordered or not) by discretizing all the continuous predictors beforehand, using the internal \code{convert_num} argument or another one, and then works with the conditional importance measures (\code{RF_condi = TRUE}) which give unbiased estimations.
#' In the spirit of a partial correlation, the conditional importance measure related to a variable \eqn{X} for the prediction of an outcome \eqn{Y}, only uses the subset of variables the most correlated to \eqn{X} for its computation. The argument \code{RF_condi_thr} that corresponds exactly to the argument \code{threshold} of the function \code{\link[party]{varimp}},
#' fixes a ratio below which a variable Z is considered sufficiently correlated to \eqn{X} to be used as an adjustment variable in the computation of the importance measure of \eqn{X} (In other words, Z is included in the conditioning for the computation, see (4) and (5) for more details). A threshold value of zero will include all variables in the computation of
#' conditional importance measure of each predictor \eqn{X}, while a threshold \eqn{< 1}, will only include a subset of variables.
#' Two remarks related to this method: firstly, notice that taking into account only subsets of predictors in the computation of the variable importance measures could lead to a relevant saving of execution time.
#' Secondly, because this approach does not take into account incomplete information, the method will only be applied to complete data (incomplete rows will be temporarily removed for the study).}
#' \item{The second possibility, always in presence of mixed types predictors, consists in the execution of two successive RF procedures. The first one will be used to select an unique candidate in each susbset of correlated predictors (detecting in the 1st section), while the second one will extract the permutation measures from the remaining subset
#' of uncorrelated predictors (\code{RF_condi = FALSE}, by default). This second possibility has the advantage to work in presence of incomplete predictors.}
#' \item{The third scenario consists in running a first time the function without RF process (\code{RF = FALSE}), and according to the presence of highly correlated predictors or not, users can choose to extract redundant predictors manually and re-runs the function with the subset of remaining non-collinear predictors to avoid potential biases introduced by the standard permutations measures.}
#' }
#' The three scenarios finally lead to a list of uncorrelated predictors of the outcome sorted in importance order. The argument \code{thresh_Y} corresponds to the minimal percent of importance required (and fixed by user) for a variable to be considered as a reliable predictor of the outcome.
#' Finally, because all random forest results are subjects to random variation, users can check whether the same importance ranking is achieved by varying the random seed parameter (\code{RF_SEED}) or by increasing the number of trees (\code{RF_ntree}).
#'
#'
#' @param databa a data.frame with a column of identifiers (of row or of database in the case of two concatened databases), an outcome, and a set of predictors. The number of columns can exceed the number of rows.
#' @param Y the label of a first target variable with quotes
#' @param Z the label of a second target variable with quotes when \code{databa} is the result of two overlayed databases.
#' @param ID the column index of the database identifier (The first column by default) in the case of two concatened databases, a row identifier otherwise
#' @param OUT a character that indicates the outcome to predict in the context of overlayed databases. By default, the outcome declared in the argument \code{Y} is predicted. Another possible outcome to predict can be set with the related argument \code{Z}.
#' @param quanti a vector of integers corresponding to the column indexes of all the numeric predictors.
#' @param nominal a vector of integers which corresponds to the column indexes of all the categorical nominal predictors.
#' @param ordinal a vector of integers which corresponds to the column indexes of all the categorical ordinal predictors.
#' @param logic a vector of integers indicating the indexes of logical predictors. No index remained by default
#' @param convert_num a vector of integers indicating the indexes of quantitative variables to convert in ordered factors. No index remained by default. Each index selected has to be defined as quantitative in the argument \code{quanti}.
#' @param convert_clss a vector of integers indicating the number of classes related to each transformation of quantitative variable in ordered factor. The length of this vector can not exceed the length of the argument \code{convert_num}. Nevertheless, if length(\code{convert_num}) > 1 and length(\code{convert_clss}) = 1,
#' all quantitative predictors selected for discretization will have by default the same number of classes.
#' @param thresh_cat a threshold associated to the Cramer's V coefficient (= 0.30 by default)
#' @param thresh_num a threshold associated to the Spearman's coefficient of correlation (= 0.70 by default)
#' @param thresh_Y a threshold linked to the RF approach, that corresponds to the minimal cumulative percent of importance measure required to be kept in the final list of predictors.
#' @param RF a boolean sets to TRUE (default) if a random forest procedure must be applied to select the best subset of predictors according to the outcome.Otherwise, only pairwise associations between predictors are used for the selection.
#' @param RF_ntree the number of bootsrap samples required from the row datasource during the random forest procedure
#' @param RF_condi a boolean specifying if the conditional importance measures must be assessed from the random forest procedure (\code{TRUE}) rather than the standard variable importance  measures (\code{FALSE} by default)
#' @param RF_condi_thr a threshold linked to (1 - pvalue) of an association test between each predictor \eqn{X} and the other variables, given that a threshold value of zero will include all variables in the computation of the conditional importance measure of \eqn{X} (0.20 is the default value).
#' Conversely, a larger threshold will only keeps the subset of variables that is strongly correlated to \eqn{X} for the computation of the variable importance measure of \eqn{X}.
#' @param RF_SEED an integer used as argument by the set.seed() for offsetting the random number generator (random integer by default). This value is only used for RF method.
#'
#'
#' @return A list of 14 (if \code{RF = TRUE}) or 11 objects (Only the first ten objects if \code{RF = FALSE}) is returned:
#' \item{seed}{the random number generator related to the study}
#' \item{outc}{the identifier of the outcome to predict}
#' \item{thresh}{a summarize of the different thresholds fixed for the study}
#' \item{convert_num}{the labels of the continuous predictors transformed in categorical form}
#' \item{DB_USED}{the final database used after potential transformations of predictors}
#' \item{vcrm_Y_cat}{a table of pairwise associations between the outcome and the categorical predictors (Cramer's V)}
#' \item{cor_Y_num}{a table of pairwise associations between the outcome and the continuous predictors (Rank correlation)}
#' \item{vcrm_X_cat}{a table of pairwise associations between the categorical predictors (Cramer's V)}
#' \item{cor_X_num}{a table of pairwise associations between the continuous predictors (Cramer's V)}
#' \item{FG_test}{the results of the Farrar and Glauber tests, with and without approximation form}
#' \item{collinear_PB}{a table of predictors with problem of collinearity according to the fixed thresholds}
#' \item{drop_var}{the labels of predictors to drop after RF process (optional output: only if \code{RF}=TRUE)}
#' \item{RF_PRED_Y}{the table of variable importance measurements, conditional or not, according to the argument \code{condi_RF} (optional output: Only if \code{RF}=TRUE)}
#' \item{best_pred}{the labels of the best predictors selected (optional output: Only if \code{RF}=TRUE) according to the value of the argument \code{thresh_Y}}
#'
#' @export
#'
#' @references
#' \enumerate{
#' \item Bergsma W. (2013). A bias-correction for Cramer's V and Tschuprow's T. Journal of the Korean Statistical Society, 42, 323–328.
#' \item Farrar D, and Glauber R. (1968). Multicolinearity in regression analysis. Review of Economics and Statistics, 49, 92–107.
#' \item Breiman L. (2001). Random Forests. Machine Learning, 45(1), 5–32.
#' \item Hothorn T, Buehlmann P, Dudoit S, Molinaro A, Van Der Laan M (2006). “Survival Ensembles.” Biostatistics, 7(3), 355–373.
#' \item Strobl C, Boulesteix A-L, Kneib T, Augustin T, Zeileis A (2008). Conditional Variable Importance for Random Forests. BMC Bioinformatics, 9, 307. \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-307}
#' }
#'
#' @author Gregory Guernec
#'
#' \email{gregory.guernec@@inserm.fr}
#'
#' @importFrom stats cor cor.test pnorm pchisq na.omit formula
#' @importFrom StatMatch pw.assoc
#' @importFrom party varimp cforest cforest_unbiased
#'
#'
#' @aliases select_pred
#'
#' @seealso \code{\link{simu_data}}, \code{\link{tab_test}}, \code{\link{OT_joint}}
#'
#' @examples
#'
#' ### Example 1
#' #-----
#' # - From two overlayed databases: using the table simu_data
#' # - Searching for the best predictors of "Yb1"
#' # - Using the row database
#' # - The RF approaches are not required
#' #-----
#'
#' data(simu_data)
#' test_DB1 = select_pred(simu_data,Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Y",
#'                        quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                        thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
#'                        RF = FALSE)
#'
#' ### Example 2
#' #-----
#' # - With same conditions as example 1
#' # - Searching for the best predictors of "Yb2"
#' #-----
#'
#' test_DB2 = select_pred(simu_data,Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Z",
#'                        quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                        thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
#'                        RF = FALSE)
#'
#' \donttest{
#' ### Example 3
#' #-----
#' # - With same conditions as example 1
#' # - Using a RF approach to estimate the standard variable importance measures
#' #   and determine the best subset of predictors
#' # - Here a seed is required
#' #-----
#'
#' test_DB3 = select_pred(simu_data,Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Y",
#'                        quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                        thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
#'                        RF = TRUE, RF_condi = FALSE, RF_SEED = 3023)
#'
#' ### Example 4
#' #-----
#' # - With same conditions as example 1
#' # - Using a RF approach to estimate the conditional variable importance measures
#' #   and determine the best subset of predictors
#' # - This approach requires to convert the numeric variables: Only "Age" here
#' #   discretized in 3 levels
#' #-----
#'
#' test_DB4 = select_pred(simu_data,Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Z",
#'                        quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                        convert_num = 8, convert_clss = 3,
#'                        thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
#'                        RF = TRUE, RF_condi = TRUE, RF_condi_thr = 0.60, RF_SEED = 3023)
#'
#' ### Example 5
#' #-----
#' # - Starting with a unique database
#' # - Same conditions as example 1
#' #-----
#' simu_A = simu_data[simu_data$DB == "A",-3]    # Base A
#'
#' test_DB5 = select_pred(simu_A, Y = "Yb1",
#'                        quanti = 7, nominal = c(1,3:4,6), ordinal = c(2,5),
#'                        thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
#'                        RF = FALSE)
#'
#' ### Example 6
#' #-----
#' # - Starting with an unique database
#' # - Using a RF approach to estimate the conditional variable importance measures
#' #   and determine the best subset of predictors
#' # - This approach requires to convert the numeric variables: Only "Age" here
#' #   discretized in 3 levels
#' #-----
#'
#' simu_B = simu_data[simu_data$DB == "B",-2]    # Base B
#'
#' test_DB6 = select_pred(simu_B, Y = "Yb2",
#'                        quanti = 7, nominal = c(1,3:4,6), ordinal = c(2,5),
#'                        convert_num = 7, convert_clss = 3,
#'                        thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
#'                        RF = TRUE, RF_condi = TRUE, RF_condi_thr = 0.60, RF_SEED = 3023)
#'
#' }
#'
select_pred = function(databa,Y = NULL, Z = NULL, ID = 1, OUT = "Y",
                       quanti = NULL, nominal = NULL, ordinal = NULL, logic = NULL,
                       convert_num = NULL, convert_clss = NULL,
                       thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
                       RF = TRUE, RF_ntree = 500, RF_condi = FALSE, RF_condi_thr = 0.20, RF_SEED = sample(1:1000000, 1)){

  if (!(OUT %in% c("Y","Z"))){

    stop("Improper argument for OUT: Y or Z only")

  } else {}

  if ((OUT == "Y") & (is.null(Y))){

    stop("When OUT = Y, Y can not be null")

  } else {}

  if ((OUT == "Z") & (is.null(Z))){

    stop("When OUT = Z, Z can not be null")

  } else {}


  # Exclude ID from arguments

  quanti  = setdiff(quanti ,ID)
  nominal = setdiff(nominal,ID)
  ordinal = setdiff(ordinal,ID)
  logic   = setdiff(logic  ,ID)


  # Convert boolean covariates

  if (length(logic) != 0){

    # Verification

    if (all(sapply(as.data.frame(databa[,logic]),is.logical))!=TRUE){

      stop("At least one logical variable is not boolean")

    } else {}

    # Conversion

    for (j in logic){

      databa[,j]   = as.factor(databa[,j])
      nominal      = unique(sort(c(nominal,logic)))
      logic        = NULL

    }

  } else {}


  ### Tests

  stopifnot(is.data.frame(databa))
  stopifnot(is.data.frame(databa))
  stopifnot(thresh_cat   <= 1)   ; stopifnot(thresh_cat   >=0)
  stopifnot(thresh_num   <= 1)   ; stopifnot(thresh_num   >=0)
  stopifnot(thresh_Y     <= 1)   ; stopifnot(thresh_Y     >=0)
  stopifnot(RF_condi_thr <= 1)   ; stopifnot(RF_condi_thr >= 0)
  stopifnot(!is.character(ordinal)); stopifnot(!is.character(quanti));
  stopifnot(!is.character(nominal)); stopifnot(!is.character(logic))
  stopifnot(!is.character(convert_num));
  stopifnot(!is.character(convert_clss));
  stopifnot(length(intersect(ordinal,quanti))  == 0);
  stopifnot(length(intersect(ordinal,nominal)) == 0);
  stopifnot(length(intersect(ordinal,quanti))  == 0);
  stopifnot(is.logical(RF_condi))
  stopifnot(is.numeric(ID))



  # Inconsistencies options

  idglob = unique(c(ID,quanti,ordinal,nominal,logic))

  if (length(setdiff(1:ncol(databa),idglob))!=0){

    stop("Invalid indexes: No intersection allowed between the indexes of the ID, quanti, ordinal, nominal, and logic arguments")

  } else {}


  if (length(convert_clss)>length(convert_num)){

    stop("Inconsistencies between convert_num and convert_clss")

  } else {}


  if ((length(convert_clss)>1)&(length(convert_clss)!=length(convert_num))){

    stop("Inconsistencies between convert_num and convert_clss")

  } else {}


  if (length(convert_clss) == 1){

    convert_clss = rep(convert_clss,length(convert_num))

  } else {}


  message("The select_pred function is running. Please wait ...","\n")


  # Exclude systematically ID and Y and Z before discretization

  indexY       = (1:ncol(databa))[colnames(databa)== Y]
  indexZ       = (1:ncol(databa))[colnames(databa)== Z]
  index_DB_Y_Z = c(ID,indexY,indexZ)

  convert_clss = convert_clss[!(convert_num %in% index_DB_Y_Z)]
  convert_num  = setdiff(convert_num,index_DB_Y_Z)


  ### Re-classmt to avoid NAs in complete levels of factors

  if (length(ordinal)!=0){

    for (k in ordinal){

      databa[,k] = ordered(databa[,k])

    }

  } else {}

  if (length(nominal)!=0){

    for (k in nominal){

      databa[,k] = factor(databa[,k])

    }

  } else {}



  ### Y and Z must be nominal or ordinal: conversion

  indY    = which(colnames(databa)== Y)
  indZ    = which(colnames(databa)== Z)


  if (OUT == "Y"){

    outc          = Y
    ordinal       = setdiff(ordinal,indZ)
    quanti        = setdiff(quanti ,indZ)
    nominal       = setdiff(nominal,indZ)

    if (indY %in% quanti){

      databa[,indY] = ordered(databa[,indY])
      ordinal       = sort(c(ordinal,indY))
      quanti        = setdiff(quanti,indY)

    } else {}

  } else if (OUT == "Z"){

    outc          = Z
    ordinal       = setdiff(ordinal,indY)
    quanti        = setdiff(quanti ,indY)
    nominal       = setdiff(nominal,indY)

    if (indY %in% quanti){

      databa[,indZ] = ordered(databa[,indZ])
      ordinal       = sort(c(ordinal,indZ))
      quanti        = setdiff(quanti,indZ)

    } else {}

  } else {

    stop("Improper argument for OUT: Y or Z only")

  }



  ### new indexes without Y and Z indexes

  indNOM  = setdiff(nominal,c(indY,indZ))
  indORD  = setdiff(ordinal,c(indY,indZ))
  indNUM  = setdiff(quanti ,c(indY,indZ))


  convert_nm = convert_num

  ### Turn continuous predictors into ordered categorical factors ?

  if (length(convert_num) != 0){

    tt = 0
    for (k in convert_num){
      tt = tt + 1
      databa[,k]  = cut(databa[,k],breaks = stats::quantile(databa[,k],
                                                            probs = seq(0,1,by = 1/convert_clss[tt]),na.rm = TRUE),
                        include.lowest = TRUE, ordered_result = TRUE)
    }

    indORD      = sort(c(indORD,convert_num))
    indNUM      = sort(setdiff(indNUM,convert_num))
    convert_num = NULL

  } else {}

  if (OUT == "Y"){

    databa = databa[!is.na(databa[,indY]),]

  } else {

    databa = databa[!is.na(databa[,indZ]),]

  }

  ### Pairwise relationships between categorical predictors ordered or not and between the outcome Y and
  ### each categorical covariates

  indCAT  = sort(c(indORD,indNOM))

  if (length(indCAT)!=0){

    bbb     = as.data.frame(databa[,indCAT])

    #for (k in 1:ncol(bbb)){
    #  bbb[,k] = factor(bbb[,k])
    #}

    colnames(bbb) = colnames(databa)[indCAT]

    if (OUT == "Y"){

      datacov = data.frame(Y = factor(databa[,indY]),bbb)

    } else {

      datacov = data.frame(Y = factor(databa[,indZ]),bbb)

    }

    datacov = datacov[!is.na(datacov$Y),]


    N_cat = ncol(datacov)


    # CRAMER'S V

    tab_cor2    = NULL
    kk = 0
    for (j in 1:(N_cat-1)){

      name1       = rep(colnames(datacov)[j],N_cat-kk)
      name2       = colnames(datacov)[j:N_cat]
      tab_cor     = data.frame(name1,name2)
      tab_cor     = tab_cor[as.character(tab_cor$name1) != as.character(tab_cor$name2),]
      kk          = kk + 1

      for (k in (j+1):(N_cat)){

        XX              = datacov[,j]
        vars            = datacov[,k]
        datacov2        = na.omit(data.frame(XX, vars))
        datacov2$vars   = as.factor(as.character(datacov2$vars))
        datacov2$XX     = as.factor(as.character(datacov2$XX))
        regression      = paste(colnames(datacov2)[1], "~", colnames(datacov2)[2],sep="")
        tab_cor[k-kk,3] = round(suppressWarnings(StatMatch::pw.assoc(formula(regression), data = datacov2)$V),4)
        tab_cor[k-kk,4] = round(suppressWarnings(StatMatch::pw.assoc(formula(regression), data = datacov2)$bcV),4)
        tab_cor[k-kk,5] = nrow(datacov2)

      }

      tab_cor2 = rbind(tab_cor2,tab_cor)

    }

    colnames(tab_cor2)[3] = "V_Cramer"
    colnames(tab_cor2)[4] = "CorrV_Cramer"
    colnames(tab_cor2)[5] = "N"


    tab_cor2   = tab_cor2[order(tab_cor2$V_Cramer,decreasing=TRUE),]
    tab_cor2_Y = tab_cor2[as.character(tab_cor2$name1) == "Y",]
    tab_cor2_X = tab_cor2[as.character(tab_cor2$name1) != "Y",]
    tab_cor3_X = tab_cor2_X[tab_cor2_X$V_Cramer> thresh_cat,]

  } else {

    tab_cor2_Y = NULL
    tab_cor3_X = NULL

  }


  ### Pairwise relationships between continuous predictors and between the outcome Y and
  ### each continuous covariates

  indNUM2  = sort(c(indORD,indNUM))

  if (length(indNUM2)!=0){

    if (OUT == "Y"){

      datacov3 = data.frame(Y = as.factor(databa[,indY]),databa[,sort(indNUM)],databa[,sort(indORD)])
      datacov3 = datacov3[!is.na(datacov3$Y),]
      colnames(datacov3) = c("Y",colnames(databa)[sort(indNUM)],colnames(databa)[sort(indORD)])

    } else {

      datacov3 = data.frame(Y = as.factor(databa[,indZ]),databa[,sort(indNUM)],databa[,sort(indORD)])
      datacov3 = datacov3[!is.na(datacov3$Y),]
      colnames(datacov3) = c("Y",colnames(databa)[sort(indNUM)],colnames(databa)[sort(indORD)])

    }

    N_num = ncol(datacov3)


    # Spearman's Rank correlation

    tab_cor4       = NULL
    mat_cor4       = matrix(ncol = N_num,nrow = N_num)
    diag(mat_cor4) = 1
    kk = 0
    for (j in 1:(N_num-1)){

      name1       = rep(colnames(datacov3)[j],N_num-kk)
      name2       = colnames(datacov3)[j:N_num]
      tab_cor     = data.frame(name1,name2)
      tab_cor     = tab_cor[as.character(tab_cor$name1) != as.character(tab_cor$name2),]
      kk          = kk+1

      for (k in (j+1):(N_num)){

        XX              = datacov3[,j]
        vars            = datacov3[,k]
        datacv          = na.omit(data.frame(XX, vars))
        rescor          = suppressWarnings(cor(rank(datacv$XX),rank(datacv$vars)))
        mat_cor4[j,k]   = rescor
        mat_cor4[k,j]   = rescor
        tab_cor[k-kk,3] = round(abs(rescor),4)
        tab_cor[k-kk,4] = round(suppressWarnings(stats::cor.test(rank(datacv$XX),rank(datacv$vars)))$p.value,4)
        tab_cor[k-kk,5] = nrow(datacv)

      }

      tab_cor4 = rbind(tab_cor4,tab_cor)

    }


    colnames(tab_cor4)[3] = "ABS_COR"
    colnames(tab_cor4)[4] = "pv_COR_test"
    colnames(tab_cor4)[5] = "N"
    colnames(mat_cor4)    = colnames(datacov3)
    row.names(mat_cor4)   = colnames(datacov3)

    mat_cor4 = mat_cor4[2:nrow(mat_cor4),2:ncol(mat_cor4)]

    tab_cor4    = tab_cor4[order(tab_cor4$ABS_COR,decreasing=TRUE),]
    tab_cor4_Y  = tab_cor4[as.character(tab_cor4$name1) == "Y",]
    tab_cor4_X  = tab_cor4[as.character(tab_cor4$name1) != "Y",]

    tab_cor5_X = tab_cor4_X[tab_cor4_X$ABS_COR> thresh_num,]

  } else {

    tab_cor4_Y = NULL
    tab_cor5_X = NULL

  }

  if (OUT == "Z"){

    tab_cor2_Y[,1] = rep("Z",nrow(tab_cor2_Y))
    tab_cor4_Y[,1] = rep("Z",nrow(tab_cor4_Y))

  } else {}

  # Farrar and Glauber test

  if (!is.null(ncol(mat_cor4))){

    ddl     = (N_num+1)*N_num/2
    FG_stat = - (nrow(datacov3)-1-(2*(N_num+1)+5)/6)*log(det(mat_cor4))
    pval_FG = 1-stats::pchisq(FG_stat,ddl)

    FG_stat_approx = sqrt(2*FG_stat) - sqrt(2*ddl-1)
    pval_FG_approx = 2*(1-stats::pnorm(abs(FG_stat_approx),0,1))

    FG_synt = c(DET_X = det(mat_cor4), pv_FG_test = pval_FG, pv_FG_test_appr = pval_FG)

  } else {

    FG_synt = NULL

  }

  ### Finding the best predictors of Y using Random Forest

  if (RF == TRUE){

    databa2 = databa[,-ID]

    if (OUT == "Y"){

      databa2 = databa2[,c(which(colnames(databa2)==Y),setdiff(1:ncol(databa2),c(which(colnames(databa2) %in% c(Y,Z)))))]
      colnames(databa2)[1] = "Y"

    } else {

      databa2 = databa2[,c(which(colnames(databa2)==Z),setdiff(1:ncol(databa2),c(which(colnames(databa2)%in% c(Y,Z)))))]
      colnames(databa2)[1] = "Y"

    }

    if (length(indNUM)!=0){

      RF_condi = FALSE

    } else {}


    # Conditional importance variable on complete or imputed databases only

    if (RF_condi == TRUE){

      databa2bis = na.omit(databa2)
      set.seed(RF_SEED); stocres = party::cforest(Y~., data=databa2bis, control = cforest_unbiased(mtry = (ncol(databa2bis)-1)^(1/2), ntree = RF_ntree))

    } else {

      set.seed(RF_SEED); stocres = party::cforest(Y~., data=databa2,control = cforest_unbiased(mtry = (ncol(databa2)-1)^(1/2), ntree = RF_ntree))

    }

    stocimp = sort(party::varimp(stocres,conditional= RF_condi, threshold = RF_condi_thr),decreasing=TRUE)
    nomlist = names(stocimp)


    ### Groups of highly correlated covariates

    candidates     = rbind(tab_cor3_X[,1:2],tab_cor5_X[,1:2])

    if (nrow(candidates)!=0){

      message("Risks of collinearity between predictors detected. Consult outputs for more details ...","\n")

      # } else {}

      candidates[,1] = as.character(candidates[,1])
      candidates[,2] = as.character(candidates[,2])

      candid = list()
      kk = 1

      while (nrow(candidates)!=0){

        candid[[kk]] = as.character(candidates[1,])
        tt = 1

        repeat{

          ncandid      = length(candid[[kk]])


          for (k in 2:nrow(candidates)){

            if (length(intersect(candid[[kk]],candidates[k,]))!=0){

              tt = unique(c(tt,k))
              candid[[kk]] = unique(c(candid[[kk]],as.character(candidates[k,])))

            } else {}

          }

          if (length(candid[[kk]])==ncandid){

            break

          }
        }

        candidates = candidates[-tt,]
        kk = kk + 1

      }

      ### Selection of uncorrelated covariates

      for (kk in 1:length(candid)){

        for (k in 1:length(nomlist)){

          if (nomlist[k] %in% candid[[kk]]){

            candid[[kk]] = setdiff(candid[[kk]],nomlist[k])
            break()

          } else {}

        }

      }

      drop_var   = unlist(candid)
      remain_var = setdiff(nomlist,drop_var)


      ### Best predictors for Y

      new_var = intersect(c("Y",remain_var),colnames(databa2))
      # databa3 = databa2[,intersect(c("Y",remain_var),colnames(databa2))]

      if (RF_condi == TRUE){

        databa3bis = databa2bis[,intersect(c("Y",remain_var),colnames(databa2))]
        set.seed(RF_SEED); stocres = party::cforest(Y~.,data=databa3bis,control = party::cforest_unbiased(mtry = (ncol(databa3bis)-1)^(1/2),ntree = RF_ntree))

      } else {

        databa3 = databa2[,intersect(c("Y",remain_var),colnames(databa2))]
        set.seed(RF_SEED); stocres = party::cforest(Y~.,data=databa3, control = party::cforest_unbiased(mtry = (ncol(databa3)-1)^(1/2),ntree = RF_ntree))

      }

      stocimp2 = sort(party::varimp(stocres,conditional= RF_condi, threshold = RF_condi_thr),decreasing=TRUE)

    } else {

      stocimp2   = stocimp
      drop_var   = NULL
      remain_var = nomlist

    }

    if (min(stocimp2)<0){

      stocimp2 = stocimp2 - min(stocimp2)

    } else {}


    resY = round(prop.table(stocimp2)*100,4)

    resYbis   = cumsum(rev(resY))
    best_pred = rev(names(resYbis)[resYbis>=(thresh_Y*100)])

    out1 = list(seed = RF_SEED, outc = outc, thresh = c(CATEG = thresh_cat, NUM = thresh_num, RF_IMP = thresh_Y), convert_num = colnames(databa)[convert_nm], DB_USED = databa,
                vcrm_OUTC_cat = tab_cor2_Y, cor_OUTC_num = tab_cor4_Y, vcrm_X_cat = tab_cor2_X, cor_X_num = tab_cor4_X, FG_test = FG_synt,
                collinear_PB = list(VCRAM = tab_cor3_X,SPEARM = tab_cor5_X),
                drop_var = drop_var, RF_PRED = resY, RF_best = best_pred)

  } else {

    out1 = list(seed = RF_SEED, outc = outc, thresh = c(CATEG = thresh_cat, NUM = thresh_num, RF_IMP = thresh_Y), convert_num = colnames(databa)[convert_nm], DB_USED = databa,
                vcrm_OUTC_cat = tab_cor2_Y, cor_OUTC_num = tab_cor4_Y, vcrm_X_cat = tab_cor2_X, cor_X_num = tab_cor4_X, FG_test = FG_synt,
                collinear_PB = list(VCRAM = tab_cor3_X,SPEARM = tab_cor5_X))

  }

  return(out1)

}








