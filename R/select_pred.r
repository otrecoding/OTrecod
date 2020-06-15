#' select_pred()
#'
#' Selection of a subset of non collinear predictors having relevant relationships with a given target outcome.
#'
#' The \code{select_pred} function provides several tools to identify, on the one hand, the relationships between predictors, by detecting especially potential problems of collinearity, and, on the other hand, proposes a parcimonious subset of relevant predictors using appropriate random forest procedures.
#' Even if, in many contexts, this function can be seen as a preliminary step of regression, this function is particularly adapted to the context of data integration of two or more datasources by pre-selecting relevant subsets of shared predictors before fusion.
#'
#'
#' A. STRUCTURE OF DATABASE REQUIRED
#'
#' The input expected database is a data.frame that especially requires a specific column of rows identifiers and a target variable (or outcome) having a finite number of values or classes (Therefore,in preference, of ordinal, nominal or discrete type). If the chosen outcome is in numeric form, it will be automatically converted in ordinal type.
#' The number of predictors is not a constraint (even if, with less than 3 variables a process of variables selection has no real sense...), can exceed the number of rows (no problem of high dimensionality here).
#' The predictors can be continuous (quantitative), nominal or ordinal with or without missing values. Nevertheless, this function provides optimal performances without or with only a limited number of continuous variables and with complete information. For this reason, this function is particularly adapted to
#' data integration using the \code{OT_joint} function available in this package, which imposes the same conditions.
#' In presence of numeric variables, it is therefore suggested to discretize them beforehand or to use the process of discretization directly integrated in the function. For doing this, \code{convert_num} and \code{convert_class} are the only two arguments that have to be strictly completed. With the \code{convert_num} argument, users choose the
#' continuous variables to convert and the \code{convert_clss} argument specify the related number of classes chosen for each discretization.
#' So, these 2 arguments must be 2 vectors of indexes of same length. The only exception allowed, when \code{convert_clss} equals a scalar S, supposed that the continuous predictors selected for conversion  will be disretized with a same number of classes S.
#' By example, if \code{convert_class} = 4, all the continuous variables selected with the \code{convert_num} argument will be discretized by quartiles. Finally, missing informations from a predictor to convert, are not discretized and stay unchanged during the transformation.
#'
#'
#' B. PAIRWISE ASSOCIATIONS BETWEEN PREDICTORS
#'
#' In a first step of process, \code{select_pred} calculates standard pairwise associations between predictors according to their types.
#' \enumerate{
#' \item{Between categorical predictors (ordinal, nominal and logical):
#' Cramer's V (and Bias-corrected Cramer's V, see Bergsma 2013, for more details) are calculated between categorical predictors and the \code{thres_cat} argument fixed the associated threshold beyond which two predictors can be considered as redundant.
#' A similar process is done between the target variable and the subset of categorical variables which provides in output, a first table ranking the top scoring predictor variates which summarizes their abilities to predict the target.}
#' \item{Between continuous predictors:
#' If the\code{ordinal} and \code{logic} arguments differ from NULL, all the corresponding predictors are beforehand converted in discrete form.
#' For numeric (quantitative), logical and ordinal predictors, pairwise correlations between ranks (Spearman) are calculated and the \code{thresh_num} arguments fixed the associated threshold beyond which two predictors can be considered as redundant.
#' A similar process is done between the outcome and the subset of discrete variables which provides in output, a table ranking the top scoring predictor variates which summarizes their abilities to predict the target.
#' In addition, the result of a Farrar and Glauber test is provided (see the corresponding reference for more details). This test is based on the determinant of the correlation matrix of covariates and the related null hypothesis of the test corresponds to an absence of collinearity between them.
#' In presence of a large number of numeric covariates and/or ordered factors, the approximate Farrar-Glauber test, based on the normal approximation of the null distribution is more adapted (Johnson, 1994) and its result is also provided in output.
#' These 2 tests are highly sensitive and, by consequence, it suggested to consider these results as simple indicators of collinearity between predictors rather than an essential condition of acceptability.}
#' }
#' If the initial number of predictors is not too important, these informations can be sufficient to users for the visualization of the potential problems of collinearity and for the selection of a subset of predictors (\code{RF} = FALSE).
#' It is nevertheless often necessary to complete this visualization by an automatical process of selection like the Random Forest approach (see Breiman 2001, for more details about RF).
#'
#'
#' C. RANDOM FOREST PROCEDURE
#'
#' As a final step of the process, a random forest approach (RF) is here prefered (to regressions) for its great flexibility: A non parametric recursive partitioning method applyable even when the number of variables exceeds the number of rows, whatever the types of covariates considered, able to deal with non linearity,
#' correlated predictors, ordered or not ordered categorical outcomes. For this, \code{select_pred} integrates in its algorithm the \code{\link[party]{cforest}} and \code{\link[party]{varimp}} of the package \pkg{party} (Hothorn, 2006) and so gives access
#' to its main arguments.
#'
#' A RF approach generally provides 2 types of measures for estimating the mean variable importance of each covariate in the prediction of an outcome: The Gini importance and the permutation importance.These measurements must be used with caution, by taking into account the following constraints:
#' \enumerate{
#' \item{The Gini importance criterion can produce bias in favor of continuous variables and variables with many categories. To avoid this problem, only the permutation criterion is available in the function.}
#' \item{The permutation importance criterion can overestimate the importance of highly correlated predictors.}
#' }
#' The function \code{select_pred} can so proceeds with 3 different scenarios according to the types of predictors:
#' \enumerate{
#' \item{The 1st scenario consists in boiling down to a set of categorical variables (ordered or not) by discretizing all the continuous predictors beforehand, using the internal \code{convert_num} argument or another one, and then works with the conditional importance measurements (\code{RF_condi} = TRUE) which give unbiased estimations.
#' In the spirit of a partial correlation, conditional importance measure of a variable X, only use a more or less important part of predictors (called Z) correlated to X for the computation of the estimation. The argument \code{RF_condi_thr} that corresponds exactly to the \code{threshold} argument of the function \code{\link[party]{varimp}},
#' fixes this ratio according to an association test between X and Z (see Strobl 2007 for more precision). It corresponds to 1 - pvalue of the test, kwnowing that a threshold value of zero will include all predictors in the computation.
#' Moreover, taking into account only subsets of predictors in the computation of the variable importance measures could lead to a relevant saving of execution time.
#' Nevertheless, this method does not take into account incomplete predictors an so have to be imputed beforehand if necessary (using the \code{imput_cov} function available in this package or another one).}
#' \item{The 2nd possibility, always in presence of mixed types predictors, consists in the execution of 2 successive RF procedures. The 1st one is used to select an unique candidate in each susbset of correlated predictors (detecting in the first section), and the 2nd one extracts the permutation importance
#' measurements from a second RF applies on the remaining subset of uncorrelated predictors (\code{RF_condi} = FALSE, default option). This 2nd possibility has the advantage to work in presence of incomplete predictors.}
#' \item{The 3rd scenario consists in running a first time the function without RF process (\code{RF} = FALSE), and according to the presence of highly correlated predictors or not, makes a decision about the type of criterion used as a 2nd execution of the function.}
#' }
#' The 3 scenarios finally provide a list of uncorrelated predictors of the outcome sorted in unbiased importance order and the argument \code{thresh_Y} corresponds to the minimal percent of importance required for a variable to be considered as a reliable predictor of the outcome.
#' Finally, because all random forest results are subject to random variation, users can check whether the same importance ranking is achieved with a different random seed (\code{RF_SEED} argument) or otherwise by increasing the number of trees (\code{RF_ntree}).
#'
#'
#' @param databa A data.frame with a column of identifiers, an outcome, and a set of predictors. The number of columns can exceed the number of rows.
#' @param Y The label of the outcome with quotes
#' @param ID The index of column of the row identifiers
#' @param quanti A vector of integers which corresponds to the indexes of columns of all the numeric predictors
#' @param nominal A vector of integers which corresponds to the indexes of columns of all the categorical nominal predictors
#' @param ordinal A vector of integers which corresponds to the indexes of columns of all the categorical ordinal predictors
#' @param logic A vector of integers indicating the indexes of logical predictors. No index remained by default
#' @param convert_num A vector of integers indicating the indexes of quantitative variables to convert in ordered factors. No index remained by default. Each index selected must has been defined as quantitative in the quanti argument.
#' @param convert_clss A vector of integers indicating the number of classes chosen for each transformation of quantitative variable in ordered factor. The length of this vector can not exceed the length of the convert_num argument. If length(convert_num) > 1 and length(convert_clss) = 1.
#' All quantitative predictors selected for discretisation will have the same number of classes.
#' @param thresh_cat A threshold associated to the Cramer's V coefficient (= 0.30 by default)
#' @param thresh_num A threshold associated to the Spearman's coefficient of correlation (= 0.70 by default)
#' @param thresh_Y A threshold indicating the minimal percent of importance required for a variable to be considered as a reliable predictor
#' @param RF A boolean that equals TRUE (default) if a random forest procedure must be used in the selecting of the best subset of predictors for the outcome
#' @param RF_ntree Number of bootsrap samples required from the row datasource during the random forest procedure
#' @param RF_condi A boolean specifying if the conditional importance measures must be assessed from the random forest procedure (\code{TRUE}) rather than the standard variable importance  measures (\code{FALSE} by default)
#' @param RF_condi_thr A threshold linked to (1 - pvalue) of an association test, kwnowing that a threshold value of zero wil include all covariates.
#' @param RF_SEED An integer used as argument by the set.seed() for offsetting the random number generator (Random integer by default)
#'
#'
#' @return A list of 13 (If \code{RF} = TRUE) or 10 objects (Only the 1st ten objects if \code{RF} = FALSE) is returned:
#' \item{seed}{The random number generator fixed or selected}
#' \item{thresh}{A summarize of the different thresholds fixed for the study}
#' \item{convert_num}{The labels of the continuous predictors transformed in categorical form}
#' \item{DB_USED}{The final database used after eventual transformations of predictors}
#' \item{vcrm_Y_cat}{Table of pairwise associations between the outcome and the categorical predictors (Cramer's V)}
#' \item{cor_Y_num}{Table of pairwise associations between the outcome and the continuous predictors (Rank correlation)}
#' \item{v_crm_X_cat}{Table of pairwise associations between the categorical predictors (Cramer's V)}
#' \item{cor_X_num}{Table of pairwise associations between the continuous predictors (Cramer's V)}
#' \item{FG_test}{Results of the Farrar and Glauber tests, with and without approximation form}
#' \item{colinear_PB}{Table of predictors with problem of collinearity}
#' \item{drop_var}{Labels of predictors to drop after RF process (optional output: Only if \code{RF}=TRUE)}
#' \item{RF_PRED_Y}{Table of variable importance measurements, conditional or not, according to the argument \code{condi_RF} (optional output: Only if \code{RF}=TRUE)}
#' \item{best_pred}{Labels of the best predictors selected (optional output: Only if \code{RF}=TRUE)}
#'
#' @export
#'
#' @references
#' # About the Cramer's V measures:
#' Bergsma W (2013). A bias-correction for Cramer's V and Tschuprow's T. Journal of the Korean Statistical Society, 42, 323–328.
#'
#' # About the Farrar and Glauber test:
#' Farrar D, and Glauber R (1968). Multicolinearity in regression analysis. Review of Economics and Statistics, 49, 92–107.
#'
#' # About RF:
#' Breiman L (2001). Random Forests. Machine Learning, 45(1), 5–32.
#'
#' About the estimation of conditional variable importances in RF:
#' Strobl C, Boulesteix A-L, Kneib T, Augustin T, Zeileis A (2008). Conditional Variable Importance for Random Forests. BMC Bioinformatics, 9, 307.
#'
#'
#' @author Gregory Guernec
#' \email{otrecod.pkg@@gmail.com}
#'
#' @importFrom stats cor cor.test pnorm pchisq na.omit
#' @importFrom StatMatch pw.assoc
#' @importFrom party varimp cforest cforest_unbiased
#'
#'
#' @aliases select_pred
#'
#' @examples
#' # Using the simu_data table, firstly, we separate the 2 bases A and B
#' # and remove the unknown outcome in each DB
#' data(simu_data)
#' simu_A = simu_data[simu_data$DB == "A",-3]    # Base A
#' simu_B = simu_data[simu_data$DB == "B",-2]    # Base B
#'
#' # In the 1st DB(A), the target variable is Yb1.
#' # The indexes of the column covariates goes from 3 to 7.
#' # By keeping the other default parameters:
#'
#' \dontrun{
#' ### Scenario 1: Discretization of age (7) in 3 classes (so no quantitative predictors),
#' # RF process chosen, and conditional importance measures assessed, all other arguments
#' # take their default values:
#' sel_A2 = select_pred(simu_A,"Yb1", quanti = 7, nominal = c(3,4,6), ordinal = c(2,5),
#'                      convert_num = 7,convert_clss = 3,
#'                      RF = TRUE, RF_condi = TRUE)
#'
#'
#' ### Scenario 2: Raw predictors, RF process chosen, and standard permutation
#' # importance measures assessed (Default values for other arguments)
#' sel_A3 = select_pred(simu_A,"Yb1", quanti = 7, nominal = c(3,4,6), ordinal = c(2,5),
#'                      RF = TRUE)
#' }
#'
#' ### Scenario 3: Discretization of age (7) in 3 classes and no random forest process:
#' sel_A1 = select_pred(simu_A,"Yb1", quanti = 7, nominal = c(3,4,6), ordinal = c(2,5),
#'                      convert_num = 7,convert_clss = 3,RF = FALSE)
#'
#'
select_pred = function(databa,Y, ID = 1, quanti = NULL, nominal = NULL,ordinal = NULL,logic = NULL,
                       convert_num = NULL, convert_clss = NULL,
                       thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
                       RF = TRUE, RF_ntree = 500, RF_condi = FALSE, RF_condi_thr = 0.2, RF_SEED = sample(1:1000000, 1)){

  # Exclude ID from arguments

  quanti  = setdiff(quanti,ID)
  nominal = setdiff(nominal,ID)
  ordinal = setdiff(ordinal,ID)
  logic   = setdiff(logic,ID)


  # Convert boolean covariates

  if (length(logic) != 0){

    # Verification

    if (all(apply(databa[,logic],2,is.logical))!=TRUE){

      stop("At least one logical variable is not boolean")

    } else {}

    # Conversion

    for (j in logic){

      databa[,j]   = as.ordered(databa[,j])
      ordinal      = unique(sort(c(ordinal,logic)))
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
  stopifnot(!is.character(ordinal)); stopifnot(!is.character(quanti)); stopifnot(!is.character(nominal)); stopifnot(!is.character(logic))
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


  cat("The select_pred function is in progress. Please wait ...","\n")


  # Exclude systematically ID and Y before discretization

  indexY    = (1:ncol(databa))[colnames(databa)== Y]
  index_DB_Y = c(ID,indexY)

  convert_clss = convert_clss[convert_num != index_DB_Y]
  convert_num  = setdiff(convert_num,index_DB_Y)


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



  ### Y must be nominal or ordinal: conversion

  indY    = which(colnames(databa)==Y)

  if (indY %in% quanti){

    databa[,indY] = ordered(databa[,indY])
    ordinal       = sort(c(ordinal,indY))
    quanti        = setdiff(quanti,indY)

  }


  ### new indexes without Y index

  indNOM  = setdiff(nominal,indY)
  indORD  = setdiff(ordinal,indY)
  indNUM  = setdiff(quanti,indY)


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


  ### Pairwise relationships between categorical predictors ordered or not and between the outcome Y and
  ### each categorical covariates

  indCAT  = sort(c(indORD,indNOM))

  if (length(indCAT)!=0){

    bbb     = as.data.frame(databa[,indCAT])

    #for (k in 1:ncol(bbb)){
    #  bbb[,k] = factor(bbb[,k])
    #}

    colnames(bbb) = colnames(databa)[indCAT]

    datacov = data.frame(Y = factor(databa[,indY]),bbb)
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
        regression      = paste(colnames(datacov2)[1], " ~ ", colnames(datacov2)[2])
        tab_cor[k-kk,3] = round(suppressWarnings(StatMatch::pw.assoc(regression, data = datacov2)$V),4)
        tab_cor[k-kk,4] = round(suppressWarnings(StatMatch::pw.assoc(regression, data = datacov2)$bcV),4)
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


    datacov3 = data.frame(Y = as.factor(databa[,indY]),databa[,sort(indNUM)],databa[,sort(indORD)])
    datacov3 = datacov3[!is.na(datacov3$Y),]
    colnames(datacov3) = c("Y",colnames(databa)[sort(indNUM)],colnames(databa)[sort(indORD)])


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

  # Farrar, D., and R. Glauber(1968) : “Multicolinearity in regression analysis,”Review of Economics and Statistics, 49, 92–107.
  # Ref package party: Hothorn, T., Hornik, K., Zeileis, A., 2012. Party: a laboratory for recursive partytioning. R package version 10-3, URL http://cranr-projectorg/package=party.

  ### Finding the best predictors of Y using Random Forest

  if (RF == TRUE){

    databa2 = databa[,-ID]
    databa2 = databa2[,c(which(colnames(databa2)==Y),setdiff(1:ncol(databa2),which(colnames(databa2)==Y)))]
    colnames(databa2)[1] = "Y"

    if (length(indNUM)!=0){

      RF_condi = FALSE

    } else {}


    # Conditional importance variable on complete or imputed databases only

    if (RF_condi == TRUE){

      databa2bis = na.omit(databa2)
      set.seed(RF_SEED); stocres = party::cforest(Y~., data=databa2bis, control = cforest_unbiased(mtry = 5, ntree = RF_ntree))

    } else {

      set.seed(RF_SEED); stocres = party::cforest(Y~., data=databa2,control = cforest_unbiased(mtry = 5, ntree = RF_ntree))

    }

    stocimp = sort(party::varimp(stocres,conditional= RF_condi, threshold = RF_condi_thr),decreasing=TRUE)



    ### Groups of highly correlated covariates

    candidates     = rbind(tab_cor3_X[,1:2],tab_cor5_X[,1:2])

    if (length(ncol(candidates)!=0)){

      cat("Risks of colinearity between predictors detected. Consult outputs for more details ...","\n")

    } else {}

    candidates[,1] = as.character(candidates[,1])
    candidates[,2] = as.character(candidates[,2])
    nomlist        = names(stocimp)

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
      set.seed(RF_SEED); stocres = party::cforest(Y~.,data=databa3bis,control = party::cforest_unbiased(mtry = 5,ntree = RF_ntree))

    } else {

      databa3 = databa2[,intersect(c("Y",remain_var),colnames(databa2))]
      set.seed(RF_SEED); stocres = party::cforest(Y~.,data=databa3, control = party::cforest_unbiased(mtry = 5,ntree = RF_ntree))

    }

    stocimp2 = sort(party::varimp(stocres,conditional= RF_condi, threshold = RF_condi_thr),decreasing=TRUE)

    if (min(stocimp2)<0){

      stocimp2 = stocimp2 - min(stocimp2)

    } else {}


    resY = round(prop.table(stocimp2)*100,4)

    resYbis   = cumsum(rev(resY))
    best_pred = rev(names(resYbis)[resYbis>=(thresh_Y*100)])

   return(list(seed = RF_SEED, thresh = c(CATEG = thresh_cat, NUM = thresh_num, RF_IMP = thresh_Y), convert_num = colnames(databa)[convert_nm], DB_USED = databa,
               vcrm_Y_cat = tab_cor2_Y, cor_Y_num = tab_cor4_Y, v_crm_X_cat = tab_cor3_X, cor_X_num = tab_cor5_X, FG_test = FG_synt,
               colinear_PB = list(VCRAM = tab_cor3_X,SPEARM = tab_cor5_X),
               drop_var = drop_var, RF_PRED_Y = resY, best_pred = best_pred))

  } else {

    return(list(seed = RF_SEED, thresh = c(CATEG = thresh_cat, NUM = thresh_num, RF_IMP = thresh_Y), convert_num = colnames(databa)[convert_nm], DB_USED = databa,
                vcrm_Y_cat = tab_cor2_Y, cor_Y_num = tab_cor4_Y, v_crm_X_cat = tab_cor3_X, cor_X_num = tab_cor5_X, FG_test = FG_synt,
                colinear_PB = list(VCRAM = tab_cor3_X,SPEARM = tab_cor5_X)))

   }

}




