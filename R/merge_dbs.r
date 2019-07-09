
#' merge_dbs()
#'
#' A function which merge vertically two databases by selecting the common covariates with a variable Y encoded
#  in a given scale in DB1 (Y1) and another scale in DB2 (Y2)
#'
#' @param DB1 A data.frame corresponding to the 1st DB to match
#' @param DB2 A data.frame corresponding to the 2nd DB to match
#' @param NAME_Y1 Name of the outcome (with quotes) in its specific scale/encoding gathered from the 1st DB
#' @param NAME_Y2 Name of the outcome (with quotes) in its specific scale/encoding gathered from the 2nd DB
#' @param order_levels_Y1 To complete only if Y is considered as an ordinal factor (scale). A vector of labels of levels (with quotes) sorted in ascending order in DB1
#' @param order_levels_Y2 To complete only if Y is considered as an ordinal factor.A vector of labels of levels sorted in ascending order in DB2 with different scale from DB1
#' @param ordinal_DB1 Vector of number of columns corresponding to ordinal variables in DB1 (No variable by default)
#' @param ordinal_DB2 Vector of number of columns corresponding to ordinal variables in DB2 (No variable by default)
#' @param impute NO imputation on covariates (By default), CC (for Complete Case), MICE (For MICE multiple imputation), MDA (for single imputation using Factorial Analysis for Mixed Data)
#' @param R_MICE Number of multiple imputations for MICE only (5 by default)
#' @param NCP_MDA Integer corresponding to the number of components used to predict NA in MDA imputation (3 by default)
#' @param seed_func Integer used as argument by the set.seed() for offsetting the random number generator (Random integer by default)
#'
#' @return A list containing 10 elements (11 if impute = "MICE"):
#'
#'         DB_READY The database matched from the 2 initial BDDs with common covariates and imputed or not accordind to the impute option.
#'         Y1_LEVELS Levels retained for the target variable in the DB1.
#'         Y2_LEVELS Levels retained for the target variable in the DB2.
#'         REMOVE1 Labels of deleted covariates because of their different types from DB1 to DB2
#'         REMOVE2 Removed factor(s) because of their different levels from DB1 to DB2
#'         REMAINING_VAR Labels of the covaraites remained for the fusion with OT algorithm
#'         IMPUTE_TYPE Method chosen to handle missing data in covariates
#'         MICE_DETAILS A list containing the details of the MICE imputation. Databases imputed for DB1 and DB2 according to the number of mutliple imputation selected (Only if impute = "MICE").
#'         DB1_RAW A data.frame corresponding to the 1st raw database
#'         DB2_RAW A data.frame corresponding to the 2nd raw database
#'         SEED An integer used as argument by the set.seed() for offsetting the random number generator (random selection by default)
#'
#' @export
#'
#' @examples
#' # Require samp.A database from the StatMatch package. c.neti and c.neti.bis coded voluntarily in 2 distinct encodings.
#' library(StatMatch)
#' data(samp.A)
#' samp.A = samp.A[,c(1:11,13,12)]
#' c.neti            = as.numeric(samp.A$c.neti)
#'
#' samp.A$c.neti.bis = as.factor(ifelse(c.neti %in% c(1,2),1,
#'                                     ifelse(c.neti %in% c(3,4),2,
#'                                            ifelse(c.neti %in% c(5,6),3,4))))
#' data1 = samp.A[1:1000,c(2:9,13)]
#' data2 =samp.A[1001:nrow(samp.A),c(5:11,12,14)]
#'
#' # Insert the variable marital in 2 different types:
#' data1$marital = as.numeric(data1$marital)
#'
#' # Insert different levels in a factor variable:
#' data2$c.age = as.character(data2$c.age)
#' data2$c.age[data2$c.age %in% c("[16,34]","(34,44]")] = "[16,44]"
#' data2$c.age = as.factor(data2$c.age)
#'
#' # Add NA in covariates:
#' add_NA = function(DB,tx){
#' DB_NA = DB
#' for (j in 1:ncol(DB)){
#'    NA_indic = sample(1:nrow(DB),round(nrow(DB)*tx/100),replace=FALSE)
#'    DB_NA[NA_indic,j] = rep(NA,length(NA_indic))
#'  }
#'  return(DB_NA)
#'}
#'
#' set.seed(4036); data3 = add_NA(data1,10); data4 = add_NA(data2,10)
#' soluc1  = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "MICE",R_MICE = 2, seed_func = 4036)
#' soluc2  = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
#' soluc3  = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "CC")
#' soluc4  = merge_dbs(data3,data4,"c.neti","c.neti.bis", order_levels_Y1 = NULL, order_levels_Y2 = NULL, ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "MDA", NCP_MDA = 3)
#' soluc5  = merge_dbs(data1,data2,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "MICE")
#' soluc6  = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",order_levels_Y2 = c("4","1","3","2"),ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
#'
#'
merge_dbs = function(DB1,DB2,NAME_Y1,NAME_Y2, order_levels_Y1 = levels(DB1[,NAME_Y1]), order_levels_Y2 = levels(DB2[,NAME_Y2]),
                     ordinal_DB1 = FALSE, ordinal_DB2 = FALSE, impute = "NO",R_MICE = 5,
                     NCP_MDA = 3, seed_func = sample(1:1000000,1)){


  cat("DBS MERGING in progress. Please wait ...","\n")

  DB1_raw= DB1
  DB2_raw= DB2


  # Constraints on the 2 databases
  #--------------------------------

  if ((!is.data.frame(DB1))|(!is.data.frame(DB2))){

    stop("At least one of your two objects is not a data.frame!")

  } else {}


  if ((is.character(NAME_Y1)!=TRUE)|(is.character(NAME_Y2)!=TRUE)){

    stop ("NAME_Y1 and NAME_Y2 must be declared as strings of characters with quotes")

  } else {}


  if ((length(NAME_Y1)!=1)|(length(NAME_Y2)!=1)){

    stop ("More than one target declared by DB !")

  } else {}




  # Remove subjects in DB1 (resp. DB2) with Y1 (resp.Y2) missing

  NB1_1 = nrow(DB1); NB2_1 = nrow(DB2);

  DB1   = DB1[!is.na(DB1[,NAME_Y1]),]
  DB2   = DB2[!is.na(DB2[,NAME_Y2]),]

  NB1_2 = nrow(DB1); NB2_2 = nrow(DB2);


  REMOVE_SUBJECT1 = NB1_1 - NB1_2
  REMOVE_SUBJECT2 = NB2_1 - NB2_2




  # Covariates selection
  #----------------------

  # Impute NA in covariates from DB1 and DB2 if necessary


  count_lev1  = apply(DB1,2,function(x){length(names(table(x)))})
  num_lev1    = sapply(DB1,is.numeric)

  typ_cov_DB1 = ifelse (((1:ncol(DB1)) %in% ordinal_DB1) & (count_lev1 > 2),"polr",
                        ifelse ((count_lev1 == 2) & (num_lev1 == FALSE),"logreg",
                                ifelse (num_lev1 == TRUE,"pmm","polyreg")))


  count_lev2  = apply(DB2,2,function(x){length(names(table(x)))})
  num_lev2    = sapply(DB2,is.numeric)

  typ_cov_DB2 = ifelse (((1:ncol(DB2)) %in% ordinal_DB2) & (count_lev2 > 2),"polr",
                        ifelse ((count_lev2 == 2) & (num_lev2 == FALSE),"logreg",
                                ifelse (num_lev2 == TRUE,"pmm","polyreg")))

  typ_cov_DB1b = typ_cov_DB1[setdiff(names(DB1),NAME_Y1)]
  typ_cov_DB2b = typ_cov_DB2[setdiff(names(DB2),NAME_Y2)]


  if ((setequal(DB1,na.omit(DB1)))&(setequal(DB2,na.omit(DB2)))&(impute!="NO")){


    cat("No missing values in covariates: No imputation methods required","\n")
    impute = "NO"

  } else {}


  if (impute == "CC"){

    DB1 = na.omit(DB1); DB2 = na.omit(DB2);
    DB1_new = na.omit(DB1[,setdiff(names(DB1),NAME_Y1)])
    DB2_new = na.omit(DB2[,setdiff(names(DB2),NAME_Y2)])

  } else if (impute == "MICE"){


    DB1bis       = imput_cov(DB1[,setdiff(names(DB1),NAME_Y1)],R_mice = R_MICE,meth = typ_cov_DB1b, missMDA = FALSE, seed_choice = seed_func)
    DB2bis       = imput_cov(DB2[,setdiff(names(DB2),NAME_Y2)],R_mice = R_MICE,meth = typ_cov_DB2b, missMDA = FALSE, seed_choice = seed_func)


  } else if (impute == "MDA"){

    DB1bis       = imput_cov(DB1[,setdiff(names(DB1),NAME_Y1)],meth = typ_cov_DB1b, missMDA = TRUE, NB_COMP = NCP_MDA, seed_choice = seed_func)
    DB2bis       = imput_cov(DB2[,setdiff(names(DB2),NAME_Y2)],meth = typ_cov_DB2b, missMDA = TRUE, NB_COMP = NCP_MDA, seed_choice = seed_func)


  } else if (impute == "NO"){

    DB1_new = DB1[,setdiff(names(DB1),NAME_Y1)]
    DB2_new = DB2[,setdiff(names(DB2),NAME_Y2)]

  } else {stop("The specification of the impute option is false: NO, CC, MICE, MDA only")}


  # Transform Y
  #-------------

  cat("Y1","\n")
  Y1            = transfo_target(DB1[,NAME_Y1], levels_order = order_levels_Y1)$NEW
  cat("Y2","\n")
  Y2            = transfo_target(DB2[,NAME_Y2], levels_order = order_levels_Y2)$NEW

  if (setequal(levels(Y1),levels(Y2))){

    stop("Your target has identical labels in the 2 databases !")

  } else {}




  # Y1 and Y2 extracted from DB1 and DB2

  if (impute %in% c("MICE","MDA")){

      DB1_new = DB1bis$DATA_IMPUTE
      DB2_new = DB2bis$DATA_IMPUTE

  } else {}


  # Covariates re-ordered by their names in each DB

  DB1_new          = DB1_new[,order(names(DB1_new))]
  DB2_new          = DB2_new[,order(names(DB2_new))]


  # Names of the common covariates between the two DB

  same_cov      = intersect(names(DB1_new),names(DB2_new))
  n_col         = length(same_cov)


  # Remaining common covariates in each DB

  DB1_4FUS     = DB1_new[,same_cov]
  DB2_4FUS     = DB2_new[,same_cov]


  # Removed covariate(s) because of their different types from DB1 to DB2

  list1   = as.list(lapply(DB1_4FUS,class))
  list1   = lapply(list1,paste,collapse=" ")
  l1      = unlist(list1)

  list2   = as.list(lapply(DB2_4FUS,class))
  list2   = lapply(list2,paste,collapse=" ")
  l2      = unlist(list2)


  ind_typ       = (l1!= l2)

  if (sum(ind_typ)!=0){

    remove_var1 = same_cov[ind_typ]

  } else {remove_var1 = NULL}



  # Removed factor(s) because of their different levels from DB1 to DB2

  modif_factor = (1:n_col)[l1 %in% c("factor","ordered factor")]

  same_cov2 = same_cov[modif_factor]

  levels_DB1 = sapply(DB1_4FUS[,modif_factor],levels)
  levels_DB2 = sapply(DB2_4FUS[,modif_factor],levels)


  ind_fac = compare_lists(levels_DB1,levels_DB2)

  if (sum(ind_fac)!=0){

    remove_var2 = same_cov2[ind_fac]

  } else {remove_var2 = NULL}


  # Names of variables removed before merging

  remove_var = c(remove_var1,remove_var2)


  # Names of variables remained before merging

  remain_var = setdiff(same_cov,remove_var)

  if (length(remain_var) == 0){

    stop("no common variable selected in the 2 DBS except the target !")

  } else {}


  DB1_4FUS = DB1_4FUS[,remain_var]
  DB2_4FUS = DB2_4FUS[,remain_var]

  Y2b = sample(Y2,length(Y1),replace=T)

  DB_COV = rbind(

    data.frame(DB = rep(1,nrow(DB1_4FUS)),Y1,Y2 = Y2b,DB1_4FUS),
    data.frame(DB = rep(2,nrow(DB2_4FUS)),Y1 = rep(NA,nrow(DB2_4FUS)),Y2,DB2_4FUS)

  )

  DB_COV$Y2[1:nrow(DB1_4FUS)] = rep(NA,nrow(DB1_4FUS))

  cat("DBS MERGING OK","\n")
  cat(rep("-",23),sep=""); cat("\n"); cat("SUMMARY OF DBS MERGING:","\n")
  cat("Nb of removed subjects because of NA on Y:",REMOVE_SUBJECT1 + REMOVE_SUBJECT2,
      "(",round((REMOVE_SUBJECT1 + REMOVE_SUBJECT2)*100/(NB1_1 + NB2_1)),"%)","\n")
  cat("Nb of removed covariates because of their different types:",length(remove_var),"\n")
  cat("Nb of remained covariates:",ncol(DB_COV)-3,"\n")
  cat("More details in output ...","\n")

  if (impute %in% c("NO","CC","MDA")){

    return(list(DB_READY = DB_COV,Y1_LEVELS = levels(DB_COV$Y1), Y2_LEVELS = levels(DB_COV$Y2),
                REMOVE1 = remove_var1, REMOVE2 = remove_var2, REMAINING_VAR = colnames(DB_COV)[4:ncol(DB_COV)],
                IMPUTE_TYPE = impute, DB1_RAW = DB1_raw, DB2_RAW = DB2_raw, SEED = seed_func))

  } else if (impute == "MICE"){

    return(list(DB_READY = DB_COV,Y1_LEVELS = levels(DB_COV$Y1), Y2_LEVELS = levels(DB_COV$Y2),
                REMOVE1 = remove_var1, REMOVE2 = remove_var2, REMAINING_VAR = colnames(DB_COV)[4:ncol(DB_COV)],
                IMPUTE_TYPE = impute,
                MICE_DETAILS = list(DB1 = list(RAW = DB1bis[[1]],LIST_IMPS = DB1bis[[4]]),
                                    DB2 = list(RAW = DB2bis[[1]],LIST_IMPS = DB2bis[[4]])),
                DB1_RAW = DB1_raw, DB2_RAW = DB2_raw, SEED = seed_func))

  }

}





