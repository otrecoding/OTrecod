#' prep_dbs()
#'
#' This function formats the variables considered as covariates for the use of the OT algorithm by presenting the nominal (unordered) variables
#' in the form of disjunctive tables and the ordinal variables as discrete numerous variables.
#'
#' @param DB A data.frame in a specific form
#' @param nominal Vector of index of columns corresponding to nominal variables
#' @param ordinal Vector of index of columns corresponding to ordinal variables
#' @param logic   Vector of index of columns corresponding to boolean variables
#'
#' @return A data.frame with only numeric values
#'
#' @export
#' @author Gregory Guernec
#' \email{gregory.guernec@@inserm.fr}
#'
#' @export
#'
#' @examples
#' # Using the previoulsy object soluc1 previously generated
#' soluc1  = merge_dbs(data3,data4,
#' NAME_Y = "c.neti",NAME_Z = "c.neti.bis",
#' ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9),
#' impute = "MICE",R_MICE = 2, seed_func = 4036)
#'
#' summary(soluc1$DB_READY)
#'
#' Discretize the continuous covariates if necessary, like "age" in this example
#' soluc1$DB_READY$age    = cut(soluc1$DB_READY$age,c(0,30,60,100))
#' summary(soluc1$DB_READY)
#' try1                   = prep_dbs(soluc1$DB_READY,nominal = 6,ordinal = 1:5)
#'
#' # A more explicit example included new variables of different types
#' new_DB        = soluc1$DB_READY
#' new_DB$age    = cut(new_DB$age,c(0,30,60,100))
#' new_DB$test   = sample(c(T,F),nrow(new_DB),replace=TRUE)
#' new_DB$multi1 = sample(c("A","B","C"),nrow(new_DB),replace=TRUE)
#' new_DB$multi2 = sample(c("FF","GG"),nrow(new_DB),replace=TRUE)
#'
#' try2                 = prep_dbs(new_DB,nominal = c(6,8:9),ordinal = 1:5,logic =7)

prep_dbs = function(DB,nominal = NULL,ordinal = NULL,logic = NULL){


  DB[,1] = as.factor(DB[,1])


  if (length(levels(DB[,1]))>2){

    stop("The fusion of DB take into account only 2 DBs")

  }

  test_NA_cov = apply(DB[,4:ncol(DB)],1,function(x){sum(is.na(x))})

  if (sum(test_NA_cov == (ncol(DB) - 3)) !=0 ){

    stop("At least one individual has only missing data among the covariates. Please, exclude him before continue")

  }

  if (length(Reduce(intersect,list(nominal,ordinal,logic))) != 0){

    stop("No possible intersection between the nominal,ordinal and logical index")

  } else {}


  if (ncol(DB) != length(c(nominal,ordinal,logic))){

    stop("The number of nominal and ordinal indexes is different from the number of column of DB")

  } else {}


  if (length(logic) != 0){

    if ((max(logic)>ncol(DB))|(min(logic) < 1)){

      stop("Incorrect index in the logic option")} else {}

  } else {}


  if (length(ordinal) != 0){

    if ((max(ordinal)>ncol(DB))|(min(ordinal) < 1)){

      stop("Incorrect index in the ordinal option")} else {}

  } else {}


  if (length(nominal) != 0){

    if ((max(nominal)>ncol(DB))|(min(nominal) < 1)){

      stop("Incorrect index in the nominal option")} else {}

  } else {}



  DB_READY2 = data.frame(base = DB[,1],Y = as.numeric(DB[,2]),Z = as.numeric(DB[,3]))

  ordinal2 = setdiff(ordinal,1:3)
  logic2   = setdiff(logic,1:3)
  nominal2 = setdiff(nominal,1:3)


  if (length(ordinal2)!=0){

    for (j in 1:length(ordinal2)){

      DB_READY2 = data.frame(DB_READY2,as.numeric(DB[,ordinal2[j]]))
      colnames(DB_READY2)[3+j] = colnames(DB)[ordinal2[j]]

    }

  } else {}


  if (length(logic2)!=0){

    for (j in 1:length(logic2)){

      DB_READY2 = data.frame(DB_READY2,as.numeric(DB[,logic2[j]]))
      colnames(DB_READY2)[3+length(ordinal2)+j] = colnames(DB)[logic2[j]]

    }

  } else {}



  if (length(nominal2)!=0){

    for (j in 1:length(nominal2)){

      blok      = transfo_quali(as.factor(DB[,nominal2[j]]),colnames(DB)[nominal2[j]])
      DB_READY2 = data.frame(DB_READY2,blok)

    }

  } else {}


  return(DB_READY2)

}







