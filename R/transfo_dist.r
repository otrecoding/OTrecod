#' transfo_dist()
#'
#' This function prepares an overlayed database for data fusion according to the distance function chosen to evaluate the proximities between units.
#'
#'
#' A. EXPECTED STRUCTURE FOR THE INPUT DATABASE
#'
#' In input of this function, the expected database is the result of an overlay between two databases A and B.
#' This structure can be guaranteed using the specific outputs of the functions \code{\link{merge_dbs}} or \code{\link{select_pred}}.
#' Nevertheless, it is also possible to apply directly the function \code{transfo_dist} on a raw database provided that a specific structure is respected in input.
#' The overlayed database (A placed on top of B) must count at least four columns (in an a unspecified order of appearance in the database):
#' \itemize{
#' \item A column indicating the database identifier (two classes or levels if factor: A and B, 1 and 2, ...)
#' \item A column dedicated to the outcome (or target variable) of the first database and denoted \eqn{Y} for example. This variable can be of categorical (nominal or ordinal factor) or continuous type. Nevertheless, in this last case, a warning will appear and the variable will be automatically converted in ordered factors as a prerequisite format of the database before using data fusion algorithms.
#' \item A column dedicated to the outcome (or target variable) of the second database and denoted \eqn{Z} for example. As before, this variable can be of categorical (nominal or ordinal factor) or continuous type, and the variable will be automatically converted in ordered factors as a prerequisite format of the database before using data fusion algorithms.
#' \item At least one shared variable (same encoding in the two databases). Incomplete information is possible on shared covariates only with more than one shared covariate in the final database.
#' }
#' In this context, the two databases are overlayed and the information related to \eqn{Y} in the second database must be missing as well as the information related to \eqn{Z} in the first one.
#' The column indexes related to the database identifier, \eqn{Y} and \eqn{Z} must be specified in this order in the argument \code{index_DB_Y_Z}.
#' Moreover, all column indexes (including those related to identifier and target variables \eqn{Y} and \eqn{Z}) of the overlayed database (DB) must be declared once (and only once), among the arguments \code{quanti}, \code{nominal}, \code{ordinal}, and \code{logic}.
#'
#'
#' B. TRANSFORMATIONS OF CONTINUOUS COVARIATES
#'
#' Because some algorithms dedicated to solving recoding problems like \code{JOINT} and \code{R-JOINT} (see (1) and/or the documentation of \code{\link{OT_joint}}) requires the use of no continuous covariates, the function \code{transfo_dist} integrates in is syntax
#' a process dedicated to the categorization of continuous variables. For this, it is necessary to rigorously fill in the arguments \code{convert_num} and \code{convert_clss}. The first one specifies the indexes of continuous variables to transform
#' in ordered factors while the second one assigns the corresponding desired number of levels.
#' Only covariates should be transformed (not outcomes) and missing informations are not taken into account for the transformations.
#' Notice that all the indexes informed in the argument \code{convert_num} must also be informed in the argument \code{quanti}.
#'
#'
#' C. TRANSFORMATIONS ON THE DATABASE ACCORDING TO THE CHOSEN DISTANCE FUNCTION
#'
#'
#' These necessary transformations are related to the type of each covariate.
#' It depends on the distance function chosen by user in the \code{prep_choice} argument.
#'
#' 1. For the Euclidean ("E") and Manhattan ("M") distances (see (2) and (3)):
#' all the remaining continuous variables are standardized.
#' The related recoding to a boolean variable is 1 for \code{TRUE} and 0 for \code{FALSE}.
#' The recoding of a nominal variable of k classes corresponds to its related disjunctive table (of (k-1) binary variables)).
#' The ordinal variables are all converted to numeric variables (please take care that the order of the classes of each of these variables is well specified at the beginning).
#'
#' 2. For the Hamming ("H") distance (see (2) and (3)):
#' all the continuous variables must be transformed beforehand in categorical forms using the internal process described in section B or via another external approach.
#' The boolean variables are all converted in ordinal forms and then turned into binaries.
#' The recoding for nominal or ordinal variable of k classes corresponds to its related disjunctive table (i.e (k-1) binary variables)).
#'
#' 3. For the Gower ("G") distance (see (4)):
#' all covariates remain unchanged
#'
#' 4. Using the principal components from a factor analysis for mixed data (FAMD (5)):
#' a factor analysis for mixed data is applied on the covariates of the database and a specific number of the related principal components is remained (depending on the minimal part of variability explained by the covariates that the user wishes to keep by varying the \code{info} option).
#' The function integrates in its syntax the function \code{\link[FactoMineR]{FAMD}} of the package \pkg{FactoMiner} (6) using default parameters.
#' After this step, the covariates are replaced by the remaining principal components of the FAMD, and each value corresponds to coordinates linked to each component.
#' Please notice that this method supposed complete covariates in input, nevertheless in presence of incomplete covariates, each corresponding rows will be dropped from the study, a warning will appear and the number of remaining rows will be indicated.
#'
#' @param DB a data.frame composed of exactly two overlayed databases with a column of database identifier, two columns corresponding to a same information
#' differently encoded in the two databases and covariates. The order of the variables have no importance.
#' @param index_DB_Y_Z a vector of exactly three integers. The first integer must correspond to the column index of the database identifier. The second integer corresponds
#' to the index of the target variable in the first database while the third integer corresponds to the index of column related to the target variable in the second database.
#' @param quanti the column indexes of all the quantitative variables (database identificatier and target variables included) stored in a vector.
#' @param nominal the column indexes of all the nominal (not ordered) variables (DB identification and target variables included) stored in a vector.
#' @param ordinal the column indexes of all the ordinal variables (DB identification and target variables included) stored in a vector.
#' @param logic the column indexes of all the boolean variables stored in a vector.
#' @param convert_num the column indexes of the continuous (quantitative) variables to convert in ordered factors. All indexes declared in this argument must have been declared in the argument \code{quanti} (no conversion by default).
#' @param convert_clss according to the argument \code{convert_num}, a vector indicating for each continuous variable to convert the corresponding desired number of levels. If the length of the argument \code{convert_num} exceeds 1 while the length of \code{convert_clss} is equal to 1 (only one integer),
#' each discretization will count the same number of levels.
#' @param prep_choice a character string corresponding to the distance function chosen between: the euclidean distance ("E", by default), the Manhattan distance ("M"),
#' the Gower distance ("G"), the Hamming (also called binary) distance ("H"), and a distance computed from principal components of a factor analysis of mixed data ("FAMD").
#' @param info a ratio (between 0 and 1, 0.8 is the default value) that corresponds to the minimal part of variability that must be taken into account by the remaining principal components of the FAMD when this approach is required.
#' This ratio will fix the number of components that will be kept with this approach. When the argument is set to 1, all the variability is considered.
#'
#' @return A data.frame whose covariates have been transformed according to the distance function or approach (for FAMD) chosen. The columns of the data.frame could have been reordered so that the database identifier, \eqn{Y} and \eqn{Z} correspond to the first three columns respectively.
#' Moreover the order of rows remains unchanged during the process.
#'
#' @export
#'
#' @importFrom stats na.omit sd
#' @importFrom FactoMineR FAMD
#'
#' @seealso \code{\link{transfo_quali}},\code{\link{merge_dbs}}
#'
#' @references
#' \enumerate{
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association. \doi{10.1080/01621459.2020.1775615}
#' \item Anderberg, M.R. (1973). Cluster analysis for applications, 359 pp., Academic Press, New York, NY, USA.
#' \item Borg, I. and Groenen, P. (1997). Modern Multidimensional Scaling. Theory and Applications. Springer.
#' \item Gower, J. C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 27, 623--637.
#' \item Pages J. (2004). Analyse factorielle de donnees mixtes. Revue Statistique Appliquee. LII (4). pp. 93-111.
#' \item Lê S, Josse J, Husson, F. (2008). FactoMineR: An R Package for Multivariate Analysis. Journal of Statistical Software. 25(1). pp. 1-18.
#' }
#'
#' @aliases transfo_dist
#'
#'
#' @author Gregory Guernec
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @examples
#'
#' ### Using the table simu_data:
#'
#' data(simu_data)
#'
#' # 1. the Euclidean distance (same output with Manhattan distance),
#' try1 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),
#'                     ordinal = c(2,6), logic = NULL, prep_choice = "E")
#' # Here Yb2 was stored in numeric: It has been automatically converted in factor
#'
#' # You can also convert beforehand Yb2 in ordered factor by example:
#' sim_data     = simu_data
#' sim_data$Yb2 = as.ordered(sim_data$Yb2)
#' try1 = transfo_dist(sim_data,quanti = 8, nominal = c(1,4:5,7),
#'                     ordinal = c(2,3,6), logic = NULL, prep_choice = "E")
#'
#' # 2. The Euclidean distance generated on principal components
#' #    by a factor analysis for mixed data (FAMD):
#' try2 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),
#'                     ordinal = c(2,6), logic = NULL, prep_choice = "FAMD")
#'
#' # Please notice that this method works only with rows that have complete
#' # information on covariates.
#'
#' # 3. The Gower distance for mixed data:
#' try3 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),
#'                     ordinal = c(2,6), logic = NULL, prep_choice = "G")
#'
#' # 4. The Hamming distance:
#' # Here the quanti option could only contain indexes related to targets.
#' # Column indexes related to potential binary covariates or covariates with
#' # finite number of values must be include in the ordinal option.
#' # So in simu_data, the discretization of the variable age is required (index=8),
#' # using the convert_num and convert_clss arguments (for tertiles = 3):
#'
#' try4 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),ordinal = c(2,6),
#' convert_num = 8, convert_clss = 3, prep_choice = "H")
#'
#'
#' ### This function works whatever the order of your columns in your database:
#' # Suppose that we re-order columns in simu_data:
#' simu_data2 = simu_data[,c(2,4:7,3,8,1)]
#'
#' # By changing the corresponding indexes in the index_DB_Y_Z argument,
#' # we observe the desired output:
#' try5 = transfo_dist(simu_data2,index_DB_Y_Z = c(8,1,6),quanti = 6:7, nominal = c(2:3,5,8),
#'                      ordinal = c(1,4), logic = NULL, prep_choice = "E")
#'
transfo_dist = function(DB,index_DB_Y_Z = 1:3,
                        quanti = NULL, nominal = NULL, ordinal = NULL, logic = NULL,
                        convert_num = NULL, convert_clss = NULL,
                        prep_choice = "E", info = 0.8){


  if (ncol(DB) < 4){

    stop("Invalid number of columns in DB: At least 4")

  } else {}



  if (!(prep_choice %in% c("E","M","FAMD","H","G"))){

    stop ("Invalid distance chosen: Please consult the possible options for prep_choice")

  } else {}


  if (length(index_DB_Y_Z)!= 3){

    stop("Invalid length for index_DB_Y_Z: This argument must contain the column indexes
          related to the identifiation of DBs and of the 2 chosen targets")

  } else {}


  if (max(index_DB_Y_Z)>ncol(DB)){

    stop("Invalid index in the index_DB_Y_Z argument")

  } else {}


  if (length(unique(c(quanti,nominal,ordinal,logic)))!= ncol(DB)){

    stop("The type of at least one variable is missing")

  } else {}

  if (!is.data.frame(DB)){

    stop("The DB must be a data.frame!")

  } else {}

  if (length(quanti)>0){

    if (max(quanti)>ncol(DB)){

      stop("Incorrect index of columns for quanti")

    } else {}

  } else {}

  if (length(nominal)>0){

    if (max(nominal)>ncol(DB)){

      stop("Incorrect index of columns for nominal")

    } else {}

  } else {}

  if (length(ordinal)>0){

    if (max(ordinal)>ncol(DB)){

      stop("Incorrect index of columns for ordinal")

    } else {}

  } else {}

  if (length(logic)>0){

    if (max(logic)>ncol(DB)){

      stop("Incorrect index of columns for logic")

    } else {}

  } else {}


  if ((length(quanti)>ncol(DB))|(length(nominal)>ncol(DB))|(length(ordinal)>ncol(DB))|(length(logic)>ncol(DB))){

    stop ("The number of at least one type of variables declared exceeds the number of columns of DB")

  } else {}

  if (length(setdiff(quanti,c(index_DB_Y_Z,convert_num))!=0)&(prep_choice == "H")){

    stop("Incompatible type(s) of covariate(s) with distance chosen: No numeric variable with Hamming distance.
         If your variable is binary or has a finite number of values, please put its corresponding index of column
         in the ordinal option, rather than in the quanti option.If not, categorize it, or change distance function.")

  } else {}

  if (all(convert_num %in% quanti) == FALSE){

    stop("Inconsistencies between convert_num and quanti arguments")

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

  # Exclude systematically Y and Z from discretization

  convert_clss = convert_clss[convert_num %in% setdiff(convert_num,index_DB_Y_Z)]
  convert_num  = setdiff(convert_num,index_DB_Y_Z)


  if (length(convert_num) != 0){

    tt = 0
    for (k in convert_num){
      tt = tt + 1
      DB[,k]  = cut(DB[,k],breaks = stats::quantile(DB[,k],
                                                    probs = seq(0,1,by = 1/convert_clss[tt]),na.rm = TRUE),
                    include.lowest = TRUE, ordered_result = TRUE)
    }

    ordinal     = sort(c(ordinal,convert_num))
    quanti      = sort(setdiff(quanti,convert_num))
    convert_num = NULL

  } else {}


  if (length(Reduce(intersect,list(quanti,nominal,ordinal,logic))) != 0){

    stop("Several types declared for at least one variable. Please consult help to complete the corresponding options")

  } else {}


  typ_var = sort(unique(c(quanti,nominal,ordinal,logic)))

  if (length(typ_var) != ncol(DB)){

    stop("The type of at least one variable is missing. Please consult help to complete the corresponding options.")

  } else {}



  colnames(DB)[index_DB_Y_Z[2]] = "Y"
  colnames(DB)[index_DB_Y_Z[3]] = "Z"

  DB$Y            = transfo_target(DB[,"Y"], levels_order = levels(DB[,"Y"]))$NEW
  DB$Z            = transfo_target(DB[,"Z"], levels_order = levels(DB[,"Z"]))$NEW



  if (index_DB_Y_Z[2] %in% nominal){

    # DB$Y = as.character(DB$Y)
    DB$Y = factor(as.character(DB$Y), levels = levels(DB$Y))

  } else {}

  if (index_DB_Y_Z[3] %in% nominal){

    # DB$Z = as.character(DB$Z)
    DB$Z = factor(as.character(DB$Z), levels = levels(DB$Z))

  } else {}

  ordinal2 = setdiff(ordinal,index_DB_Y_Z)
  nominal2 = setdiff(nominal,index_DB_Y_Z)
  quanti2  = setdiff(quanti ,index_DB_Y_Z)

  if (prep_choice == "H"){

    nominal2 = sort(unique(c(ordinal2,nominal2)))
    ordinal2 = NULL

  } else {}


  for (j in sort(unique(c(ordinal2,nominal2)))){

    DB[,j] = as.factor(DB[,j])

  }


  if (prep_choice %in% c("M","E","H")){


    if (length(logic) != 0){

      for (j in logic){

        DB[,j]   = as.numeric(DB[,j])
        ordinal2 = unique(sort(c(ordinal2,logic)))
        logic    = NULL

      }

    } else {}


    if (length(union(ordinal2,quanti2))!=0){

      for (j in sort(unique(c(ordinal2,quanti2)))){

        DB[,j] = as.numeric(DB[,j])

      }

    } else {}


    ### Standardization

    if (length(quanti2)!=0){

      for (j in quanti2){

        DB[,j] = round((DB[,j] - mean(DB[,j], na.rm = TRUE))/stats::sd(DB[,j], na.rm = TRUE),4)

      }

    } else {}


    if (length(nominal2)!=0){

      name_quali = colnames(DB)[nominal2]


      # Transformation des variables qualis a k modalites en (k-1) binaires:
      # The 1st level is taken as reference

      bin_quali     = transfo_quali(DB[,nominal2[1]])
      nbmod         = length(levels(DB[,nominal2[1]]))


      # Transforming column names

      if (nbmod >=2){

        name_col_quali = paste(name_quali[1],2:nbmod,sep="_")

      } else {

        name_col_quali = name_quali[1]

      }

      # Generalization to T variables (T>1)

      if (length(nominal2)>1){

        for (j in 2:length(nominal2)){

          bin_quali     = cbind(bin_quali,transfo_quali(DB[,nominal2[j]]))
          nbmod         = length(levels(DB[,nominal2[j]]))

          if (nbmod >=2){

            name_col_quali = c(name_col_quali,paste(name_quali[j],2:nbmod,sep="_"))

          } else {

            name_col_quali = c(name_col_quali,name_quali[j])

          }

        }

      } else {}


      if((length(quanti2)!=0)|(length(ordinal2)!=0)){

        DB_NEW           = data.frame(DB[,index_DB_Y_Z],bin_quali,DB[,sort(unique(c(quanti2,ordinal2)))])
        colnames(DB_NEW) = c(colnames(DB)[index_DB_Y_Z],name_col_quali,colnames(DB)[sort(unique(c(quanti2,ordinal2)))])

      } else {

        DB_NEW           = data.frame(DB[,index_DB_Y_Z],bin_quali)
        colnames(DB_NEW) = c(colnames(DB)[index_DB_Y_Z],name_col_quali)

      }

    } else {

      bin_quali      = NULL
      name_col_quali = NULL
      DB_NEW         = DB

    }



  } else if (prep_choice == "G"){

    ### Standardization

    # if (length(quanti2)!=0){

    #  for (j in quanti2){

    #    DB[,j] = as.numeric(DB[,j])
    #    DB[,j] = round((DB[,j] - mean(DB[,j], na.rm = TRUE))/sd(DB[,j], na.rm = TRUE),4)

    # }

    #} else {}

    DB_NEW = DB

    for (j in ordinal2){

      DB_NEW[,j]   = as.ordered(DB[,j])

    }

    col_nameDB                       = colnames(DB_NEW)[sort(unique(c(nominal2,quanti2,ordinal2,logic)))]
    DB_NEW                           = data.frame(DB_NEW[,index_DB_Y_Z],DB_NEW[,sort(unique(c(nominal2,quanti2,ordinal2,logic)))])
    colnames(DB_NEW)[4:ncol(DB_NEW)] = col_nameDB


    } else if (prep_choice == "FAMD"){


    DB_NEW  = DB
    DB_NEW2 = DB_NEW[,-index_DB_Y_Z]

    if (nrow(stats::na.omit(DB_NEW2))!= nrow(DB_NEW)){

      warning("Presence of NA on covariates. You should work with complete or imputed DB before using this method.
              By default, only rows with no NA on covariates have been kept.")


      countNA = apply(DB_NEW2,1,function(x){sum(is.na(x))})
      DB_NEW  = DB_NEW[countNA==0,]
      DB_NEW2 = DB_NEW[,-index_DB_Y_Z]

      message(
        "Only ",
        nrow(DB_NEW)," rows ",
        "(",
        round(nrow(DB_NEW) * 100 / nrow(DB),0),
        "%)","are kept here corresponding to complete cases",
        "\n"
      )

    } else {}


    res     = FactoMineR::FAMD(DB_NEW2, ncp = 11, graph = F)

    nb_CP   = min((1:nrow(res$eig))[res$eig[,3] > info*100])
    DB_NEW  = data.frame(DB_NEW[,index_DB_Y_Z],res$ind$coord[,1:nb_CP])

  } else {

    stop("Bad specification for prep_choice option: Please choose between E,M,G or FAMD")

  }

  DB_NEW[,1] = as.factor(DB_NEW[,1])
  DB_NEW[DB_NEW[,1] == levels(DB_NEW[,1])[1],3] = NA
  DB_NEW[DB_NEW[,1] == levels(DB_NEW[,1])[2],2] = NA

  if (unique(DB_NEW[,1])[1] != levels(DB_NEW[,1])[1]){

    stop("Please change the name of your databases in the ID column so that the names of the 2 databases will be alphanumerically ranked in ascending order.
          By example, usee A for data1 and B for data2 or simply 1 and 2.")

  } else {}


  return(DB_NEW)

}
