#' transfo_dist()
#'
#' A function that prepares the database according to the distance function chosen to evaluate the proximities between rows of databases in the OT algorithm
#'
#' The database must have at least 4 columns (stored in a nonspecific order) : An identifier column (of 2 classes) to distinguish the 2 superimposed databases,
#' a (nominal or ordinal) variable corresponding to the information to merge with its specific encoding in the 1st database (called Y by example),
#' the corresponding target variable in the 2nd database (called Z by example), and at least a covariate common to both databases (that is to say with same encoding in the 2 bases).
#'
#' In this context, the information related to Y in the 2nd database is missing as the information related to Z in the 1st database.
#' We remind users that the purpose of this package is to predict the missing informations in Y and Z.
#'
#' All column indexes (including those related to identifier and target variables Y and Z) of the used database must be declare once, among the \code{quanti}, \code{nominal}, \code{ordinal}, or \code{logic} options.
#'
#' TRANSFORMATIONS ON THE DATABASE ACCORDING TO THE CHOICE OF THE DISTANCE
#' -----------------------------------------------------------------------
#'
#' These necessary transformations are related to the type of each of the covariates.
#' Obviously it depends on the choice of the distance function chooses by user in the \code{prep_choice} option.
#'
#' 1. For the Euclidean (E) and Manhattan (M) distances:
#' The numeric variables are unchanged and the potential transformations only concerned all other types of covariates (if there is any).
#' The related recoding to a boolean variable is 1 for \code{TRUE} and 0 for \code{FALSE}.
#' The recoding for a nominal variable of k classes corresponds to its related disjunctive table (i.e (k-1) binary variables))
#' The ordinal variables are all converted to numeric variables (Please take care that the order of the classes of each of these variables is well specified at the beginning).
#'
#' 2. For the Hamming (H) distance:
#' All the numeric variables must be categorize beforehand as described in the 4th example.
#' The related recoding to a boolean variable is 1 for \code{TRUE} and 0 for \code{FALSE}
#' The recoding for nominal or ordinal variable of k classes corresponds to its related disjunctive table (i.e (k-1) binary variables))
#'
#' 3. For the Gower (G) distance:
#' All covariates remain unchanged
#'
#' 4. Using the principal components from a factor analysis for mixed data (FAMD):
#' A factor analysis for mixed data is done on the covariates of the database and a specific number of the related principal components is remained (depending on the variability part explained by the covariates that the user wishes to keep by varying the \code{info} option).
#' After this step, all new covariates will be obviously in numeric forms.
#'
#' @param DB A data.frame composed of exactly 2 superimposed databases with a column of database identification, 2 columns corresponding to a same information
#' differently encoded in the 2 databases and covariates. The order of the variables have no importance.
#' @param index_DB_Y_Z A vector of exactly 3 integers. The 1st integer must correspond to the index of the database identifiaction column. The 2nd integer corresponds
#' to the index of the target variable in the 1st database while the 3rd integer corresponds to the index of column related to the target varaible in the 2nd database.
#' @param quanti A vector of integers that corresponds to the indexes of columns of all the quantitative variables (DB identification and target variables included)
#' @param nominal A vector of integers that corresponds to the indexes of columns of all the nominal (not ordered) variables (DB identification and target variables included)
#' @param ordinal A vector of integers that corresponds to the indexes of columns of all the ordinal variables (DB identification and target variables included)
#' @param logic A vector of integers that corresponds to the indexes of columns of all the boolean variables.
#' @param prep_choice A character (with quotes) corresponding to the distance function chosen between: The euclidean distance ("E", by default), The Manhattan distance ("M"),
#' the Gower distance ("G"), the Hamming (also called binary) distance and the Euclidean or Manhattan distance, calculated from principal components of a factor analysis of mixed data ("FAMD").
#' @param info A percent value (between 0 and 100) that corresponds to the part of variability taken into account by the principal components of the FAMD when this option is required.
#'
#' @return A data.frame which covariates have been transformed according to the distance function chosen. The columns of the data.frame could have been reordered so that the identifier, Y and Z correspond to the 1st three columns respectively.
#' @export
#'
#' @importFrom stats na.omit
#' @importFrom FactoMineR FAMD
#'
#' @references
#' # For Factor Analysis with mixed data:
#' Pages J. (2004). Analyse factorielle de donnees mixtes. Revue Statistique Appliquee. LII (4). pp. 93-111.
#'
#' # About the Gower distance:
#' Gower, J. C. (1971), “A general coefficient of similarity and some of its properties”. Biometrics, 27, 623--637.
#'
#' # About the other distance measurements:
#' Anderberg, M.R. (1973), Cluster analysis for applications, 359 pp., Academic Press, New York, NY, USA.
#'
#' Borg, I. and Groenen, P. (1997) Modern Multidimensional Scaling. Theory and Applications. Springer.
#'
#' @aliases transfo_dist
#'
#'
#' @author Gregory Guernec
#' \email{gregory.guernec@@inserm.fr}
#'
#' @examples
#'
#' ### Using the simu_data example, suppose that you want to prepare
#' ### your data before applying the OT algorithm with:
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
#' # 2. The Euclidean distance on principal components generated
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
#' # Indexes columns related to potential binary covariates or covariates with
#' # finite number of integer values must be include in the ordinal option.
#' # So in simu_data, the Age covariate can not be directly use without prior
#' # categorization.
#'
#' simu_dat = simu_data[,-8]    # categorization
#' simu_dat$AgeC = cut(simu_data$Age,breaks = c(34,50,54,66))
#' try4 = transfo_dist(simu_dat,quanti = 3, nominal = c(1,4:5,7),ordinal = c(2,6,8),
#' prep_choice = "H")
#'
#' # Other situation with a new numeric covariate with a finite number of values
#' simu_dat = simu_data[,-8]
#' simu_dat$new = sample(1:3,nrow(simu_data),replace = TRUE)
#' try5 = transfo_dist(simu_dat,quanti = 3, nominal = c(1,4:5,7),ordinal = c(2,6,8),
#' prep_choice = "H")
#'
#'
#' ### This function works whatever the order of your columns in your database:
#' # Suppose that we re-order columns in simu_data:
#' simu_data2 = simu_data[,c(2,4:7,3,8,1)]
#'
#' # By changing the corresponding indexes in the index_DB_Y_Z option,
#' # we observe the desired output:
#' try5 = transfo_dist(simu_data2,index_DB_Y_Z = c(8,1,6),quanti = 6:7, nominal = c(2:3,5,8),
#'                      ordinal = c(1,4), logic = NULL, prep_choice = "E")
#'
transfo_dist = function(DB,index_DB_Y_Z = 1:3,quanti = NULL,nominal = NULL,ordinal = NULL,logic = NULL,
                        prep_choice = "E",info = 80){


  if (ncol(DB) < 4){

    stop("Invalid number of columns in your DB: At least 4")

  } else {}



  if (!(prep_choice %in% c("E","M","FAMD","H","G"))){

    stop ("Invalid distance chosen: Please consult the possible options for prep_choice")

  } else {}


  if (length(index_DB_Y_Z)!= 3){

    stop("Invalid length for index_DB_Y_Z: This option must contain the column indexes
          related to the identifiation of DBs and of the 2 chosen targets")

  } else {}


  if (max(index_DB_Y_Z)>ncol(DB)){

    stop("Invalid index in the index_DB_Y_Z option")

  } else {}


  if (length(unique(c(quanti,nominal,ordinal,logic)))!= ncol(DB)){

    stop("The type of at least one variable is missing")

  } else {}

  if (!is.data.frame(DB)){

    stop("The DB must be a data.frame!")

  } else {}

  if (!is.null(quanti)){

    if (max(quanti)>ncol(DB)){

      stop("Incorrect index of columns for quanti")

    } else {}

  } else {}

  if (!is.null(nominal)){

    if (max(nominal)>ncol(DB)){

      stop("Incorrect index of columns for nominal")

    } else {}

  } else {}

  if (!is.null(ordinal)){

    if (max(ordinal)>ncol(DB)){

      stop("Incorrect index of columns for ordinal")

    } else {}

  } else {}

  if (!is.null(logic)){

    if (max(logic)>ncol(DB)){

      stop("Incorrect index of columns for logic")

    } else {}

  } else {}


  if ((length(quanti)>ncol(DB))|(length(nominal)>ncol(DB))|(length(ordinal)>ncol(DB))|(length(logic)>ncol(DB))){

    stop ("The number of at least one type of variables declared exceeds the number of columns of DB")

  } else {}

  if ((length(setdiff(quanti,index_DB_Y_Z))!=0)&(prep_choice == "H")){

    stop("Incompatible type(s) of covariate(s) with distance chosen: No numeric variable with Hamming distance.
         If your variable is binary or has a finite number of values, please put its corresponding index of column
         in the ordinal option, rather than in the quanti option.If not, categorize it, or change distance function.")

  } else {}


  if (length(logic) != 0){

    for (j in logic){

      DB[,j]   = as.numeric(DB[,j])
      quanti   = unique(sort(c(quanti,logic)))
      logic    = NULL

    }

  } else {}

  if (length(Reduce(intersect,list(quanti,nominal,ordinal))) != 0){

    stop("Several types declared for at least one variable. Please consult help to complete the corresponding options")

  } else {}


  typ_var = sort(unique(c(quanti,nominal,ordinal)))

  if (length(typ_var) != ncol(DB)){

    stop("The type of at least one variable is missing. Please consult help to complete the corresponding options.")

  } else {}



  colnames(DB)[index_DB_Y_Z[2]] = "Y"
  colnames(DB)[index_DB_Y_Z[3]] = "Z"

  DB$Y            = transfo_target(DB[,"Y"], levels_order = levels(DB[,"Y"]))$NEW
  DB$Z            = transfo_target(DB[,"Z"], levels_order = levels(DB[,"Z"]))$NEW



  if (index_DB_Y_Z[2] %in% nominal){

    DB$Y = as.character(DB$Y)

  } else if (index_DB_Y_Z[3] %in% nominal){

    DB$Z = as.character(DB$Z)

  } else {}

  ordinal2 = setdiff(ordinal,index_DB_Y_Z)
  nominal2 = setdiff(nominal,index_DB_Y_Z)
  quanti2  = setdiff(quanti,index_DB_Y_Z)

  if (prep_choice == "H"){

    nominal2 = sort(unique(c(ordinal2,nominal2)))
    ordinal2 = NULL

  } else {}


  for (j in sort(unique(c(ordinal2,nominal2)))){

    DB[,j] = as.factor(DB[,j])

  }


  if (prep_choice %in% c("M","E","H")){


    if (length(union(ordinal2,quanti2))!=0){

      for (j in sort(unique(c(ordinal2,quanti2)))){

        DB[,j] = as.numeric(DB[,j])

      }

   } else {}


    if (length(nominal2)!=0){

      name_quali = colnames(DB)[nominal2]


      # Transformation des variables qualis à k modalités en (k-1) binaires:
      # The 1st level is taken as reference

      bin_quali     = transfo_quali(DB[,nominal2[1]])
      nbmod         = length(levels(DB[,nominal2[1]]))


      # Transforming column names

      if (nbmod >=2){

        name_col_quali = paste(name_quali[1],2:nbmod,sep="_")

      } else {

        name_col_quali = name_quali[1]

      }

      # Generalization à T variables (T>1)

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

    DB_NEW = DB

    for (j in ordinal2){

      DB_NEW[,j]   = as.ordered(DB[,j])

    }

    col_nameDB                       = colnames(DB_NEW)[sort(unique(c(nominal2,quanti2,ordinal2)))]
    DB_NEW                           = data.frame(DB_NEW[,index_DB_Y_Z],DB_NEW[,sort(unique(c(nominal2,quanti2,ordinal2)))])
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

      cat(
        "Only",
        nrow(DB_NEW),"rows",
        "(",
        round(nrow(DB_NEW) * 100 / nrow(DB)),
        "%)","are kept here corresponding to complete cases",
        "\n"
      )

    } else {}


    res     = FactoMineR::FAMD(DB_NEW2, ncp = 9, graph = F)



    nb_CP = min((1:nrow(res$eig))[res$eig[,3]> info])

    DB_NEW     = data.frame(DB_NEW[,index_DB_Y_Z],res$ind$coord[,1:nb_CP])

  } else {

    stop("Bad specification for prep_choice option: Please choose between E,M,G or FAMD")

  }

  DB_NEW[,1] = as.factor(DB_NEW[,1])
  DB_NEW[DB_NEW[,1] == levels(DB_NEW[,1])[1],3] = NA
  DB_NEW[DB_NEW[,1] == levels(DB_NEW[,1])[2],2] = NA

  if (unique(DB_NEW[,1])[1] != levels(DB_NEW[,1])[1]){

    stop("Please change the name of your databases in the ID column so that the names of the 2 databases will be alphanumerically ranked in ascending order. You can use A for data1 and B for data2 or simply 1 and 2 by example.")

  } else {}


  return(DB_NEW)

}
