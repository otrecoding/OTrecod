
#' transfo_dist()
#'
#' A function that prepares the database according to the distance function chosen to evaluate the proximity between the lines in the OT algorithm
#'
#' @param DB A data.frame
#' @param quanti A vector of integers
#' @param nominal A vector of integers
#' @param ordinal A vector of integers
#' @param logical A vector of integers
#' @param prep_choice A string of characters corresponding to the distance function chosen.
#'                    The possibilities are: "G"(Default) for the Gower distance, "M" for the Manhattan distance,
#'                    "E" for the Euclidean distance, and "FAMD" for the use of the principal components of a factor analysis of mixed data
#' @param info A double corresponding to the desired percentage of information to take into account in the FAMD (if this option is selected)
#'
#' @return A data.frame which covariates are adapted to the distance function chosen
#' @export
#'
#' @examples
#' soluc2 = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
#' prep2  = transfo_dist(soluc2[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")
#' prep3  = transfo_dist(soluc2[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "E")
#'
#' soluc1 = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "MICE",R_MICE = 2, seed_func = 4036)
#' prep1  = transfo_dist(soluc1[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "FAMD")
#'
#'

transfo_dist = function(DB,quanti = NULL,nominal = NULL,ordinal = NULL,logical = NULL, prep_choice = "G",info = 90){

  if (!is.data.frame(DB)){

    stop("The DB must be a data.frame!")

  } else {}


  if (info > 100){

    stop ("info is a percentage and so can not exceed 100")

  } else {}

  if ((length(quanti)>ncol(DB))|(length(nominal)>ncol(DB))|(length(ordinal)>ncol(DB))|(length(logical)>ncol(DB))){

    stop ("The number of at least one type of variables declared exceeds the number of columns of DB")

  } else {}


  if (length(logical) != 0){

    for (j in logical){

      DB[,j]   = as.numeric(DB[,j])
      quanti   = unique(sort(c(quanti,logical)))
      logical  = NULL

    }

  } else {}

  if (length(Reduce(intersect,list(quanti,nominal,ordinal,logical))) != 0){

    stop("Several types declared for at least one variable. Please consult help to complete the corresponding options")

  } else {}


  typ_var = sort(unique(c(quanti,nominal,ordinal)))

  if (length(typ_var) != ncol(DB)){

    stop("The type of at least one variable is missing. Please consult help to complete the corresponding options.")

  } else {}



  colnames(DB)[2] = "Y1"
  colnames(DB)[3] = "Y2"

  Y1            = transfo_target(DB[,"Y1"], levels_order = levels(DB[,"Y1"]))$NEW
  Y2            = transfo_target(DB[,"Y2"], levels_order = levels(DB[,"Y2"]))$NEW


  # if (2 %in% nominal){

  #  Y1 = transfo_target(DB[,"Y1"])$NEW

  #} else if (3 %in% nominal){

  #  Y2 = transfo_target(DB[,"Y1"])$NEW

  # } else {}

  ordinal2 = setdiff(ordinal,1:3)
  nominal2 = setdiff(nominal,1:3)
  quanti2  = setdiff(quanti,1:3)


  for (j in sort(unique(c(ordinal2,nominal2)))){

    DB[,j] = as.factor(DB[,j])

  }


  if (prep_choice %in% c("M","E")){


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


      DB_NEW = data.frame(DB[,1],Y1,Y2,bin_quali,DB[,sort(unique(c(quanti2,ordinal2)))])
      colnames(DB_NEW)[1] = "DB"
      colnames(DB_NEW) = c(colnames(DB)[1:3],name_col_quali,colnames(DB)[sort(unique(c(quanti2,ordinal2)))])

    } else {

      bin_quali      = NULL
      name_col_quali = NULL
      DB_NEW         = DB

    }




 } else if (prep_choice == "G"){

    DB_NEW = data.frame(DB[,1],Y1,Y2,DB[,sort(unique(c(quanti2,ordinal2,nominal2)))])
    colnames(DB_NEW)[1] = "DB"

    for (j in ordinal2){

        DB_NEW[,j]   = as.ordered(DB[,j])

    }




#    indic_gow            = sort(unique(c(1,quanti,ordinal)))

#    for (j in indic_gow){

#      DB_NEW[,j]   = as.numeric(DB[,j])

#    }


#    for (j in nominal){

#      DB_NEW[,j]   = as.character(DB[,j])

#    }


  } else if (prep_choice == "FAMD"){


    if ((nrow(na.omit(DB[,-(2:3)]))!=(nrow(DB[,-(2:3)])))){

      stop("Presence of NA in the database. Please work with complete or imputed DB when using this method")

    } else {}

    DB_NEW        = DB
    # DB_NEW[,"Y1"] = as.numeric(DB[,"Y1"])
    # DB_NEW[,"Y2"] = as.numeric(DB[,"Y2"])


    DB_NEW2 = DB_NEW[,-(2:3)]
    res     = FAMD(DB_NEW2, ncp = 9, graph = F)



    nb_CP = min((1:nrow(res$eig))[res$eig[,3]> info])

    DB_NEW = data.frame(DB[,1],Y1,Y2,res$ind$coord[,1:nb_CP])
    colnames(DB_NEW)[1] = "DB"


  } else {

    stop("Bad specification for prep_choice option: Please choose between E,M,G or FAMD")

  }

  return(DB_NEW)

}


soluc2 = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
prep2  = transfo_dist(soluc2[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")


prep3  = transfo_dist(soluc2[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "E")





