
#' imput_cov()
#'
#' A function which performs imputations on variables using MICE (Van Buuren's Multiple Imputation) or missMDA (Simple Imputation wiht Multivariate data analysis)
#' (missMDA and mice packages required)
#'
#' @param dat1 A data.frame. containing the variables to be imputed and the one participating in the imputations
#' @param indcol A vector of integers. The corresponding column numbers corresponding to the variables to be imputed and those which participate in the imputations.
#' @param R_mice An integer. The number of imputed database to generate with MICE method (5 by default).
#' @param meth A vector of characters which specify the imputation method to be used for each column in dat1.
#'             "pmm" for a continuous covariate, "logreg" for a binary covariate, "polr" for an ordinal covariate, "polyreg" for an ordinal covariates.
#' @param missMDA A boolean. If TRUE, missing values are imputed using the factoral analysis for mixed data (FAMD) from the missMDA package (Simple imputation)
#' @param NB_COMP An integer corresponding to the number of components used to predict the missing entries (3 by default)
#' @param seed_choice An integer used as argument by the set.seed() for offsetting the random number generator (Random integer by default)
#'
#' @return A list of 3 objects:
#'         RAW A data.frame corresponding to the raw database
#'         IMPUTE A character corresponding to the type of imputation selected
#'         DATA_IMPUTE A dta.frame The imputed (consensus if multiple imputations) database
#'         MICE_IMPS
#'
#' @export
#'
#' @examples
#' # From the example of the merge.dbs() function ...
#' soluc2  = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
#'
#' imput_mice = imput_cov(soluc2$DB_READY,indcol = 4:6,R_mice = 3,meth = c("pmm","polr","logreg"))
#' summary(imput_mice)
#' imput_famd = imput_cov(soluc2$DB_READY,indcol = 4:6,meth = c("pmm","polr","logreg"),missMDA = TRUE)
#' summary(imput_famd)
#'
imput_cov = function(dat1,indcol = 1:ncol(dat1),R_mice = 5,meth = rep("pmm",ncol(dat1)), missMDA = FALSE,NB_COMP = 3,
                     seed_choice = sample(1:1000000,1)){


  if (is.data.frame(dat1) == FALSE){

    stop ("Your objet must be a data.frame")

  } else {}


  if (length(indcol)>ncol(dat1)){

    stop ("The length of indcol can not be greater than the number of columns of the declared data.frame")

  } else {}


  if (length(meth)>length(indcol)){

    stop ("The length of meth must be equal to the length of indcol")

  } else {}


  datcov = dat1[,indcol]

  typ_fact   = c("polr","polyreg","logreg")

  # Numbers of columns corresponding to missing covariates:

  indic_NA   = (1:ncol(datcov))[apply(datcov,2,function(x){sum(is.na(x))!=0}) == TRUE]
  datcov_IMP = datcov


  # Converting categorical variables to factor before imputation

  if (TRUE %in% is.element(typ_fact,meth)){

    indic_typ = (1:(length(meth)))[meth %in% typ_fact]

    for (j in 1:length(indic_typ)){

      datcov[,indic_typ[j]] = as.factor(datcov[,indic_typ[j]])

    }

  }


  if (missMDA == FALSE){


    # MICE imputation

    stoc_mice = mice(as.data.frame(datcov),method = meth, m = R_mice,print=FALSE,remove.collinear=FALSE,seed = seed_choice)
    list_mice = complete(stoc_mice, "all")

    # Storage of the multiple imputations of the covariate in a data.frame

    for (k in 1:length(indic_NA)){

      base_mice_IMP = as.data.frame(datcov[,indic_NA[k]])


      for (u in 1:R_mice){

        base_mice_IMP = data.frame(base_mice_IMP,complete(stoc_mice,u)[,indic_NA[k]])

      }


      # Imputation de la catégorie la plus fréquente (dans le cas de données
      # catégorielles)

      if (meth[indic_NA[k]]!="pmm"){

        col_imp = as.factor(apply(base_mice_IMP[,2:ncol(base_mice_IMP)],1,
                                                    function(x){as.character(names(table(x))[which.max(table(x))])}))

        datcov_IMP[, indic_NA[k]] = ordered(col_imp,levels = levels(datcov[,indic_NA[k]]))




        # Imputation de la moyenne (dans le cas de données continues)

      } else {

        datcov_IMP[, indic_NA[k]] = apply(base_mice_IMP[,2:ncol(base_mice_IMP)],1,mean)

      }

    }


  } else if (missMDA == TRUE){

    FAMD_imp = imputeFAMD(datcov,ncp = NB_COMP,seed = seed_choice)$completeObs

    fact_var = sapply(datcov,is.factor)
    typ_var  = sapply(datcov,is.integer)

    for (k in 1:length(typ_var)){


      # dat1[,typvar] = apply(dat1[,typvar],2,function(x){sort(levels(x))})
      #datcov_IMP[,k] = ifelse(typ_var[k]==TRUE,as.integer(FAMD_imp[,k]),FAMD_imp[,k])
      #datcov_IMP[,k] = ifelse(fact_var[k]==TRUE,mapvalues(FAMD_imp[,k],from = levels(FAMD_imp[,k]),to = levels(dat1[,k])),FAMD_imp[,k])



      if (fact_var[k]==TRUE){

        datcov_IMP[,k]   = mapvalues(FAMD_imp[,k],from = levels(FAMD_imp[,k]),to = sort(levels(dat1[,indcol[k]])))
        # datcov_IMP[,k] = mapvalues(FAMD_imp[,k],from = levels(FAMD_imp[,k]),to = sort(levels(datcov[,k])))
        datcov_IMP[,k]   = ordered(datcov_IMP[,k],levels = levels(dat1[,indcol[k]]))

      } else {

        datcov_IMP[,k] = FAMD_imp[,k]

      }


      if (typ_var[k]==TRUE){

        datcov_IMP[,k] = as.integer(FAMD_imp[,k])

      } else {

        datcov_IMP[,k] = datcov_IMP[,k]

      }

    }

  }


  # Return the imputed database

  if (missMDA == FALSE){

    return(list(RAW = dat1,IMPUTE = "MICE", DATA_IMPUTE = datcov_IMP,MICE_IMPS = list_mice))

  } else {

    return(list(RAW = dat1,IMPUTE = "MDA", DATA_IMPUTE = datcov_IMP))

    }

}


imput_cov(soluc2$DB_READY,indcol = 4:6,R_mice = 3,meth = c("pmm","polr","logreg"))

dat1 = soluc2$DB_READY
indcol = 4:6
R_mice = 2
meth = c("pmm","polr","logreg")
missMDA = FALSE










