
#' imput_cov()
#'
#' This function performs imputations on incomplete covariates, whatever their types, using functions from the package \pkg{MICE} (Van Buuren's Multiple Imputation) or functions from the package \pkg{missMDA} (Simple Imputation with Multivariate data analysis)
#'
#' By default, the function \code{impute_cov} handles missing information using multivariate imputation by chained equations (MICE, see (1) for more details about the method) by integrating in its syntax the function \code{\link[mice]{mice}}.
#' All values of this last function are taken by default, excepted the required number of multiple imputations, which can be fixed by using the argument \code{R_mice}, and the chosen imputation method for each variable (\code{meth} argument),
#' that corresponds to the argument \code{defaultMethod} of the function \code{\link[mice]{mice}}.
#' When multiple imputations are required (for MICE only), each missing information is imputed by a consensus value:
#' the average of the candidate values will be retained for numerical variables, while the most frequent class will be retained for categorical variables (ordinal or not).
#' The output \code{MICE_IMPS} stores the imputed databases to allow users to build their own consensus values by themselves and(or) to eventually assess the variabilities related to the proposed imputed values if necessary.
#' For this method, a random number generator must be fixed or sampled using the argument \code{seed_choice}.
#'
#' When the argument \code{missMDA} is equalled to \code{TRUE}, incomplete values are replaced (single imputation) using a method based on dimensionality reduction called factor analysis for mixed data (FAMD) using the the \code{\link[missMDA]{imputeFAMD}} function of the \pkg{missMDA} package (2).
#' Using this approach, the function \code{imput_cov} keeps all the default values integrated in the function \code{imputeFAMD} excepted the number of dimensions used for FAMD which can be fixed by users (3 by default).
#'
#' @param dat1 A data.frame containing the variables to be imputed and the one participating in the imputations
#' @param indcol A vector of integers. The corresponding column numbers corresponding to the variables to be imputed and those which participate in the imputations.
#' @param R_mice An integer. The number of imputed database generated with MICE method (5 by default).
#' @param meth A vector of characters which specify the imputation method to be used for each column in \code{dat1}.
#' "pmm" for continuous covariates only, "logreg" for binary covariates, "polr" for ordinal covariates, "polyreg" for categorical covariates (no order), (cf \code{\link[mice]{mice}} for more details).
#' @param missMDA A boolean. If \code{TRUE}, missing values are imputed using the factoral analysis for mixed data (\code{\link[missMDA]{imputeFAMD}}) from the \pkg{missMDA} package (Josse, 2016).
#' @param NB_COMP An integer corresponding to the number of components used in FAMD to predict the missing entries (3 by default) when the \code{missMDA} option is TRUE.
#' @param seed_choice An integer used as argument by the set.seed() for offsetting the random number generator (Random integer by default)
#'
#' @return A list of 3 or 4 objects (depending on the missMDA argument). The first three following objects if \code{missMDA} = TRUE, otherwise 4 objects are returned:
#' \item{RAW}{A data.frame corresponding to the initial (or raw) database}
#' \item{IMPUTE}{A character corresponding to the type of selected imputation}
#' \item{DATA_IMPUTE}{A data.frame corresponding to the completed (consensus if multiple imputations) database}
#' \item{MICE_IMPS}{Only if missMDA = FALSE. A list object containing the R imputed databases generated by MICE}
#'
#' @import mice missMDA
#'
#' @importFrom plyr mapvalues
#'
#' @export
#' @author Gregory Guernec
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @aliases imput_cov
#'
#' @references
#' \enumerate{
#' \item van Buuren, S., Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1–67. url{https://www.jstatsoft.org/v45/i03/}
#' \item Josse J, Husson F. (2016). missMDA: A Package for Handling Missing Values in Multivariate Data Analysis. Journal of Statistical Software, 70(1), 1–31. DOI: 10.18637/jss.v070.i01
#' }
#'
#' @examples
#' # Imputation of all incomplete covariates in the table simu_data:
#' data(simu_data)
#'
#' # Here we keep the complete variable "Gender" in the imputation model.
#' # Using MICE (REP = 3):
#' imput_mice = imput_cov(simu_data,indcol = 4:8,R_mice = 3,
#'                        meth = c("logreg","polyreg","polr","logreg","pmm"))
#' summary(imput_mice)
#'
#'
#' # Using FAMD (NB_COMP = 3):
#' imput_famd = imput_cov(simu_data,indcol = 4:8,
#'                        meth = c("logreg","polyreg","polr","logreg","pmm"),
#'                        missMDA = TRUE)
#' summary(imput_famd)
#'
imput_cov = function(dat1,
                     indcol  = 1:ncol(dat1),
                     R_mice  = 5,
                     meth    = rep("pmm", ncol(dat1)),
                     missMDA = FALSE,
                     NB_COMP = 3,
                     seed_choice = sample(1:1000000, 1)){


  if (is.data.frame(dat1) == FALSE) {
    stop ("Your objet must be a data.frame")

  } else {
  }


  if (length(indcol) > ncol(dat1)) {
    stop (
      "The length of indcol can not be greater than the number of columns of the declared data.frame"
    )

  } else {
  }


  if (length(meth) > length(indcol)) {
    stop ("The length of meth must be equal to the length of indcol")

  } else {
  }


  datcov = dat1[, indcol]

  typ_fact   = c("polr", "polyreg", "logreg")

  # Numbers of columns corresponding to missing covariates:

  indic_NA   = (1:ncol(datcov))[apply(datcov, 2, function(x) {
    sum(is.na(x)) != 0
  }) == TRUE]
  datcov_IMP = datcov


  # Converting categorical variables to factor before imputation

  if (TRUE %in% is.element(typ_fact, meth)) {
    indic_typ = (1:(length(meth)))[meth %in% typ_fact]

    for (j in 1:length(indic_typ)) {
      datcov[, indic_typ[j]] = as.factor(datcov[, indic_typ[j]])

    }

  }


  if (missMDA == FALSE) {

    # MICE imputation

    stoc_mice = mice::mice(
      as.data.frame(datcov),
      method = meth,
      m = R_mice,
      print = FALSE,
      remove.collinear = FALSE,
      # printFlag = FALSE,
      seed = seed_choice
    )
    list_mice = mice::complete(stoc_mice, "all")


    # Storage of the multiple imputations of the covariate in a data.frame

    for (k in 1:length(indic_NA)) {
      base_mice_IMP = as.data.frame(datcov[, indic_NA[k]])


      for (u in 1:R_mice) {
        base_mice_IMP = data.frame(base_mice_IMP, mice::complete(stoc_mice, u)[, indic_NA[k]])

      }


      # Imputation by the most frequent class (for categorical variables)

      if (meth[indic_NA[k]] != "pmm") {
        col_imp = as.factor(apply(base_mice_IMP[, 2:ncol(base_mice_IMP)], 1,
                                  function(x) {
                                    as.character(names(table(x))[which.max(table(x))])
                                  }))

         # NEW

         if (meth[indic_NA[k]] == "polr"){

          datcov_IMP[, indic_NA[k]] = ordered(col_imp, levels = levels(datcov[, indic_NA[k]]))

        } else {

          datcov_IMP[, indic_NA[k]] = col_imp

        }


        # Imputation by the mean (for continuous vaiables)

      } else {
        datcov_IMP[, indic_NA[k]] = apply(base_mice_IMP[, 2:ncol(base_mice_IMP)], 1, mean)

      }

    }

  } else if (missMDA == TRUE) {

    FAMD_imp = missMDA::imputeFAMD(datcov, ncp = NB_COMP, seed = seed_choice)$completeObs

    fact_var = sapply(datcov, is.factor)
    typ_var  = sapply(datcov, is.integer)

    for (k in 1:length(typ_var)) {

      if (fact_var[k] == TRUE) {
        datcov_IMP[, k]   = plyr::mapvalues(FAMD_imp[, k],
                                            from = levels(FAMD_imp[, k]),
                                            to = sort(levels(datcov[,k])))

        if (meth[k] == "polr"){
            datcov_IMP[, k]   = ordered(datcov_IMP[, k], levels = levels(datcov[, k]))
        } else {}

      } else {
        datcov_IMP[, k] = FAMD_imp[, k]

      }


      if (typ_var[k] == TRUE) {
        datcov_IMP[, k] = as.integer(FAMD_imp[, k])

      } else {
        datcov_IMP[, k] = datcov_IMP[, k]

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
