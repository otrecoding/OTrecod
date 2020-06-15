#' verif_OT()
#'
#' This function proposes post-process verifications after data fusion by Optimal Transportation algorithms
#'
#' In a context of data fusion, where a same latent information is summarized via 2 specific variables Y and Z (two ordinal or nominal factors with different number of levels nY and nZ), stored in two distinct databases A and B respectively,
#' Optimal Transportation (OT) algorithms (\code{OUTCOME}, \code{R_OUTCOME}, \code{JOINT}, and \code{R_JOINT}, see the reference (2) for more details)
#' propose an original method for the recoding of Y in B and Z in A. The functions \code{OT_outcome} and \code{OT_joint} available in this package, make the use of these algorithms accessible,
#' and the function \code{verif_OT} uses their outputs by providing a set of tools (optional or not, depending on the choices done by user in input) to estimate:
#' \enumerate{
#' \item The association between Y and Z after recoding
#' \item The stability of the predictions proposed by the algorithm
#' }
#'
#' A. PAIRWISE ASSOCIATION BETWEEN Y AND Z
#'
#' The first step uses standard criterions (Cramer's V, chi square test of independence, Spearman coefficient of correlation) to evaluate associations between two ordinal variables in both databases or only in one.
#' When the argument \code{group.clss = TRUE}, these informations can be completed by those provided by the function \code{error_group} (available separately in the package), which is directly integrate in the function \code{verif_OT}.
#' Assuming that nY > nZ, this function studies all possible groupings of classes of the variable Y so that nY = nZ, in order to determine the grouping which will optimize the proximity between Y and Z.
#'
#' B.STABILITY OF THE PREDICTIONS
#'
#' These optional results are available by filling in the argument \code{stab.prob} (Provided by default).
#' Assuming that the missing information of Z in base A was predicted from an OT algorithm (the reasoning will be identical with the prediction of Y in B), the function \code{verif_OT} uses the conditional probabilities stored in the
#' object \code{estimatorZA} (see outputs of the functions \code{OT_outcome} and \code{OT_joint}) which contains the estimates of all the conditional probabilities of Z in A, kwowing a given profile of covariates x and kwowing a level of Y = y.
#' Each individual (or row) from A, can be so associated with a conditional probability \code{p_z,xy} and this latter is used to assess the stability of the algorithm in A.
#'
#' Assuming that the stability of the algorithm in A is defined as its average ability to assign a same prediction of Z to individuals that have a same given profile of covariates x and a same given level of Y, for each individual we can draw one value from a Bernoulli distribution of parameter \code{p_z,xy} and
#' the average of these individual draws gives an estimation of the result.
#' Nevertheless, because estimates of conditional probabilities depends on the number of neighbors related to each profiles of covariates, some probabilities calculated from too small a number of individuals can be assumed to be unreliable and must be exclude of the final computations.
#' User can fix this minimal number of neighbor required in the estimation of a conditional probability to be remained in the final computation by using the argument \code{min.neigb} (By default all conditional probabilities are used).
#' In the case of the predictions of Z in A and Y in B are available, the function \code{verif_OT} provides in output, global and by databases results.
#'
#'
#' @param ot_out An output object of the function \code{OT_outcome} or \code{OT_joint}
#' @param group.clss A boolean indicating if proximity between outcomes is assessed by grouping levels of factor (FALSE by default)
#' @param ordinal A boolean that indicates if Y and Z are ordinal (\code{TRUE} by default) or not. This argument is only used in the context of grouping of level (\code{group.clss}=TRUE).
#' @param stab.prob A boolean indicating if the results related to the stability of the predictions are requested.
#' @param min.neigb A value indicating the minimal number of neighbors required to consider a conditional reliability as reliable
#' @param R A positive integer that indicates the number of repetitions for the bernoulli simulations.
#' @param seed.stab An integer used as argument by the set.seed() for offsetting the random number generator (Random integer by default). Only use if the study about the stability of predictions is required.
#'
#' @return A list of 7 objects is returned:
#' \item{seed}{The list of random number generator used. The first one is fixed by user or randomly chosen}
#' \item{nb.profil}{The number of profiles of covariates}
#' \item{conf.mat}{Global confusion matrix}
#' \item{res.prox}{A summary table related to the association measurements between Y and Z}
#' \item{res.grp}{A summary table related to the study of the proximity of Y and Z using groupings of levels}
#' \item{eff.neig}{A table which corresponds to a count of conditional probabilities according to the number of neighbors used in their computation (Only the first ten values)}
#' \item{res.stab}{A summary table related to the stability of the algorithm}
#'
#'
#' @importFrom stats cor rbinom sd addmargins
#' @importFrom StatMatch pw.assoc
#'
#' @author Gregory Guernec
#' \email{otrecod.pkg@@gmail.com}
#'
#' @aliases verif_OT
#'
#' @seealso \code{\link{OT_outcome}}, \code{\link{OT_joint}}, \code{\link{proxim_dist}}, \code{\link{error_group}}
#'
#'
#' @references
#' \enumerate{
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679 | \url{https://doi.org/10.1515/ijb-2018-0106}
#' \item Gares V, Omer J. Regularized optimal transport of covariates and outcomes in datarecoding(2019).hal-02123109 \url{https://hal.archives-ouvertes.fr/hal-02123109/document}
#' }
#'
#' @export
#'
#' @examples
#'
#' ### Example 1
#' #-----
#' # - Using the data simu_data
#' # - Studying the proximity between Y and Z using standard criterions
#' # - When Y and Z are predicted in B and A respectively
#' # - Using an outcome model (individual assignment with knn)
#' #-----
#' data(simu_data)
#' try1 = OT_outcome(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                  dist.choice = "G",percent.knn = 0.90, maxrelax = 0,
#'                  convert.num = 8, convert.clss = 3,
#'                  indiv.method = "sequential",which.DB = "BOTH",prox.dist = 0.30)
#'
#' ver1 = verif_OT(try1); ver1
#'
#'
#' \dontrun{
#'
#' ### Example 2
#' #-----
#' # - Using the data simu_data
#' # - Studying the proximity between Y and Z using standard criterions and studying
#' #   associations by grouping levels of Z
#' # - When only Y is predicted in B
#' # - Using an outcome model (individual assignment with knn)
#' #-----
#'
#' data(simu_data)
#' try2 = OT_outcome(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                  dist.choice = "G",percent.knn = 0.90, maxrelax = 0,
#'                  convert.num = 8, convert.clss = 3,
#'                  indiv.method = "sequential",which.DB = "B",prox.dist = 0.30)
#'
#' ver2 = verif_OT(try2, group.clss = TRUE, ordinal = TRUE); ver2
#'
#'
#' ### Example 3
#' #-----
#' # - Using the data simu_data
#' # - Studying the proximity between Y and Z using standard criterions and studying
#' #   associations by grouping levels of Z
#' # - Studying the stability of the conditional probabilities
#' # - When Y and Z are predicted in B and A respectively
#' # - Using an outcome model (individual assignment with knn)
#' #-----
#'
#' ver3 = verif_OT(try2, group.clss = TRUE, ordinal = TRUE, stab.prob = TRUE, min.neigb = 5,
#'        seed.stab = 991891); ver3
#'
#' }
#'

verif_OT = function(ot_out, group.clss = FALSE, ordinal = TRUE, stab.prob = FALSE,
                    min.neigb = 1, R = 10, seed.stab = sample(1:1000000, 1)){

  stopifnot(is.list(ot_out))
  stopifnot(is.logical(group.clss))
  stopifnot(is.logical(ordinal))
  stopifnot(is.logical(stab.prob))
  stopifnot(min.neigb >= 1)


  ### Test 1: Evaluation of the proximity between the distributions of Y and Z

  DATA1_OT = ot_out$DATA1_OT
  DATA2_OT = ot_out$DATA2_OT
  inst     = ot_out$res_prox

  n1 = nrow(DATA1_OT)

  if (!("OTpred" %in% colnames(ot_out$DATA2_OT))){

    DATA2_OT = NULL

  } else if (!("OTpred" %in% colnames(ot_out$DATA1_OT))){

    DATA1_OT = NULL

  } else {}


  ID.DB1 = unique(as.character(DATA1_OT[,1]))
  ID.DB2 = unique(as.character(DATA2_OT[,1]))

  predZ = as.factor(c(as.character(DATA1_OT$OTpred),as.character(DATA2_OT[,3])))
  predY = as.factor(c(as.character(DATA1_OT[,2])   ,as.character(DATA2_OT$OTpred)))
  ID.DB = c(as.character(DATA1_OT[,1])             ,as.character(DATA2_OT[,1]))


  if ((!is.null(DATA1_OT))&(!is.null(DATA2_OT))){

  # Using standard criterions: Global

  stoc      = data.frame(predY,predZ)
  N         = nrow(stoc)
  vcram     = round(suppressWarnings(StatMatch::pw.assoc(predY~predZ,data = stoc)$V),2)
  # chisqT    = suppressWarnings(stats::chisq.test(table(predY,predZ))$p.value)
  rankor    = round(stats::cor(rank(predY),rank(predZ)),3)


  # Standard criterions: 1st DB

  stoc1     = stoc[ID.DB == ID.DB1,]
  N1        = nrow(stoc1)
  vcram1    = round(suppressWarnings(StatMatch::pw.assoc(predY~predZ,data = stoc1)$V),2)
  # chisqT1   = suppressWarnings(stats::chisq.test(table(stoc1$predY,stoc1$predZ))$p.value)
  rankor1   = round(stats::cor(rank(stoc1$predY),rank(stoc1$predZ)),3)


  # Standard criterions: 2nd DB

  stoc2     = stoc[ID.DB == ID.DB2,]
  N2        = nrow(stoc2)
  vcram2    = round(suppressWarnings(StatMatch::pw.assoc(predY~predZ,data = stoc2)$V),2)
  # chisqT2   = suppressWarnings(chisq.test(table(stoc2$predY,stoc2$predZ))$p.value)
  rankor2   = round(stats::cor(rank(stoc2$predY),rank(stoc2$predZ)),3)


  restand = rbind(c(N = N , vcram = vcram , rank_cor = rankor),
                  c(N = N1, vcram = vcram1, rank_cor = rankor1),
                  c(N = N2, vcram = vcram2, rank_cor = rankor2))

  colnames(restand)[2] = "V_cram"
  row.names(restand)   = c("Global","1st DB","2nd DB")

  } else {

    stoc      = data.frame(predY,predZ)
    N         = nrow(stoc)
    vcram     = round(suppressWarnings(StatMatch::pw.assoc(predY~predZ,data = stoc)$V),2)
    # chisqT    = suppressWarnings(stats::chisq.test(table(predY,predZ))$p.value)
    rankor    = round(stats::cor(rank(predY),rank(predZ)),3)

    restand           = c(N = N , vcram = vcram, rank_cor = rankor)
    names(restand)[2] = "V_cram"

  }


  # Using comparisons by grouping levels: error_group

  if (group.clss == TRUE){

    n1l = length(levels(as.factor(predZ)))
    n2l = length(levels(as.factor(predY)))

    if (n1l > n2l){

      resgrp = error_group(predY,predZ, ord = ordinal)
      colnames(resgrp)[1] = paste("combi","Y",paste ="_")

    } else {

      resgrp = error_group(predZ,predY, ord = ordinal)
      colnames(resgrp)[1] = paste("combi","Z",paste ="_")

      }


  } else {

    resgrp = NULL

  }

  if (stab.prob == TRUE){


 ### Test 2: Average ability of OT to give a same prediction according to the profile of covariates
 ###         and according to the reliability of the conditional prediction which depends on the
 ###         number of neigbhors for each profile


    inst   = ot_out$res_prox

    # assignment of a profile of covariates to each individual

    aaa    = duplicated(inst$Xobserv)
    bbb    = inst$Xobserv[aaa == TRUE, ]
    ccc    = inst$Xobserv[aaa == FALSE,]
    profil = vector(length = nrow(inst$Xobserv))

    # A profile for each individual
    profil[as.numeric(row.names(ccc))] = 1:nrow(ccc)


    xo     = inst$Xobserv
    incpro = which(profil == 0)


    for (k in 1:length(incpro)){

      dup                    =  duplicated(rbind(xo[incpro[k],],xo[1:(incpro[k]-1),]))
      dup[incpro[0:(k-1)]+1] = FALSE

      profil[incpro[k]] = profil[min(which(dup))-1]

    }


  # assignment of conditional probabilities to each individual

    seed.stb  = seed.stab
    simu1_avg = simu2_avg = simuglb_avg = NULL

    # DATABASE A

    simu1 = list()

    if (!is.null(DATA1_OT)){

      estimatorZA = ot_out$estimatorZA


      for (i in 1:nrow(DATA1_OT)){

        coord1 = which(row.names(estimatorZA[profil[i],,]) == as.character(DATA1_OT$Y)[i])
        coord2 = as.numeric(DATA1_OT$OTpred)[i]

        DATA1_OT$prob[i] = estimatorZA[profil[i],coord1,coord2]

      }

      # Simulations

      for (i in 1:nrow(DATA1_OT)){

        profil1   = profil[1:nrow(DATA1_OT)]
        freq_prof = tapply(rep(1,nrow(DATA1_OT[inst$indXA[[profil1[i]]],])),DATA1_OT[inst$indXA[[profil1[i]]],2],sum)
        coord1    = as.numeric(DATA1_OT[i,2])
        DATA1_OT$eff[i] = freq_prof[coord1]

      }

      for (k in 1:R){

        simu_vf1 = vector(length=nrow(DATA1_OT))

        for (i in 1:nrow(DATA1_OT)){

          if (DATA1_OT$prob[i]>1){DATA1_OT$prob[i] = 1} else {}
          set.seed(seed.stb); simu_vf1[i] = rbinom(1,1,DATA1_OT$prob[i])

        }

        simu1[[k]] = simu_vf1[DATA1_OT$eff >= min.neigb]
        simu1_avg  = c(simu1_avg,mean(simu1[[k]]))
        seed.stb  = seed.stb + 1

      }

    } else {

      simu1 = NULL

    }



  ### DATABASE B

    simu2 = list()

    if (!is.null(DATA2_OT)){

      seed.stb   = seed.stab
      estimatorYB = ot_out$estimatorYB

      for (i in 1:nrow(DATA2_OT)){

        coord1 = which(row.names(estimatorYB[profil[i+n1],,]) == as.character(DATA2_OT$Z)[i])
        coord2 = as.numeric(DATA2_OT$OTpred)[i]

        DATA2_OT$prob[i] = estimatorYB[profil[i+n1],coord1,coord2]

      }

      for (i in 1:nrow(DATA2_OT)){


          profil2   = profil[(n1+1):(n1 + nrow(DATA2_OT))]

          freq_prof       = tapply(rep(1,nrow(DATA2_OT[inst$indXB[[profil2[i]]],])),DATA2_OT[inst$indXB[[profil2[i]]],3],sum)
          coord1          = as.numeric(DATA2_OT[i,3])
          DATA2_OT$eff[i] = freq_prof[coord1]

      }

    # Simulations

    for (k in 1:R){

      simu_vf2  = vector(length=nrow(DATA2_OT))

      for (i in 1:nrow(DATA2_OT)){

        if (DATA2_OT$prob[i]>1){DATA2_OT$prob[i] = 1} else {}
        set.seed(seed.stb); simu_vf2[i] = stats::rbinom(1,1,DATA2_OT$prob[i])

      }

      simu2[[k]] = simu_vf2[DATA2_OT$eff >= min.neigb]
      simu2_avg  = c(simu2_avg, mean(simu2[[k]]))
      seed.stb  = seed.stb + 1


    }

    } else {

      simu2 = NULL
  }

  for (k in 1:R){

    simuglb     = c(simu1[[k]],simu2[[k]])
    simuglb_avg = c(simuglb_avg,mean(simuglb))

  }

  if (is.null(DATA1_OT)){

    restand3 = data.frame(N = length(simu2[[1]])  , min.N = min.neigb, R = R, mean = round(mean(simu2_avg)  ,3), sd = round(stats::sd(simu2_avg),3))
    row.names(restand3) = "2nd DB"

  } else if (is.null(DATA2_OT)){

    restand3 = data.frame(N = length(simu1[[1]])  , min.N = min.neigb, R = R, mean = round(mean(simu1_avg)  ,3), sd = round(stats::sd(simu1_avg),3  ))
    row.names(restand3) = "1st DB"

  } else {

    restand3 = rbind(c(N = length(simuglb), min.N = min.neigb, R = R, mean = round(mean(simuglb_avg),3), sd = round(stats::sd(simuglb_avg),3)),
                     c(N = length(simu1[[1]])  , min.N = min.neigb, R = R, mean = round(mean(simu1_avg)  ,3), sd = round(stats::sd(simu1_avg),3  )),
                     c(N = length(simu2[[1]])  , min.N = min.neigb, R = R, mean = round(mean(simu2_avg)  ,3), sd = round(stats::sd(simu2_avg),3  )))

    row.names(restand3)   = c("Global","1st DB","2nd DB")

  }

  eff              = c(DATA1_OT$eff,DATA2_OT$eff)
  eft              = table(eff)
  eff_tb           = data.frame(as.numeric(names(eft)),Nb.Prob = eft)
  eff_tb           = eff_tb[,-2]
  colnames(eff_tb) = c("Nb.neighbor","Nb.Prob")
  eff_tb           = eff_tb[1:10,]


 } else {

    restand3 = NULL
    eff_tb   = NULL

 }

  out_verif = list(seed = seed.stab, nb.profil = length(inst$indXA), conf.mat = stats::addmargins(table(predY,predZ)),
                   res.prox = restand, res.grp = resgrp, eff.neig = eff_tb, res.stab = restand3)

}







