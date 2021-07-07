#' verif_OT()
#'
#' This function proposes post-process verifications after data fusion by optimal transportation algorithms.
#'
#' In a context of data fusion, where information from a same target population is summarized via two specific variables \eqn{Y} and \eqn{Z} (two ordinal or nominal factors with different number of levels \eqn{n_Y} and \eqn{n_Z}), never jointly observed and respectively stored in two distinct databases A and B,
#' Optimal Transportation (OT) algorithms (see the models \code{OUTCOME}, \code{R_OUTCOME}, \code{JOINT}, and \code{R_JOINT} of the reference (2) for more details)
#' propose methods for the recoding of \eqn{Y} in B and/or \eqn{Z} in A. Outputs from the functions \code{OT_outcome} and \code{OT_joint} so provides the related predictions to \eqn{Y} in B and/or \eqn{Z} in A,
#' and from these results, the function \code{verif_OT} provides a set of tools (optional or not, depending on the choices done by user in input) to estimate:
#' \enumerate{
#' \item the association between \eqn{Y} and \eqn{Z} after recoding
#' \item the similarities between observed and predicted distributions
#' \item the stability of the predictions proposed by the algorithm
#' }
#'
#' A. PAIRWISE ASSOCIATION BETWEEN \eqn{Y} AND \eqn{Z}
#'
#' The first step uses standard criterions (Cramer's V, and Spearman's rank correlation coefficient) to evaluate associations between two ordinal variables in both databases or in only one database.
#' When the argument \code{group.clss = TRUE}, these informations can be completed by those provided by the function \code{\link{error_group}}, which is directly integrate in the function \code{verif_OT}.
#' Assuming that \eqn{n_Y > n_Z}, and that one of the two scales of \eqn{Y} or \eqn{Z} is unknown, this function gives additional informations about the potential link between the levels of the unknown scale.
#' The function proceeds to this result in two steps. Firsty, \code{\link{error_group}} groups combinations of modalities of \eqn{Y} to build all possible variables \eqn{Y'} verifying \eqn{n_{Y'} = n_Z}.
#' Secondly, the function studies the fluctuations in the association of \eqn{Z} with each new variable \eqn{Y'} by using adapted comparisons criterions (see the documentation of \code{\link{error_group}} for more details).
#' If grouping successive classes of \eqn{Y} leads to an improvement in the initial association between \eqn{Y} and \eqn{Z} then it is possible to conclude in favor of an ordinal coding for \eqn{Y} (rather than nominal)
#' but also to emphasize the consistency in the predictions proposed by the algorithm of fusion.
#'
#' B. SIMILARITIES BETWEEN OBSERVED AND PREDICTED DISTRIBUTIONS
#'
#' When the predictions of \eqn{Y} in B and/or \eqn{Z} in A are available in the \code{datab} argument, the similarities between the observed and predicted probabilistic distributions of \eqn{Y} and/or \eqn{Z} are quantified from the Hellinger distance (see (1)).
#' This measure varies between 0 and 1: a value of 0 corresponds to a perfect similarity while a value close to 1 (the maximum) indicates a great dissimilarity.
#' Using this distance, two distributions will be considered as close as soon as the observed measure will be less than 0.05.
#'
#' C. STABILITY OF THE PREDICTIONS
#'
#' These results are based on the decision rule which defines the stability of an algorithm in A (or B) as its average ability to assign a same prediction
#' of \eqn{Z} (or \eqn{Y}) to individuals that have a same given profile of covariates \eqn{X} and a same given level of \eqn{Y} (or \eqn{Z} respectively).
#'
#' Assuming that the missing information of \eqn{Z} in base A was predicted from an OT algorithm (the reasoning will be identical with the prediction of \eqn{Y} in B, see (2) and (3) for more details), the function \code{verif_OT} uses the conditional probabilities stored in the
#' object \code{estimatorZA} (see outputs of the functions \code{\link{OT_outcome}} and \code{\link{OT_joint}}) which contains the estimates of all the conditional probabilities of \eqn{Z} in A, given a profile of covariates \eqn{x} and given a level of \eqn{Y = y}.
#' Indeed, each individual (or row) from A, is associated with a conditional probability \eqn{P(Z= z|Y= y, X= x)} and averaging all the corresponding estimates can provide an indicator of the predictions stability.
#'
#' The function \code{\link{OT_joint}} provides the individual predictions for subject \eqn{i}: \eqn{\widehat{z}_i}, \eqn{i=1,\ldots,n_A} according to the the maximum a posteriori rule:
#' \deqn{\widehat{z}_i= \mbox{argmax}_{z\in \mathcal{Z}} P(Z= z| Y= y_i, X= x_i)}
#' The function \code{\link{OT_outcome}} directly deduces the individual predictions from the probablities \eqn{P(Z= z|Y= y, X= x)} computed in the second part of the algorithm (see (3)).
#'
#' It is nevertheless common that conditional probabilities are estimated from too rare covariates profiles to be considered as a reliable estimate of the reality.
#' In this context, the use of trimmed means and standard deviances is suggested by removing the corresponding probabilities from the final computation.
#' In this way, the function provides in output a table (\code{eff.neig} object) that provides the frequency of these critical probabilities that must help the user to choose.
#' According to this table, a minimal number of profiles can be imposed for a conditional probability to be part of the final computation by filling in the \code{min.neigb} argument.
#'
#' Notice that these results are optional and available only if the argument \code{stab.prob = TRUE}.
#' When the predictions of \eqn{Z} in A and \eqn{Y} in B are available, the function \code{verif_OT} provides in output, global results and results by database.
#' The \code{res.stab} table can produce NA with \code{OT_outcome} output in presence of incomplete shared variables: this problem appears when the \code{prox.dist} argument is set to 0 and can
#' be simply solved by increasing this value.
#'
#' @param ot_out an otres object from \code{\link{OT_outcome}} or \code{\link{OT_joint}}
#' @param group.clss a boolean indicating if the results related to the proximity between outcomes by grouping levels are requested in output (\code{FALSE} by default).
#' @param ordinal a boolean that indicates if \eqn{Y} and \eqn{Z} are ordinal (\code{TRUE} by default) or not. This argument is only useful in the context of groups of levels (\code{group.clss}=TRUE).
#' @param stab.prob a boolean indicating if the results related to the stability of the algorithm are requested in output (\code{FALSE} by default).
#' @param min.neigb a value indicating the minimal required number of neighbors to consider in the estimation of stability (1 by default).
#'
#' @return A list of 7 objects is returned:
#' \item{nb.profil}{the number of profiles of covariates}
#' \item{conf.mat}{the global confusion matrix between \eqn{Y} and \eqn{Z}}
#' \item{res.prox}{a summary table related to the association measures between \eqn{Y} and \eqn{Z}}
#' \item{res.grp}{a summary table related to the study of the proximity of \eqn{Y} and \eqn{Z} using group of levels. Only if the \code{group.clss} argument is set to TRUE.}
#' \item{hell}{Hellinger distances between observed and predicted distributions}
#' \item{eff.neig}{a table which corresponds to a count of conditional probabilities according to the number of neighbors used in their computation (only the first ten values)}
#' \item{res.stab}{a summary table related to the stability of the algorithm}
#'
#'
#' @importFrom stats cor rbinom sd addmargins
#' @importFrom StatMatch pw.assoc
#'
#' @author Gregory Guernec
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @aliases verif_OT
#'
#' @seealso \code{\link{OT_outcome}}, \code{\link{OT_joint}}, \code{\link{proxim_dist}}, \code{\link{error_group}}
#'
#'
#' @references
#' \enumerate{
#' \item Liese F, Miescke K-J. (2008). Statistical Decision Theory: Estimation, Testing, and Selection. Springer
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679. doi:10.1515/ijb-2018-0106
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association. \doi{10.1080/01621459.2020.1775615}
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
#' \donttest{
#'
#' ### Example 2
#' #-----
#' # - Using the data simu_data
#' # - Studying the proximity between Y and Z using standard criterions and studying
#' #   associations by grouping levels of Z
#' # - When only Y is predicted in B
#' # - Tolerated distance between a subject and a profile: 0.30 * distance max
#' # - Using an outcome model (individual assignment with knn)
#' #-----
#'
#' data(simu_data)
#' try2 = OT_outcome(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
#'                  dist.choice = "G",percent.knn = 0.90, maxrelax = 0, prox.dist = 0.3,
#'                  convert.num = 8, convert.clss = 3,
#'                  indiv.method = "sequential",which.DB = "B")
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
#' ver3 = verif_OT(try2, group.clss = TRUE, ordinal = TRUE, stab.prob = TRUE, min.neigb = 5); ver3
#'
#' }
#'

verif_OT = function(ot_out, group.clss = FALSE, ordinal = TRUE, stab.prob = FALSE, min.neigb = 1){

  if (class(ot_out) != "otres"){

    stop("ot_out must be an otres object: output from OT_outcome or OT_joint")

  } else {}


  stopifnot(is.logical(group.clss))
  stopifnot(is.logical(ordinal))
  stopifnot(is.logical(stab.prob))


  ### Test 1: Evaluation of the proximity between the distributions of Y and Z

  DATA1_OT = ot_out$DATA1_OT
  lev1     = levels(DATA1_OT[,2])
  DATA2_OT = ot_out$DATA2_OT
  lev2     = levels(DATA2_OT[,3])
  inst     = ot_out$res_prox

  yy       = factor(c(as.character(DATA1_OT$Y)     , as.character(DATA2_OT$OTpred)) ,levels = lev1)
  zz       = factor(c(as.character(DATA1_OT$OTpred), as.character(DATA2_OT$Z))      ,levels = lev2)

  # Hellinger distance
  YA = prop.table(table(yy[1:nrow(DATA1_OT)]))
  YB = prop.table(table(yy[(nrow(DATA1_OT)+1):(nrow(DATA1_OT)+nrow(DATA2_OT))]))

  ZA = prop.table(table(zz[1:nrow(DATA1_OT)]))
  ZB = prop.table(table(zz[(nrow(DATA1_OT)+1):(nrow(DATA1_OT)+nrow(DATA2_OT))]))

  if (length(yy) == length(zz)){

    hh = data.frame(
      YA_YB = round(sqrt(0.5 * sum((sqrt(YA) - sqrt(YB))^2)),3),
      ZA_ZB = round(sqrt(0.5 * sum((sqrt(ZA) - sqrt(ZB))^2)),3)
    )
    row.names(hh) = "Hellinger dist."

  } else if (length(yy) > length(zz)){

    hh            = data.frame(YA_YB  = round(sqrt(0.5 * sum((sqrt(YA) - sqrt(YB))^2)),3),
                               ZA_ZB  = NA)
    row.names(hh) = "Hellinger dist."

  } else {

    hh = data.frame(
      YA_YB = NA,
      ZA_ZB = round(sqrt(0.5 * sum((sqrt(ZA) - sqrt(ZB))^2)),3)
    )
    row.names(hh) = "Hellinger dist."

  }

  n1 = n2 = nrow(DATA1_OT)

  if (!("OTpred" %in% colnames(ot_out$DATA2_OT))){

    DATA2_OT = NULL

  } else if (!("OTpred" %in% colnames(ot_out$DATA1_OT))){

    DATA1_OT = NULL
    n2       = 0

  } else {}


  if (is.null(DATA1_OT)){

    predZ = DATA2_OT[,3]

    if (is.ordered(DATA1_OT[,2])){

      predY = ordered(c(as.character(DATA1_OT[,2]), as.character(DATA2_OT$OTpred)), levels = lev1)

    } else {

      predY = factor(c(as.character(DATA1_OT[,2]), as.character(DATA2_OT$OTpred)), levels = lev1)

    }

  } else {}


  if (is.null(DATA2_OT)){

    predY = DATA1_OT[,2]

    if (is.ordered(DATA1_OT[,2])){

      predZ = ordered(c(as.character(DATA1_OT$OTpred), as.character(DATA2_OT[,3])), levels = lev2)

    } else {

      predZ = factor(c(as.character(DATA1_OT$OTpred), as.character(DATA2_OT[,3])), levels = lev2)

    }
  }

  if ((!is.null(DATA1_OT))&(!is.null(DATA2_OT))){

    if (is.ordered(DATA1_OT[,2])){

      predY = ordered(c(as.character(DATA1_OT[,2]), as.character(DATA2_OT$OTpred)), levels = lev1)

    } else {

      predY = factor(c(as.character(DATA1_OT[,2]), as.character(DATA2_OT$OTpred)), levels = lev1)

    }

    if (is.ordered(DATA1_OT[,2])){

      predZ = ordered(c(as.character(DATA1_OT$OTpred), as.character(DATA2_OT[,3])), levels = lev2)

    } else {

      predZ = factor(c(as.character(DATA1_OT$OTpred), as.character(DATA2_OT[,3])), levels = lev2)

    }

  }

  ID.DB1 = unique(as.character(DATA1_OT[,1]))
  ID.DB2 = unique(as.character(DATA2_OT[,1]))

  # predZ = ordered(as.factor(c(as.character(DATA1_OT$OTpred), as.character(DATA2_OT[,3])))   , levels = lev2)
  # predY = ordered(as.factor(c(as.character(DATA1_OT[,2])   , as.character(DATA2_OT$OTpred))), levels = lev1)
  ID.DB = c(as.character(DATA1_OT[,1]), as.character(DATA2_OT[,1]))


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
      colnames(resgrp)[1] = "grp levels Z to Y"

    } else {

      resgrp = error_group(predZ,predY, ord = ordinal)
      colnames(resgrp)[1] = paste("combi","Z",paste ="_")
      colnames(resgrp)[1] = "grp levels Y to Z"

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

    row.names(inst$Xobserv) = 1:nrow(inst$Xobserv)
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

    # seed.stb  = seed.stab
    simu1_avg = simu2_avg = simuglb_avg = NULL

    # DATABASE A

    if (!is.null(DATA1_OT)){

      estimatorZA = ot_out$estimatorZA


      for (i in 1:nrow(DATA1_OT)){

        coord1 = which(row.names(estimatorZA[profil[i],,]) == as.character(DATA1_OT$Y)[i])
        coord2 = as.numeric(predZ)[i]

        DATA1_OT$prob[i] = estimatorZA[profil[i],coord1,coord2]

      }

      for (i in 1:nrow(DATA1_OT)){

        profil1   = profil[1:nrow(DATA1_OT)]
        freq_prof = tapply(rep(1,nrow(DATA1_OT[inst$indXA[[profil1[i]]],])),DATA1_OT[inst$indXA[[profil1[i]]],2],sum)
        coord1    = as.numeric(DATA1_OT[i,2])
        DATA1_OT$eff[i] = freq_prof[coord1]

      }

      N1         = length(DATA1_OT$prob[DATA1_OT$eff >= min.neigb])
      simu1_avg  = mean(DATA1_OT$prob[DATA1_OT$eff >= min.neigb])
      simu1_sd   = stats::sd(DATA1_OT$prob[DATA1_OT$eff >= min.neigb])

    }


    ### DATABASE B


    if (!is.null(DATA2_OT)){

      estimatorYB = ot_out$estimatorYB

      for (i in 1:nrow(DATA2_OT)){

        coord1 = which(row.names(estimatorYB[profil[i+n1],,]) == as.character(DATA2_OT$Z)[i])
        #coord2 = as.numeric(DATA2_OT$OTpred)[i]
        coord2 = as.numeric(predY)[i+n2]

        DATA2_OT$prob[i] = estimatorYB[profil[i+n1],coord1,coord2]

      }

      for (i in 1:nrow(DATA2_OT)){


        profil2   = profil[(n1+1):(n1 + nrow(DATA2_OT))]

        freq_prof       = tapply(rep(1,nrow(DATA2_OT[inst$indXB[[profil2[i]]],])),DATA2_OT[inst$indXB[[profil2[i]]],3],sum)
        coord1          = as.numeric(DATA2_OT[i,3])
        DATA2_OT$eff[i] = freq_prof[coord1]

      }


       N2         = length(DATA2_OT$prob[DATA2_OT$eff >= min.neigb])
       simu2_avg  = mean(DATA2_OT$prob[DATA2_OT$eff >= min.neigb])
       simu2_sd   = stats::sd(DATA2_OT$prob[DATA2_OT$eff >= min.neigb])

      }


      simuglb     = c(DATA1_OT$prob[DATA1_OT$eff >= min.neigb], DATA2_OT$prob[DATA2_OT$eff >= min.neigb])
      Nglob       = length(simuglb)
      simuglb_avg = mean(simuglb)
      simuglb_sd  = stats::sd(simuglb)


    if (is.null(DATA1_OT)){

      restand3 = data.frame(N = N2, min.N = min.neigb, mean = round(simu2_avg, 3), sd = round(simu2_sd,3))
      row.names(restand3) = "2nd DB"

    } else if (is.null(DATA2_OT)){

      restand3 = data.frame(N = N1, min.N = min.neigb, mean = round(simu1_avg, 3), sd = round(simu1_sd,3))
      row.names(restand3) = "1st DB"

    } else {

      restand3 = rbind(c(N = Nglob, min.N = min.neigb,  mean = round(simuglb_avg, 3), sd = round(simuglb_sd,3)),
                       c(N = N1   , min.N = min.neigb,  mean = round(simu1_avg, 3)  , sd = round(simu1_sd,3)),
                       c(N = N2   , min.N = min.neigb,  mean = round(simu2_avg, 3)  , sd = round(simu2_sd,3)))

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

  if (is.null(DATA1_OT)){

    Z    = predZ
    conf = stats::addmargins(table(predY,Z))

  } else if (is.null(DATA2_OT)){

    Y    = predY
    conf = stats::addmargins(table(Y,predZ))

  } else {

    conf = stats::addmargins(table(predY,predZ))

  }

  out_verif = list(nb.profil = length(inst$indXA), conf.mat = conf,
                   res.prox = restand, res.grp = resgrp, hell = hh, eff.neig = eff_tb, res.stab = restand3)

}
