#' OT()
#'
#' Implementation of the Optimal Transportation (OT) algorithm for data integration with or without relaxation on the constraints of marginal
#' distributions
#'
#' Assuming that:
#'
#' \enumerate{
#' \item{Y and Z summarize a same information of interest encoded in two distinct forms stored in two independent databases A and B (no individuals or rows in common)}
#' \item{The two databases have a set of common covariates X (i.e the same variables stored in the same forms in the two bases)}
#' \item{You would like to merge vertically A and B keeping this variable of interest in at least one of its two encodings}
#' }
#' So, this function gives you a possible answer to this recoding problem by using the Optimal Transportation Theory.
#' The models implemented in this function correspond to the OUTCOME and R-OUTCOME models described in the referenced article 2.
#'
#' The OUTCOME model performs data fusion of the two databases assuming that the transport of the target variable is equally distributed in the two databases of interest (ie L(Y) is equal in A and B, L(Z) is equal in A and B).
#' By default, the OUTCOME model is implemented by fixing the \code{maxrelax} option equals to 0.
#'
#' When the \code{maxrelax} option > 0, a relaxation of the constraints of marginal is possible and the R-OUTCOME model is performed.
#' This model allows to override the distributional assumption previously described (Here Y and Z can be differently distributed in A and B).
#'
#' In this function, these models are implemented in two parts:
#' 1. The 1st part partially resolves the recoding problem using a linear programming model: It gives user an estimator of the joint distributions of Y and Z in database A (and B).
#' 2. The 2nd part uses the result of the 1st part and an independent procedure to gives you the individuals predictions of Y in B and Z in A.
#'
#' 2 procedures are implemented in this function and can be chosen using the \code{indiv_method} option.
#' The 1st procedure (\code{indiv_method} = "sequential") corresponds to a nearest neighbor algorithm and can make some arbitrary decisions.
#' The 2nd procedure (\code{indiv_method} = "optimal") searchs for the individual predictions (of Z in A by example) that minimize the sum of the indivual distances in A with the modalities of Z in B, by solving an optimization problem (simplex algorithm).
#' Be a little patient please, the time execution of this last procedure may be quite long.
#'
#'
#' REQUIRED COHERENCES BETWEEN OPTIONS
#' ------------------------------------
#' The users must respect some coherences between the \code{prep_choice} and \code{norm} options for the choice of the distance measurements.
#' The \code{prep_choice} option prevails on the \code{norm} option that is why:
#' \enumerate{
#' \item If the Gower distance is selected in the \code{prep_choice} option ("G"), the \code{norm} option will be forced to 3.
#' \item If the FAMD approach is retained (\code{prep_choice} = "FAMD"), the \code{norm} option can never be equal to 0 and such a selection will return a message of error.
#' \item If the Hamming distance is selected in the \code{prep_choice} option ("H"), all covariates in your databases must be binaries covariates.
#' }
#'
#' For more details, please consult the reference article 2.
#'
#' @aliases OT OT_group ot
#'
#' @param datab A data.frame that must have at least 4 columns sorted in a non-specific order. One column must be a key column of 2 classes for the databases identification, where the names of the two databases
#' must be alphanumerically ranked in ascending order (By examples: 1 for the top database and 2 for the database from below, or more logically here A and B  ...But NOT B and A!).One column (Y here but other names are allowed)
#' must correspond to the target variable related to the information of interest to merge with its specific encoding in the database A (corresponding encoding should be so missing in the database B). In the same way,
#' one column (Z here) corresponds to the target variable that summarizes the same information as Y but with its specific encoding in the database B (corresponding encoding should be so missing in the database A).
#' Finally, your database must have at least one covariate with same encoding in A and B. Please not that, if your data.frame has only 4 columns, that is to say, only one covariate, if this latter has NA, and unless user
#' has previously imputed the corresponding missing information, the OT algorithm will only run with complete cases.
#' @param index_DB_Y_Z A vector of exactly 3 integers. The 1st integer must correspond to the index of the databases identification column (DB identifier). The 2nd integer corresponds
#' to the index of the target variable in the 1st database (A) while the 3rd integer corresponds to the index of column related to the target variable in the 2nd database (B).
#' @param quanti A vector of integers that corresponds to the indexes of columns of all the quantitative variables (DB identification and target variables included if it is the case for them).
#' @param nominal A vector of integers that corresponds to the indexes of columns of all the nominal (not ordered) variables (DB identification and target variables included if it is the case for them).
#' @param ordinal A vector of integers that corresponds to the indexes of columns of all the ordinal variables (DB identification and target variables included if it is the case for them).
#' @param logic A vector of integers that corresponds to the indexes of columns of all the boolean variables of the data.frame.
#' @param prep_choice A character (with quotes) corresponding to the distance function chosen between: The euclidean distance ("E", by default), The Manhattan distance ("M"),
#' the Gower distance ("G"), the Hamming (also called binary) distance and the Euclidean or Manhattan distance, calculated from principal components of a factor analysis of mixed data ("FAMD").
#' @param infoFAMD A percent value (between 0 and 100) that corresponds to the part of variability taken into account by the principal components of the FAMD when this option is required.
#' @param percent_c Percent of closest neighbors taken in the computation of the costs.
#' @param maxrelax Maximum percentage of deviation from expected probability masses (0 for the OUTCOME model, a non-zero value otherwise). Please consult the reference article for more details.
#' @param norm A value (with no quotes) corresponding to the distance function chosen to calculate the distances between covariates. Equals to 0 for the Hamming distance, 1 for the Manhattan distance (by default), 2 for the euclidean distance and 3 for Gower distance.
#' @param prox_dist A value between 0 and 1 that corresponds to a threshold (0.1 by default) below which an individual is considered significantly close (neighbor) to a given profile of covariates.
#' @param indiv_method A string of characters that specifies the method used to get individual transport from group joint probabilities, "sequential" by default, or "optimal". Please consult the referenced article 2 for more details.
#'
#' @return A list containing 7 elements:
#' \describe{
#'     \item{TIME_EXE}{Running time of the function}
#'     \item{TRANSPORT_A}{Cost matrix corresponding to an estimation (gamma, see reference for more details) to the joint distribution of (YA,ZA)}
#'     \item{TRANSPORT_B}{Cost matrix corresponding to an estimation to the joint distribution of (YB,ZB)}
#'     \item{profile}{A data.frame that gives all details about the remaining P profiles of covariates. These informations can be linked to the \code{estimatorZA} and the \code{estimatorYB} objects for a better interpretation of the results}
#'     \item{estimatorZA}{Estimates of the probability distribution of Z conditional to X and Y in database A from predictions stored in an array. The number of rows of each table corresponds to the total number of the P profiles of covariates.
#'     The number of columns of each table corresponds to the number of levels of Y. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of Z}
#'     \item{estimatorYB}{Estimates of the probability distribution of Y conditional to X and Z in database B from predictions stored in an array. The number of rows of each table corresponds to the total number of profiles of covariates.
#'     The number of columns of each table corresponds to the number of levels of Z. The row names of each table corresponds to the values of the covariates sorted by order of appearance in the merged database. The third element of the array is the possible level of Y}
#'     \item{DATA1_OT}{database A with imputed individual prediction on Z using OT}
#'     \item{DATA2_OT}{database B with imputed individual prediction on Y using OT}
#' }
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#' \email{gregory.guernec@@inserm.fr}
#'
#' @references
#' # Article 1:
#' Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' 0, 20180106 (2019) | \url{https://doi.org/10.1515/ijb-2018-0106}
#'
#' # Article 2:
#' Gares V, Omer J. Regularized optimal transport of covariates and outcomes in datarecoding(2019).hal-02123109 \url{https://hal.archives-ouvertes.fr/hal-02123109/document}
#'
#'
#' @seealso \code{\link{transfo_dist}},\code{\link{proxim_dist}}, \code{\link{avg_dist_closest}}, \code{\link{indiv_grp_closest}}, \code{\link{indiv_grp_optimal}}
#'
#' @import ROI ROI.plugin.glpk
#' @importFrom ompr sum_expr
#' @importFrom ompr.roi with_ROI
#' @importFrom plyr mapvalues
#' @export
#'
#'@examples
#'
#' #' ### Using the \code{simu_data} object
#' ### Y and Z are a same variable encoded in 2 different forms:
#' ### (3 levels for Y and 5 levels for Z)
#' #--------
#' data(simu_data)
#'
#' ### An example of OUTCOME model with nearest neighbor procedure and Manhattan distance
#' ### (original OT algorithm from article 1)
#' try1 = OT(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6), prep_choice = "M",
#'           percent_c = 0.90, maxrelax = 0, norm = 1,indiv_method = "sequential")
#'
#'
#' head(try1$estimatorZA[,,1])
#' # ... Corresponds to P[Z = 1|Y,P1] when P1 corresponds to the 1st profile of covariates (P_1)
#' # detailed in the 1st row of the profile object:
#' try1$profile[1,]   # Details of P_1
#'
#' # So estimatorZA[1,1,1]= 0.2 corresponds to an estimation of:
#' # P[Z = 1|Y=[20-40],Gender_2=0,Treatment_2=1,Treatment_3=0,Smoking_2=1,Dosage=3,Age=65.44]
#' # Thus, we can conclude that all individuals with the P_1 profile of covariates have
#' # 20 percents of chance to be affect to the 1st level of Z in database A.
#' # ... And so on, the reasoning will be the same for the estimatorYB object
#'
#' \dontrun{
#' ### An example of OUTCOME model using an optimization procedure for the individual predictions
#' ### and Manhattan distance
#' try2 = OT(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6), prep_choice = "M",
#'           percent_c = 0.90, maxrelax = 0, norm = 1, indiv_method = "optimal")
#'
#'
#' ### An example of R-OUTCOME model with nearest neighbor procedure and FAMD approach
#' try3 = OT(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6), prep_choice = "FAMD",
#'           percent_c = 0.90, maxrelax = 0.5, norm = 2, indiv_method = "sequential")
#'
#'
#' ### An example of R-OUTCOME model with nearest neighbor procedure and Gower distance
#' simu_data2 = simu_data[,c(4,1,3,5:7,2,8)]
#'
#' try3G  = OT(simu_data2, index_DB_Y_Z = c(2,7,3),quanti = c(3,8), nominal = c(1:2,4,6),
#'             ordinal = c(5,7), prep_choice = "G",percent_c = 0.90, maxrelax = 0.5,
#'             indiv_method = "sequential")
#'
#'
#' ### An example of R-OUTCOME model using an optimization procedure for the individual predictions
#' # With Manhattan distance
#' try4 = OT(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6), prep_choice = "M",
#'           percent_c = 0.90, maxrelax = 0.4, norm = 1, indiv_method = "optimal")
#'
#' # With Gower distance
#' try4G = OT(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6), prep_choice = "G",
#'            percent_c = 0.90, maxrelax = 0.4, indiv_method = "optimal")
#'
#'
#'  ### An example of R-OUTCOME model using an optimization procedure for the individual predictions
#'
#' data(tab_test)
#' # Adding NAs in Y1 and Y2 and keeping only a part of tab_test
#' tab_test2 = tab_test[c(1:100,5001:5100),]
#' tab_test2[tab_test2$ident == 2,2] = NA
#' tab_test2[tab_test2$ident == 1,3] = NA
#'
#' try5 = OT(tab_test2, nominal = c(1,4:5), ordinal = c(2,3,6), prep_choice = "M",
#'           percent_c = 0.90, maxrelax = 0.7, norm = 1, indiv_method = "optimal")
#' }

OT = function(datab, index_DB_Y_Z = 1:3, quanti = NULL, nominal = NULL, ordinal = NULL,logic = NULL, prep_choice = "E",
              infoFAMD = 80, percent_c = 1.0, maxrelax = 0, norm = 1, prox_dist = 1, indiv_method = "sequential"){



  tstart = Sys.time()

  norm = ifelse(prep_choice == "G",3,norm)

  cat("---------------------------------------","\n")
  cat("OT PROCEDURE in progress ..."           ,"\n")
  cat("---------------------------------------","\n")
  cat("Type                  = ", ifelse(maxrelax == 0,"OUTCOME","R-OUTCOME"),"\n")
  cat("Distance              = ", ifelse(norm == 0,"Hamming",ifelse(norm == 1,"Manhattan",ifelse(norm == 2,"Euclidean","Gower"))),"\n")
  cat("Percent closest       = ", 100.0*percent_c, "%","\n")
  cat("Relaxation term       = ", ifelse(maxrelax == 0,"NO","YES"),"\n")
  cat("Relaxation weight     = ", maxrelax,"\n")
  cat("Recoding method       = ", ifelse(indiv_method == "sequential","Sequential","Optimal"),"\n")
  cat("---------------------------------------","\n")



  if (!(indiv_method %in% c("sequential","optimal"))){

    stop("Incorrect syntax int the indiv_method option")

  } else {}


  if ((prep_choice == "FAMD") & (norm == 3)){

    stop("Invalid norm option with prep_choice option")

  } else {}


  dataB = transfo_dist(datab,index_DB_Y_Z = index_DB_Y_Z,quanti = quanti, nominal = nominal, ordinal = ordinal,
                       logic = logic, prep_choice = prep_choice, info = infoFAMD)


  if (prep_choice == "H"){

    datac  = dataB[,-index_DB_Y_Z]

    test_H = apply(as.data.frame(datac),2,function(x){length(table(x))})

    if (max(test_H)>2){stop("With Hamming distance, all your covariates must be binaries !")}else{}

  } else {}


  inst = proxim_dist(dataB,norme = norm, prox = prox_dist)


  # Redefine A and B for the model
  nA      = inst$nA
  nB      = inst$nB
  A       = 1:inst$nA
  B       = 1:inst$nB
  Y       = inst$Y
  Z       = inst$Z
  indY    = inst$indY
  indZ    = inst$indZ
  indXA   = inst$indXA
  indXB   = inst$indXB
  #prof  = do.call(paste0,unique(inst$Xobserv))
  prof    = as.data.frame(unique(inst$Xobserv))
  ID_prof = paste(rep("P",nrow(prof)),1:nrow(prof),sep="_")

  nbindY = numeric(0)

  for (y in Y){

    nbindY  = c(nbindY,length(inst$indY[[y]]))}
    nbindZ  = numeric(0)

  for (z in Z){

    nbindZ  = c(nbindZ,length(inst$indZ[[z]]))}
    freqY   = numeric(0)

  for (y in Y){

    freqY   = c(freqY,nbindY[y] / length(A))}
    freqZ   = numeric(0)

  for (z in Z){

    freqZ = c(freqZ,nbindZ[z] / length(B))

  }


  ###########################################################################
  # Compute data for aggregation of the individuals
  ###########################################################################

  nbX = length(indXA);

  # Computation of the cost matrix as average distances between the
  # individuals of two groups
  C = avg_dist_closest(inst, percent_closest = percent_c)[[1]]

  # absdevA = absdevB = deviationA = deviationB = transport = transportA = transportB = NULL

  if (maxrelax == 0){

    # Objective: minimize the distance between individuals of A and B

    Min = c(as.numeric(C),rep(0,2*length(Y)+2*length(Z)))


    result          = ompr::MIPModel() %>%

      ompr::add_variable(transport[y,z], y = Y, z = Z, type = "continuous") %>%
      ompr::add_variable(deviationA[y] , y = Y       , type = "continuous") %>%
      ompr::add_variable(absdevA[y]    , y = Y       , type = "continuous") %>%
      ompr::add_variable(deviationB[z]        , z = Z, type = "continuous") %>%
      ompr::add_variable(absdevB[z]           , z = Z, type = "continuous") %>%

      ompr::set_objective(sum_expr(C[y,z]*transport[y,z], y = Y, z = Z) , "min") %>%

      ompr::add_constraint(sum_expr(transport[y,z]       , z = Z) == freqY[y] + deviationA[y], y = Y) %>%
      ompr::add_constraint(sum_expr(transport[y,z], y = Y       ) == freqZ[z] + deviationB[z], z = Z) %>%
      ompr::add_constraint(sum_expr(deviationA[y] , y = Y) == 0 ) %>%
      ompr::add_constraint(sum_expr(deviationB[z]        ,z = Z )== 0) %>%
      ompr::add_constraint(deviationB[z]<= absdevB[z]    , z = Z) %>%
      ompr::add_constraint(deviationB[z]>= -absdevB[z]   , z = Z) %>%
      ompr::add_constraint(sum_expr(absdevB[z],z=Z)<= maxrelax/2.0) %>%
      ompr::add_constraint(deviationA[y] <= absdevA[y],y =Y) %>%
      ompr::add_constraint(deviationA[y] >= -absdevA[y],y =Y) %>%
      ompr::add_constraint(sum_expr(absdevA[y],y = Y)<= maxrelax/2.0) %>%

      ompr::solve_model(with_ROI(solver = "glpk"))

    # Solve the problem

    solution       = ompr::get_solution(result, transport[y,z])

    # Extract the values of the solution

    transportA_val = matrix(solution[,4], length(Y),length(Z))
    transportB_val = transportA_val


  } else {

    b2 = numeric(length(Z))

    for (z in Z){

      stc_sum2 = vector(length = nbX)

      for (i in 1:nbX){

        stc_sum2[i] = ifelse(length(indXB[[i]])==0,1/length(Z),length(indXB[[i]][inst$Zobserv[indXB[[i]]+nA] == z])/ length(indXB[[i]]))*length(indXA[[i]])/nA

      }
      b2[z]= sum(stc_sum2)
    }
    b2   = b2/sum(b2)


    result <-  MIPModel() %>%

      ompr::add_variable(transportA[y,z],y = Y, z = Z, type = "continuous",lb=0) %>%
      ompr::add_variable(deviationB[z]  ,       z = Z, type = "continuous") %>%
      ompr::add_variable(absdevB[z]     ,       z = Z, type = "continuous",lb=0) %>%

      ompr::set_objective(sum_expr(C[y,z]*transportA[y,z],y = Y,z=Z) , "min") %>%

      ompr::add_constraint(sum_expr(transportA[y,z], z = Z) == freqY[y], y =Y) %>%
      ompr::add_constraint(sum_expr(transportA[y,z],y = Y) - deviationB[z] == b2[z], z = Z) %>%
      ompr::add_constraint(sum_expr(deviationB[z],z = Z)== 0) %>%
      ompr::add_constraint(deviationB[z]<= absdevB[z], z = Z) %>%
      ompr::add_constraint(deviationB[z]>= -absdevB[z], z = Z) %>%
      ompr::add_constraint(sum_expr(absdevB[z],z=Z)<= maxrelax/2.0) %>%

      ompr::solve_model(with_ROI(solver = "glpk"))

    solution       = ompr::get_solution(result, transportA[y,z])

    transportA_val = matrix(solution$value, length(Y),length(Z))


    b1 = numeric(length(Y))

    for (y in Y){

      stc_sum2 = vector(length = nbX)

      for (i in 1:nbX){

        stc_sum2[i] = ifelse(length(indXA[[i]])==0,1/length(Y),length(indXA[[i]][inst$Yobserv[indXA[[i]]] == y])/ length(indXA[[i]]))*length(indXB[[i]])/nB

      }

      b1[y]= sum(stc_sum2)}
      b1   = b1/sum(b1)

      result <-  MIPModel() %>%

        ompr::add_variable(transportB[y,z],y = Y, z= Z, type = "continuous",lb=0) %>%
        ompr::add_variable(deviationA[y]  , y = Y     , type = "continuous") %>%
        ompr::add_variable(absdevA[y]     , y = Y     , type = "continuous",lb=0) %>%

        ompr::set_objective(sum_expr(C[y,z]*transportB[y,z],y = Y,z=Z) , "min") %>%

        ompr::add_constraint(sum_expr(transportB[y,z], z = Z) - deviationA[y] == b1[y], y =Y) %>%
        ompr::add_constraint(sum_expr(transportB[y,z],y = Y) == freqZ[z], z = Z) %>%
        ompr::add_constraint(sum_expr(deviationA[y],y = Y)== 0) %>%
        ompr::add_constraint(deviationA[y] <= absdevA[y],y =Y) %>%
        ompr::add_constraint(deviationA[y] >= -absdevA[y],y =Y) %>%
        ompr::add_constraint(sum_expr(absdevA[y],y = Y)<= maxrelax/2.0) %>%

        ompr::solve_model(with_ROI(solver = "glpk"))


      solution       = ompr::get_solution(result, transportB[y,z])
      transportB_val = matrix(solution$value, length(Y),length(Z))


  }

  ####
  # Get the individual transport from the group transport

  if (indiv_method == "sequential"){

    YApred = indiv_grp_closest(inst, transportA_val, transportB_val, percent_closest = percent_c)$YAtrans
    YBpred = indiv_grp_closest(inst, transportA_val, transportB_val, percent_closest = percent_c)$ZBtrans

  } else if (indiv_method == "optimal"){

    YApred = indiv_grp_optimal(inst, transportA_val, transportB_val, percent_closest = percent_c)$YAtrans
    YBpred = indiv_grp_optimal(inst, transportA_val, transportB_val, percent_closest = percent_c)$ZBtrans
  }

  # Compute the estimated probability distributions from predictions
  indXA = inst$indXA; indXB = inst$indXB;
  nbX = length(indXA);

  estimatorZA = array(rep(0,nbX*length(Y)*length(Z)),dim = c(nbX,length(Y),length(Z)))
  # estimatorYB = estimatorZA
  estimatorYB = array(rep(0,nbX*length(Y)*length(Z)),dim = c(nbX,length(Z),length(Y)))


  for (x in 1:nbX){

    for (i in indXA[[x]]){

      estimatorZA[x,inst$Yobserv[i],YBpred[i]] = estimatorZA[x,inst$Yobserv[i],YBpred[i]] +  1/sum(inst$Yobserv[indXA[[x]]] == inst$Yobserv[i])
    }


    for (y in Y){

      if (sum(inst$Yobserv[indXA[[x]]] == y) == 0){

        estimatorZA[x,y,] = 1/length(Z)*rep(1,length(Z))
      }

    }

    for (i in indXB[[x]]){

      # estimatorYB[x,YApred[i],inst$Zobserv[i+nA]] = estimatorYB[x,YApred[i],inst$Zobserv[i+nA]] + 1/ sum(inst$Zobserv[indXB[[x]] + nA] == inst$Zobserv[i + nA])
      estimatorYB[x,inst$Zobserv[i+nA],YApred[i]] = estimatorYB[x,inst$Zobserv[i+nA],YApred[i]] + 1/ sum(inst$Zobserv[indXB[[x]] + nA] == inst$Zobserv[i + nA])
    }

    for (z in Z){

      if (sum(inst$Zobserv[indXB[[x]]+inst$nA] == z) == 0){

        # estimatorYB[x,,z] = 1/length(Y)*rep(1,length(Y))
        estimatorYB[x,z,] = 1/length(Y)*rep(1,length(Y))

      }
    }
  }

  row.names(estimatorZA) = row.names(estimatorYB) = ID_prof
  colnames(estimatorZA)  = as.character(levels(dataB[,2]))
  colnames(estimatorYB)  = as.character(levels(dataB[,3]))

  DATA1_OT         = dataB[dataB[,1] == unique(dataB[,1])[1],]
  DATA1_OT$OTpred  = as.factor(plyr::mapvalues(YBpred,from = sort(unique(YBpred)), to = levels(dataB[,3])[sort(unique(YBpred))]))

  if (is.ordered(dataB[,3])){

    DATA1_OT$OTpred = as.ordered(DATA1_OT$OTpred)

  } else {}


  DATA2_OT         = dataB[dataB[,1] == unique(dataB[,1])[2],]
  DATA2_OT$OTpred  = as.factor(plyr::mapvalues(YApred,from = sort(unique(YApred)), to = levels(dataB[,2])[sort(unique(YApred))]))

  if (is.ordered(dataB[,2])){

    DATA2_OT$OTpred = as.ordered(DATA2_OT$OTpred)

  } else {}


  tend = Sys.time()

  return(list(TIME_EXE = difftime(tend,tstart),TRANSPORT_A = transportA_val,TRANSPORT_B = transportB_val,
              profile = data.frame(ID = ID_prof,prof),estimatorZA= estimatorZA,estimatorYB = estimatorYB,DATA1_OT = DATA1_OT,DATA2_OT  = DATA2_OT))
}




