
#' avg_dist_closest()
#'
#' \code{avg_dist_closest} is an intermediate function of the Optimal Transportation algorithm (OT) that especially computes average distances between each individual of 2 databases (A and B) and levels of target factors that synthetise a same information in 2 different encodings.
#'
#' \code{avg_dist_closest} is directly implemented in the \code{OT} and \code{OT_JOINT} functions but can also be used separately.
#' Nevertheless, this function can not be used without the \code{\link{proxim_dist}} function because it uses the objects from the latter.
#' So do not use this function directly on your data, and, if needed, do not hesitate to consult the example provided for a better understanding.
#'
#'
#' Definition of the cost matrix
#' ------------------------------
#'
#' \code{avg_dist_closest} uses, in particular, the D distance matrix (that stores distances between individuals of A and B) from \code{\link{proxim_dist}} to produce 3 distinct matrices stored in a list.
#' This matrix is used to calculate average distances between individuals in A (resp. B) and individuals of B (resp. A) related to each level of Z (resp. Y) in B (resp. A). This information is stored in \code{DinvidA} and \code{DinvidB} respectively.
#'
#' A cost matrix called \code{Davg} which dimensions correspond to the number of levels of Y (rows) and Z (columns), is also generated to calculate average distances between covariations of a given part of individuals (the percent closest neighbors)  with identical levels of Y in A and identical levels of Z in B.
#' Davg[i,j] corresponds to average distances between covariations of a given part of individuals (in A) which related values for Y are i, and another given part of  individuals (in B) which related values for Z are j.
#'
#' Thus, this matrix can be seen as the ability for an individual to move from a given class of a target variable (Y) in the 1st database (A) to another given class of the corresponding target variable (Z) in the 2nd database (B).
#'
#'
#' @param proxim An object corresponding to the output of the \code{proxim_dist} function
#' @param percent_closest A value between 0 and 1 (by default) corresponding to the desired percent of indivisuals that will be taken in the computation of the distances
#'
#' @return A list of 3 matrix is returned:
#' \item{Davg}{The cost matrix which number of rows corresponds to the number of levels of the target variable in database A, and which number of columns corresponds to the number of levels of the target variable in database B.
#' This cost matrix can be interpreted as the cost of moving from one level of the target factor in the 1st database A to one level of the target factor in the 2nd database B.
#' Davg[P,Q] so refers to the distance between covariates vectors of individuals in database A having modality P and individuals in database B having modality Q}
#' \item{DindivA}{A matrix which number of rows corresponds to the number of rows of the 1st database A and the number of columns corresponds to the number of levels of the target variable in the 2nd database B.
#' DindivA[i,Q] refers to the average distance of the i^th individual of the 1st database and a given proportion of individuals (\code{percent_closest} set by the user) of the 2nd database having the modality Q as value of the target variable}
#' \item{DindivB}{A matrix which number of rows corresponds to the number of rows of the 2nd database B and the number of columns corresponds to the number of levels of the target variable in the 1st database A.
#' DindivB[k,P] refers to the average distance of the k^th individual of the 1st database and a given proportion of individuals (\code{percent_closest} set by the user) of the 2nd database having the modality P as value of the target variable}
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#' \email{gregory.guernec@@inserm.fr}
#'
#' @references
#' Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' 0, 20180106 (2019) | \url{https://doi.org/10.1515/ijb-2018-0106}
#'
#' @seealso \code{\link{proxim_dist}}
#'
#' @aliases avg_dist_closest
#'
#' @example
#' data(simu_data)
#' ### The covariates of the data are prepared according to the distance chosen
#' ### using the transfo_dist function
#'
#' ### Example with The Manhattan distance
#'
#' try1 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),
#'                     ordinal = c(2,6), logic = NULL, prep_choice = "M")
#' res1 = proxim_dist(try1,norme = 1)  # norme = 2 for Euclidean
#'
#' # proxim_dist() fixes the distance measurement chosen,
#' # and defines neighborhoods between profiles and individuals
#'
#' res_new  = avg_dist_closest(res1,percent_closest = 0.80)
#'
#'
avg_dist_closest = function(proxim, percent_closest = 1){

  if (!is.list(proxim)){

    stop("This object must be a list returned by the proxim_dist function")

  } else {}

  if ((percent_closest>1)|(percent_closest <= 0)){

    stop("Incorrect value for the percent_closest option")

  } else {}


  # Redefine A and B for the model
  A = (1:proxim$nA)
  B = (1:proxim$nB)
  Y = proxim$Y
  Z = proxim$Z
  indY = proxim$indY
  indZ = proxim$indZ

  # Compute average distances
  Davg     = matrix(0,length(Y),length(Z));
  DindivA  = matrix(0,proxim$nA,length(Z));
  DindivB  = matrix(0,proxim$nB,length(Y));

  for (y in Y){
    for (i in indY[[y]]){
      for (z in Z){
        nbclose = max(round(percent_closest*length(indZ[[z]])),1);

        distance     = sort(proxim$D[i,indZ[[z]]])
        DindivA[i,z] = sum(distance[1:nbclose])/nbclose;
        Davg[y,z]    = Davg[y,z] + sum(distance[1:nbclose])/nbclose/length(indY[[y]])/2.0;
      }
    }
  }
  for (z in Z){
    for (j in indZ[[z]]){
      for (y in Y){
        nbclose      = max(round(percent_closest*length(indY[[y]])),1);

        distance     = sort(proxim$D[indY[[y]],j])
        DindivB[j,y] = sum(distance[1:nbclose])/nbclose;
        Davg[y,z]    = Davg[y,z] + sum(distance[1:nbclose])/nbclose/length(indZ[[z]])/2.0;
      }
    }
  }

  return(list(Davg=Davg, DindivA=DindivA, DindivB=DindivB))
}



