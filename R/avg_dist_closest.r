
#' avg_dist_closest()
#'
#'
#' Assuming that A and B are 2 databases with a common subset of variables and a same target factor (called Y in A and Z in B, such that Y is unknown in B and Z is unknown in A) whose encoding depends on the database (nA levels in A and nB levels in B). \code{avg_dist_closest} is an intermediate function of the Optimal Transportation algorithm (OT) that especially computes average distances between each row of one of the two databases and subsets of rows having a same level of target factor in the other database.
#'
#' \code{avg_dist_closest} is directly implemented in the \code{OT} and \code{OT_JOINT} functions but can also be used separately.
#' Nevertheless, this function can not be used without the \code{\link{proxim_dist}} function because it uses the objects from the latter.
#' So do not use this function directly on your data, and, if needed, do not hesitate to consult the example provided for a better understanding.
#'
#'
#' DEFINITION OF THE COST MATRIX
#'
#'
#' \code{avg_dist_closest} uses, in particular, the D distance matrix (that stores distances between rows of A and B) from \code{\link{proxim_dist}} to produce 3 distinct matrices stored in a list.
#' This matrix is used to calculate average distances between rows of A (resp. B) and rows of B (resp. A) related to each level of Z (resp. Y) in B (resp. A). This information is stored in \code{DinvidA} and \code{DinvidB} respectively.
#'
#' A cost matrix called \code{Davg} whose dimensions (nA*nB) correspond to the number of levels of Y (rows) and Z (columns), is also generated to calculate average distances between subsets of rows (whose size depends of the percent of nearest neighbours tolerated) with identical levels of Y in A and identical levels of Z in B.
#' Thus, Davg[i,j] corresponds to average distance between subsets of rows (from A) whose related level of Y is i, and a subset of rows (from B) whose related level of Z is j.
#'
#' Thus, this matrix can be seen as the ability for an individual (a row) to move from a given level of the target variable (Y) in A to given level of Z in the database B.
#'
#'
#' @param proxim An object corresponding to the output of the \code{proxim_dist} function
#' @param percent_closest A value between 0 and 1 (by default) corresponding to the desired percent of rows (or individuals) that will be taken into acccount for the computation of the distances
#'
#' @return A list of 3 matrix is returned:
#' \item{Davg}{The cost matrix whose number of rows corresponds to nA, the number of levels of the target variable in database A, and whose number of columns corresponds to nB: the number of levels of the target variable in B.
#' This cost matrix can be interpreted as the ability to move from one level of Y in A to one level of Z in the 2nd database B.
#' Davg[P,Q] so refers to the average distance between the modality P in Y (only known in A) and modality Q of Z (only known in B).}
#' \item{DindivA}{A matrix whose number of rows corresponds to the number of rows of the 1st database A and number of columns corresponds to nB, the number of levels of the target variable in the 2nd database B.
#' DindivA[i,Q] refers to the average distance between the i_th individual (or row) of the 1st database and a chosen proportion of individuals (\code{percent_closest} set by the user) of the 2nd database having the modality Q of Z.}
#' \item{DindivB}{A matrix whose number of rows corresponds to the number of rows of the 2nd database B and number of columns corresponds to nA, the number of levels of the target variable in the 1st database A.
#' DindivB[k,P] refers to the average distance between the k_th individual (or row) of the 2nd database and a chosen proportion of individuals (see \code{percent_closest}) of the 1st database having the modality P of Y.}
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#'
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
#' @export
#'
#' @examples
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



