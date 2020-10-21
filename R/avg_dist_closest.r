
#' avg_dist_closest()
#'
#' This function computes average distances between levels of two categorical variables located in two distinct databases.
#'
#'
#' The function \code{avg_dist_closest} is an intermediate function for the implementation of original algorithms dedicated to the solving of recoding problems in data fusion using Optimal Transportation theory (for more details, consult the corresponding algorithms called
#' \code{OUTCOME}, \code{R_OUTCOME}, \code{JOINT} and \code{R_JOINT}, in the reference (2)). The function \code{avg_dist_closest} is so directly implemented in the \code{OT_outcome} and \code{OT_joint} functions but can also be used separately.
#' The function \code{avg_dist_closest} uses, in particular, the distance matrix D (that stores distances between rows of A and B) from the function \code{\link{proxim_dist}} to produce three distinct matrices saved in a list object.
#' Therefore, the function requires in input, the specific output of the function \code{\link{proxim_dist}} which is available in the package and so must be used beforehand.
#' In consequence, do not use this function directly on your database, and do not hesitate to consult the provided examples provided for a better understanding.
#'
#'
#' DEFINITION OF THE COST MATRIX
#'
#' Assuming that A and B are two databases with a set of shared variables and that a same information (referred to a same target population) is stored as a variable \eqn{Y} in A and \eqn{Z} in B, such that \eqn{Y} is unknown in B and \eqn{Z} is unknown in A, whose encoding depends on the database (\eqn{n_Y} levels in A and \eqn{n_Z} levels in B).
#' A distance between one given level y of \eqn{Y} and one given level z of Z is estimated by averaging the distances between the two subsets of individuals (units or rows) assigned to y in A and z in B, characterized by their vectors of covariates.
#' The distance between two individuals depends on the variations between the shared covariates, and so depends on the chosen distance function using the function \code{proxim_dist}.
#' For these computations, all the individuals concerned by these two levels can be taken into account, or only a part of them, depending on the argument \code{percent_closest}.
#' When \code{percent_closest} < 1, the average distance between an individual i and a given level of factor z only uses the corresponding part of individuals related to z that are the closest to i.
#' Therefore, this choice influences the estimations of average distances between levels of factors but also permits to reduce time computation when necessary.
#'
#' The average distance between each individual of \eqn{Y} (resp. \eqn{Z}) and each levels of \eqn{Z} (resp. \eqn{Y}) are returned in output, in the object \code{DindivA} (\code{DindivB} respectively).
#' The average distance between each levels of \eqn{Y} and each levels of \eqn{Z} are returned in a matrix saved in output (the object \code{Davg}).
#' \code{Davg} returns the computation of the cost matrix D, whose dimensions (\eqn{n_Y \times n_Z}) correspond to the number of levels of \eqn{Y} (rows) and \eqn{Z} (columns).
#' This matrix can be seen as the ability for an individual (row) to move from a given level of the target variable (\eqn{Y}) in A to a given level of \eqn{Z} in the database B (or vice versa).
#'
#'
#' @param proxim A \code{proxim_dist} object
#' @param percent_closest A value between 0 and 1 corresponding to the desired percent of rows (or statistical units, or individuals) that will participate in the computation of the average distances between levels of factors or
#' between individuals and levels of only one factor. The default value 1 means that all rows (individuals) with a same level of \eqn{Y} or a same level of \eqn{Z} will be kept for the average computations. See 'Details'.
#'
#' @return A list of 3 matrices is returned:
#' \item{Davg}{The cost matrix whose number of rows corresponds to \eqn{n_Y}, the number of levels of the target variable \eqn{Y} in the database A, and whose number of columns corresponds to \eqn{n_Z}: the number of levels of the target variable in B.
#' In this case, the related cost matrix can be interpreted as the ability to move from one level of \eqn{Y} in A to one level of \eqn{Z} in B.
#' Davg[P,Q] refers to the average distance between the modality P of \eqn{Y} (only known in A) and modality \eqn{Q} of \eqn{Z} (only known in B).}
#' \item{DindivA}{A matrix whose number of rows corresponds to the number of rows of the first database A and number of columns corresponds to \eqn{n_Z}, the number of levels of the target variable \eqn{Z} in the second database B.
#' DindivA[i,Q] refers to the average distance between the \eqn{i^{th}} individual (or row) of the first database and a chosen proportion of individuals (\code{percent_closest} set by the user) of the second database having the modality \eqn{Q} of \eqn{Z}.}
#' \item{DindivB}{A matrix whose number of rows corresponds to the number of rows of the second database B and number of columns corresponds to nA, the number of levels of the target variable in the first database A.
#' DindivB[k,P] refers to the average distance between the \eqn{k^{th}} individual (or row) of the second database and a chosen proportion of individuals (depending on \code{percent_closest}) of the first database having the modality P of \eqn{Y}.}
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @references
#' \enumerate{
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679 | \url{https://doi.org/10.1515/ijb-2018-0106}
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association, doi: 10.1080/01621459.2020.1775615 | \url{https://doi.org/10.1080/01621459.2020.1775615}
#' }
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
#' res1 = proxim_dist(try1,norm = "M")
#'
#' # proxim_dist() fixes the chosen distance function,
#' # and defines neighborhoods between profiles and individuals
#'
#' # The following row uses only 80 percents of individuals of each level
#' # of factors for the computation of the average distances:
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



