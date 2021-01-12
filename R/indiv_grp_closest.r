
#' indiv_grp_closest()
#'
#' This function sequentially assigns individual predictions using a nearest neighbors procedure to solve recoding problems of data fusion.
#'
#' A. THE RECODING PROBLEM IN DATA FUSION
#'
#' Assuming that \eqn{Y} and \eqn{Z} are two variables which refered to the same target population in two separate databases A and B respectively (no overlapping rows),
#' so that \eqn{Y} and \eqn{Z} are never jointly observed. Assuming also that A and B share a subset of common covariates \eqn{X} of any types (same encodings in A and B)
#' completed or not. Integrating these two databases often requires to solve the recoding problem by creating an unique database where
#' the missing information of \eqn{Y} and \eqn{Z} is fully completed.
#'
#'
#' B. DESCRIPTION OF THE FUNCTION
#'
#' The function \code{indiv_grp_closest} is an intermediate function used in the implementation of an algorithm called {OUTCOME} (and its enrichment {R-OUTCOME}, see the reference (2) for more details) dedicated to the solving of recoding problems in data fusion using Optimal Transportation theory.
#' The model is implemented in the function \code{\link{OT_outcome}} which integrates the function \code{indiv_grp_closest} in its syntax as a possible second step of the algorithm.
#' The function \code{indiv_grp_closest} can also be used separately provided that the argument \code{proxim} receives an output object of the function \code{\link{proxim_dist}}.
#' This latter is available in the package and is so directly usable beforehand.
#'
#' The algorithms \code{OUTCOME} (and \code{R-OUTCOME}) are made of two independent parts. Assuming that the objective consists in the prediction of \eqn{Z} in the database A:
#' \itemize{
#' \item The first part of the algorithm solves the optimization problem by providing a solution called \eqn{\gamma} that corresponds here to an estimation of the joint distribution \eqn{(Y,Z)} in A.
#' \item From the first part, a nearest neighbor procedure is carried out as a second part to provide the individual predictions of \eqn{Z} in A: this procedure is implemented in the function \code{indiv_group_closest}.
#' In other words, this function sequentially assigns to each individual of A the modality of \eqn{Z} that is closest.
#' }
#' Obviously, this algorithm  runs in the same way for the prediction of \eqn{Y} in the database B.
#' The function \code{indiv_grp_closest} integrates in its syntax the function \code{\link{avg_dist_closest}}. Therefore, the related argument \code{percent_closest} is identical in the two functions.
#' Thus, when computing average distances between an individual \eqn{i} and a subset of individuals assigned to a same level of \eqn{Y} or \eqn{Z} is required, user can decide if all individuals from the subset of interest can participate to the computation (\code{percent_closest}=1) or only a fixed part p (<1) corresponding to the closest neighbors of \eqn{i} (in this case \code{percent_closest} = p).
#'
#' The arguments \code{jointprobaA} and \code{jointprobaB} correspond to the estimations of \eqn{\gamma} (sum of cells must be equal to 1) in A and/or B respectively, according to the \code{which.DB} argument.
#' For example, assuming that \eqn{n_{Y_1}} individuals are assigned to the first modality of \eqn{Y} in A, the objective consists in the individual predictions of \eqn{Z} in A. Then, if \code{jointprobaA}[1,2] = 0.10,
#' the maximum number of individuals that can be assigned to the second modality of \eqn{Z} in A, can not exceed \eqn{0.10 \times n_A}.
#' If \eqn{n_{Y_1} \leq 0.10 \times n_A} then all individuals assigned to the first modality of \eqn{Y} will be assigned to the second modality of \eqn{Z}.
#' At the end of the process, each individual with still no affectation  will receive the same modality of \eqn{Z} as those of his nearest neighbor in B.
#'
#' @param proxim a \code{\link{proxim_dist}} object or an object of similar structure
#' @param jointprobaA a matrix whose number of columns corresponds to the number of modalities of the target variable \eqn{Y} in database A, and which number of rows corresponds to the number of modalities of Z in database B. It gives an estimation of the joint probability of \eqn{(Y,Z)} in A.
#' The sum of cells of this matrix must be equal to 1
#' @param jointprobaB a matrix whose number of columns equals to the number of modalities of the target variable \eqn{Y} in database A, and which number of rows corresponds to the number of modalities of \eqn{Z} in database B. It gives an estimation of the joint probability of \eqn{(Y,Z)} in B.
#' The sum of cells of this matrix must be equal to 1
#' @param percent_closest a value between 0 and 1 (by default) corresponding to the fixed \code{percent closest} of individuals remained in the computation of the average distances
#' @param which.DB a character string (with quotes) that indicates which individual predictions need to be computed: only the individual predictions of \eqn{Y} in B ("B"), only those of \eqn{Z} in A ("A") or the both ("BOTH" by default)
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @seealso \code{\link{proxim_dist}},\code{\link{avg_dist_closest}}, ,\code{\link{OT_outcome}}
#'
#' @return A list of two vectors of numeric values:
#' \item{YAtrans}{a vector corresponding to the individual predictions of \eqn{Y} (numeric form) in the database B using the Optimal Transportation algorithm}
#' \item{ZBtrans}{a vector corresponding to the individual predictions of \eqn{Z} (numeric form) in  the database A using the Optimal Transportation algorithm}
#'
#' @references
#' \enumerate{
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679. doi:10.1515/ijb-2018-0106
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association. \doi{10.1080/01621459.2020.1775615}
#' }
#'
#' @export
#'
#' @examples
#' data(simu_data)
#'
#' ### Example with the Manhattan distance
#'
#' try1 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),
#'                     ordinal = c(2,6), logic = NULL, prep_choice = "M")
#' res1 = proxim_dist(try1,norm = "M")
#'
#' ### Y(Yb1) and Z(Yb2) are a same information encoded in 2 different forms:
#' ### (3 levels for Y and 5 levels for Z)
#' ### ... Stored in two distinct databases, A and B, respectively
#' ### The marginal distribution of Y in B is unknown,
#' ### as the marginal distribution of Z in A ...
#'
#' # Empirical distribution of Y in database A:
#' freqY = prop.table(table(try1$Y)); freqY
#'
#' # Empirical distribution of Z in database B
#' freqZ = prop.table(table(try1$Z)); freqZ
#'
#' # By supposing that the following matrix called transport symbolizes
#' # an estimation of the joint distribution L(Y,Z) ...
#' # Note that, in reality this distribution is UNKNOWN and is
#' # estimated in the OT function by resolving an optimisation problem.
#'
#' transport1 = matrix(c(0,0.35285714,0,0.09142857,0,0.03571429,
#'                       0,0,0.08285714,0,0.07857143,0.03142857,
#'                       0.32714286,0,0),ncol = 5,byrow = FALSE)
#'
#' # ... So that the marginal distributions of this object corresponds to freqY and freqZ:
#' apply(transport1,1,sum)  # = freqY
#' apply(transport1,2,sum)  # = freqZ
#'
#' # The affectation of the predicted values of Y in database B and Z in database A
#' # are stored in the following object:
#'
#' res2      = indiv_grp_closest(res1,jointprobaA = transport1, jointprobaB = transport1,
#'                               percent_closest= 0.90)
#' summary(res2)
#'
#' # For the prediction of Z in A only, add the corresponding argument:
#' res3 = indiv_grp_closest(res1,jointprobaA = transport1, jointprobaB = transport1,
#' percent_closest= 0.90,which.DB="A")
#'
indiv_grp_closest=function(proxim, jointprobaA = NULL, jointprobaB = NULL, percent_closest = 1.0, which.DB = "BOTH"){

  if (!is.list(proxim)){

    stop("This object must be a list returned by the proxim_dist function")

  } else {}

  if ((!is.matrix(jointprobaA))|(!is.matrix(jointprobaA))){

    stop("The joint distributions must be store in matrix objects")

  } else {}

  if ((ncol(jointprobaA) != ncol(jointprobaB))|(nrow(jointprobaA) != nrow(jointprobaB))){

    stop("The joint distributions must be store in matrix of same size")

  } else {}


  if ((format(sum(jointprobaA)) != "1")| (format(sum(jointprobaB)) != "1")){

    stop("The sum of the jointprobaA matrix or the sum of the jointprobaB matrix differs from 1 !")

  } else {}


  if ((percent_closest>1)|(percent_closest <= 0)){

    stop("Improper value for the percent_closest argument")

  } else {}

  if (!(which.DB %in% c("A","B","BOTH"))){

    stop("Improper value for the which.DB argument")

  } else {}

  if ((which.DB == "BOTH") & (is.null(jointprobaA)|is.null(jointprobaB))){

    stop("The jointprobaA and B arguments must be filled")

  } else {}

  if ((which.DB == "A") & (is.null(jointprobaA))){

    stop("The jointprobaA argument is missing")

  } else {}

  if ((which.DB == "B") & (is.null(jointprobaB))){

    stop("The jointprobaB argument is missing")

  } else {}


  # Redefine A and B for the model

  A = 1:proxim$nA
  B = 1:proxim$nB
  Y = proxim$Y
  Z = proxim$Z
  Yobserv = proxim$Yobserv
  Zobserv = proxim$Zobserv
  indY = proxim$indY
  indZ = proxim$indZ
  nbindY = numeric(0)
  for (y in Y){
    nbindY = c(nbindY,length(proxim$indY[[y]]))}
  nbindZ= numeric(0)
  for (z in Z){
    nbindZ = c(nbindZ,length(proxim$indZ[[z]]))}
  freqY= numeric(0)
  for (y in Y){
    freqY = c(freqY,nbindY[y] / length(A))}
  freqZ= numeric(0)
  for (z in Z){
    freqZ = c(freqZ,nbindZ[z] / length(B))}

  # In essence, assign to each individual the modality that is closest,
  # where the distance from an individual to a modality is computed as the
  # average distance to the individuals having this modality (in the other base):

  YAtrans          = NULL
  YBtrans          = NULL
  average_distance = avg_dist_closest(proxim, percent_closest);
  Davg             = average_distance$Davg
  DindivA          = average_distance$DindivA
  DindivB          = average_distance$DindivB


  DA = DB =  list()

  for (z in Z){

    DA[[z]] = vector(length = length(A))

    for (i in A){

      DA[[z]][i] = DindivA[i,z]

    }
  }

  for (y in Y){

    DB[[y]] = vector(length = length(B))

    for (i in B){

      DB[[y]][i] = DindivB[i,y]

    }
  }

  if (which.DB %in% c("BOTH","A")){

    YBtrans = numeric(proxim$nA);

    for (y in Y){

      indtrans = indY[[y]]

      for (z in Z){

        nbtrans = min(round(jointprobaA[y,z]/freqY[y] * nbindY[y]),length(indtrans));
        distance = numeric(0)

        for (i in indtrans){

          distance = c(distance,DA[[z]][i])

        }

        names(distance)          = indtrans
        distance                 = sort(distance)
        ident_dist               = as.numeric(names(distance))

        for (k in 1:nbtrans){

          YBtrans[ident_dist[k]] = z;
          indtrans               = setdiff(indtrans,ident_dist[k])

        }

      }


      # affect potential individuals that have not been transported due to
      # rounding
      for (i in indtrans){

        YBtrans[i] = Zobserv[which.min(proxim$D[i,])]

      }
    }

  } else {}

  if (which.DB %in% c("BOTH","B")){

    YAtrans = numeric(proxim$nB)

    for (z in Z){

      indtrans = indZ[[z]]

      for (y in Y){

        nbtrans = min(round(jointprobaB[y,z]/freqZ[z] * nbindZ[z]), length(indtrans));
        distance = numeric(0)

        for (j in indtrans){

          distance = c(distance,DB[[y]][j])

        }

        names(distance)          = indtrans
        distance                 = sort(distance)
        ident_dist               = as.numeric(names(distance))

        for (k in 1:nbtrans){
          YAtrans[ident_dist[k]] = y
          indtrans = indtrans    = setdiff(indtrans,ident_dist[k])
        }
      }

      # affect potential individuals that have not been transported due to
      # rounding:

      for (j in indtrans){

        YAtrans[j] = Yobserv[which.min(proxim$D[,j])]

      }
    }

  } else {}

  return(list(YAtrans = YAtrans, ZBtrans = YBtrans))

}
