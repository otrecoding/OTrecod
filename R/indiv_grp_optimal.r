
#' indiv_grp_optimal()
#'
#' This function assigns individual predictions related to recoding problems of data fusion by solving a linear optimization problem
#'
#' A. THE RECODING PROBLEM IN DATA FUSION
#'
#' Assuming that Y and Z are two variables that summarize a same latent information in two separate (no overlapping rows) databases A and B respectively,
#' so that Y and Z are never jointly observed in A and B. Assuming also that A and B share a subset of common covariates X of any types (but same encodings in A and B)
#' completed or not. Integrating these two databases often requires to solve the recoding problem observed between Y and Z by creating an unique database where
#' the missing information of Y and Z is completed.
#'
#'
#' B. DESCRIPTION OF THE FUNCTION
#'
#' The function \code{indiv_grp_optimal} is an intermediate function used in the implementation of a model called \code{OUTCOME} (and its enrichment \code{R-OUTCOME}, see the reference (2) for more details) dedicated to the solving of recoding problems in data fusion using Optimal Transportation theory.
#' The model is implemented in the function \code{\link{OT_outcome}} which integrates the function \code{indiv_grp_optimal} in its syntax as a possible second step of the algorithm.
#' The function \code{indiv_grp_optimal} can nevertheless be used separately provided that the argument \code{proxim} receives an output object of the function \code{\link{proxim_dist}}.
#' This latter is available in the package and is so directly usable beforehand.
#'
#' The function \code{indiv_grp_optimal} constitutes an alternative method to the nearest neighbor procedure implemented in the function \code{\link{indiv_grp_closest}}.
#' As for the function \code{\link{indiv_grp_closest}}, assuming that the objective consists in the prediction of Z in the database A, the first step of the algorithm related to \code{OUTCOME} provides an estimate of \eqn{\gamma} which can is an estimation of the joint distribution \code{(Y,Z)} in A.
#' Rather than using a nearest neighbor approach to provide individual predictions,the function \code{indiv_grp_optimal} solves an optimization problem using the simplex algorithm which searches for the individual predictions of Z that minimize the computed total distance satisfying the joint probability distribution estimated in the first part.
#' More details about the theory related to the solving of this optimization problem is described in the section 5.3 of (2).
#'
#' Obviously, this algorithm  runs in the same way for the prediction of Y in the database B.
#' The function \code{indiv_grp_optimal} integrates in its syntax the function \code{\link{avg_dist_closest}} and the related argument \code{percent_closest} is identical in the two functions.
#' Thus, when computing average distances between an individual i and a subset of individuals assigned to a same level of Y or Z is required, user can decide if all individuals from the subset of interest can participate to the computation (\code{percent_closest = 1}) or only a fixed part p (<1) corresponding to the closest neighbors of i (in this case \code{percent_closest} = p).
#'
#' The arguments \code{jointprobaA} and \code{jointprobaB} can be seen as estimations of \eqn{\gamma} (sum of cells must be equal to 1) that correspond to estimations of the joint distributions of \code{(Y,Z)} in A and B respectively.
#'
#' @param proxim An object corresponding to the output of the function \code{\link{proxim_dist}}
#' @param jointprobaA A matrix whose number of columns is equal to the number of modalities of the target variable Y in database A, and which number of rows is equal to the number of modalities of Z in database B. It gives an estimation of the joint probability (Y,Z) in the database A.
#' The sum of cells of this matrix must be equal to 1
#' @param jointprobaB A matrix whose number of columns is equal to the number of modalities of the target variable Y in database A, and which number of rows is equal to the number of modalities of Z in database B. It gives an estimation of the joint probability (Y,Z) in the database B.
#' The sum of cells of this matrix must be equal to 1
#' @param percent_closest A value between 0 and 1 (by default) corresponding to the fixed \code{percent closest} of individuals used in the computation of the average distances
#' @param which.DB A character string (with quotes) that indicates which individual predictions compute: only the individual predictions of Y in B ("B"), only those of Z in A ("A") or the both ("BOTH" by default)
#'
#' @return A list of two vectors of numeric values:
#' \item{YAtrans}{A vector corresponding to the predicted values of Y in database B (numeric form) using the Optimal Transportation theory}
#' \item{ZBtrans}{A vector corresponding to the predicted values of Z in database A (numeric form) using the Optimal Transportation theory}
#'
#'
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer
#'
#' \email{otrecod.pkg@@gmail.com}
#'
#' @seealso \code{\link{proxim_dist}}, \code{\link{avg_dist_closest}}, \code{\link{indiv_grp_closest}}
#'
#' @import ompr ROI ROI.plugin.glpk
#' @importFrom ompr.roi with_ROI
#' @importFrom dplyr %>%
#'
#' @references
#' \enumerate{
#' \item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' Volume 16, Issue 1, 20180106, eISSN 1557-4679 | \url{https://doi.org/10.1515/ijb-2018-0106}
#' \item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association, DOI: 10.1080/01621459.2020.1775615
#' }
#'
#' @aliases indiv_grp_optimal
#'
#' @export
#'
#' @examples
#'
#' ### Example using The Euclidean distance on a complete database
#' # For this example we keep only 200 rows:
#'
#' data(tab_test)
#' tab_test2 = tab_test[c(1:80,5001:5080),]; dim(tab_test2)
#'
#' # Adding NAs in Y1 and Y2
#' tab_test2[tab_test2$ident == 2, 2] = NA
#' tab_test2[tab_test2$ident == 1, 3] = NA
#'
#' # Because all covariates are ordered in numeric form,
#' # the transfo_dist function is not required here
#'
#' res3      = proxim_dist(tab_test2,norm = "M")
#'
#' #' ### Y(Y1) and Z(Y2) are a same variable encoded in 2 different forms:
#' ### 4 levels for Y1 and 3 levels for Y2
#' ### ... Stored in two distinct databases, A and B, respectively
#' ### The marginal distribution of Y in B is unknown,
#' ### as the marginal distribution of Z in A ...
#'
#' # Assuming that the following matrix called transport symbolizes
#' # an estimation of the joint distribution L(Y,Z) ...
#' # Note that, in reality this distribution is UNKNOWN and is
#' # estimated in the OT function by resolving the optimization problem.
#'
#' # By supposing:
#'
#' val_trans  = c(0.275,0.115,0,0,0,0.085,0.165,0,0,0,0.095,0.265)
#' transport2 = matrix(val_trans,ncol = 3,byrow = FALSE)
#'
#' # Getting the individual predictions of Z in A (only)
#' # by computing average distances on 90% of the nearest neighbors of
#' # each modality of Z in B
#' res4      = indiv_grp_optimal(res3,jointprobaA = transport2,
#'             jointprobaB = transport2, percent_closest= 0.90,
#'             which.DB = "A")
#'
#' \dontrun{
#' ### Example 2 using The Manhattan distance with incomplete covariates
#' data(simu_data)
#'
#' try1 = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),
#'                     ordinal = c(2,6), logic = NULL, prep_choice = "M")
#' res1 = proxim_dist(try1,norm = "M")
#'
#'
#' ### Y and Z are a same variable encoded in 2 different forms:
#' ### (3 levels for Y and 5 levels for Z)
#' ### ... Stored in two distinct databases, A and B, respectively
#' ### The marginal distribution of Y in B is unknown,
#' ### as the marginal distribution of Z in A ...
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
#'
#' # The predicted values of Y in database B and Z in
#' database A are stored in the following object:
#'
#' res2 = indiv_grp_optimal(res1,jointprobaA = transport1,
#'                          jointprobaB = transport1,
#'                          percent_closest= 0.90)                               )
#' summary(res2)
#' }
#'

indiv_grp_optimal =function(proxim, jointprobaA, jointprobaB, percent_closest = 1.0, which.DB = "BOTH"){


  if (!is.list(proxim)){

    stop("This object must be a list returned by the proxim_dist function")

  } else {}

  if ((!is.matrix(jointprobaA))|(!is.matrix(jointprobaA))){

    stop("The joint distributions must be store in matrix objects")

  } else {}

  if ((ncol(jointprobaA) != ncol(jointprobaB))|(nrow(jointprobaA) != nrow(jointprobaB))){

    stop("The joint distributions must be store in matrix of same size")

  } else {}

  # if (!(solvR %in% c("glpk","cbc","GLPK","CBC"))){
  #
  #  stop("Improper argument for solvR")

  #} else {}


  if ((format(sum(jointprobaA)) != "1")| (format(sum(jointprobaB)) != "1")){

    stop("The sum of the jointprobaA matrix or the sum of the jointprobaB matrix differs from 1 !")

  } else {}


  if ((percent_closest>1)|(percent_closest <= 0)){

    stop("Incorrect value for the percent_closest option")

  } else {}

  # Redefine A and B for the model
  A      = 1:proxim$nA
  B      = 1:proxim$nB
  Y      = proxim$Y
  Z      = proxim$Z
  indY   = proxim$indY
  indZ   = proxim$indZ
  nbindY = numeric(0)

  for (y in Y){

    nbindY = c(nbindY,length(proxim$indY[[y]]))}
  nbindZ= numeric(0)

  for (z in Z){

    nbindZ = c(nbindZ,length(proxim$indZ[[z]]))}

  # Create a model for the optimal transport of individuals
  # indiv = Model(solver=IpoptSolver(print_level=4))
  # indiv = Model(with_optimizer(Clp.Optimizer,LogLevel=0))
  # indiv = Model(solver=GurobiSolver(Method=2,LogToConsole=0)); #Presolve=0,Method=2,Crossover=0))


  # Variables
  # - assignA[i][z] : fraction of individual i assigned to modality z
  # @variable(indiv, assignA[i in A, z in Z] >= 0, base_name="assignA")
  # - assignB[j][y] : fraction of individual j assigned to modality y
  # @variable(indiv, assignB[j in B, y in Y] >= 0, base_name="assignB")

  # compute the average distance between the individuals and the modalities of
  # the other base
  CA = matrix(rep(0,proxim$nA * length(Z)), nrow = proxim$nA,ncol = length(Z))
  CB = matrix(rep(0,proxim$nB * length(Y)), nrow = proxim$nB,ncol = length(Y))

  for (i in A){

    for (z in Z){

      nbclose      = round(percent_closest*nbindZ[z])
      distance     = sort(proxim$D[i,indZ[[z]]])
      CA[i,z]      = sum(distance[1:nbclose])/nbclose

    }
  }

  for (j in B){

    for (y in Y){

      nbclose      = round(percent_closest*nbindY[y])
      distance     = sort(proxim$D[indY[[y]],j])
      CB[j,y]      = CB[j,y] + sum(distance[1:nbclose])/nbclose

    }
  }

  # Objective: minimize the distance between individuals of A and B
  #@objective(indiv, Min, sum(CA[i,z]*assignA[i,z] for i in A, z in Z)
  #                        + sum(CB[j,y]*assignB[j,y] for j in B, y in Y))

  CAf <- function(i,z){

    CA[i,z]

  }

  CBf <- function(j,y){

    CB[j,y]

  }

  result <-  MIPModel() %>%
    add_variable(assignA[i,z],  i = A, z = Z,type = "continuous",lb=0) %>%
    add_variable(assignB[j,y],  j = B, y = Y,type = "continuous",lb=0) %>%
    set_objective(sum_expr(CAf(i,z)*assignA[i,z], i = A,z=Z) + sum_expr(CBf(j,y)*assignB[j,y],j = B,y=Y), "min") %>%
    add_constraint(sum_expr(assignA[i,z], i=indY[[y]])   == jointprobaA[y,z], z = Z, y = Y) %>%
    add_constraint(sum_expr(assignB[j,y], j = indZ[[z]]) == jointprobaB[y,z], z = Z, y = Y) %>%
    add_constraint(sum_expr(assignA[i,z], z = Z) == 1/(length(A)),i = A) %>%
    add_constraint(sum_expr(assignB[j,y], y = Y) == 1/(length(B)),j = B) %>%
    solve_model(with_ROI(solver = "glpk"))

  solution = get_solution(result, assignA[i,z])
  assignA = matrix(solution$value, length(A),length(Z))

  solution = get_solution(result, assignB[j,y])
  assignB=  matrix(solution$value, length(B),length(Y))


  # Extract the values of the solution
  # assignA_val = [value(assignA[i,z]) for i in A, z in Z]
  # assignB_val = [value(assignB[j,y]) for j in B, y in Y]



  # Transport the modality that maximizes frequency
  # YBtrans = [findmax([assignA_val[i,z]  for z in Z])[2] for i in A]
  # YAtrans = [findmax([assignB_val[j,y]  for y in Y])[2] for j in B]
  YBtrans = apply(assignA,1,which.max)
  YAtrans = apply(assignB,1,which.max)

  if (which.DB == "A"){

    YAtrans = NULL

  } else if (which.DB == "B"){

    YBtrans = NULL

  } else {}

  return(list(YAtrans = YAtrans, ZBtrans = YBtrans))
}
