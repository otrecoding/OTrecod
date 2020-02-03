
<<<<<<< HEAD
#' individual_from_group_optimal()
#' 
#' Assuming that Y and Z summarize a same information encoded in two distinct forms stored in two distinct databases A and B.
#' This function is an alternative to the \code{\link{individual_from_group_closest}} approach.
#' It solves an optimization problem to get the individual transport that minimizes total distance while satisfying the joint probability computed by the model.
#' See section 5.3 (p24) of the reference article for more details about this method.
#' 
#'
#' @param inst An object corresponding to the output of the \code{Instance} function
#' @param jointprobaA A matrix which number of columns equals the number of modalities of the target variable Y in database A, and which number of rows equals the number of modalities of Z in database B. It gives an estimation of the joint probability (Y,Z) in DB A.
#' The sum of cells of this matrix must be equal to 1 
#' @param jointprobaB A matrix which number of columns equals the number of modalities of the target variable Y in database A, and which number of rows equals the number of modalities of Z in database B. It gives an estimation of the joint probability (Y,Z) in DB B.
#' The sum of cells of this matrix must be equal to 1 
#' @param percent_closest A value between 0 and 1 (by default) corresponding to the desired \code{percent closest} indivisuals taken in the computation of the distances
#'
#' @return A list of two vectors of numeric values:
#' \item{YAtrans}{A vector corresponding to the predicted values of Y in database B using the Optimal Transportation theory}
#' \item{YBtrans}{A vector corresponding to the predicted values of Z in database A using the Optimal Transportation theory}
#' 
#' 
#' @author Gregory Guernec, Valerie Gares, Jeremy Omer 
#' \email{gregory.guernec@@inserm.fr}
#' 
#' @seealso \code{\link{Instance}}, \code{\link{average_distance_to_closest}}, \code{\link{individual_from_group_closest}}
#' 
#' @import ompr ROI ROI.plugin.glpk 
#' @importFrom ompr.roi with_ROI
#' 
#' @references
#' Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
#' 0, 20180106 (2019) | \url{https://doi.org/10.1515/ijb-2018-0106}
#'
#'@example
#' # See the example of the merge_dbs() function to obtain the soluc1 object
#' try1 = prep_dbs(soluc1$DB_READY,nominal = 6,ordinal = 1:5) 
#' res1 = Instance(try1)                          # Using the Manhattan distance
#' 
#' ### Y and Z are the same variable encoded in 2 different forms (7 levels for Y and 4 levels for Z) stored in two distinct databases, A and B, respectively
#' # By supposing that the following matrix called transport symbolizes an estimation of the joint distribution L(Y,Z) ...
#' 
#' transport1 = matrix(c(0.19222222,0,0,0,0.20555556,0,0,0,0.01704523,0.1718437,0,0,0,0.1655457,0.01000983,0,0,0,0.10333333,0,0,0,0.08190020,0.005877581,0,0,0,0.046666667),ncol = 4,byrow = TRUE)
#' 
#' # The affectation of the predicted values of Y in database B and Z in database A are stored in the following object:
#' 
#' res4      = individual_from_group_optimal(res1,jointprobaA = transport1, jointprobaB = transport1, percent_closest= 0.90)
#' ### Be patient ... This execution can take a while (About 10 or 15 minutes)
#' summary(res4)
#' 
#' res2            = Instance(data_test,norme = 1)
#' res5            = individual_from_group_optimal(res2,jointprobaA = transport1, jointprobaB = transport1, percent_closest= 0.90)

individual_from_group_optimal=function(inst, jointprobaA, jointprobaB, percent_closest = 1.0){

  
  if (!is.list(inst)){
    
    stop("This object must be a list returned by the Instance function")
    
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
    
    stop("Incorrect value for the percent_closest option")
    
  } else {}
  
=======
#' individual_from_group_optimal(inst, jointprobaA, jointprobaB, percent_closest=1.0)
#'
#' @param inst todo list
#' @param jointprobaA todo list
#' @param jointprobaB todo list
#' @param percent_closest todo list
#'
#' @return todo list
#' 
#' @importFrom ompr get_solution
#' @importFrom ompr.roi with_ROI
#'
# @examples
individual_from_group_optimal=function(inst, jointprobaA, jointprobaB, percent_closest=1.0){
  
  
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
  # Redefine A and B for the model
  A = 1:inst$nA
  B = 1:inst$nB
  Y = inst$Y
  Z = inst$Z
  indY = inst$indY
  indZ = inst$indZ
  nbindY = numeric(0)
  for (y in Y){
    nbindY = c(nbindY,length(inst$indY[[y]]))}
  nbindZ= numeric(0)
  for (z in Z){
    nbindZ = c(nbindZ,length(inst$indZ[[z]]))}
  
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
  CA = matrix(rep(0,inst$nA * length(Z)), nrow = inst$nA,ncol = length(Z))
  CB = matrix(rep(0,inst$nB * length(Y)), nrow = inst$nB,ncol = length(Y))
<<<<<<< HEAD
  
=======
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
  for (i in A){
    for (z in Z){
      nbclose      = round(percent_closest*nbindZ[z])
      distance     = sort(inst$D[i,indZ[[z]]])
      CA[i,z]      = sum(distance[1:nbclose])/nbclose
    }
  }
<<<<<<< HEAD
  
=======
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
  for (j in B){
    for (y in Y){
      nbclose      = round(percent_closest*nbindY[y])
      distance     = sort(inst$D[indY[[y]],j])
      CB[j,y]      = CB[j,y] + sum(distance[1:nbclose])/nbclose
    }
  }
  
  # Objective: minimize the distance between individuals of A and B
  #@objective(indiv, Min, sum(CA[i,z]*assignA[i,z] for i in A, z in Z)
  #                        + sum(CB[j,y]*assignB[j,y] for j in B, y in Y))
  
  CAf <- function(i,z) {
    CA[i,z]
  }
  
  CBf <- function(j,y) {
    CB[j,y]
  }
  
<<<<<<< HEAD
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
=======
  result <-  ompr::MIPModel() %>%
    ompr::add_variable(assignA[i,z],  i = A, z = Z,type = "continuous",lb=0) %>%
    ompr::add_variable(assignB[j,y],  j = B, y = Y,type = "continuous",lb=0) %>%
    ompr::set_objective(ompr::sum_expr(CAf(i,z)*assignA[i,z], i = A,z=Z) + ompr::sum_expr(CBf(j,y)*assignB[j,y],j = B,y=Y), "min") %>%
    ompr::add_constraint(ompr::sum_expr(assignA[i,z], i=indY[[y]])   == jointprobaA[y,z], z = Z, y = Y) %>%
    ompr::add_constraint(ompr::sum_expr(assignB[j,y], j = indZ[[z]]) == jointprobaB[y,z], z = Z, y = Y) %>%
    ompr::add_constraint(ompr::sum_expr(assignA[i,z], z = Z) == 1/(length(A)),i = A) %>%
    ompr::add_constraint(ompr::sum_expr(assignB[j,y], y = Y) == 1/(length(B)),j = B) %>%
    ompr::solve_model(with_ROI(solver = "glpk"))
  
  solution = ompr::get_solution(result, assignA[i,z]) 
  assignA = matrix(solution$value, length(A),length(Z))
  solution = ompr::get_solution(result, assignB[j,y])  
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
  assignB=  matrix(solution$value, length(B),length(Y))
  
  
  # Extract the values of the solution
  # assignA_val = [value(assignA[i,z]) for i in A, z in Z]
  # assignB_val = [value(assignB[j,y]) for j in B, y in Y]
  
  
  
  # Transport the modality that maximizes frequency
  # YBtrans = [findmax([assignA_val[i,z]  for z in Z])[2] for i in A]             
  # YAtrans = [findmax([assignB_val[j,y]  for y in Y])[2] for j in B]
  YBtrans = apply(assignA,1,which.max)
  YAtrans = apply(assignB,1,which.max)
  
<<<<<<< HEAD
  return(list(YAtrans = YAtrans, ZBtrans = YBtrans))
=======
  return(list(YAtrans = YAtrans, YBtrans = YBtrans))
>>>>>>> 6f4cbad61224e0f392a0d80f4240360f0a0367ae
}
