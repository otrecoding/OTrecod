
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
  for (i in A){
    for (z in Z){
      nbclose      = round(percent_closest*nbindZ[z])
      distance     = sort(inst$D[i,indZ[[z]]])
      CA[i,z]      = sum(distance[1:nbclose])/nbclose
    }
  }
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
  assignB=  matrix(solution$value, length(B),length(Y))
  
  
  # Extract the values of the solution
  # assignA_val = [value(assignA[i,z]) for i in A, z in Z]
  # assignB_val = [value(assignB[j,y]) for j in B, y in Y]
  
  
  
  # Transport the modality that maximizes frequency
  # YBtrans = [findmax([assignA_val[i,z]  for z in Z])[2] for i in A]             
  # YAtrans = [findmax([assignB_val[j,y]  for y in Y])[2] for j in B]
  YBtrans = apply(assignA,1,which.max)
  YAtrans = apply(assignB,1,which.max)
  
  return(list(YAtrans = YAtrans, YBtrans = YBtrans))
}
