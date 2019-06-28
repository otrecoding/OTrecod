result <-  MIPModel() %>%
  add_variable(gammaA[x,y,z],x = 1:nbX, y=Y,z= Z,type = "continuous") %>%
  add_variable(errorA_XY[x,y],x=1:nbX, y = Y, type = "continuous") %>%
  add_variable(abserrorA_XY[x,y],x=1:nbX, y = Y, type = "continuous") %>%
  add_variable(errorA_XZ[x,z],x=1:nbX, z = Z, type = "continuous") %>%
  add_variable(abserrorA_XZ[x,z],x=1:nbX, z = Z, type = "continuous") %>%
  add_variable(reg_absA[x1, x2,y,z],x1=1:nbX, x2= voisins_X[x1,],y=Y,z= Z, type = "continuous") %>%
  set_objective(sum_expr(C[y,z] * gammaA[x,y,z],y = Y,z=Z, x = 1:nbX) + lambda_reg * sum_expr(1/length(voisins_X[x1,]) *reg_absA[x1,x2,y,z]) , "min") %>%
  add_constraint(sum_expr(gammaA[x,y,z], z = Z) == estim_XA_YA[x,y] + errorA_XY[x,y], x = 1:nbX,y =Y) %>%
  add_constraint(estim_XB[x]*sum_expr(gammaA[x,y,z],y = Y) == estim_XB_ZB[x,z] * estim_XA[x] + estim_XB[x]*errorA_XZ[x,z], x = 1:nbX, z = Z) %>%
  add_constraint(errorA_XY[x,y] <= abserrorA_XY[x,y], x = 1:nbX,y =Y) %>%
  add_constraint(-errorA_XY[x,y] <= abserrorA_XY[x,y], x = 1:nbX,y =Y) %>%
  add_constraint(sum_expr(abserrorA_XY[x,y], x = 1:nbX,y = Y)<= maxrelax/2.0) %>%
  add_constraint(sum_expr(errorA_XY[x,y], x = 1:nbX,y = Y)== 0.0) %>%
  
  add_constraint(errorA_XZ[x,z] <= abserrorA_XZ[x,z], x = 1:nbX,z =Z) %>%
  add_constraint -errorA_XZ[x,z] <= abserrorA_XZ[x,z], x = 1:nbX,z =Z) %>%
  add_constraint(sum_expr(abserrorA_XZ[x,z], x = 1:nbX,z = Z)<= maxrelax/2.0) %>%
  add_constraint(sum_expr(errorA_XZ[x,z], x = 1:nbX,z = Z)<= maxrelax/2.0) %>%
  solve_model(with_ROI(solver = "glpk"))

solution <- get_solution(result, gammaA[x,y,z]) 


result <-  MIPModel() %>%
  add_variable(gammaB[x,y,z],x = 1:nbX, y=Y,z= Z,type = "continuous") %>%
  add_variable(errorB_XY[x,y],x=1:nbX, y = Y, type = "continuous") %>%
  add_variable(abserrorB_XY[x,y],x=1:nbX, y = Y, type = "continuous") %>%
  add_variable(errorB_XZ[x,z],x=1:nbX, z = Z, type = "continuous") %>%
  add_variable(abserrorB_XZ[x,z],x=1:nbX, z = Z, type = "continuous") %>%
  add_variable(reg_absB[x1, x2,y,z],x1=1:nbX, x2= voisins_X[x1,],y=Y,z= Z, type = "continuous") %>%
  set_objective(sum_expr(C[y,z] * gammaB[x,y,z],y = Y,z=Z, x = 1:nbX) + lambda_reg *  sum(1/length(voisins_X[x1,]) *reg_absB[x1,x2,y,z] ) , "min") %>%
  add_constraint(sum_expr(gammaB[x,y,z], y = Y) == estim_XB_ZB[x,z] + errorB_XZ[x,z], x = 1:nbX,z =Z) %>%
  add_constraint(estim_XA[x]*sum_expr(gammaB[x,y,z] ,z = Z) == estim_XA_YA[x,y] * estim_XB[x] + estim_XA[x] * errorB_XY[x,y], x = 1:nbX, y = Y) %>%
  add_constraint(errorB_XY[x,y] <= abserrorB_XY[x,y], x = 1:nbX,y =Y) %>%
  add_constraint(-errorB_XY[x,y] <= abserrorB_XY[x,y], x = 1:nbX,y =Y) %>%
  add_constraint(sum_expr( abserrorB_XY[x,y], x = 1:nbX,y = Y)<= maxrelax/2.0) %>%
  add_constraint(sum_expr(errorB_XY[x,y], x = 1:nbX,y = Y)== 0.0) %>%
  
  add_constraint(errorB_XZ[x,z] <= abserrorB_XZ[x,z], x = 1:nbX,z =Z) %>%
  add_constraint -errorB_XZ[x,z] <= abserrorB_XZ[x,z], x = 1:nbX,z =Z) %>%
  add_constraint(sum_expr(abserrorB_XZ[x,z], x = 1:nbX,z = Z)<= maxrelax/2.0) %>%
  add_constraint(sum_expr(errorB_XZ[x,z], x = 1:nbX,z = Z)<= maxrelax/2.0) %>%
  solve_model(with_ROI(solver = "glpk"))

solution <- get_solution(result, gammaB[x,y,z]) 