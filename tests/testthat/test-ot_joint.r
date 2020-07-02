test_that("ot_joint works", {

  data(tab_test)

  tab_test2 = tab_test[c(1:5000,5001:10000),1:5]

  try = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
                   dist.choice = "M", which.DB = "BOTH")

  print(try$GAMMA_A)
  ref_B = matrix(data=c(1.406767833, -0.072824738, -0.02409854, 0.260466308,
                        0.006403946,  0.90557653 ,  0.003684786,  1.036314488,
                        0.04408662,   0.073081074,  0.536106304,  0.41723539),  byrow = TRUE, nrow=4)

  ref_A = matrix(data=c(1.321600000,  0.000000000,  0.000000000,  0.464276415,
                        0.708279912,  0.009843674, -0.013584190,  0.809853194,
                        0.296130997, -0.012644716, 0.000000000,  1.045844716), byrow = TRUE, nrow=4)

  expect_equal(all(mapply(all.equal, try$gamma_A, ref_A, tolerance=1e-6)), TRUE)
  expect_equal(all(mapply(all.equal, try$gamma_B, ref_B, tolerance=1e-6)), TRUE)
})
