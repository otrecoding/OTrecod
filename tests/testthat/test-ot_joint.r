test_that("ot_joint works", {

  data(tab_test)

  tab_test2 = tab_test[c(1:5000,5001:10000),1:5]

  try = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
                   dist.choice = "M", which.DB = "BOTH")

  print(try$GAMMA_A)
  ref_B = matrix(data=c( 1.309845, 0.4341554, 0.0000000, 0.000000,
                          0.000000, 0.7382913, 0.7677087, 0.000000,
                          0.000000, 0.0000000, 0.3163772, 1.026423), nrow=4)

  ref_A = matrix(data=c(1.3216, 0.4380475, 0.0000000, 0.0000,
                         0.0000, 0.7443525, 0.7737806, 0.0000,
                         0.0000, 0.0000000, 0.3186194, 1.0332), nrow=4)

  expect_equal(all(mapply(all.equal, try$gamma_A, ref_A, tolerance=1e-6)), TRUE)
  expect_equal(all(mapply(all.equal, try$gamma_B, ref_B, tolerance=1e-6)), TRUE)
})
