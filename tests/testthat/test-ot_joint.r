test_that("ot_joint works", {

  data(tab_test)

  n <- 75
  tab_test2 = tab_test[c(1:n,5001:(5000+n)),1:5]

  try = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
                   dist.choice = "M", which.DB = "B")

  print(try$GAMMA_A)
  ref = matrix(data=c(1.4133405, 0.4799928, 0.0000000,
                      0.0000000, 0.0000000, 0.6266879,
                      0.5199788, 0.0000000, 0.0000000,
                      0.0000000, 0.4758268, 1.0175065), nrow=4)

  expect_equal(all(mapply(all.equal, try$gamma_B, ref, tolerance=1e-6)), TRUE)
})
