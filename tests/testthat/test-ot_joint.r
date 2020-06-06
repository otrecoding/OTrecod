test_that("ot_joint works", {

  data(tab_test)
  tab_test2 = tab_test[c(1:75,5001:5075),1:5]

  try1J = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
                   dist.choice = "M", which.DB = "B")

  expect_equal(2 * 2, 4)
})
