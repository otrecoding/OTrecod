test_that("compare_lists works", {

  data1 = data.frame(Gender = rep(c("m", "f"), 5), Age = rnorm(5, 20, 4))
  data2 = data.frame(Gender = rep(c("m", "f"), 5), Age = rnorm(5, 21, 5))

  list1 = list(A = 1:4,
               B = as.factor(c("A", "B", "C")),
               C = matrix(1:6, ncol = 3))
  list2 = list(A = 1:4,
               B = as.factor(c("A", "B")),
               C = matrix(1:6, ncol = 3))
  list3 = list(A = 1:4,
               B = as.factor(c("A", "B", "C")),
               C = matrix(c(1:5, 7), ncol = 3))
  list4 = list(A = 1:4,
               B = as.factor(c("A", "B", "C")),
               C = matrix(1:6, ncol = 2))
  list5 = list(A = 1:4,
               B = as.factor(c("A", "B")),
               C = matrix(1:6, ncol = 2))
  list6 = list(A = 1:4, B = as.factor(c("A", "B")), C = data1)
  list7 = list(A = 1:4, B = as.factor(c("A", "B")), C = data2)

  expect_equal(compare_lists(list1, list2),
               c(FALSE, TRUE, FALSE))
  expect_equal(compare_lists(list1, list3),
               c(FALSE, FALSE, TRUE))
  expect_equal(compare_lists(list1, list4),
               c(FALSE, FALSE, TRUE))
  expect_equal(compare_lists(list1, list5),
               c(FALSE, TRUE, TRUE))
  expect_equal(compare_lists(list6, list7),
               c(FALSE, FALSE, TRUE))

})
