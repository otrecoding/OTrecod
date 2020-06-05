test_that("otjoin works", {

### Y and Z are a same variable encoded in 2 different forms in DB A and B:
### (3 levels for Y and 5 levels for Z)

data(simu_data)

### using a sample of the tab_test object (3 complete covariates)
### Y1 and Y2 are a same variable encoded in 2 different forms in DB 1 and 2:
### (4 levels for Y1 and 3 levels for Y2)

data(tab_test)
# Example with n1 = n2 = 75 and only X1 and X2 as covariates

tab_test2 = tab_test[c(1:75,5001:5075),1:5]

### An example of JOINT model (Manhattan distance)
# Suppose we want to impute the missing parts of Y1 in DB2 only ...
sol = OT_joint(tab_test2, nominal = c(1,4:5), ordinal = c(2,3),
                  prep_choice = "M", norm = 1, which.DB = 2)


ID = c("P_1", "P_2", "P_3", "P_4", "P_5", "P_6")
X1_2 = c(0, 0, 1, 1, 1, 0)
X2_2 = c(0, 1, 0, 0, 1, 0)
X2_3 = c(0, 0, 0, 1, 0, 1)
ref = data.frame(ID, X1_2, X2_2, X2_3)
expect_equal(all(ref == sol), TRUE)
})
