% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_dbs.R
\name{merge_dbs}
\alias{merge_dbs}
\title{merge_dbs()}
\usage{
merge_dbs(
  DB1,
  DB2,
  row_ID1 = NULL,
  row_ID2 = NULL,
  NAME_Y,
  NAME_Z,
  order_levels_Y = levels(DB1[, NAME_Y]),
  order_levels_Z = levels(DB2[, NAME_Z]),
  ordinal_DB1 = NULL,
  ordinal_DB2 = NULL,
  impute = "NO",
  R_MICE = 5,
  NCP_FAMD = 3,
  seed_choice = sample(1:1e+06, 1)
)
}
\arguments{
\item{DB1}{a data.frame corresponding to the 1st database to merge (top database)}

\item{DB2}{a data.frame corresponding to the 2nd database to merge (bottom database)}

\item{row_ID1}{the column index of the row identifier of DB1 if it exists (no identifier by default)}

\item{row_ID2}{the column index of the row identifier of DB2 if it exists (no identifier by default)}

\item{NAME_Y}{the name of the outcome (with quotes) in its specific scale/encoding from the 1st database (DB1)}

\item{NAME_Z}{the name of the outcome (with quotes) in its specific scale/encoding from the 2nd database (DB2)}

\item{order_levels_Y}{the levels of \eqn{Y} stored in a vector and sorted in ascending order in the case of ordered factors. This option permits to reorder the levels in the 1st database (DB1) if necessary.}

\item{order_levels_Z}{the levels of \eqn{Z} stored in a vector and sorted in ascending order in the case of ordered factors. This option permits to reorder the levels in the 2nd database (DB2) if necessary.}

\item{ordinal_DB1}{a vector of column indexes corresponding to ordinal variables in the 1st database (no ordinal variable by default)}

\item{ordinal_DB2}{a vector of column indexes corresponding to ordinal variables in the 2nd database (no ordinal variable by default)}

\item{impute}{a character equals to "NO" when missing data on covariates are kept (Default option), "CC" for Complete Case by keeping only covariates with no missing information , "MICE" for MICE multiple imputation approach, "FAMD" for single imputation approach using Factorial Analysis for Mixed Data}

\item{R_MICE}{the chosen number of multiple imputations required for the  MICE approach (5 by default)}

\item{NCP_FAMD}{an integer corresponding to the number of components used to predict missing values in FAMD imputation (3 by default)}

\item{seed_choice}{an integer used as argument by the set.seed() for offsetting the random number generator (Random integer by default, only useful with MICE)}
}
\value{
A list containing 12 elements (13 when \code{impute} equals "MICE"):
\item{DB_READY}{the database matched from the two initial databases with common covariates and imputed or not according to the impute option}
\item{ID1_drop}{the row numbers or row identifiers excluded of the data merging because of the presence of missing values in the target variable of DB1. NULL otherwise}
\item{ID2_drop}{the row numbers or row identifiers excluded of the data merging because of the presence of missing values in the target variable of DB2. NULL otherwise}
\item{Y_LEVELS}{the remaining levels of the target variable \eqn{Y} in the DB1}
\item{Z_LEVELS}{the remaining Levels of the target variable \eqn{Z} in the DB2}
\item{REMOVE1}{the labels of the deleted covariates because of type incompatibilies of type from DB1 to DB2}
\item{REMOVE2}{the removed factor(s) because of levels incompatibilities from DB1 to DB2}
\item{REMAINING_VAR}{labels of the remained covariates for data fusion}
\item{IMPUTE_TYPE}{a character with quotes that specify the method eventually chosen to handle missing data in covariates}
\item{MICE_DETAILS}{a list containing the details of the imputed datasets using \code{MICE} when this option is chosen. Raw and imputed databases imputed for DB1 and DB2 according to the number of multiple imputation selected (Only if impute = "MICE")}
\item{DB1_raw}{a data.frame corresponding to DB1 after merging}
\item{DB2_raw}{a data.frame corresponding to DB2 after merging}
\item{SEED}{an integer used as argument by the \code{set.seed} function for offsetting the random number generator (random selection by default)}
}
\description{
Harmonization and merging before data fusion of two databases with specific outcome variables and shared covariates.
}
\details{
Assuming that DB1 and DB2 are two databases (two separate data.frames with no overlapping rows) to be merged vertically before data fusion, the function \code{merge_dbs} performs this merging and checks the harmonization of the shared variables.
Firslty, the two databases declared as input to the function (via the argument \code{DB1} and \code{DB2}) must have the same specific structure.
Each database must contain a target variable (whose label must be filled in the argument \code{Y} for DB1 and in \code{Z} for DB2 respectively, so that the final synthetic database in output will contain an incomplete variable \code{Y} whose corresponding values will be missing in DB2 and another incomplete target \code{Z} whose values will be missing in DB1), a subset of shared covariates (by example, the best predictors of \eqn{Y} in DB1, and \eqn{Z} in DB2).
Each database can have a row identifier whose label must be assigned in the argument \code{row_ID1} for DB1 and \code{row_ID2} for DB2. Nevertheless, by default DB1 and DB2 are supposed with no row identifiers. The merging keeps unchanged the order of rows in the two databases provided that \eqn{Y} and \eqn{Z} have no missing values.
By building, the first declared database (in the argument \code{DB1}) will be placed automatically above the second one (declared in the argument \code{DB2}) in the final database.

Firstly, by default, a variable with the same name in the two databases is abusively considered as shared. This condition is obviously insufficient to be kept in the final subset of shared variables,
and the function \code{merge_dbs} so performs checks before merging described below.

A. Discrepancies between shared variables
\itemize{
\item Shared variables with discrepancies of types between the two databases (for example, a variable with a common name in the two databases but stored as numeric in DB1, and stored as character in DB2) will be removed from the merging and the variable name will be saved in output (\code{REMOVE1}).
\item Shared factors with discrepancies of levels (or number of levels) will be also removed from the merging and the variable name will be saved in output (\code{REMOVE2}).
\item covariates whose names are specific to each database will be also deleted from the merging.
\item If some important predictors have been improperly excluded from the merging due to the above-mentioned checks, it is possible for user to transform these variables a posteriori, and re-run the function.
}

B. Rules for the two outcomes (target variables)

The types of \code{Y} and \code{Z} must be suitable:
\itemize{
\item Categorical (ordered or not) factors are allowed.
\item Numeric and discrete outcomes with a finite number of values are allowed but will be automatically converted as ordered factors using the function \code{transfo_target} integrated in the function \code{merge_dbs}.
}

C. The function \code{merge_dbs} handles incomplete information of shared variables, by respecting the following rules:
\itemize{
\item If \code{Y} or \code{Z} have missing values in DB1 or DB2, corresponding rows are excluded from the database before merging. Moreover, in the case of incomplete outcomes,
if A and B have row identifiers, the corresponding identifiers are removed and these latters are stored in the objects \code{DB1_ID} and \code{DB2_ID} of the output.
\item Before overlay, the function deals with incomplete covariates according to the argument \code{impute}.
Users can decide to work with complete case only ("CC"), to keep ("NO") or impute incomplete information ("MICE","FAMD").
\item The function \code{imput_cov}, integrated in the syntax of \code{merge_dbs} deals with imputations. Two approaches are actually available:
the multivariate imputation by chained equation approach (MICE, see (3) for more details about the approach or the corresponding package \pkg{mice}),
and an imputation approach from the package \pkg{missMDA} that uses a dimensionality reduction method (here a factor analysis for mixed data called FAMD (4)), to provide single imputations.
If multiple imputation is required (\code{impute} = "MICE"), the default imputation methods are applied according to the type of the variables. The average of the plausible values will be kept for a continuous variable, while the most frequent candidate will be kept as a consensus value for a categorical variable or factor (ordinal or not).
}

As a finally step, the function checks that all values related to \eqn{Y} in B are missing and inversely for \eqn{Z} in A.
}
\examples{

### Assuming two distinct databases from simu_data: data_A and data_B
### Some transformations will be made beforehand on variables to generate
### heterogeneities between the two bases.
data(simu_data)
data_A <- simu_data[simu_data$DB == "A", c(2, 4:8)]
data_B <- simu_data[simu_data$DB == "B", c(3, 4:8)]

# For the example, a covariate is added (Weight) only in data_A
data_A$Weight <- rnorm(300, 70, 5)

# Be careful: the target variables must be in factor (or ordered) in the 2 databases
# Because it is not the case for Yb2 in data_B, the function will convert it.
data_B$Yb2 <- as.factor(data_B$Yb2)

# Moreover, the Dosage covariate is stored in 3 classes in data_B (instead of 4 classes in data_B)
# to make the encoding of this covariate specific to each database.
data_B$Dosage <- as.character(data_B$Dosage)
data_B$Dosage <- as.factor(ifelse(data_B$Dosage \%in\% c("Dos 1", "Dos 2"), "D1",
  ifelse(data_B$Dosage == "Dos 3", "D3", "D4")
))

# For more diversity, this covariate iis placed at the last column of the data_B
data_B <- data_B[, c(1:3, 5, 6, 4)]

# Ex 1: The two databases are merged and incomplete covariates are imputed using MICE
merged_ex1 <- merge_dbs(data_A, data_B,
  NAME_Y = "Yb1", NAME_Z = "Yb2",
  ordinal_DB1 = c(1, 4), ordinal_DB2 = c(1, 6),
  impute = "MICE", R_MICE = 2, seed_choice = 3011)

summary(merged_ex1$DB_READY)


# Ex 2: The two databases are merged and missing values are kept
merged_ex2 <- merge_dbs(data_A, data_B,
  NAME_Y = "Yb1", NAME_Z = "Yb2",
  ordinal_DB1 = c(1, 4), ordinal_DB2 = c(1, 6),
  impute = "NO", seed_choice = 3011
)

# Ex 3: The two databases are merged by only keeping the complete cases
merged_ex3 <- merge_dbs(data_A, data_B,
  NAME_Y = "Yb1", NAME_Z = "Yb2",
  ordinal_DB1 = c(1, 4), ordinal_DB2 = c(1, 6),
  impute = "CC", seed_choice = 3011
)

# Ex 4: The two databases are merged and incomplete covariates are imputed using FAMD
merged_ex4 <- merge_dbs(data_A, data_B,
  NAME_Y = "Yb1", NAME_Z = "Yb2",
  ordinal_DB1 = c(1, 4), ordinal_DB2 = c(1, 6),
  impute = "FAMD", NCP_FAMD = 4, seed_choice = 2096
)

# Conclusion:
# The data fusion is successful in each situation.
# The Dosage and Weight covariates have been normally excluded from the fusion.
# The covariates have been imputed when required.

}
\references{
\enumerate{
\item Gares V, Dimeglio C, Guernec G, Fantin F, Lepage B, Korosok MR, savy N (2019). On the use of optimal transportation theory to recode variables and application to database merging. The International Journal of Biostatistics.
Volume 16, Issue 1, 20180106, eISSN 1557-4679. doi:10.1515/ijb-2018-0106
\item Gares V, Omer J (2020) Regularized optimal transport of covariates and outcomes in data recoding. Journal of the American Statistical Association. \doi{10.1080/01621459.2020.1775615}
\item van Buuren, S., Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1–67. url{https://www.jstatsoft.org/v45/i03/}
\item Josse J, Husson F (2016). missMDA: A Package for Handling Missing Values in Multivariate Data Analysis. Journal of Statistical Software, 70(1), 1–31. \doi{10.18637/jss.v070.i01}
}
}
\seealso{
\code{\link{imput_cov}}, \code{\link{transfo_target}}, \code{\link{select_pred}}
}
\author{
Gregory Guernec

\email{otrecod.pkg@gmail.com}
}
