
#------------------ Fonction AFFECTATION ---------------------------------------
simplx      = x2$solution
dat         = prep2
dist_choice = "E"
which.DB    = "DB2"

simplx      = x2$solution
n2 = 772
val1 = levels(dat$Y1)
val2 = levels(dat$Y2)
mat_simplx   = matrix(round(simplx*n2,4),nrow=length(val2),ncol=length(val1))
mat_simplx2  = matrix(round(simplx*n1,4),nrow=length(val2),ncol=length(val1))

table(dat$Y2)
table(dat$Y1)





#' affectation()
#'
#' Function which realizes the second (and last) step of the OT algorithm
#'
#' @param simplx A vector of double. It corresponds to the solution of the simplex from the cost function
#' @param dat A data.frame which first column corresponds to the number of the database (1 or 2), the second column to the target of Y1 and the third column to the target of Y2
#'            The following columns corresponds to common covariates at DB1 and DB2 in an unspecific order.
#' @param dist_choice A character with quotes. The distance function selected among the Gower distance ("G", by default), the Manhattan distance ("M"), or the Euclidean distance ("E")
#' @param which.DB A character string. If "BOTH" (Default), the OT algorithm is used on the targets variable of the 2 stacked DBs. If "DB1" the imputation is only done on the target of the 1st DB.
#'                 If "DB2" the imputation is only done on the target of the 2nd DB
#
#' @return A list of 2 data.frames:
#'         DB1_NEW A data.frame with the imputed target (OT column) in DB1 with the encoding from DB2 (NULL if which.DB = "DB2")
#'         DB2_NEW A data.frame with the imputed target (OT column) in DB2 with the encoding from DB1 (NULL if which.DB = "DB1")
#'
#' @export
#'
#' @examples
#'
#' soluc   = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
#' prep    = transfo_dist(soluc[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")
#' xx      = cost(prep,"G")
#' yy_DB1  = affectation(xx$solution,prep,dist_choice = "G",which.DB = "DB1")
#' yy_DB2  = affectation(xx$solution,prep,dist_choice = "G",which.DB = "DB2")
#' yy      = affectation(xx$solution,prep,dist_choice = "G",which.DB = "BOTH")
#'
#' soluc2  = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "CC")
#' prep2   = transfo_dist(soluc2[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "FAMD")
#' x2      = cost(prep2,"E")
#' yy2     = affectation(x2$solution,prep2,dist_choice = "E",which.DB = "BOTH")
#' yy2_DB1 = affectation(x2$solution,prep2,dist_choice = "E",which.DB = "DB1")
#' yy2_DB2 = affectation(x2$solution,prep2,dist_choice = "E",which.DB = "DB2")
#'
#'
affectation = function(simplx,dat,dist_choice = "E",which.DB = "BOTH"){



  if (!(which.DB %in% c("BOTH","DB1","DB2"))){

    stop("Bad specification for which.DB option: Please choose between BOTH,DB1, or DB2")

  } else {}


  # if (setequal(names(table(dat[,1])),1:2)==FALSE){

  #  stop("Encoding problem with the 1st column specifying the DB. Values 1 and 2 only ...")

  # } else {}

  dat[,1] = as.numeric(as.factor(dat[,1]))

    if (length(names(table(dat[,1]))) != 2){

      stop("Your 1st column counts an inappropriate number of DBs (only 2 DBs required)")

    } else {}


  if (all(is.na(dat[dat[,1]==1,3])) == FALSE){

    stop("All encodings from DB2 are not null in DB1. Problem of NA in col3")

  } else {}

  if (all(is.na(dat[dat[,1]==2,2])) == FALSE){

    stop("All encodings from DB1 are not null in DB2. Problem of NA in col2")

  } else {}


  n1      = nrow(dat[dat[,1]==1,])               # Nbr of subjects in DB1
  n2      = nrow(dat[dat[,1]==2,])               # Nbr of subjects in DB2

  # Add a topic identification column
  dat$suj = 1:nrow(dat)

  # Add an imputation column
  dat$OT  = NA

  val1    = sort(unique(dat[dat[,1]==1,2]))   # Levels of Y in DB1
  val2    = sort(unique(dat[dat[,1]==2,3]))   # Levels of Y in DB2


  if (which.DB %in% c("BOTH","DB1")){



    # I. CONVERSION OF Y[DB1] WITH ENCODING FROM DB2
    #-----------------------------------------------

    # 1) The indices of the maximum from the COST matrix are retrieved

    mat_simplx = matrix(round(simplx*n1,4),nrow=length(val2),ncol=length(val1))

    # Matrix of integers only
    mat_simplx = t(simpond(t(mat_simplx),dat[,2]))


    print(paste("Number of remaining NA to encode in DB1:",sum(mat_simplx)))

    ind = which(mat_simplx==max(mat_simplx), arr.ind = T)

    # Handling the ex-aecquos: The index of the first met is conserved
    if (is.matrix(ind)== TRUE){ind=ind[1,]}


    # 2) We create a new base from which we will progressively remove the subjects for which a value of Y has been imputed in DB1

    dat2 = dat


    # As long as the transition matrix is not empty, the following conditions are repeated:

    while (length(dat$OT[(dat[,1] == 1) & (is.na(dat$OT)==TRUE)])!=0){


      # We extract from DB1 and DB2, the subgroups of subjects which corresponds to the indices of the maximum retained

      sujetsB1 = subset(dat2,dat2[,1] == 1 & dat2[,2] == val1[ind[2]])
      sujetsB2 = subset(dat2,dat2[,1] == 2 & dat2[,3] == val2[ind[1]])
      datnew   = rbind(sujetsB1,sujetsB2)

      if (dist_choice == "G"){

        mat_gov  = gower.dist(datnew[,4:ncol(datnew)])

      } else {

        mat_gov = NULL

      }


      # SUM OF DISTANCES BETWEEN A SUBJECT BELONGING TO A GROUP IN DB1 AND ALL OTHER SUBJECTS OF A GROUP IN DB2
      #--------------

      dista = numeric(nrow(sujetsB1))

      for (m in (1:nrow(sujetsB1))){

        V =0

        for (n in (1:nrow(sujetsB2))){


          if (dist_choice == "G"){

            # GOWER DISTANCE
            V = sum(c(V,mat_gov[m,nrow(sujetsB1)+n]),na.rm = TRUE)

          } else {

            for (k in (4:(ncol(dat2)-2))){

              if (dist_choice == "H"){

                # MANHATTAN DISTANCE
                V = sum(c(V,as.numeric(sujetsB1[m,k]!=sujetsB2[n,k])),na.rm = TRUE)

              } else {

                # EUCLIDEAN DISTANCE
                # V = V + (sujetsB1[m,k]-sujetsB2[n,k])^2
                V = sum(c(V,(sujetsB1[m,k]-sujetsB2[n,k])^2),na.rm = TRUE)

              }

            }
         }
        }

        dista[m] = V / (nrow(sujetsB2)*(ncol(dat)-3))

      }


      indo = sujetsB1$suj[order(dista)]
       suj = numeric(0)


      # Patients belonging to the subgroup extracted from base 1 are classified by increasing mean distance from each of these patients
      # to other subjects in the group extracted from base 2.
    #---------------------

      for (k in (1:mat_simplx[ind[1],ind[2]])){

        dat[indo[k],ncol(dat)] = val2[ind[1]]
        suj = c(suj,indo[k])

      }


      # we remove the subjects of DB1 for which Y has just been recoded in DB2:

      dat2=dat2[!(dat2$suj %in% suj),]


      # Une fois le max. de la matrice de passage exploité, on le remplace
      # par un 0:


      mat_simplx[ind[1],ind[2]] = 0
      print(paste("Number of remaining NA to encode in DB1:",sum(mat_simplx)))


      # We are looking for a new max. and we select the first met in case of ex-aecquos

      ind = which(mat_simplx==max(mat_simplx), arr.ind = T)

      if (is.matrix(ind)== TRUE){

        ind=ind[1,]

      } else {}

    }

    DB1_NEW     = dat[dat[,1] == 1,]
    DB1_NEW$OT  = as.factor(DB1_NEW$OT)
    DB1_NEW$OT  = mapvalues(DB1_NEW$OT,from = levels(DB1_NEW$OT), to = levels(DB1_NEW$Y2))


  } else {

    DB1_NEW = NULL
    dat$OT[dat[,1]==1] = rep(9,n1)

  }


  if (which.DB %in% c("BOTH","DB2")){



    # II. CONVERSION OF Y[DB2] WITH ENCODING FROM DB1
    #------------------------------------------------

    # Same reasoning as before ...

    mat_simplx = matrix(round(simplx*n2,4),nrow=length(val2),ncol=length(val1))
    mat_simplx = simpond(mat_simplx,dat[,3])

    print(paste("Number of remaining NA to encode in DB2:",sum(mat_simplx)))

    ind = which(mat_simplx==max(mat_simplx), arr.ind = T)


    if (is.matrix(ind)== TRUE){ind=ind[1,]}

      dat2 = dat


      while (length(dat$OT[is.na(dat$OT)==TRUE])!=0){

        sujetsB1 = subset(dat2,dat2[,1] == 1 & dat2[,2] == val1[ind[2]])
        sujetsB2 = subset(dat2,dat2[,1] == 2 & dat2[,3] == val2[ind[1]] )
        datnew   = rbind(sujetsB2,sujetsB1)

        if (dist_choice == "G"){

          mat_gov  = gower.dist(datnew[,4:ncol(datnew)])

        } else {

          mat_gov = NULL

        }

        distb = numeric(nrow(sujetsB2))

        for (m in (1:nrow(sujetsB2))){

          V =0

          for (n in (1:nrow(sujetsB1))){


            if (dist_choice == "G"){

          # GOWER DISTANCE
          V = sum(c(V,mat_gov[m,nrow(sujetsB2)+n]),na.rm = TRUE)

            } else {

              for (k in (4:(ncol(dat2)-2))){

                if (dist_choice == "H"){

                  # MANHATTAN DISTANCE
                  V = sum(c(V,as.numeric(sujetsB1[n,k]!=sujetsB2[m,k])),na.rm = TRUE)

                } else {

                  # EUCLIDEAN DISTANCE
                  V = sum(c(V,(sujetsB1[n,k]-sujetsB2[m,k])^2),na.rm = TRUE)


                }

              }
            }
          }

          distb[m] = V / (nrow(sujetsB1)*(ncol(dat)-3))

        }

        indo = sujetsB2$suj[order(distb)]

        suj = numeric(0)


        for (k in (1:mat_simplx[ind[1],ind[2]])){

          dat[indo[k],ncol(dat)] = val1[ind[2]]
          suj=c(suj,indo[k])

        }

        dat2=dat2[!(dat2$suj %in% suj),]


        mat_simplx[ind[1],ind[2]] = 0
        print(paste("Number of remaining NA to encode in DB2:",sum(mat_simplx)))

        ind = which(mat_simplx==max(mat_simplx), arr.ind = T)

        if (is.matrix(ind)== TRUE){

          ind=ind[1,]

        } else {}


    }

    DB2_NEW     = dat[dat[,1] == 2,]
    DB2_NEW$OT  = as.factor(DB2_NEW$OT)
    DB2_NEW$OT  = mapvalues(DB2_NEW$OT,from = levels(DB2_NEW$OT), to = levels(DB2_NEW$Y1))

  } else {

    DB2_NEW = NULL

  }

  return(list(DB1_NEW = DB1_NEW[,-(ncol(DB1_NEW)-1)],DB2_NEW = DB2_NEW[,-(ncol(DB2_NEW)-1)]))

}

