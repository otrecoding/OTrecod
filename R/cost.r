
#' cost()
#'
#'    Function which realizes the first step (on 2) of the OT algorithm
#'
#' @param dat         A data.frame that overlays the two initial databases in a specific form (after using transfo_dist() function by example): The first column corresponds to the number of the databases
#'                    the second column corresponds to Y1, the values of the target Y encoded in the first database (unknown in the second one), and the third column corresponds to Y2
#'                    the values of the target Y encoded in the second database (unknown in the first one). The other columns corresponds to covariates ranked in an undifferentiated order.
#' @param dist_choice Choice of the distance mesure to use for the OT fusion. It depends on the type of the covariates.
#'                    The Gower's distance ("G") is used by default. Otherwise, this can be "E", or "H" for euclidean or Manhattan distance respectively.
#'
#' @return A list of 4 objects:
#'         Risk A vector corresponding to the risk rates of coupling modalities from Y1 with modalities from Y2. The risk is defined from the difference between the entropies of the covariable distributions.
#'         Cval A vector corresponding to the marginal distributions of Y1 and Y2
#'         A The matrix of constraints corresponding to the possible passages from Y1 to Y2
#'         solution The solution of the simplex corresponding to the joint distribution estimation of Y1 and Y2
#'
#'
#' @export
#'
#' @examples
#' soluc  = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")
#' prep   = transfo_dist(soluc[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")
#' xx     = cost(prep,"G"); xx$solution
#'
#' prep1  = transfo_dist(soluc[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "E")
#' x1     = cost(prep1,"E"); x1$solution
#'
#' soluc2 = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "CC")
#' prep2  = transfo_dist(soluc2[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "FAMD")
#' x2     = cost(prep2,"E"); x2$solution
#'
#'
#' dat = prep3
#' dist_choice = "E"
#'
#'
cost = function(dat, dist_choice = "G"){

  val1  = sort(unique(dat[dat[,1]==1,2]))
  val2  = sort(unique(dat[dat[,1]==2,3]))


  Risk=0
  V=0

  #  SUM OF INDIVIDUAL DISTANCES BETWEEN GROUPS

   for (i in 1:length(val1)){


    # print(paste(i,"/",length(val1)))
    cc = round(i*100/length(val1))
    progress(cc,max.value = 100)
    Sys.sleep(0.01)
    if (cc == 100) message("COST ESTIMATIONS OK")

    for (j in 1:length(val2)){

      V = 0

      sujetsB1 = subset(dat,dat[,1] == 1 & dat[,2] == val1[i])
      sujetsB2 = subset(dat,dat[,1] == 2 & dat[,3] == val2[j])
      datnew   = rbind(sujetsB1,sujetsB2)

      if (dist_choice == "G"){mat_gov  = gower.dist(datnew[,4:ncol(datnew)])} else {mat_gov = NULL}

      for (m in (1:nrow(sujetsB1))){

        for (n in (1:nrow(sujetsB2))){

          if (dist_choice == "G"){

            # GOWER DISTANCE

            V = sum(c(V,mat_gov[m,nrow(sujetsB1)+n]),na.rm = TRUE)

          } else {

            for (k in (4:ncol(dat))){

              if (dist_choice == "H"){

                # Hamming distance

                V = sum(c(V,as.numeric(sujetsB1[m,k]!=sujetsB2[n,k])),na.rm = TRUE)


              } else {

                # Euclidian distance

                V = sum(c(V,(sujetsB1[m,k]-sujetsB2[n,k])^2),na.rm = TRUE)

              }

            }

          }

        }
      }

      V=V/((ncol(dat)-3)*nrow(sujetsB1)*nrow(sujetsB2))
      Risk = c(Risk,V)

    }

  }

  Risk = Risk[-1]


  # PREPARING A,x and B FOR SIMPLEX

  maxi = c(length(val1),length(val2))


  S = matrix(ncol = maxi[1]*maxi[2],nrow = maxi[1]+maxi[2])


  seq1 = rep(1,maxi[2])
  seq0 = rep(0,maxi[2])
  seq2 = c(1,rep(0,maxi[2]-1))


  S[1,]         = L1 = c(seq1,rep(seq0,maxi[1]-1))
  S[maxi[1]+1,] = L2 = rep(seq2,maxi[1])


  for (i in setdiff(2:nrow(S),maxi[1]+1)){

    if (i <= maxi[1]){

      L1    = c(seq0,L1)
      S[i,] = L1[1:ncol(S)]
    }

    else {

      L2    = c(0,L2)
      S[i,] = L2[1:ncol(S)]

    }
  }


  CVal = as.numeric(c(table(dat[dat[,1]==1,2])/nrow(dat[dat[,1]==1,]),
                      table(dat[dat[,1]==2,3])/nrow(dat[dat[,1]==2,])))



  # SIMPLEX FUNCTION

  res = solveLP(cvec=Risk, bvec=CVal, Amat=S, const.dir = rep( "==",
                                                               length(CVal)),lpSolve=TRUE)


  x = res$solution

  return(list(Risk=Risk,CVal=CVal,A=S,solution=x))

}




