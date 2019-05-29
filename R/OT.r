
soluc   = merge_dbs(data3,data4,NAME_Y1 = "c.neti",NAME_Y2 = "c.neti.bis",ordinal_DB1 = c(2,4,6,9), ordinal_DB2 = c(1,3,9), impute = "NO")

prep    = transfo_dist(soluc[[1]],quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")
xx      = cost(prep,"G")
yy_DB1  = affectation(xx$solution,prep,dist_choice = "G",which.DB = "DB1")
yy_DB2  = affectation(xx$solution,prep,dist_choice = "G",which.DB = "DB2")
yy      = affectation(xx$solution,prep,dist_choice = "G",which.DB = "BOTH")


DB1_OT  = yy[[1]]
DB2_OT  = yy[[2]]

# DB_choice = "DB1"

tx       = 0.10

indicbd = c(2,1)
indicot = c(3,2)

OTDB  = list()
Ytheo = list()

for (k in 1:2){

  DB_OT  = yy[[k]]
  new_size = round(tx*nrow(DB_OT))

  set.seed(3026); samp = sample(1:nrow(DB_OT),new_size,replace = FALSE)

  DBOT_new = DB_OT[samp,]


  Ytheo[[k]] = as.character(DB_OT[samp,ncol(DB_OT)])

  DBOT_new[,k] = rep(k,nrow(DBOT_new))


  DBOT_new[,c(1:3,ncol(DBOT_new))] = apply(DBOT_new[,c(1:3,ncol(DBOT_new))],2,as.character)

  DBOT_new[,1] = rep(indicbd[k],nrow(DBOT_new))

  DBOT_new[,indicot[k]] = DBOT_new[,ncol(DBOT_new)]
  DBOT_new[,k+1] = rep(NA,nrow(DBOT_new))

  DBOT_new = DBOT_new[,-ncol(DBOT_new)]

  OTDB[[indicbd[k]]] = DBOT_new


}

OT_ready       = rbind(OTDB[[1]],OTDB[[2]])

for (j in 1:2){

#----> PBBBBBBB !

  if (is.ordered(DB1_OT[,j])){

    OT_ready[,j] = as.factor(OT_ready[,j])

  } else {

    OT_ready[,j] = as.character(OT_ready[,j])

  }

}

prep_new       = transfo_dist(OT_ready,quanti = c(1,4),nominal = 6,ordinal = c(2,3,5),prep_choice = "G")



for (j in 1:3){

  OT_ready[,j] = as.factor(OT_ready[,j])

}












