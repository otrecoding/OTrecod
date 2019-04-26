#-------------------------------------------------------------------------------
#
#   PROJET BIG DATA [Equipe 5]: FUSION DE BASES
#                               V. Gares / C. DiMeglio / G. Guernec / N. Savy
#-------------------------------------------------------------------------------


### Packages utiles au programme
library(gdata); library(FNN); library(RColorBrewer);
library(gplots); library(boot); library(linprog);

library(mvtnorm); library(mice)
library(nnet)
library(linprog)      #  ... Pour le Simplexe
library(mvtnorm)      #  ... Pour la simulation de vecteurs norm?s
library(mice)         #  ... Pour la m?thode d'imputation MICE
library(FactoMineR)   #  ... Pour l'Analyse Factorielle sur Donn?es Mixtes
library(nnet)
library(linprog);

### Quelques jeux de donn?es initiaux

# Solution evidente simulee (10022017):

#datadd =  read.table("D:\\statistiques\\4_Dossier partage\\Big data_VG_GG\\Program\\transport1502\\soluce_evidente.txt", sep="\t",header=TRUE)


# Augmentation de sa taille, conservation de ses distributions

# dataf = datadd
# for (i in 1:100){
#   
#   dataf = rbind(dataf,datadd)
#   
# }
# 
# dataf = dataf[sort.list(dataf$base),]
# dataf=dataf[,c(1,2,3,5,4)]
# 
# 
# # Inversion des colonnes Yb1 et Yb2
# 
# datae=datadd[,c(1,2,3,5,4)]
# colnames(datae)[4:5] = colnames(datadd)[4:5]
# 
# 
# # R?-arrangement des colonnes de "datae" (avec les covariables ? la fin)
# 
# datae_new = datae[,c(1,4,5,2,3)]

#---------

