
#' IV.  FONCTIONS D'ESTIMATION DES ERREURS
#' Cette fonction retourne sur une ligne, les erreurs moyennes de pr?diction
#' des m?thodes FOT et MICE appliqu?es sur donn?es simul?es par approche
#' d?terministe ou non-d?terministe
#'
#' @param  M          : Nombre de r?p?titions
#' @param  n_1        : Taille de la base 1
#' @param  coeff      : Coefficient de proportionnalit? taille base 2
#' @param  p1cc       : Proportion des cat?gories pour la covariable 1
#' @param  p2cc       : Proportion des cat?gories pour la covariable 1
#' @param  qual       : Les n? de colonnes des variables qualitatives
#' @param  quan       : Les n? de colonnes des variables quantitatives
#' @param  choose_dist: Choix de la distance pour la m?thode FOT ("E" ou "H")
#' @param  rep_MICE   : Nbre de BDD ? imputer pour MICE
#' @param  cov_MICE   : Sp?cifier si Y1, Y2, X1, X2 et X3 doivent ?tre consid?r?es dans
#'                      MICE au format num?rique ("pmm") ou facteur ("polyreg")
#' @param  calcu      : "FREQ" ou "MOY" selon que les variables imput?es soient quantis
#'                       ou qualis
#' @param  typ_simu   : D?terministe "D" ou Non D?terministe "ND"
#'
#' @param  Les coefficients suivants ne sont utilis?s que dans le cas non-d?terministe:
#' @param  r2t        : R^2 de la r?gression lin?aire Ycontinue = f(X1,X2,X3)
#' @param  py1cc      : Proportion de patients dans chaque classe de Y en base 1
#' @param  py2cc      : Proportion de patients dans chaque classe de Y en base 2
#' @param  vlY1       : Les labels de chaque classe pour Y cat?gorielle en base 1
#' @param  vlY2       : Les labels de chaque classe pour Y cat?gorielle en base 2
#' @param  bdd        : Jeu de donn?es non simul?
#' @param  seedsim    : G?n?rateur al?atoire (tir? au sort par d?faut)
#'
#' @return
#' @export
#'
#' @examples
#' n.sim = 1
#' n_1         = 1000
#' coeff       = 1
#' p1cc        = c(0.5,0.5)
#' p2cc        = c(0.3,0.3,0.4)
#' qual        = 4:5
#' quan        = 6
#' choose_dist = "E"
#' rho_cov     = 0.2
#' rep_MICE    = 5
#' cov_MICE    = c("polyreg","polyreg","polyreg","polyreg","pmm")
#' calcu       = "FREQ"
#' typ_simu    = "ND"
#' r2t         = 0.5
#' py1cc       = c(0.25,0.25,0.25,0.25)
#' py2cc       = c(0.333,0.334,0.333)
#' vlY1        = c(1,2,3,4)
#' vlY2        = c(1,2,3)
#' n1,k,cormat,px1c,px2c,r2,py1c,py2c,valY1,valY2
#' simul_glob (n_1, coeff, p1cc, p2cc, qual , quan ,
#'             choose_dist="E", rho_cov, rep_MICE=5, cov_MICE = "pmm", calcu="FREQ",
#'             typ_simu,r2t , py1cc, py2cc,
#'             vlY1, vlY2)
simul_glob = function(n_1=NULL, coeff=NULL, p1cc=NULL, p2cc=NULL, qual = NULL, quan = NULL,
                      choose_dist="E", rho_cov=NULL, rep_MICE=5, cov_MICE = "pmm", calcu="FREQ",
                      typ_simu=NULL,r2t = NULL, py1cc = NULL, py2cc = NULL,
                      vlY1 = NULL, vlY2 = NULL){
  # Stockage dans un vecteur des taux d'erreur de pr?diction et des temps
  # d'ex?cution
  cor_cov = c(corX1_X2 = rho_cov,corX1_X3 = rho_cov,
              corX2_X3 = rho_cov)
  peraffec_FOT  = time_FOT  = per_affec_COMP1 = numeric(1)
  peraffec_MICE = time_MICE = per_affec_COMP2 = numeric(1)
  peraffec_PR = time_PR = per_affec_COMP3 = numeric(1)
  per_affec_COMP3 = per_affec_COMP4 = per_affec_COMP5= per_affec_COMP6 = numeric(1)
  pbFOT = pbMICE = 0




  # I] PREPARATION DE LA BDD



  #set.seed(seedsim + i);
  datgen = gendata_ND(n_1,coeff,cor_cov,p1cc,p2cc,r2t,py1cc,py2cc,vlY1,vlY2)



  # Transformation ?ventuelle des outcomes

  if (class(datgen$Y1) %in% c("factor","character")){

    datgen$Y1         = as.factor(as.character(datgen$Y1))
    levels(datgen$Y1) = 1:length(levels(datgen$Y1))
    datgen$Y1         = as.numeric(datgen$Y1)

  }

  if (class(datgen$Y2) %in% c("factor","character")){

    datgen$Y2         = as.factor(as.character(datgen$Y2))
    levels(datgen$Y2) = 1:length(levels(datgen$Y2))
    datgen$Y2         = as.numeric(datgen$Y2)

  }


  # On sp?cifie les Y manquants en base 1 et 2:

  datgen2 = datgen
  datgen2[datgen[,1]==1,"Y2"] = NA
  datgen2[datgen[,1]==2,"Y1"] = NA



  # Transformation des covariables

  datgen3 = prepar_cov(datgen2, quali = qual, quanti = quan)




  # II] Imputation des Y selon les m?thodes FOT et MICE


  #--------------#
  # METHODE FOT  #
  #--------------#


  t1 = Sys.time()

  # datnew = affectation(cout(datgen3,choose_dist)$solution,datgen3)
  datnew = try(affectation(cout(datgen3,choose_dist)$solution,datgen3),TRUE)

  t2 = Sys.time()


  # R?cup?ration du tps d'ex?cution de FOT:

  time_FOT[1] = as.numeric(round(difftime(t2,t1,units="secs"),2))


  # R?cup?ration de l'erreur de pr?diction FOT:

  if(is.character(datnew)==FALSE){

    datnew$old = c(datgen[1:nrow(datgen[datgen[,1]==1,]),"Y2"],
                   datgen[(nrow(datgen[datgen[,1]==1,])+1):nrow(datgen),
                          "Y1"])
    peraffec_FOT[1] = mean(as.numeric(datnew$FOT != datnew$old))

  } else {

    datnew = datgen3
    datnew$old = c(datgen[1:nrow(datgen[datgen[,1]==1,]),"Y2"],
                   datgen[(nrow(datgen[datgen[,1]==1,])+1):nrow(datgen),
                          "Y1"])
    datnew$FOT = NA
    peraffec_FOT[1] = 1
    pbFOT = pbFOT + 1

  }

  #---------------#
  # METHODE MICE  #
  #---------------#

  # Conversion de variables en facteur avant d'appliquer MICE

  if ("polyreg" %in% cov_MICE){

    indic_fac = (2:(length(cov_MICE) + 1))[cov_MICE != "pmm"]

    for (j in 1:length(indic_fac)){

      datgen2[,indic_fac[j]] = as.factor(datgen2[,indic_fac[j]])

    }

  }


  # Ex?cution de MICE

  t1 = Sys.time()
  #set.seed(seedsim + i)
  datnew$MICE = try(as.numeric(imp_MICE(datgen2,rep_MICE,cov_MICE,
                                        calcu)),TRUE)
  t2 = Sys.time()


  if ((is.numeric(datnew$MICE)==TRUE)){

    # R?cup?ration du tps d'ex?cution de MICE:

    time_MICE[1] = as.numeric(round(difftime(t2,t1,units="secs"),2))


    # R?cup?ration de l'erreur de pr?diction MICE:

    peraffec_MICE[1] = 1-(sum(as.numeric(datnew$MICE == datnew$old),
                              na.rm=TRUE))/(nrow(datnew))


  } else {

    peraffec_MICE[1] = 1
    pbMICE = pbMICE + 1
    datnew$MICE = NA

  }


  if ((sum(is.na(datnew$MICE))!=nrow(datnew))&(sum(is.na(datnew$FOT))!=nrow(datnew))){

    per_affec_COMP1[1] = mean(as.numeric(datnew$MICE == datnew$FOT),na.rm=TRUE)
    count1 = ifelse(as.numeric(datnew$FOT  == datnew$old),1,0)
    count2 = ifelse(as.numeric(datnew$MICE == datnew$old),1,0)
    # per_affec_COMP2[1] = mean(count1 == count2)
    per_affec_COMP2[1] = mean(ifelse((count1==1)&(count2==1),1,0))

  } else {

    per_affec_COMP1[1] = per_affec_COMP2[1] = NA

  }



  # METHODE PR:
  #-------------
  t1 = Sys.time()
  data1= subset(datgen,ident==1)
  data2= subset(datgen,ident==2)
  reg1 = nnet::multinom(Y1  ~ X1 + X2 + X3, data = data1)
  data2$predY1 = predict(reg1, newdata = data2[,c(4,5,6)] )
  reg2 = nnet::multinom(Y2  ~ X1 + X2 + X3, data = data2)
  data1$predY2 = predict(reg2, newdata = data1[,c(4,5,6)] )
  datnew$PR = c(data1$predY2,data2$predY1)
  peraffec_PR[1] = mean(as.numeric(datnew$PR != datnew$old))

  if ((sum(is.na(datnew$PR))!=nrow(datnew))&(sum(is.na(datnew$FOT))!=nrow(datnew))){

    per_affec_COMP3[1] = mean(as.numeric(datnew$PR == datnew$FOT),na.rm=TRUE)
    count1 = ifelse(as.numeric(datnew$FOT  == datnew$old),1,0)
    count2 = ifelse(as.numeric(datnew$PR == datnew$old),1,0)
    # per_affec_COMP2[i] = mean(count1 == count2)
    per_affec_COMP4[1] = mean(ifelse((count1==1)&(count2==1),1,0))

  } else {

    per_affec_COMP3[1] = per_affec_COMP4[1] = NA

  }

  if ((sum(is.na(datnew$MICE))!=nrow(datnew))&(sum(is.na(datnew$PR))!=nrow(datnew))){

    per_affec_COMP5[1] = mean(as.numeric(datnew$MICE == datnew$PR),na.rm=TRUE)
    count1 = ifelse(as.numeric(datnew$FOT  == datnew$old),1,0)
    count2 = ifelse(as.numeric(datnew$PR == datnew$old),1,0)
    # per_affec_COMP2[i] = mean(count1 == count2)
    per_affec_COMP6[1] = mean(ifelse((count1==1)&(count2==1),1,0))

  } else {

    per_affec_COMP5[1] = per_affec_COMP6[1] = NA

  }



  results=c(peraffec_FOT,peraffec_MICE,peraffec_PR,per_affec_COMP1,per_affec_COMP2,per_affec_COMP3,
            per_affec_COMP4,per_affec_COMP5,per_affec_COMP6)
  return(results)

  # if (is.data.frame(bdd) == FALSE){
  #
  #   return(round(c(n1 = n_1, K = coeff, prob1c = p1cc[2] , prob2c = p2cc[2],
  #                  cor_cov, R2 = ifelse(typ_simu=="D",NA,r2t),
  #
  #                  meanError_FOT   = mean(peraffec_FOT),
  #                  sdError_FOT     = sd(peraffec_FOT),
  #                  meanTime_FOT    = mean(time_FOT),sdTime_FOT = sd(time_FOT),
  #
  #                  meanError_MICE  = mean(peraffec_MICE),
  #                  sdError_MICE    = sd(peraffec_MICE),
  #                  meanTime_MICE   = mean(time_MICE),
  #                  sdTime_MICE     = sd(time_MICE),
  #                  seed            = seedsim,
  #                  pb_FOT          = pbFOT,
  #                  pb_MICE         = pbMICE,
  #
  #                  meanOK_COMP1 = ifelse(sum(is.na(per_affec_COMP1))!=length(per_affec_COMP1),mean(per_affec_COMP1,na.rm=TRUE),NA),
  #                  sdOK_COMP1   = ifelse(sum(is.na(per_affec_COMP1))!=length(per_affec_COMP1),sd(per_affec_COMP1,na.rm=TRUE),NA),
  #                  meanOK_COMP2 = ifelse(sum(is.na(per_affec_COMP2))!=length(per_affec_COMP2),mean(per_affec_COMP2,na.rm=TRUE),NA),
  #                  sdOK_COMP2   = ifelse(sum(is.na(per_affec_COMP2))!=length(per_affec_COMP2),sd(per_affec_COMP2,na.rm=TRUE),NA)),6))
  #
  # } else {
  #
  #   return(round(c(n1 = nrow(datnew[datnew$base==1,]), K = nrow(datnew[datnew$base==2,])/nrow(datnew[datnew$base==1,]),
  #
  #                  Error_FOT       = peraffec_FOT,
  #                  Time_FOT        = time_FOT,
  #
  #                  Error_MICE      = peraffec_MICE,
  #                  Time_MICE       = time_MICE,
  #                  seed            = seedsim,
  #                  pb_FOT          = pbFOT,
  #                  pb_MICE         = pbMICE,
  #
  #                  OK_COMP1        = ifelse(sum(is.na(per_affec_COMP1))!=length(per_affec_COMP1),mean(per_affec_COMP1,na.rm=TRUE),NA),
  #                  OK_COMP2        = ifelse(sum(is.na(per_affec_COMP2))!=length(per_affec_COMP2),mean(per_affec_COMP2,na.rm=TRUE),NA)),6))
  #
  #
  #
  # }

}
