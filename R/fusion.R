
#     III.  FONCTIONS DE FUSION DE BASES

#' Fonction COUT
#'
#' Cette fonction est la partie I/II de la méthode FOT
#' Elle estime les co?ts de passage entre les modalit?s Y1 et Y2
#'
#' @param dat BDD de la forme [base|Y1|Y2|???], ??? ?tant des covariables
#'            binaires ou quantis (apr?s transformation par fct prepar-cov()
#'            par exemple)
#' @param choix_dist (FACULTATIF) Distance euclidienne "E" par d?faut, sinon "H" pour
#'             distance de Hamming
#'
#' @return
#' @export
#'
#' @examples
#'
#' summary(datae_new)
#' t1 = Sys.time(); xx = cout(datae_new);     t2 = Sys.time(); difftime(t2,t1); xx
#' t1 = Sys.time(); zz = cout(datae_new,"H"); t2 = Sys.time(); difftime(t2,t1); zz
#'
#' summary(newBDD)
#' t1 = Sys.time(); c1 = cout(newBDD);      t2 = Sys.time(); difftime(t2,t1); c1
#'
#' summary(newBDD2)
#' t1 = Sys.time(); c2 = cout(newBDD2);     t2 = Sys.time(); difftime(t2,t1); c2
#' t1 = Sys.time(); c3 = cout(newBDD2,"H"); t2 = Sys.time(); difftime(t2,t1); c3
#'
#'
#' # AJOUT DES NA: On sp?cifie les Y manquants en base 1 et 2: --> OK aussi
#' datae_new_NA = datae_new
#' datae_new_NA[datae_new_NA[,1]==1,"Ybase2"] = NA
#' datae_new_NA[datae_new_NA[,1]==2,"Ybase1"] = NA
#' t1 = Sys.time(); xx_NA = cout(datae_new_NA);
#' t2 = Sys.time(); difftime(t2,t1); xx_NA
cout = function(dat, choix_dist = "E"){


  val1  = sort(unique(dat[dat[,1]==1,2]))    # Les modalit?s de Y dans la base 1
  val2  = sort(unique(dat[dat[,1]==2,3]))    # Les modalit?s de Y dans la base 2


  Risk=0  # Initialisation du vecteur Risk
  V=0


  #  CALCUL DES SOMMES DES DISTANCES INDIVIDUELLES ENTRE 2 GROUPES

  for (i in 1:length(val1)){

    for (j in 1:length(val2)){

      V = 0

      sujetsB1 = subset(dat,dat[,1] == 1 & dat[,2] == val1[i])
      sujetsB2 = subset(dat,dat[,1] == 2 & dat[,3] == val2[j])

      for (m in (1:nrow(sujetsB1))){

        for (n in (1:nrow(sujetsB2))){

          for (k in (4:ncol(dat))){

            if (choix_dist == "H"){

              # Distance de Hamming
              V = V + as.numeric(sujetsB1[m,k]!=sujetsB2[n,k])

            } else {

              # Distance Euclidienne
              V = V + (sujetsB1[m,k]-sujetsB2[n,k])^2

            }

          }

        }

      }

      V=V/((ncol(dat)-3)*nrow(sujetsB1)*nrow(sujetsB2))
      Risk = c(Risk,V)

    }

  }


  # Suppression du 0 initial de la fonction de risque
  Risk = Risk[-1]


  # On pr?pare A, x et b de la contrainte Ax = b en vue d'appliquer le simplexe

  # Stockage des modalit?s max de la base 1 et 2

  maxi = c(length(val1),length(val2))


  # Contruction de la matrice compl?te A(S):

  S = matrix(ncol = maxi[1]*maxi[2],nrow = maxi[1]+maxi[2])


  # S?quences utiles de 1 et de 0 avant affectation dans matrice

  seq1 = rep(1,maxi[2])
  seq0 = rep(0,maxi[2])
  seq2 = c(1,rep(0,maxi[2]-1))


  # Remplissage des lignes de S correspondant aux 1?res modalit?s de Yb1 et Yb2

  S[1,]         = L1 = c(seq1,rep(seq0,maxi[1]-1))
  S[maxi[1]+1,] = L2 = rep(seq2,maxi[1])


  # On compl?te les lignes restantes ? partir des 2 lignes pr?c?dentes:

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


  # le "b" de la contrainte:

  CVal = as.numeric(c(table(dat[dat[,1]==1,2])/nrow(dat[dat[,1]==1,]),
                      table(dat[dat[,1]==2,3])/nrow(dat[dat[,1]==2,])))



  # Simplexe avec R: On utilise la fonction solveLP(linprog)

  res = linprog::solveLP(cvec=Risk, bvec=CVal, Amat=S, const.dir = rep( "==",
                                                               length(CVal)),lpSolve=TRUE)


  x = res$solution

  # Cette fonction retourne la matrice S, le vecteur de Risque et la solution
  # du simplexe:

  return(list(Risk=Risk,CVal=CVal,A=S,solution=x))

}









#' simpond(x)
#' Cette fonction propose des arrondis aux effectifs de la matrice de simplexe
#' lorsque cette dernière n'est pas uniquement composée d'entiers
#' En cas d'ex-aecquos, le premier pr?sent? est le premier servi
#'
#' @param x un vecteur
#' @return un vecteur aussi
#' @export
#' @examples
#' x = c(0.5,0,0.5,0.5,3.5); x; sum(x)
#' simplex_new = simpond(x); simplex_new; sum(simplex_new)
# #[1] 0.5 0.0 0.5 0.5 3.5
# #[1] 5
# 'x = c(127,18,0,14); x; sum(x)
# 'simplex_new = simpond(x); simplex_new; sum(simplex_new)
# # [1] 127  18  0  14
# # [1] 159
# 'x = c(0.5,3.5,0.25,0.75,0.5,0.5); sum(x)
#' simplex_new = simpond(x); simplex_new; sum(simplex_new)
# #[1] 1 4 0 1 0 0
# #[1] 6
simpond = function(x){

  y = floor(x)

  if (sum(x)==sum(y)){xnew = y} else {

    nb_malC    = sum(x) - sum(y)
    candidat   = ceiling(x)-y
    ind_cand   = (1:length(candidat))[candidat == 1]
    val_cand   = x[ind_cand] - floor(x[ind_cand])
    choix_cand = ind_cand[order(val_cand,decreasing=T)][1:nb_malC]
    ajt        = rep(0,length(x))
    ajt[choix_cand] = rep(1,length(choix_cand))

    xnew = y + ajt

  }

  return(xnew)

}



#------------------ Fonction AFFECTATION ------------------------------------------------------------------------------------------------

#' Cette fonction est la partie II/II de la methode FOT
#'
#' @param simplx Vecteur des couts de passage estime par la fonction cout()
#' @param dat BDD de la forme [base|Y1|Y2|???], ??? etant des covariables binaires
#' ou quantis (apres transformation par fct prepar-cov() par exemple)
#' @param choix_dist (FACULTATIF) Distance euclidienne par defaut sinon "H" pour distance
#' de Hamming
#'
#' @export
#'
#' @examples
#'
#' t1 = Sys.time(); yy = affectation(xx$solution,datae_new);
#' t2 = Sys.time(); difftime(t2,t1); yy
#'
#' t1 = Sys.time(); yy = affectation(xx_NA$solution,datae_new_NA);
#' t2 = Sys.time(); difftime(t2,t1); yy
#'
#' t1 = Sys.time(); yy = affectation(zz$solution,datae_new,"H");
#' t2 = Sys.time(); difftime(t2,t1); yy
#'
#' t1 = Sys.time(); yy = affectation(c1$solution,newBDD);
#' t2 = Sys.time(); difftime(t2,t1); yy
#'
#' t1 = Sys.time(); yy = affectation(c2$solution,newBDD2);
#' t2 = Sys.time(); difftime(t2,t1); yy
#'
#' t1 = Sys.time(); yy = affectation(c3$solution,newBDD2,"H");
#' t2 = Sys.time(); difftime(t2,t1); yy
#'
#' set.seed(410747+4);
#' datgen = gendata(10,0.5,0,c(0.90,0.10),c(0.90,0.10))
#' datgen2 = datgen
#' datgen2[datgen2[,1]==1,"Y2"] = NA
#' datgen2[datgen2[,1]==2,"Y1"] = NA
#' affectation(cout(datgen2,"E")$solution,datgen2)
affectation = function(simplx,dat,choix_dist = "E"){


  n1 = nrow(dat[dat[,1]==1,])               # nbr de sujets ds base 1
  n2 = nrow(dat[dat[,1]==2,])               # nbr de sujets ds base 2

  # Ajout d'une colonne d'identification des sujets
  dat$suj=1:nrow(dat)

  # Ajout d'une colonne d'imputation
  dat$FOT = NA

  val1  = sort(unique(dat[dat[,1]==1,2]))   # Les modalit?s de Y dans la base 1
  val2  = sort(unique(dat[dat[,1]==2,3]))   # Les modalit?s de Y dans la base 2



  # I. CONVERSION DES Y[base 1] EN Y[base 2]

  # 1) On r?cup?re les indices du max. de la matrice de COUT:

  mat_simplx = matrix(round(simplx*n1,4),nrow=length(val2),
                      ncol=length(val1))

  if (sum(mat_simplx == floor(mat_simplx))==
      ncol(mat_simplx)*nrow(mat_simplx)){

    mat_simplx = mat_simplx} else {
      mat_simplx = apply(mat_simplx,2,simpond)}


  ind = which(mat_simplx==max(mat_simplx), arr.ind = TRUE)

  # Gestion des max. ex-aecquos: On conserve l'indice du 1er rencontr?
  if (is.matrix(ind)== TRUE){ind=ind[1,]}


  # 2) On cr?e une nouvelle base dont on enlevera au fur et a mesure les sujets
  #    pour lesquels une valeur de Y a ?t? imput?e en base 1

  dat2 = dat


  # Tant que la matrice de passage n'est pas vide, on r?p?te les conditions
  # suivantes:

  while (length(dat$FOT[(dat[,1] == 1) & (is.na(dat$FOT)==TRUE)])!=0){


    # On extrait des bases 1 et 2, les 2 sous-groupes de sujets
    # correspondant aux indices du max retenu pr?c?demment:

    sujetsB1 = subset(dat2,dat2[,1] == 1 & dat2[,2] == val1[ind[2]])
    sujetsB2 = subset(dat2,dat2[,1] == 2 & dat2[,3] == val2[ind[1]])



    #  CALCUL DE LA SOMME DES DISTANCES ENTRE UN SUJET APPARTENANT A UN
    #  GROUPE DE LA BASE 1 ET TOUS LES AUTRES SUJETS D'UN GROUPE DE
    #  LA BASE 2

    # https://johanndejong.wordpress.com/2015/09/23/fast-hamming-distance-in-r/

    dist = numeric(nrow(sujetsB1))

    for (m in (1:nrow(sujetsB1))){

      V =0

      for (n in (1:nrow(sujetsB2))){

        for (k in (4:(ncol(dat2)-2))){

          if (choix_dist == "H"){

            # Distance de Hamming
            V = V + as.numeric(sujetsB1[m,k]!=sujetsB2[n,k])

          } else {

            # Distance Euclidienne
            V = V + (sujetsB1[m,k]-sujetsB2[n,k])^2

          }

        }
      }

      dist[m] = V / (nrow(sujetsB2)*(ncol(dat)-3))

    }


    indo = sujetsB1$suj[order(dist)]

    suj = numeric(0)


    # On classe les patients appartenant au sous-groupe extrait de la
    # base 1 par distance moyenne croissante de chacun de ces patients
    # aux autres sujets du groupe extrait de la base 2.
    # On code les Y de ces sujets en base 2 en respectant la marge
    # indiqu?e par mat_simplx[ind]


    for (k in (1:mat_simplx[ind[1],ind[2]])){

      dat[indo[k],ncol(dat)] = val2[ind[1]]
      suj = c(suj,indo[k])

    }

    # on enl?ve les sujets de la base 1 pour lesquels Y vient d'?tre
    # recod? en base 2:

    dat2=dat2[!(dat2$suj %in% suj),]


    # Une fois le max. de la matrice de passage exploit?, on le remplace
    # par un 0:

    mat_simplx[ind[1],ind[2]] = 0


    # On recherche un nouveau max. et on s?lectionne le 1er rencontr?
    # en cas d'ex-aecquos

    ind = which(mat_simplx==max(mat_simplx), arr.ind = TRUE)
    if (is.matrix(ind)== TRUE){ind=ind[1,]}

  }


  # II. CONVERSION DES Y[base 2] EN Y[base 1]

  # 1) On r?cup?re les indices du max. de la matrice de COUT:


  mat_simplx = matrix(round(simplx*n2,4),nrow=length(val2),ncol=length(val1))

  if (sum(mat_simplx == floor(mat_simplx))==
      ncol(mat_simplx)*nrow(mat_simplx)){

    mat_simplx = mat_simplx} else {
      mat_simplx = t(apply(mat_simplx,1,simpond))}


  ind = which(mat_simplx==max(mat_simplx), arr.ind = TRUE)

  # Gestion des max. ex-aecquos: On conserve l'indice du 1er rencontr?
  if (is.matrix(ind)== TRUE){ind=ind[1,]}


  # 2) On cr?e une nouvelle base dont on enlevera au fur et a mesure les sujets
  #    pour lesquels une valeur de Y a ?t? imput?e en base 2

  dat2 = dat

  # Tant que la matrice de passage n'est pas vide, on r?p?te les conditions
  # suivantes:


  while (length(dat$FOT[is.na(dat$FOT)==TRUE])!=0){

    # On extrait des bases 1 et 2, les 2 sous-groupes de sujets
    # correspondant aux indices du max retenu pr?c?demment:

    sujetsB1 = subset(dat2,dat2[,1] == 1 & dat2[,2] == val1[ind[2]])
    sujetsB2 = subset(dat2,dat2[,1] == 2 & dat2[,3] == val2[ind[1]] )


    #  CALCUL DE LA SOMME DES DISTANCES ENTRE UN SUJET APPARTENANT A UN
    #  GROUPE DE LA BASE 2 ET TOUS LES AUTRES SUJETS D'UN GROUPE DE
    #  LA BASE 1


    dist = numeric(nrow(sujetsB2))

    for (m in (1:nrow(sujetsB2))){

      V =0

      for (n in (1:nrow(sujetsB1))){


        for (k in (4:(ncol(dat2)-2))){

          if (choix_dist == "H"){

            # Distance de Hamming
            V = V + as.numeric(sujetsB1[m,k]!=sujetsB2[n,k])

          } else {

            # Distance Euclidienne
            V = V + (sujetsB1[n,k]-sujetsB2[m,k])^2

          }

        }
      }

      dist[m] = V / (nrow(sujetsB1)*(ncol(dat)-3))

    }

    indo = sujetsB2$suj[order(dist)]

    suj = numeric(0)


    # On classe les patients appartenant au sous-groupe extrait de la
    # base 2 par distance moyenne croissante de chacun de ces patients
    # aux autres sujets du groupe extrait de la base 1.
    # On code les Y de ces sujets en base 1 en respectant la marge
    # indiqu?e par mat_simplx[ind]


    for (k in (1:mat_simplx[ind[1],ind[2]])){

      dat[indo[k],ncol(dat)] = val1[ind[2]]
      suj=c(suj,indo[k])

    }

    # on enl?ve les sujets de la base 2 pour lesquels Y vient d'?tre
    # recod? en base 1:

    dat2=dat2[!(dat2$suj %in% suj),]


    # Une fois le max. de la matrice de passage exploit?, on le remplace
    # par un 0:

    mat_simplx[ind[1],ind[2]] = 0


    # On recherche un nouveau max. et on s?lectionne le 1er rencontr?
    # en cas d'ex-aecquos:

    ind = which(mat_simplx==max(mat_simplx), arr.ind = TRUE)
    if (is.matrix(ind)== TRUE){ind=ind[1,]}


  }


  return(dat)

}

#' MICE - Liste des methodes possible pour l'imputation:
#'
#' pmm         :  Predictive mean matching (any)
#' norm        : Bayesian linear regression (numeric)
#' norm.nob    : Linear regression ignoring model error (numeric)
#' norm.boot   : Linear regression using bootstrap (numeric)
#' norm.predict: Linear regression, predicted values (numeric)
#' mean        : Unconditional mean imputation (numeric)
#' 2l.norm     : Two-level normal imputation (numeric)
#' 2l.pan      : Two-level normal imputation using pan (numeric
#' 2lonly.mean : Imputation at level-2 of the class mean (numeric)
#' 2lonly.norm : Imputation at level-2 by Bayesian linear regression (numeric)
#' 2lonly.pmm  : Imputation at level-2 by Predictive mean matching (any)
#' quadratic   : Imputation of quadratic terms (numeric)
#' logreg      : Logistic regression (factor, 2 levels)
#' logreg.boot : Logistic regression with bootstrap
#' polyreg     : Polytomous logistic regression (factor, >= 2 levels)
#' polr        : Proportional odds model (ordered, >=2 levels)
#' lda         : Linear discriminant analysis (factor, >= 2 categories)
#' cart        : Classification and regression trees (any)
#' rf          : Random forest imputations (any)
#' ri          : Random indicator method for nonignorable data (numeric)
#' sample      : Random sample from the observed values (any)
#' fastpmm     : Experimental: Fast predictive mean matching using C++ (any)
#'
#' @param datessai Jeu de donn?es AVEC les NA sur Y1 et Y2
#'            requiert que les outcomes s'appellent "Y1" et "Y2"
#' @param R_mice Nbr de r?p?titions pour MICE
#' @param meth Vecteur des m?thodes d'imputation ? appliquer ?
#' chaque variable et covariable selon son type: quanti ou facteur "pmm" si toutes quantis
#' "polyreg" si qualis (facteur) ? >= 2 modalit?s
#' @param calc Doit-on moyenner les valeurs imput?es ou conserver les plus
#'            fr?quentes ? ="MOY" ou ="FREQ"
#'
#' @return un resultat
#' @export
#'
#' @examples
#'
#' colnames(datae_new)[2:3] = c("Y1","Y2")
#' datae_new[datae_new[,1]==1,"Y2"] = NA
#' datae_new[datae_new[,1]==2,"Y1"] = NA
#'
#' imp_MICE(datae_new,5,"pmm","FREQ")
#' #[1] 4 3 4 4 3 2 1 4 4 1 3 1 3 3 3 2 2 2 3 2
#' set.seed(1236); imp_MICE(datae_new,5,"pmm","FREQ")
#' #[1] 4 4 4 4 1 1 1 4 3 1 2 1 3 3 3 1 1 2 3 1
#'
#' datae_nw  = apply(datae_new[,2:5],2,as.factor)
#' datae_nw2 = data.frame(base = datae_new[,1],datae_nw)
#' summary(datae_nw2)
#' set.seed(1236); imp_MICE(datae_nw2,5,c("polyreg","polyreg","pmm","pmm"),"FREQ")
#' #[1] 4 3 4 3 2 1 1 4 3 1 3 1 3 3 3 2 1 1 3 2
imp_MICE = function(datessai,R_mice,meth,calc){

  # Imputation MICE:

  stoc_mice = mice::mice(as.data.frame(datessai[,2:ncol(datessai)]),method = meth,
                   m = R_mice, print=FALSE)


  # Stockage des imputations r?p?t?es dans un data.frame

  base_mice_Y1 = as.data.frame(datessai[,"Y1"])
  base_mice_Y2 = as.data.frame(datessai[,"Y2"])

  for (u in 1:R_mice){

    base_mice_Y1 = data.frame(base_mice_Y1,mice::complete(stoc_mice,u)[,"Y1"])
    base_mice_Y2 = data.frame(base_mice_Y2,mice::complete(stoc_mice,u)[,"Y2"])

  }


  # Imputation de la cat?gorie la plus fr?quente (dans le cas de donn?es
  # cat?gorielles) - La covariable "base" n'est pas incluse dans le mod?le

  if (calc=="FREQ"){

    newcol_Y1 = apply(base_mice_Y1[,2:ncol(base_mice_Y1)],1,
                      function(x){as.integer(names(table(x))[which.max(table(x))])})

    newcol_Y2 = apply(base_mice_Y2[,2:ncol(base_mice_Y2)],1,
                      function(x){as.integer(names(table(x))[which.max(table(x))])})

  }


  # Imputation de la moyenne (dans le cas de donn?es continues)

  if (calc=="MOY"){

    newcol_Y1 = apply(base_mice_Y1[,2:ncol(base_mice_Y1)],1,mean)
    newcol_Y2 = apply(base_mice_Y2[,2:ncol(base_mice_Y2)],1,mean)

  }

  imput_mice = c(newcol_Y2[datessai[,1]==1],newcol_Y1[datessai[,1]==2])

  return(imput_mice)

}




