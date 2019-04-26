#-------------------------------------------------------------------------------
#     II.  FONCTIONS DE PREPARATION DES DONNEES
#-------------------------------------------------------------------------------



#' fonction TRANSFO
#'   Cette fonction transforme une variables qualitative (format facteur) ? k
#'   classes en (k-1) variables binaires
#'  
#'
#' @param x Variable facteur
#' @param nom Possibilite de donner un nom a votre variable pour mieux reperer les
#'        binaires correspondantes en sortie
#'        La modalite de reference est la 1ere (cf levels(x) au prealable)
#'
#' @return
#' @export
#'
transfo = function(x,nom=NULL){

  # Il faut que la variable en entr?e soit au format facteur:

  if (is.factor(x)==FALSE){
    stop("Convertissez votre variable en facteur svp !!!")
  }

  lev_x = levels(x)


  # Gestion des covariables ? une seule modalit?:

  if (length(lev_x) == 1){

    # Distinction n?cessaire lorsque l'unique modalit? est 0
    #(x ?tant en facteur)

    if (lev_x == "0"){

      x_bin = as.matrix(as.numeric(x)-1)
      colnames(x_bin)[1] = nom

    } else if (lev_x != "0"){

      x_bin = as.matrix(rep(1,length(x)))
      colnames(x_bin)[1] = nom


    }

  } else {

    # Nb modalit?s >= 2
    # Stockage de toutes les binaires dans une matrice:

    x_bin = matrix(nrow = length(x), ncol = length(lev_x) - 1)

    for (j in 1:(length(lev_x)-1)){

      x_bin[,j] = ifelse(x == lev_x[j+1],1,0)

    }


    # On nomme les binaires g?n?r?es:

    colnames(x_bin) = paste(nom,2:(length(lev_x)),sep="_")

  }

  return(x_bin)

}




#' CAS 1: TOUTES LES COVARIABLES SONT QUALITATIVES
#'
#' Cette fonction effectue un travail pr?paratoire sur les covariables
#'
#' Si toutes les covariables sont qualis --> Tableau Disjonctif Complet(binaires)
#' Si les covariables sont qualis et quantis  --> Coordonn?es sur CP d'une AFDM
#' Si toutes les covariables sont quantis --> Simple standardisation
#'
#' @param dat   : BDD compos?e de 3 1?res colonnes fixes et ordonn?es (base,Y1,Y2)
#'         la partie variable est compos?e de covariables souhait?es qualis
#'         (?ventuellement) puis quantis (?ventuellement)
#' @param quali : A REMPLIR UNIQUEMENT si certaines covariables sont souhait?es qualis,
#'         indiquez les n? de colonnes correspondant
#' @param quanti: A REMPLIR UNIQUEMENT si certaines covariables sont souhait?es quantis,
#'         indiquez les n? de colonnes correspondant
#' @param info  : A REMPLIR UNIQUEMENT si vous avez un m?lange de covariables quantis et
#'         qualis ... % d'information mini prise en compte par les CP de l'AFDM
#'-------------------------------------------------------------------------------
#'
#' @return
#' @export
#'
prepar_cov = function(dat, quali = NULL,quanti = NULL, info = 90){
  #------------------------------------------------

  if ((length(quanti)==0)&(length(quali)!=0)){


    # On convertit les covariables souhait?es qualitatives en facteur:
    # (simple ?tape de pr?caution)

    dat2 = apply(dat[,c(1,quali)],2,as.factor)
    dat3 = data.frame(dat2[,1],dat[,2:3],dat2[,-1],dat[,quanti])
    colnames(dat3) = colnames(dat)


    # R?cup?ration des noms des covariables qualitatives

    nom_col = colnames(dat3)[quali]


    # Transformation des variables qualis ? k modalit?s en (k-1) binaires:
    # La 1?re modalit? est la r?f?rence

    bin_quali     = transfo(dat3[,quali[1]])
    nbmod         = length(levels(dat3[,quali[1]]))


    # Transformation des noms de colonnes

    if (nbmod >=2){

      nom_col_quali = paste(nom_col[1],2:nbmod,sep="_")

    } else {

      nom_col_quali = nom_col[1]

    }

    # G?n?ralisation ? T variables (T>1)

    for (i in 2:length(quali)){

      bin_quali     = cbind(bin_quali,transfo(dat3[,quali[i]]))
      nbmod         = length(levels(dat3[,quali[i]]))

      if (nbmod >=2){

        nom_col_quali = c(nom_col_quali,paste(nom_col[i],2:nbmod,sep="_"))

      } else {

        nom_col_quali = c(nom_col_quali,nom_col[i])

      }

    }

    # R?cup?ration du jeu de donn?es transform?

    tabnew = data.frame(dat[,1:3],bin_quali,dat[,quanti])
    colnames(tabnew) = c(colnames(dat)[1:3],nom_col_quali,
                         colnames(dat3)[quanti])

  }


  # CAS 2: LES COVARIABLES SONT QUALITATIVES ET QUANTITATIVES
  #----------------------------------------------------------

  if ((length(quanti)!=0)&(length(quali)!=0)){

    #  AFDM: Ref: Pages J. (2004). Analyse factorielle de donnees mixtes.
    #        Revue Statistique Appliquee. LII (4). pp. 93-111.


    # ... par d?faut, l'AFDM ne conserve que les coordonn?es des 5 premi?res CP.
    #     --> Option ncp ? modifier si n?cessaire, graph = F pour ne pas tracer
    #         les graphes de l'AFDM
    # ... La variable indiquant la base est ici incluse en facteur dans l'AFDM:

    dat2 = apply(dat[,c(1,quali)],2,as.factor)
    dat3 = data.frame(dat2[,1],dat[,2:3],dat2[,-1],dat[,quanti])
    colnames(dat3) = colnames(dat)

    dat4 = dat3[,-(2:3)]
    res = FactoMineR::FAMD(dat4, ncp = 9, graph = F)


    #  On conserve les coordonn?es des k?mes premi?res CPs qui prennent en
    #  compte plus de 90% de l'information:

    nb_CP = min((1:nrow(res$eig))[res$eig[,3]> info])


    # Les coordonn?es individuelles sont d?j? centr?es et r?duites d'une CP ?
    # l'autre:

    tabnew = data.frame(dat[,1:3],res$ind$coord[,1:nb_CP])


  }


  # CAS 3: LES COVARIABLES SONT QUANTITATIVES
  #------------------------------------------

  if ((length(quanti)!=0)&(length(quali)==0)){

    # On standardise juste les variables quantitatives:

    dat2 = scale(dat[,quanti])
    tabnew = dat
    tabnew[,quanti] = dat2


  }

  return(tabnew)

}


