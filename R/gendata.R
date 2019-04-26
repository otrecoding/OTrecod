#-------------------------------------------------------------------------------
#     I.  FONCTIONS DE SIMULATION DE DONNEES
#-------------------------------------------------------------------------------

 #j'essai de comprendre github
#--------- fonction GENDATA  [DETERMINISTE]

#' SIMULATION DE DONNEES POUR APPROCHE DETERMINISTE
#'  
#'   -> Y1 et Y2 sont enti?rement d?termin?es par les valeurs de X1 et X2
#'   -> Cette fonction g?n?re une BDD de la forme [Base|Y1|Y2|X1|X2] si X1 et X2
#'      sont cat?gorielles ? 2 modalit?s en m?me temps.
#'   -> Cette fonction g?n?re une BDD de la forme [Base|X1|X2] si X1 et X2
#'      sont cat?gorielles avec X1 ou X2 comptent plus de 2 modalit?s
#'  
#' @param n1 Nombre de sujets dans la base 1
#' @param k Coefficient de proportionnalit? entre le nbr de sujets dans la Base 1 et le nombre de sujets dans la Base 2 
#' @param rho Coefficient de corr?lation entre les covariables X1 et X2 
#' @param p1c Distribution de la covariable X1 cat?gorielle (format: c(,))c
#' @param p2c Distribution de la covariable X2 cat?gorielle (format: c(,)) 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' essai1 = gendata(20,1,0.3,p1c = c(0.2,0.4,0.1,0.3),p2c = c(0.6,0.4))
#' essai2 = gendata(20,1,0.3,p1c = c(0.6,0.4),p2c = c(0.7,0.3))
gendata = function(n1,k,rho,p1c,p2c){

  # Si la somme des proportions des cat?gories des chacune des covariables
  # ne vaut pas 1--> Erreur

  if ((sum(p1c)!=1)| (sum(p2c)!=1)) stop("Distribution incorrecte pour au
                                         moins 1 covariable!")


  # Si k*n1 n'est pas un nbr entier, on prend le 1er entier sup?rieur

  n2 = ceiling(k*n1)
  n  = n1 + n2
  corrcov = matrix(c(1,rho,rho,1),ncol=2)


  # Distributions cumul?es

  p1cc = cumsum(p1c[1:(length(p1c)-1)])
  p2cc = cumsum(p2c[1:(length(p2c)-1)])


  # Simulation du vecteur norm? centr?(X1,X2) et r?cup?ration des quantiles
  # correspondants:

  Y_glob = mvtnorm::rmvnorm(n,mean=c(0,0),sigma = corrcov)
  q1c = qnorm(p1cc,mean = 0, sd = 1)
  q2c = qnorm(p2cc,mean = 0, sd = 1)


  # Discr?tisation des covariables X1 et X2

  X1 = cut(Y_glob[,1],breaks= c(min(Y_glob[,1])-1,q1c,max(Y_glob[,1])+1))
  X2 = cut(Y_glob[,2],breaks= c(min(Y_glob[,2])-1,q2c,max(Y_glob[,2])+1))
  levels(X1) = 1:length(p1c)
  levels(X2) = 1:length(p2c)


  # Variable d'identification de la Base

  ident = c(rep(1,n1),rep(2,n2))
  dat = data.frame(ident)


  # On g?n?re Y1 et Y2 en fonction des r?gles d?terministes choisies:

  if ((length(levels(X1))==2)&(length(levels(X2))==2)){

    X1 = as.numeric(X1)-1
    X2 = as.numeric(X2)-1

    dat$Y1[X1 == 0 & X2 == 0] = 3
    dat$Y2[X1 == 0 & X2 == 0] = 4
    dat$Y1[X1 == 0 & X2 == 1] = 2
    dat$Y2[X1 == 0 & X2 == 1] = 3
    dat$Y1[X1 == 1 & X2 == 0] = 3
    dat$Y2[X1 == 1 & X2 == 0] = 2
    dat$Y1[X1 == 1 & X2 == 1] = 1
    dat$Y2[X1 == 1 & X2 == 1] = 1

  }

  # R?cup?ration de la BDD simul?e:
  # Si X1 et X2 ne sont pas binaires, elles restent en facteur dans dat.

  dat = data.frame(dat,X1,X2)

  return(dat)

}

