#--------- fonction GENDATA_ND  [NON DETERMINISTE]

#' SIMULATION DE DONNEES POUR APPROCHE NON DETERMINISTE
#' -> Cette fonction genere une BDD de la forme [Base|Y1|Y2|X1|X2|X3] avec les
#'    covariables X1, X2 categorielles ? k1 et K2 modalites resp., et X3 continue
#'
#' @param n1     : Nombre de patients dans la base 1
#' @param k      : Coefficient de proportionnalite entre les 2 bases
#' @param cormat : vecteur contenant les correlations: rho = cor(X1,X2),
#'                 delta = cor(X1,X3), mu = cor(X2,X3)
#' @param rho    : Correlation entre X1 et X2 (continues)
#' @param delta  : Correlation entre X1 et X3 (continues)
#' @param mu     : Correlation entre X2 et X3 (continues)
#' @param px1c   : Proportion de patients dans chaque classe de X1
#' @param px2c   : Proportion de patients dans chaque classe de X2
#' @param  r2    : R^2 de la regression lineaire Ycontinue = f(X1,X2,X3)
#' @param py1c   : Proportion de patients dans chaque classe de Y en base 1
#' @param py2c   : Proportion de patients dans chaque classe de Y en base 2
#' @param valY1  : Les labels de chaque classe pour Y cat?gorielle en base 1
#' @param valY2  : Les labels de chaque classe pour Y cat?gorielle en base 2
#'
#' @return
#' @export
#'
#' # Quel que soit le nombre de classes de X1 et X2, elles 
#' # sont au format numerique en sortie
gendata_ND =function(n1,k,cormat,px1c,px2c,r2,py1c,py2c,valY1,valY2){

  # Si la somme des proportions des cat?gories des chacune des covariables
  # ne vaut pas 1--> Erreur

  if ((sum(px1c)!=1)|(sum(px2c)!=1)|(sum(py1c)!=1)|(sum(py2c)!=1)){
    stop("Distribution incorrecte pour au moins 1 covariable!")
  }

  # Taille de la base 2

  n2 = ceiling(k*n1)
  n  = n1 + n2

  # Matrice de corr?lation entre X1, X2 et X3

  rho   = cormat[1]
  delta = cormat[2]
  mu    = cormat[3]
  corrcov = matrix(c(1,rho,delta,rho,1,mu,delta,mu,1),ncol=3)

  # Les distributions cumul?es souhait?es pour X1 et X2

  px1cc = cumsum(px1c[1:(length(px1c)-1)])
  px2cc = cumsum(px2c[1:(length(px2c)-1)])


  # Simulation du vecteur norm? centr? (X1,X2,X3) et r?cup?ration des quantiles
  # correspondants pour X1 et X2:

  X_glob = mvtnorm::rmvnorm(n,mean=c(0,0,0),sigma = corrcov)

  qx1c = qnorm(px1cc,mean = 0, sd = 1)
  qx2c = qnorm(px2cc,mean = 0, sd = 1)


  # Discr?tisation des covariables X1 et X2 en fonction des dsitributions
  # souhait?es

  X1 = cut(X_glob[,1],breaks= c(min(X_glob[,1])-1,qx1c,max(X_glob[,1])+1))
  X2 = cut(X_glob[,2],breaks= c(min(X_glob[,2])-1,qx2c,max(X_glob[,2])+1))
  X3 =  X_glob[,3]
  levels(X1) = 1:length(px1c)
  levels(X2) = 1:length(px2c)


  # Conversion en covariables discr?tes

  X1 = as.numeric(X1)-1
  X2 = as.numeric(X2)-1


  # Simulation de Y quanti sachant que Y ~ X1 + X2 + X3 et ? R2 fix? (cf pdf)

  sigma2 = ((3 + 2*rho + 2*delta + 2*mu)*(1 - r2))/r2

  Y = X_glob[,1] + X_glob[,2] + X_glob[,3] + rnorm(n,0,sqrt(sigma2))


  # --> sd(Y) ?tant diff?rent de 1, on la standardise:

  #Ynorm = (Y-mean(Y))/sqrt(var(Y))

  #Ynorm = (Y_0)/sqrt(3 + 2*rho + 2*delta + 2*mu + sigma2)
  #Ynorm= (Y-mean(Y))/sqrt(var(Y))
  # On discr?tise Y en 2 types de cat?gories suivant que l'on se place en
  # base1 ou en base2
  Ynorm= Y
  Y1 = Y2 = vector(length=n)

  #py1cc = cumsum(py1c[1:(length(py1c)-1)])
  #py2cc = cumsum(py2c[1:(length(py2c)-1)])

  #qy1c = qnorm(py1cc,mean = 0, sd = 1)

  Y1 = cut(Ynorm, breaks = c(min(Ynorm)-1,quantile(Ynorm, c( 0.25, 0.5, 0.75)),max(Ynorm)+1), include.lowest = TRUE,labels=valY1)
  Y2 = cut(Ynorm, breaks = c(min(Ynorm)-1,quantile(Ynorm, c(1/3,  2/3)),max(Ynorm)+1), include.lowest = TRUE,labels=valY2)

  # Conversion en outcomes discrets

  Y1 = as.numeric(factor(Y1,levels(Y1)[valY1]))
  Y2 = as.numeric(factor(Y2,levels(Y2)[valY2]))


  # Identification de la base et r?cup?ration de la BDD simul?e

  ident = c(rep(1,n1),rep(2,n2))
  dat = data.frame(ident,Y1,Y2,X1,X2,X3)

  return(dat)

}

