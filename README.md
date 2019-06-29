
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OTrecod

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/otrecoding/OTrecod.svg?branch=master)](https://travis-ci.org/otrecoding/OTrecod)
[![CRAN
status](https://www.r-pkg.org/badges/version/OTrecod)](https://cran.r-project.org/package=OTrecod)
[![codecov](https://codecov.io/gh/otrecoding/OTrecod/branch/master/graph/badge.svg)](https://codecov.io/gh/otrecoding/OTrecod)
[![Launch
binder](http://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/otrecoding/OTrecod/master)
<!-- badges: end -->

The goal of OTrecod is to …

## Installation

You can install the released version of OTrecod from
[CRAN](https://CRAN.R-project.org) with:

``` r
# Install release version from CRAN
install.packages("OTrecod")
# Install development version from GitHub
devtools::install_github("otrecoding/OTrecod")
```

## Example

This is a basic example which shows you how to solve a common problem:

todo: insert a working {r example}

``` r
n_1         = 1000
coeff       = 1
p1cc        = c(0.5,0.5)
p2cc        = c(0.3,0.3,0.4)
qual        = 4:5
quan        = 6
choose_dist = "E"
rho_cov     = 0.2
rep_MICE    = 5
cov_MICE    = c("polyreg","polyreg","polyreg","polyreg","pmm")
calcu       = "FREQ"
typ_simu    = "ND"
r2t         = 0.5
py1cc       = c(0.25,0.25,0.25,0.25)
py2cc       = c(0.333,0.334,0.333)
vlY1        = c(1,2,3,4)
vlY2        = c(1,2,3)
simul_glob(n_1, coeff, p1cc, p2cc, qual,
                        quan ,choose_dist="E", rho_cov, rep_MICE=5,
                        cov_MICE , calcu="FREQ",typ_simu,r2t,py1cc,py2cc,vlY1 , vlY2 )
```

    base    X1  X2  Ybase1  Ybase2
    1   1   1   4   3
    1   1   0   3   2
    1   1   1   4   3
    1   1   0   3   2
    1   0   1   2   3
    1   0   0   1   1
    1   0   0   1   1
    1   1   1   4   3
    1   1   0   3   2
    1   0   0   1   1
    2   0   1   2   3
    2   0   0   1   1
    2   1   1   4   3
    2   1   1   4   3
    2   1   1   4   3
    2   1   0   3   2
    2   0   0   1   1
    2   0   0   1   1
    2   1   1   4   3
    2   1   0   3   2

## Generate data

    # Taille des bases
    
    taille_n1    = c(50,100,500,1000,5000,10000)
    ratio_k      = c(0.25,0.5,0.75,1)
    
    
    # Propriétés des covariables
    
    rho_cov      = c(0.2,0.4,0.6,0.8)
    r2_reg       = c(0.2,0.4,0.6,0.8)
    
    vlis1 <- function(nsim=100) {
      vList <- simsalapar::varlist(
        n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
        n_1         = list(value =c(50,100,500,1000,5000,10000),type="grid"),
        coeff       = list(value =1),
        p1cc        = list(value =c(0.5,0.5)),
        p2cc        = list(value =c(0.3,0.3,0.4)),
        qual        = list(value =4:5),
        quan        = list(value =6),
        choose_dist = list(value ="E"),
        rho_cov     = list(value =0.2),
        rep_MICE    = list(value =5),
        cov_MICE    = list(value =c("polyreg","polyreg","polyreg","polyreg","pmm")),
        calcu       = list(value ="FREQ"),
        typ_simu    = list(value ="ND"),
        r2t         = list(value =0.5),
        py1cc       = list(value =c(0.25,0.25,0.25,0.25)),
        py2cc       = list(value =c(0.333,0.334,0.333)),
        vlY1        = list(value =c(1,2,3,4)),
        vlY2        = list(value =c(1,2,3)))
      return(vList)
    }
    
    vlis2 <- function(nsim=100) {
      vList <- simsalapar::varlist(
        n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
        n_1         = list(value =1000),
        coeff       = list(value =c(0.25,0.5,0.75,1),type="grid"),
        p1cc        = list(value =c(0.5,0.5)),
        p2cc        = list(value =c(0.3,0.3,0.4)),
        qual        = list(value =4:5),
        quan        = list(value =6),
        choose_dist = list(value ="E"),
        rho_cov     = list(value =0.2),
        rep_MICE    = list(value =5),
        cov_MICE    = list(value =c("polyreg","polyreg","polyreg","polyreg","pmm")),
        calcu       = list(value ="FREQ"),
        typ_simu    = list(value ="ND"),
        r2t         = list(value =0.5),
        py1cc       = list(value =c(0.25,0.25,0.25,0.25)),
        py2cc       = list(value =c(0.333,0.334,0.333)),
        vlY1        = list(value =c(1,2,3,4)),
        vlY2        = list(value =c(1,2,3)))
      return(vList)
    }
    
    vlis3 <- function(nsim=100) {
      vList <- simsalapar::varlist(
        n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
        n_1         = list(value =1000),
        coeff       = list(value =1),
        p1cc        = list(value =c(0.5,0.5)),
        p2cc        = list(value =c(0.3,0.3,0.4)),
        qual        = list(value =4:5),
        quan        = list(value =6),
        choose_dist = list(value ="E"),
        rho_cov     = list(value =0.2),
        rep_MICE    = list(value =5),
        cov_MICE    = list(value =c("polyreg","polyreg","polyreg","polyreg","pmm")),
        calcu       = list(value ="FREQ"),
        typ_simu    = list(value ="ND"),
        r2t         = list(value =c(0.2,0.4,0.6,0.8),type="grid"),
        py1cc       = list(value =c(0.25,0.25,0.25,0.25)),
        py2cc       = list(value =c(0.333,0.334,0.333)),
        vlY1        = list(value =c(1,2,3,4)),
        vlY2        = list(value =c(1,2,3)))
      return(vList)
    }
    
    vlis4 <- function(nsim=100) {
      vList <- simsalapar::varlist(
        n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
        n_1         = list(value =1000),
        coeff       = list(value =1),
        p1cc        = list(value =c(0.5,0.5)),
        p2cc        = list(value =c(0.3,0.3,0.4)),
        qual        = list(value =4:5),
        quan        = list(value =6),
        choose_dist = list(value ="E"),
        rho_cov     = list(value =c(0.2,0.4,0.6,0.8),type="grid"),
        rep_MICE    = list(value =5),
        cov_MICE    = list(value =c("polyreg","polyreg","polyreg","polyreg","pmm")),
        calcu       = list(value ="FREQ"),
        typ_simu    = list(value ="ND"),
        r2t         = list(value =0.5),
        py1cc       = list(value =c(0.25,0.25,0.25,0.25)),
        py2cc       = list(value =c(0.333,0.334,0.333)),
        vlY1        = list(value =c(1,2,3,4)),
        vlY2        = list(value =c(1,2,3)))
      return(vList)
    }
    
    
    
    runSims <- function(vList=vlis1(nsim=100), parallel=TRUE,
                        seedList=NULL, # set_seed(vList=vList)
                        sfile = NULL # .rds save simulation
    ) {
      #if(is.null(seedList)){seedList <- set_seed(vList=vList)}
      if(parallel){
        res <- simsalapar::doMclapply(vList, seed = seedList, doOne = simul_glob, 
                                      sfile = sfile, #"res42_n1000_lapply_LEc.rds"
                                      monitor = interactive())
      } else{
        res <- doLapply(vList, seed = seedList, doOne = simul, 
                        sfile = sfile, #"res42_n1000_lapply_LEc.rds"
                        monitor = interactive()) 
      }
      return(res)
    }
    
    
    runSims2 =runSims(vList=vlis2(nsim=100), parallel=TRUE, seedList=NULL,sfile = NULL)
    runSims3 =runSims(vList=vlis3(nsim=100), parallel=TRUE, seedList=NULL,sfile = NULL)
    runSims4 =runSims(vList=vlis4(nsim=100), parallel=TRUE, seedList=NULL,sfile = NULL)
    
      
    summry1 <- function(
      res,   # output of simsalapar doApply commands, an array of lists
      ipar=5, 
      data_stats,
      caption.nam = list("peraffec_FOT","peraffec_MICE","peraffec_PR","per_affec_COMP1","per_affec_COMP2","per_affec_COMP3",
                         "per_affec_COMP4","per_affec_COMP5","per_affec_COMP6")[[ipar]],
      sub="median, n=100",
      #rv = c("cs", "rho"),
      cv = c("n"), # "n" if more than one simulated; TODO add a dimension to code in this case
      fun=median, 
      digits=2
    ){
      val <- getArray(res)
      err <- getArray(res, "error")
      warn <- getArray(res, "warning")
      time <- getArray(res, "time")
      df <- array2df(val)
      # str(val)
      print(dimnames(val))
      ft.val <- ftable(round(apply(val[ipar,,,,],1:3,fun),digits), row.vars = rv, col.vars = cv)
      countNA <- apply(val[ipar,,,,],1:3, FUN=function(x) sum(is.na(x)))
      ft.NA <- ftable(round(countNA,digits=0), row.vars = rv, col.vars = cv)
      tp.val <- tablePrint(val[ipar,,,,],row.vars=c("cs","rho"),col.vars=c("beta2"),
                           caption.nam=caption.nam,sub=sub)
      ft.time <- ftable(time, row.vars = rv, col.vars = cv) # run times in ms 
      time.df <- array2df(time)
      errorstat <- array2df(err)
      numerrors <- with(errorstat,sum(value))
      warnstat <- array2df(warn)
      numwarn <- with(warnstat,sum(value))
      list(
        dmn=dimnames(val)[-1], # [ipar,,,,]
        ft.val=ft.val, tp.val=tp.val, 
        ft.NA = ft.NA,
        ft.time=ft.time,
        time.df=time.df,
        numerrors=numerrors,numwarn=numwarn
      )
    }
