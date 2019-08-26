##---------------------------------------------------------------##
## Bias correction for measures of association estimated in      ##
## EHR data using probabilistic phenotypes                       ##
##                                                               ##
## Last updated: 5/1/2019                                        ##
##---------------------------------------------------------------##

library(gee)

## Inputs: 
##  p - probabilistic phenotype (vector)
##  X - covariate of interest (vector)
##  W - confounder variable
##
##  If mean phenotype probability in cases and controls is unknown:
##  pstar - proposed cutpoint for dichtomizing phenotype
##  C - specificity of phenotype dichtomoized at proposed cutpoint
##  S - sensitivity of phenotype dichtomoized at proposed cutpoint
##
##  If mean phenotype probability in cases and controls is known:
##  mu0 - mean phenotype probability in controls
##  mu1 - mean phenotype probability in cases


## Function for bias correction with unknown values for mu0 and mu1
## link can take values "ident", "log", or "logit"
bias.adjust.unk <- function(p,X,W,S,C,pstar, link = "ident"){
  
  # dichotomized phenotype
  hatY <- 1*(p > pstar)
  
  # sample size
  N <- length(hatY)
  
  # regress probabilistic phenotype on predictors
  fitp = gee(p~X + W,id = seq(1,N), fam = "gaussian")
  
  # estimate proportion of controls misclassified (rank0cut) and cases misclassified (rank1cut)
  trueND = (sum(hatY==0)-C*N)/(1-S-C)
  rank0cut = (1-C)*(N-trueND)
  rank1cut = (1-S)*trueND
  
  # estimate mu0 and mu1, removing upper tail of probabilistic phenotype for controls
  # and lower tail of probabilistic phenotype for cases
  muhat1star = mean(p[hatY == 1][rank(p[hatY==1]) < (sum(hatY==1)-rank1cut)],na.rm = T)
  muhat0star = mean(p[hatY == 0][rank(p[hatY==0]) > rank0cut])
  
  betastar = fitp$coef/(muhat1star - muhat0star)
  
  # estimated prevalence
  p0 = sum(p)/N
  
  if (link == "ident"){
    betastar <- betastar
  } else if (link == "log"){
    betastar = betastar/p0
  } else if (link == "logit"){
    betastar = betastar/(p0*(1-p0))
  } else return("unsupported link function")
  
  # return association parameters (drop intercept)
  return(betastar[-1])
}

## Function for bias correction with known values for mu0 and mu1
## link can take values "ident", "log", or "logit"
bias.adjust.known <- function(p,X,W,mu0,mu1, link = "ident"){
  
  # sample size
  N <- length(p)

  # regress probabilistic phenotype on predictors
  fitp = gee(p~X + W,id = seq(1,N), fam = "gaussian")
  
  # make bias correctioon
  betastar = fitp$coef/(mu1 - mu0)
  
  # estimated prevalence
  p0 = sum(p)/N
  
  if (link == "ident"){
    betastar = betastar
  } else if (link == "log"){
    betastar = betastar/p0
  } else if (link == "logit"){
    betastar <- betastar/(p0*(1-p0))
  } else return("unsupported link function")
  
  # return association parameters (drop intercept)
  return(betastar[-1])
}


