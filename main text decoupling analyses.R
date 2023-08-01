

library(phytools)
library(car)
library(MASS)
library(nlme)
library(MuMIn)
library(phylosignal)
library(geiger)

traits <- read.csv('drags.decoup.traits.csv')
na.drags <- read.tree('drags.decoup.phylo.tre')

rownames(traits) <- traits$binom

#### between-stage analyses!


### body size first

size00a <- gls(log(adult.size) ~ log(larval.size), data = traits, cor = corPagel(value = 1, phy = na.drags, form = ~binom, fixed = TRUE))
size00b <- gls(log(adult.size) ~ log(larval.size), data = traits, cor = corPagel(value = 1, phy = na.drags, form = ~binom, fixed = FALSE))
size00c <- gls(log(adult.size) ~ log(larval.size), data = traits, cor = corPagel(value = 0.5, phy = na.drags, form = ~binom, fixed = FALSE))
size00d <- gls(log(adult.size) ~ log(larval.size), data = traits, cor = corPagel(value = 0, phy = na.drags, form = ~binom, fixed = TRUE))
AICc(size00a, size00b, size00c, size00d)

        # df      AICc
# size00a  3 -128.2837
# size00b  4 -137.3412
# size00c  4 -137.3412
# size00d  3  -89.8952

summary(size00b)

# Generalized least squares fit by REML
  # Model: log(adult.size) ~ log(larval.size) 
  # Data: traits 
        # AIC       BIC   logLik
  # -137.8676 -128.3898 72.93378

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.9436966 

# Coefficients:
                     # Value  Std.Error   t-value p-value
# (Intercept)      1.0114908 0.21327907  4.742569       0
# log(larval.size) 0.9151839 0.06149487 14.882279       0

 # Correlation: 
                 # (Intr)
# log(larval.size) -0.967

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -1.7969749 -0.5005463  0.1098194  0.7115810  2.6612987 

# Residual standard error: 0.1552036 
# Degrees of freedom: 81 total; 79 residual

confint(size00b)
                     # 2.5 %   97.5 %
# (Intercept)      0.5934715 1.429510
# log(larval.size) 0.7946561 1.035712

log.as <- log(traits$adult.size)
names(log.as) <- traits$binom
log.ls <- log(traits$larval.size)
names(log.ls) <- traits$binom

pic.a.size <- pic(log.as, phy = na.drags)
pic.l.size <- pic(log.ls, phy = na.drags)
size.pic <- lm(pic.a.size ~ pic.l.size -1)
summary(size.pic)

# Call:
# lm(formula = pic.a.size ~ pic.l.size - 1)

# Residuals:
      # Min        1Q    Median        3Q       Max 
# -0.024456 -0.007497  0.002381  0.007743  0.046269 

# Coefficients:
           # Estimate Std. Error t value Pr(>|t|)    
# pic.l.size  0.91180    0.06792   13.43   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.01254 on 79 degrees of freedom
# Multiple R-squared:  0.6952,	Adjusted R-squared:  0.6914 
# F-statistic: 180.2 on 1 and 79 DF,  p-value: < 2.2e-16

# size Rpic
sqrt(summary(size.pic)$r.squared) # 0.8338101

# essentially no decoupling of body size


## shape
shape00a <- gls(log(adult.elong) ~ log(larval.elong), data = traits, cor = corPagel(value = 1, phy = na.drags, form = ~binom, fixed = TRUE))
shape00b <- gls(log(adult.elong) ~ log(larval.elong), data = traits, cor = corPagel(value = 1, phy = na.drags, form = ~binom, fixed = FALSE))
shape00c <- gls(log(adult.elong) ~ log(larval.elong), data = traits, cor = corPagel(value = 0.5, phy = na.drags, form = ~binom, fixed = FALSE))
shape00d <- gls(log(adult.elong) ~ log(larval.elong), data = traits, cor = corPagel(value = 0, phy = na.drags, form = ~binom, fixed = TRUE))
AICc(shape00a, shape00b, shape00c, shape00d)

         # df      AICc
# shape00a  3 -187.1485
# shape00b  4 -206.4948
# shape00c  4 -206.4948
# shape00d  3 -136.9914

summary(shape00b)

# Generalized least squares fit by REML
  # Model: log(adult.elong) ~ log(larval.elong) 
  # Data: traits 
        # AIC       BIC   logLik
  # -207.0211 -197.5433 107.5106

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.9074821 

# Coefficients:
                       # Value Std.Error   t-value p-value
# (Intercept)        0.7109994 0.0660831 10.759171  0.0000
# log(larval.elong) -0.1477301 0.1225311 -1.205654  0.2315

 # Correlation: 
                  # (Intr)
# log(larval.elong) -0.867

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -2.5238789 -1.3957499 -0.5110539  0.1694671  1.9837989 

# Residual standard error: 0.09599704 
# Degrees of freedom: 81 total; 79 residual

confint(shape00b)
                       # 2.5 %     97.5 %
# (Intercept)        0.5814789 0.84051986
# log(larval.elong) -0.3878865 0.09242633

#shape pics
log.a.elong <- log(traits$adult.elong)
names(log.a.elong) <- traits$binom

log.l.elong <- log(traits$larval.elong)
names(log.l.elong) <- traits$binom

pic.a.elong <- pic(log.a.elong, phy = na.drags)
pic.l.elong <- pic(log.l.elong, phy = na.drags)
elong01 <- lm(pic.a.elong ~ pic.l.elong - 1)
summary(elong01)

# Call:
# lm(formula = pic.a.elong ~ pic.l.elong - 1)

# Residuals:
      # Min        1Q    Median        3Q       Max 
# -0.018532 -0.003720 -0.000451  0.005978  0.036670 

# Coefficients:
            # Estimate Std. Error t value Pr(>|t|)
# pic.l.elong  -0.1529     0.1318   -1.16     0.25

# Residual standard error: 0.008753 on 79 degrees of freedom
# Multiple R-squared:  0.01675,	Adjusted R-squared:  0.004304 
# F-statistic: 1.346 on 1 and 79 DF,  p-value: 0.2495

#shape Rpic
sqrt(summary(elong01)$r.squared) # 0.1294207

### body shape evolves quite independently between life-cycle stages, but body size does not. Trait-dependent decoupling


##### Within-stage analyses


# phylogenetic signal
a.size.kappa <- fitContinuous(na.drags, dat = log.as, model = 'kappa', control = list(hessian = TRUE))
a.shape.kappa <- fitContinuous(na.drags, dat = log.a.elong, model = 'kappa', control = list(hessian = TRUE))
l.size.kappa <- fitContinuous(na.drags, dat = log.ls, model = 'kappa', control = list(hessian = TRUE))
l.shape.kappa <- fitContinuous(na.drags, dat = log.l.elong, model = 'kappa', control = list(hessian = TRUE))


l.shape.kappa$opt$kappa # 0.8240745
l.shape.kappa$opt$CI
       # kappa        sigsq
# lb 0.6278404 4.243483e-05
# ub 1.0816425 2.192245e-04

l.size.kappa$opt$kappa # 0.9991374
l.size.kappa$opt$CI
       # kappa        sigsq
# lb 0.7312208 0.0001383856
# ub 1.3652176 0.0012877568


a.shape.kappa$opt$kappa # 0.396754
a.shape.kappa$opt$CI
       # kappa        sigsq
# lb 0.1944265 0.0001716758
# ub 0.8096309 0.0011844638

a.size.kappa$opt$kappa #0.9560179
a.size.kappa$opt$CI
       # kappa        sigsq
# lb 0.7203309 0.0002197911
# ub 1.2688198 0.0015555510


## body shape a little lower phylogenetic signal in larval stage and quite a bit lower phylogenetic signal in adult stage, but confidence intervals overlap a lot



#### compare BM rates of evo (sensu Ackerly 2009)

## first estimate the mean BM rate +/- 95% CIs

## larval traits
l.l.size.rate <- fitContinuous(na.drags, dat = log.ls, model = 'BM', control = list(hessian = TRUE))
l.l.shape.rate <- fitContinuous(na.drags, dat = log.l.elong, model = 'BM', control = list(hessian = TRUE))

l.l.shape.rate$opt$sigsq # 5.444421e-05
l.l.shape.rate$opt$CI
          # sigsq
# lb 4.001274e-05
# ub 7.408070e-05

l.l.size.rate$opt$sigsq # 0.0004208978
l.l.size.rate$opt$CI
          # sigsq
# lb 0.0003093309
# ub 0.0005727038



## adult traits
l.a.size.rate <- fitContinuous(na.drags, dat = log.as, model = 'BM', control = list(hessian = TRUE))
l.a.shape.rate <- fitContinuous(na.drags, dat = log.a.elong, model = 'BM', control = list(hessian = TRUE))

l.a.shape.rate$opt$sigsq #7.599301e-05
l.a.shape.rate$opt$CI
          # sigsq
# lb 5.584962e-05
# ub 1.034015e-04

l.a.size.rate$opt$sigsq #0.0005033137
l.a.size.rate$opt$CI
          # sigsq
# lb 0.0003699009
# ub 0.0006848448

## point estimates for BM evo rates are much faster for body size than body shape in both stages.

#### let's compare those estimated rates statistically using the LRT approach advocated in Adams 2013


## first need to create the functions using the custom scripts from Adams 2013

#rm(list=ls())

CompareRates.multTrait <- function(phy, x, TraitCov=T, ms.err=NULL, ms.cov=NULL){
  #Compares LLik of R-matrix vs. LLik of R-matrix with constrained diagonal
  
  #TraitCov = TRUE assumes covariation among traits (default)
  #ms.err allows the incorporation of within-species measurement error. Input is a matrix of species (rows) by within-species variation for each trait (columns).
  #ms.cov allows the incorporation of within-species covariation between traits. Input is a matrix of species (rows) by within-species covariation for each pair of traits (columns). These must be provided in a specific order, beginning with covariation between trait 1 and the rest, then trait 2 and the rest, etc. For instance, for 4 traits, the columns are: cov_12, cov_13, cov_14, cov_23, cov_24 cov_34.
  
  #Some calculations adapted from 'evol.vcv' in phytools (Revell, 2012)
  
  library(MASS)
  x<-as.matrix(x)
  N<-nrow(x)
  p<-ncol(x)
  C<-vcv.phylo(phy)
  C<-C[rownames(x),rownames(x)]
  if (is.matrix(ms.err)){    
    ms.err<-as.matrix(ms.err[rownames(x),])}
  if (is.matrix(ms.cov)){    
    ms.cov<-as.matrix(ms.cov[rownames(x),])}
  
  #Cholesky decomposition function for diagonal-constrained VCV
  build.chol<-function(b){
    c.mat<-matrix(0,nrow=p,ncol=p)
    c.mat[lower.tri(c.mat)] <- b[-1]  
    c.mat[p,p]<-exp(b[1])
    c.mat[1,1]<-sqrt(sum((c.mat[p,])^2))
    if(p>2){
      for (i in 2:(p-1)){
        c.mat[i,i]<-ifelse( (c.mat[1,1]^2-sum((c.mat[i,])^2) )>0,
                            sqrt(c.mat[1,1]^2-sum((c.mat[i,])^2)), 0)
      }}
    return(c.mat) 
  }
  
  #Fit Rate matrix for all traits: follows code of L. Revell (evol.vcv)
  a.obs<-colSums(solve(C))%*%x/sum(solve(C))   
  D<-matrix(0,N*p,p)
  for(i in 1:(N*p)) for(j in 1:p) if((j-1)*N<i&&i<=j*N) D[i,j]=1.0
  y<-as.matrix(as.vector(x))
  one<-matrix(1,N,1)
  R.obs<-t(x-one%*%a.obs)%*%solve(C)%*%(x-one%*%a.obs)/N
  if (TraitCov==F)    #for TraitCov = F
  { R.obs<-diag(diag(R.obs),p)  }
  #Calculate observed likelihood with or without measurement error
  LLik.obs<-ifelse(is.matrix(ms.err)==TRUE, 
                   -t(y-D%*%t(a.obs))%*%ginv((kronecker(R.obs,C)+ diag(as.vector(ms.err))))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                     determinant((kronecker(R.obs,C)+ diag(as.vector(ms.err))))$modulus[1]/2 , 
                   -t(y-D%*%t(a.obs))%*%ginv(kronecker(R.obs,C))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                     determinant(kronecker(R.obs,C))$modulus[1]/2
  ) 
  
  #Fit common rate for all traits; search over parameter space   
  sigma.mn<-0   #reasonable start value for diagonal
  
  #Within-species measurement error matrix
  if(is.matrix(ms.err)){m.e<-diag(as.vector(ms.err))}
  
  #Within-species measurement error and trait covariation matrix
  if (is.matrix(ms.err) && is.matrix(ms.cov)){	
    within.spp<-cbind(ms.err,ms.cov)
    rc.label<-NULL
    for (i in 1:p){ rc.label<-rbind(rc.label,c(i,i)) }
    for (i in 1:p){
      for (j in 2:p){ if (i!=j && i<j){rc.label<-rbind(rc.label,c(i,j))} }}
    m.e<-NULL
    for (i in 1:p){
      tmp<-NULL
      for (j in 1:p){
        for (k in 1:nrow(rc.label)){
          if(setequal(c(i,j),rc.label[k,])==T) {tmp<-cbind(tmp,diag(within.spp[,k]))}
        }
      }
      m.e<-rbind(m.e,tmp)
    }
  }
  
  #likelihood optimizer for no trait covariation
  lik.covF<-function(sigma){  
    R<-R.obs
    diag(R)<-sigma
    LLik<-ifelse(is.matrix(ms.err)==TRUE, 
                 -t(y-D%*%t(a.obs))%*%ginv((kronecker(R,C)+ m.e))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant((kronecker(R,C)+ m.e))$modulus[1]/2 , 
                 -t(y-D%*%t(a.obs))%*%ginv(kronecker(R,C))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant(kronecker(R,C))$modulus[1]/2
    ) 
    if (LLik == -Inf) { LLikk <- -1e+10  }
    return(-LLik)
  }
  
  #likelihood optimizer with trait covariation
  lik.covT<-function(sigma){  
    low.chol<-build.chol(sigma)
    R<-low.chol%*%t(low.chol)
    
    LLik<-ifelse(is.matrix(ms.err)==TRUE, 
                 -t(y-D%*%t(a.obs))%*%ginv((kronecker(R,C)+ m.e))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant((kronecker(R,C)+ m.e))$modulus[1]/2 , 
                 -t(y-D%*%t(a.obs))%*%ginv(kronecker(R,C))%*%(y-D%*%t(a.obs))/2-N*p*log(2*pi)/2-  
                   determinant(kronecker(R,C))$modulus[1]/2
    ) 
    if (LLik == -Inf)  {LLikk <- -1e+10  }
    return(-LLik)
  }
  
  ##Optimize for no trait covariation
  if (TraitCov==F)    
  { model1<-optim(sigma.mn,fn=lik.covF,method="Nelder-Mead", control = list(maxit = 3000, reltol = 0.000001))}
  ##Optimize with trait covariation
  R.offd<-rep(0,(p*(p-1)/2))
  if (TraitCov==T)  
  {model1<-optim(par=c(sigma.mn,R.offd),fn=lik.covT, method="Nelder-Mead", control = list(maxit = 3000, reltol = 0.000001))}
  
  #### Assemble R.constrained
  if (TraitCov==F){R.constr<-diag(model1$par,p)}
  if (TraitCov==T){  
    chol.mat<-build.chol(model1$par)
    R.constr<-chol.mat%*%t(chol.mat)}
  
  if(model1$convergence==0)
    message<-"Optimization has converged."
  else
    message<-"Optim may not have converged.  Consider changing start value or lower/upper limits."
  LRT<- (-2*((-model1$value-LLik.obs)))
  LRT.prob<-pchisq(LRT, (p-1),lower.tail=FALSE) #df = Nvar-1
  AIC.obs<- -2*LLik.obs+2*p+2*p #(2p twice: 1x for rates, 1x for anc. states)
  AIC.common<- -2*(-model1$value)+2+2*p #(2*1: for 1 rate 2p for anc. states)
  return(list(Robs=R.obs, Rconstrained=R.constr,Lobs=LLik.obs,Lconstrained=(-model1$value),LRTest=LRT,Prob=LRT.prob,
              AICc.obs=AIC.obs,AICc.constrained=AIC.common,optimmessage=message))   
 }

Rate.Correlation <- function(Tree, Data, Trait1, Trait2) {
 DataM<-as.matrix(Data)

  Rates<-CompareRates.multTrait(Tree, DataM[,c(which(colnames(DataM)==Trait1), which(colnames(DataM)==Trait2))])
  Rates

}


Trait.Correlation <- function(Tree, Data, Trait1, Trait2) {
T1 <- Data[,Trait1]
names(T1) <- rownames(Data)
T1.pic <- pic(T1, Tree)

T2 <- Data[,Trait2]
names(T2) <- rownames(Data)
T2.pic<-pic(T2, Tree)

lm.pic <- lm(T2.pic ~ T1.pic -1) 

plot(T1.pic, T2.pic, pch=19)
abline(lm.pic, col="red")

#summary(lm.pic)
Pvalue.pic <- anova(lm.pic)$'Pr(>F)'[1]
Pvalue.pic
}

### end of function

#### now directly compare the rates of evolution between traits within each stage

## the function takes data in matrix form, so the first step is to set up separate matricies for each life-cycle stage.
adult.columns <- data.frame(log.as, log.a.elong)
ac <- adult.columns[order(match(rownames(adult.columns), na.drags$tip.label)), ]
adult.traits <- as.matrix(ac)

larval.columns <- data.frame(log.ls, log.l.elong)
lc <- larval.columns[order(match(rownames(larval.columns), na.drags$tip.label)), ]
larval.traits <- as.matrix(lc)


## compare the rates using the function. Be sure to set the argument 'TraitCov = T' so that the model explicitly accounts for the covariance between the traits

# test within adult stage
a.rate.test <- CompareRates.multTrait(na.drags, adult.traits, TraitCov = T)
a.rate.test

# $Robs
                  # log.as  log.a.elong
# log.as      5.033137e-04 7.480833e-05
# log.a.elong 7.480833e-05 7.599301e-05

# $Rconstrained
             # [,1]         [,2]
# [1,] 2.907150e-04 7.497933e-05
# [2,] 7.497933e-05 2.907150e-04

# $Lobs
# [1] 128.4177

# $Lconstrained
# [1] 92.99277

# $LRTest
# [1] 70.8499

# $Prob
# [1] 3.854627e-17

# $AICc.obs
# [1] -248.8354

# $AICc.constrained
# [1] -179.9855

# $optimmessage
# [1] "Optimization has converged."

5.033137e-04/7.599301e-05 # 6.623158
7.599301e-05/5.033137e-04 # 0.1509854

### body size evolves 6.62 times faster than body shape in the adult stage. body shape evolves 15.1% as fast in the adult stage as body size



### test within larval stage
l.rate.test <- CompareRates.multTrait(na.drags, larval.traits, TraitCov = T)
l.rate.test

# $Robs
                   # log.ls   log.l.elong
# log.ls       4.208978e-04 -6.272469e-05
# log.l.elong -6.272469e-05  5.444421e-05

# $Rconstrained
              # [,1]          [,2]
# [1,]  2.380139e-04 -6.257582e-05
# [2,] -6.257582e-05  2.380139e-04

# $Lobs
# [1] 150.3876

# $Lconstrained
# [1] 109.1429

# $LRTest
# [1] 82.48934

# $Prob
# [1] 1.062415e-19

# $AICc.obs
# [1] -292.7752

# $AICc.constrained
# [1] -212.2858

# $optimmessage
# [1] "Optimization has converged."

4.208978e-04/5.444421e-05 # 7.730809
5.444421e-05/4.208978e-04 # 0.1293526

### body size evolves 7.73 times faster than body shape in the larval stage. Body shape evolves only 12.9% as fast in the larval stage as body size

### overall, these analyses indicate that some traits evolve more independently than others between stages, but the decoupling doesn't guarantee greater diversification