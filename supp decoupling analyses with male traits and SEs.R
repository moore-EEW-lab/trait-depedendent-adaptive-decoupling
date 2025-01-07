

library(phytools)
library(car)
library(MASS)
library(nlme)
library(MuMIn)
library(phylosignal)
library(geiger)


### load dataset and phylogeny. Subsetted to only those species for which there were at least 3 observations of males. Because of significant sexual dimorphism in body shape, re-run all analyses using just male body shape data. Additionally, use this dataset to test the effect of incorporating standard errors into within-stage analyses. 

m.traits <- read.csv('drags.decoup.m.traits.csv')
m.drags <- read.tree('drags.decoup.m.phylo.tre')



### between stages analyses for body size
m.size00a <- gls(log(adult.size) ~ log(larval.size), data = m.traits, cor = corPagel(value = 1, phy = m.drags, form = ~binom, fixed = TRUE))
m.size00b <- gls(log(adult.size) ~ log(larval.size), data = m.traits, cor = corPagel(value = 1, phy = m.drags, form = ~binom, fixed = FALSE))
m.size00c <- gls(log(adult.size) ~ log(larval.size), data = m.traits, cor = corPagel(value = 0.5, phy = m.drags, form = ~binom, fixed = FALSE))
m.size00d <- gls(log(adult.size) ~ log(larval.size), data = m.traits, cor = corPagel(value = 0, phy = m.drags, form = ~binom, fixed = TRUE))
AICc(m.size00a, m.size00b, m.size00c, m.size00d)

          # df      AICc
# m.size00a  3 -62.98937
# m.size00b  4 -70.08296
# m.size00c  4 -70.08296
# m.size00d  3 -50.75104

summary(m.size00b)

# Generalized least squares fit by REML
  # Model: log(adult.size) ~ log(larval.size) 
  # Data: m.traits 
        # AIC       BIC  logLik
  # -71.01319 -63.69863 39.5066

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.9209464 

# Coefficients:
                     # Value  Std.Error   t-value p-value
# (Intercept)      1.0644329 0.31270659  3.403935  0.0014
# log(larval.size) 0.9033262 0.09024763 10.009418  0.0000

 # Correlation: 
                 # (Intr)
# log(larval.size) -0.984

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -1.8691220 -0.4164479  0.1580847  0.7375552  1.4564202 

# Residual standard error: 0.1521538 
# Degrees of freedom: 48 total; 46 residual

confint(m.size00b)
                     # 2.5 %   97.5 %
# (Intercept)      0.4515392 1.677327
# log(larval.size) 0.7264441 1.080208

## between stage analysees for body shape -- here looking at relationship between body shape of adult males and larval body shape
m.shape00a <- gls(log(m.elong) ~ log(larval.elong), data = m.traits, cor = corPagel(value = 1, phy = m.drags, form = ~binom, fixed = TRUE))
m.shape00b <- gls(log(m.elong) ~ log(larval.elong), data = m.traits, cor = corPagel(value = 1, phy = m.drags, form = ~binom, fixed = FALSE))
m.shape00c <- gls(log(m.elong) ~ log(larval.elong), data = m.traits, cor = corPagel(value = 0.5, phy = m.drags, form = ~binom, fixed = FALSE))
m.shape00d <- gls(log(m.elong) ~ log(larval.elong), data = m.traits, cor = corPagel(value = 0, phy = m.drags, form = ~binom, fixed = TRUE))
AICc(m.shape00a, m.shape00b, m.shape00c, m.shape00d)

           # df      AICc
# m.shape00a  3 -128.7038
# m.shape00b  4 -126.4332
# m.shape00c  4 -126.4332
# m.shape00d  3  -73.7949

summary(m.shape00a)

# Generalized least squares fit by REML
  # Model: log(m.elong) ~ log(larval.elong) 
  # Data: m.traits 
        # AIC       BIC   logLik
  # -129.2493 -123.7633 67.62463

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 1 

# Coefficients:
                       # Value  Std.Error   t-value p-value
# (Intercept)        0.7045980 0.06666575 10.569114  0.0000
# log(larval.elong) -0.0633214 0.12384008 -0.511316  0.6116

 # Correlation: 
                  # (Intr)
# log(larval.elong) -0.832

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -2.7589020 -1.6012638 -0.9087847 -0.3084145  1.5659440 

# Residual standard error: 0.09849903 
# Degrees of freedom: 48 total; 46 residual

confint(m.shape00a)
                       # 2.5 %    97.5 %
# (Intercept)        0.5739355 0.8352604
# log(larval.elong) -0.3060435 0.1794007

### as in main text, body shape evolves much more independently between stages than does body size.


### within-stage analyses. First just w/o SEs to compare to main text

## create vectors to give to fitContinous()
m.log.as <- log(m.traits$adult.size)
names(m.log.as) <- m.traits$binom
m.log.ls <- log(m.traits$larval.size)
names(m.log.ls) <- m.traits$binom

m.log.m.elong <- log(m.traits$m.elong)
names(m.log.m.elong) <- m.traits$binom
m.log.l.elong <- log(m.traits$larval.elong)
names(m.log.l.elong) <- m.traits$binom




### compare Blomberg's K between traits within each life cycle stage
m.a.size.kappa <- fitContinuous(m.drags, dat = m.log.as, model = 'kappa', control = list(hessian = TRUE))
m.m.shape.kappa <- fitContinuous(m.drags, dat = m.log.m.elong, model = 'kappa', control = list(hessian = TRUE))
m.l.size.kappa <- fitContinuous(m.drags, dat = m.log.ls, model = 'kappa', control = list(hessian = TRUE))
m.l.shape.kappa <- fitContinuous(m.drags, dat = m.log.l.elong, model = 'kappa', control = list(hessian = TRUE))


m.a.size.kappa$opt$kappa #0.8200123
m.a.size.kappa$opt$CI
       # kappa        sigsq
# lb 0.5390914 0.0002124645
# ub 1.2473213 0.0028605727

m.m.shape.kappa$opt$kappa # 0.7636811
m.m.shape.kappa$opt$CI

       # kappa        sigsq
# lb 0.4437389 1.953698e-05
# ub 1.3143063 4.207565e-04


m.l.size.kappa$opt$kappa # 0.8214607
m.l.size.kappa$opt$CI

      # kappa        sigsq
# lb 0.479354 0.0001153442
# ub 1.407723 0.0030716520


m.l.shape.kappa$opt$kappa # 0.9699156
m.l.shape.kappa$opt$CI
       # kappa        sigsq
# lb 0.7134961 1.955839e-05
# ub 1.3184882 1.981309e-04



### similar phylogenetic signal between traits within each life cycle stage. Similar to main text



#### compare BM rates of evo within each life-cycle stage

## adult traits
m.l.a.size.rate <- fitContinuous(m.drags, dat = m.log.as, model = 'BM', control = list(hessian = TRUE))
m.l.a.shape.rate <- fitContinuous(m.drags, dat = m.log.m.elong, model = 'BM', control = list(hessian = TRUE))



m.l.a.size.rate$opt$sigsq #0.0004150072
m.l.a.size.rate$opt$CI
          # sigsq
# lb 0.0002781666
# ub 0.0006191651

m.l.a.shape.rate$opt$sigsq #3.960163e-05
m.l.a.shape.rate$opt$CI
          # sigsq
# lb 2.654375e-05
# ub 5.908318e-05


## larval traits
m.l.l.size.rate <- fitContinuous(m.drags, dat = m.log.ls, model = 'BM', control = list(hessian = TRUE))
m.l.l.shape.rate <- fitContinuous(m.drags, dat = m.log.l.elong, model = 'BM', control = list(hessian = TRUE))

m.l.l.size.rate$opt$sigsq # 0.000315575
m.l.l.size.rate$opt$CI
          # sigsq
# lb 0.0002115202
# ub 0.0004708184

m.l.l.shape.rate$opt$sigsq # 5.581764e-05
m.l.l.shape.rate$opt$CI
          # sigsq
# lb 3.741285e-05
# ub 8.327646e-05

## looks like mean evo rates for body size are faster than those for body shape within both life-cycle stages. 

## next, statistically compare the mean rates using the likelihood ratio test advocated by Adams 2013

#### load in script for the custom function developed by Adams for the likelihood ratio test

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

## end function

### The scripts developed by Adams need to be given the data in matrix format, so create a matrix of trait values for each life-cycle stage
adult.columns <- data.frame(m.log.as, m.log.m.elong)
adult.traits <- as.matrix(adult.columns)
larv.columns <- data.frame(m.log.ls, m.log.l.elong)
larv.traits <- as.matrix(larv.columns)



### compare rates in adult stage
m.a.rate.test <- CompareRates.multTrait(m.drags, adult.traits, TraitCov = T, ms.cov = NULL)
m.a.rate.test

# $Robs
                  # m.log.as m.log.m.elong
# m.log.as      0.0004150072  1.479490e-05
# m.log.m.elong 0.0000147949  3.960163e-05

# $Rconstrained
             # [,1]         [,2]
# [1,] 2.269775e-04 1.441069e-05
# [2,] 1.441069e-05 2.269775e-04

# $Lobs
# [1] 86.05395

# $Lconstrained
# [1] 58.34401

# $LRTest
# [1] 55.41989

# $Prob
# [1] 9.734653e-14

# $AICc.obs
# [1] -164.1079

# $AICc.constrained
# [1] -110.688

# $optimmessage
# [1] "Optimization has converged."

## compare rates in larval stage
m.l.rate.test <- CompareRates.multTrait(m.drags, larv.traits, TraitCov = T, ms.cov = NULL)
m.l.rate.test

# $Robs
                   # m.log.ls m.log.l.elong
# m.log.ls       3.155750e-04 -4.630664e-05
# m.log.l.elong -4.630664e-05  5.581764e-05

# $Rconstrained
              # [,1]          [,2]
# [1,]  1.852477e-04 -4.600898e-05
# [2,] -4.600898e-05  1.852477e-04

# $Lobs
# [1] 87.18382

# $Lconstrained
# [1] 69.48745

# $LRTest
# [1] 35.39275

# $Prob
# [1] 2.694881e-09

# $AICc.obs
# [1] -166.3676

# $AICc.constrained
# [1] -132.9749

# $optimmessage
# [1] "Optimization has converged."

## body shape evolves SLOWER than body size within both stages. All magnitudes are similar to analyses in main text


########
######## now let's do within stage analyses on this dataset but directly incorporating standard errors into the analyses.
########



### first, calculate SEs for each trait on log-scale for within-stage analyses
m.log.as.se <- (sqrt(log((m.traits$adult.size.sd/(m.log.as^2)) + 1)))/sqrt(15)
names(m.log.as.se) <- m.traits$binom

m.log.ls.se <- (sqrt(log((m.traits$larval.size.sd/(m.log.ls^2)) + 1)))/sqrt(15)
names(m.log.ls.se) <- m.traits$binom

m.log.m.elong.se <- (sqrt(log((as.numeric(m.traits$m.elong.sd)/(m.log.m.elong^2)) + 1)))/sqrt(m.traits$m.elong.n)
names(m.log.m.elong.se) <- m.traits$binom

m.log.l.elong.se <- (sqrt(log((m.traits$larval.elong.sd/(m.log.l.elong^2)) + 1)))/sqrt(15)
names(m.log.l.elong.se) <- m.traits$binom


## kappas
adult.size.kappa <- fitContinuous(m.drags, dat = m.log.as, SE = m.log.as.se, model = 'kappa', control = list(hessian = TRUE))
adult.shape.kappa <- fitContinuous(m.drags, dat = m.log.m.elong, SE = m.log.m.elong.se, model = 'kappa', control = list(hessian = TRUE))
larval.size.kappa <- fitContinuous(m.drags, dat = m.log.ls, SE = m.log.ls.se, model = 'kappa', control = list(hessian = TRUE))
larval.shape.kappa <- fitContinuous(m.drags, dat = m.log.l.elong, SE = m.log.l.elong.se, model = 'kappa', control = list(hessian = TRUE), bounds = list(kappa=c(0.0001,10)))


adult.size.kappa$opt$kappa #0.8196959
adult.size.kappa$opt$CI
       # kappa        sigsq
# lb 0.5375141 0.0002100688
# ub 1.2500162 0.0028712284

adult.shape.kappa$opt$kappa # 0.4344986
adult.shape.kappa$opt$CI

        # kappa        sigsq
# lb 0.04793287 7.280523e-06
# ub 3.93861426 4.672899e-03


larval.size.kappa$opt$kappa #0.8259835
larval.size.kappa$opt$CI
       # kappa        sigsq
# lb 0.4778359 0.0001081339
# ub 1.4277888 0.0031208617

larval.shape.kappa$opt$kappa # 3.877565
larval.shape.kappa$opt$CI
      # kappa        sigsq
# lb 1.740298 7.190601e-18
# ub 8.639617 2.198323e-04


## kappa for larval shape is very high, but doesn't change what we conclude from these analyses

#### compare BM rates of evo

## adult traits
m.l.a.size.rate <- fitContinuous(m.drags, dat = m.log.as, SE = m.log.as.se, model = 'BM', control = list(hessian = TRUE))
m.l.a.shape.rate <- fitContinuous(m.drags, dat = m.log.m.elong, SE = m.log.m.elong.se, model = 'BM', control = list(hessian = TRUE))

m.l.a.size.rate$opt$sigsq #0.0004126531
m.l.a.size.rate$opt$CI
          # sigsq
# lb 0.0002761222
# ub 0.0006166928

m.l.a.shape.rate$opt$sigsq # 2.464137e-05
m.l.a.shape.rate$opt$CI
          # sigsq
# lb 1.095111e-05
# ub 5.544616e-05


## larval traits
m.l.l.size.rate <- fitContinuous(m.drags, dat = m.log.ls, SE = m.log.ls.se, model = 'BM', control = list(hessian = TRUE))
m.l.l.shape.rate <- fitContinuous(m.drags, dat = m.log.l.elong, SE = m.log.l.elong.se, model = 'BM', control = list(hessian = TRUE))

m.l.l.size.rate$opt$sigsq # 0.0003117696
m.l.l.size.rate$opt$CI
          # sigsq
# lb 0.0002079590
# ub 0.0004674013

m.l.l.shape.rate$opt$sigsq # 3.517519e-05
m.l.l.shape.rate$opt$CI
          # sigsq
# lb 2.063176e-05
# ub 5.997035e-05



#### statistically compare evo rates within stages using the Adams 2013 test

## generate matricies for the standard errors
a.ses.frame <- data.frame(m.log.as.se, m.log.m.elong.se)
a.ses <- as.matrix(a.ses.frame)
l.ses.frame <- data.frame(m.log.ls.se, m.log.l.elong.se)
l.ses <- as.matrix(l.ses.frame)

## compare rates in adult stage
m.se.a.rate.test <- CompareRates.multTrait(m.drags, adult.traits, TraitCov = T, ms.err = a.ses, ms.cov = NULL)
m.se.a.rate.test

# $Robs
                  # m.log.as m.log.m.elong
# m.log.as      0.0004150072  1.479490e-05
# m.log.m.elong 0.0000147949  3.960163e-05

# $Rconstrained
             # [,1]         [,2]
# [1,] 1.914292e-04 3.360059e-05
# [2,] 3.360059e-05 1.914292e-04

# $Lobs
# [1] 28.8025

# $Lconstrained
# [1] 25.10164

# $LRTest
# [1] 7.401709

# $Prob
# [1] 0.006516196

# $AICc.obs
# [1] -49.60499

# $AICc.constrained
# [1] -44.20329

# $optimmessage
# [1] "Optimization has converged."

## compare rates in larval stage
m.se.l.rate.test <- CompareRates.multTrait(m.drags, larv.traits, TraitCov = T, ms.err = l.ses, ms.cov = NULL)
m.se.l.rate.test

# $Robs
                   # m.log.ls m.log.l.elong
# m.log.ls       3.155750e-04 -4.630664e-05
# m.log.l.elong -4.630664e-05  5.581764e-05

# $Rconstrained
              # [,1]          [,2]
# [1,]  1.533203e-04 -7.246432e-06
# [2,] -7.246432e-06  1.533203e-04

# $Lobs
# [1] 45.28344

# $Lconstrained
# [1] 39.72557

# $LRTest
# [1] 11.11572

# $Prob
# [1] 0.0008559899

# $AICc.obs
# [1] -82.56687

# $AICc.constrained
# [1] -73.45115

# $optimmessage
# [1] "Optimization has converged."


### in both stages, body size evolves a lot faster than body shape. magnitudes of all rates similar to male-only shape analyses AND to the analyses in the main text.
