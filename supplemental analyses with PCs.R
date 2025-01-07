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

# create columns for larval and adult width
traits$a.width <- traits$adult.size/traits$adult.elong
traits$l.width <- traits$larval.size/traits$larval.elong

# principle components for each life stage
a.pcs <- princomp(traits[,c('adult.size', 'a.width')], cor = TRUE, scores = TRUE)
l.pcs <- princomp(traits[,c('larval.size', 'l.width')], cor = TRUE, scores = TRUE)


## add pc scores for each species. PC1 is size, PC2 is shape in both stages
a.pcs$scores[,1]
traits$a.PC1 <- a.pcs$scores[,1]
traits$a.PC2 <- a.pcs$scores[,2]
traits$l.PC1 <- l.pcs$scores[,1]
traits$l.PC2 <- l.pcs$scores[,2]

# size correlation analysis
size00a <- gls(a.PC1 ~ l.PC1, data = traits, cor = corPagel(value = 1, phy = na.drags, form = ~binom, fixed = TRUE))
size00b <- gls(a.PC1 ~ l.PC1, data = traits, cor = corPagel(value = 1, phy = na.drags, form = ~binom, fixed = FALSE))
size00c <- gls(a.PC1 ~ l.PC1, data = traits, cor = corPagel(value = 0.5, phy = na.drags, form = ~binom, fixed = FALSE))
size00d <- gls(a.PC1 ~ l.PC1, data = traits, cor = corPagel(value = 0, phy = na.drags, form = ~binom, fixed = TRUE))
AICc(size00a, size00b, size00c, size00d)

        # df     AICc
# size00a  3 125.1117
# size00b  4 110.3129
# size00c  4 110.3129
# size00d  3 155.3830

summary(size00b)
# Generalized least squares fit by REML
  # Model: a.PC1 ~ l.PC1 
  # Data: traits 
       # AIC      BIC    logLik
  # 109.7866 119.2644 -50.89331

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.8995762 

# Coefficients:
                 # Value  Std.Error   t-value p-value
# (Intercept) -0.0282352 0.24383730 -0.115795  0.9081
# l.PC1        0.9368362 0.07111349 13.173818  0.0000

 # Correlation: 
      # (Intr)
# l.PC1 -0.281

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -2.09025934 -0.61845210 -0.02741702  0.74911280  2.20089964 

# Residual standard error: 0.6832164 
# Degrees of freedom: 81 total; 79 residual

confint(size00b)

                 # 2.5 %    97.5 %
# (Intercept) -0.5061475 0.4496772
# l.PC1        0.7974563 1.0762161



#### shape correlation analysis
shape00a <- gls(a.PC2 ~ l.PC2, data = traits, cor = corPagel(value = 1, phy = na.drags, form = ~binom, fixed = TRUE))
shape00b <- gls(a.PC2 ~ l.PC2, data = traits, cor = corPagel(value = 1, phy = na.drags, form = ~binom, fixed = FALSE))
shape00c <- gls(a.PC2 ~ l.PC2, data = traits, cor = corPagel(value = 0.5, phy = na.drags, form = ~binom, fixed = FALSE))
shape00d <- gls(a.PC2 ~ l.PC2, data = traits, cor = corPagel(value = 0, phy = na.drags, form = ~binom, fixed = TRUE))
AICc(shape00a, shape00b, shape00c, shape00d)

         # df       AICc
# shape00a  3 -63.296692
# shape00b  4 -74.163852
# shape00c  4 -74.163852
# shape00d  3   5.813607

summary(shape00b)
# Generalized least squares fit by REML
  # Model: a.PC2 ~ l.PC2 
  # Data: traits 
        # AIC       BIC   logLik
  # -74.69017 -65.21238 41.34508

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
   # lambda 
# 0.9365436 

# Coefficients:
                 # Value  Std.Error   t-value p-value
# (Intercept)  0.1378260 0.07971158  1.729059  0.0877
# l.PC2       -0.1572952 0.11587400 -1.357468  0.1785

 # Correlation: 
      # (Intr)
# l.PC2 0.052 

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -2.8400424 -1.4068716 -0.5400489  0.1222823  1.9948472 

# Residual standard error: 0.2293136 
# Degrees of freedom: 81 total; 79 residual


confint(shape00b)

                  # 2.5 %     97.5 %
# (Intercept) -0.01840584 0.29405781
# l.PC2       -0.38440405 0.06981367


### create vectors for PC scores to do diversification analysis
a.PC1 <- a.pcs$scores[,1]
a.PC2 <- a.pcs$scores[,2]
l.PC1 <- l.pcs$scores[,1]
l.PC2 <- l.pcs$scores[,2]

## phylogenetic signal (K) for each trait
a.size.kappa <- fitContinuous(na.drags, dat = a.PC1, model = 'kappa', control = list(hessian = TRUE))
a.shape.kappa <- fitContinuous(na.drags, dat = a.PC2, model = 'kappa', control = list(hessian = TRUE))
l.size.kappa <- fitContinuous(na.drags, dat = l.PC1, model = 'kappa', control = list(hessian = TRUE))
l.shape.kappa <- fitContinuous(na.drags, dat = l.PC2, model = 'kappa', control = list(hessian = TRUE), bounds = list(kappa = c(0, 15)))



l.shape.kappa$opt$kappa # 1.073016
l.shape.kappa$opt$CI

      # kappa        sigsq
# lb 0.873014 9.431785e-05
# ub 1.318838 4.914120e-04


l.size.kappa$opt$kappa # 0.9887713
l.size.kappa$opt$CI
       # kappa      sigsq
# lb 0.7153072 0.00225489
# ub 1.3667816 0.02219767


a.shape.kappa$opt$kappa # 0.4920153
a.shape.kappa$opt$CI
       # kappa        sigsq
# lb 0.2675359 0.0006040679
# ub 0.9048469 0.0047351190

a.size.kappa$opt$kappa #0.8920649
a.size.kappa$opt$CI
       # kappa       sigsq
# lb 0.6520005 0.004816264
# ub 1.2205203 0.035837211

## evolutionary rate (from BM model) for each trait
l.l.size.rate <- fitContinuous(na.drags, dat = l.PC1, model = 'BM', control = list(hessian = TRUE))
l.l.shape.rate <- fitContinuous(na.drags, dat = l.PC2, model = 'BM', control = list(hessian = TRUE))
l.a.size.rate <- fitContinuous(na.drags, dat = a.PC1, model = 'BM', control = list(hessian = TRUE))
l.a.shape.rate <- fitContinuous(na.drags, dat = a.PC2, model = 'BM', control = list(hessian = TRUE))


l.l.shape.rate$opt$sigsq # 0.000278406
l.l.shape.rate$opt$CI
# lb 0.0002046092
# ub 0.0003788192

l.l.size.rate$opt$sigsq # 0.006807064
l.l.size.rate$opt$CI
          # sigsq
# lb 0.005002723
# ub 0.009262180


l.a.shape.rate$opt$sigsq #0.0003560232
l.a.shape.rate$opt$CI
          # sigsq
# lb 0.0002616525
# ub 0.0004844308

l.a.size.rate$opt$sigsq #0.009138721
l.a.size.rate$opt$CI
          # sigsq
# lb 0.00671633
# ub 0.01243480



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
adult.columns <- data.frame(a.PC1, a.PC2)
ac <- adult.columns[order(match(rownames(adult.columns), na.drags$tip.label)), ]
adult.traits <- as.matrix(ac)

larval.columns <- data.frame(l.PC1, l.PC2)
lc <- larval.columns[order(match(rownames(larval.columns), na.drags$tip.label)), ]
larval.traits <- as.matrix(lc)

# test within adult stage
a.rate.test <- CompareRates.multTrait(na.drags, adult.traits, TraitCov = T)
a.rate.test

# $Robs
             # a.PC1        a.PC2
# a.PC1 0.0091387210 0.0000325971
# a.PC2 0.0000325971 0.0003560232

# $Rconstrained
             # [,1]         [,2]
# [1,] 4.748735e-03 4.653786e-05
# [2,] 4.653786e-05 4.748735e-03

# $Lobs
# [1] -57.93422

# $Lconstrained
# [1] -136.3305

# $LRTest
# [1] 156.7925

# $Prob
# [1] 5.681957e-36

# $AICc.obs
# [1] 123.8684

# $AICc.constrained
# [1] 278.6609

# $optimmessage
# [1] "Optimization has converged."

l.rate.test <- CompareRates.multTrait(na.drags, larval.traits, TraitCov = T)
l.rate.test
# $Robs
              # l.PC1         l.PC2
# l.PC1  0.0068070637 -0.0008291978
# l.PC2 -0.0008291978  0.0002784060

# $Rconstrained
              # [,1]          [,2]
# [1,]  0.0035415237 -0.0008244681
# [2,] -0.0008244681  0.0035415237

# $Lobs
# [1] -17.80543

# $Lconstrained
# [1] -110.3422

# $LRTest
# [1] 185.0736

# $Prob
# [1] 3.78202e-42

# $AICc.obs
# [1] 43.61086

# $AICc.constrained
# [1] 226.6845

# $optimmessage
# [1] "Optimization has converged."