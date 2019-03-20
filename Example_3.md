### 1-Data

We load the wheat data set in the BGLR package, extract one phenotype, scale and center genotypes and split the data into a training and a testing set.

```r
 ##### DATA #############################################
 library(BGLR)
 data(wheat); X=scale(wheat.X); Y=wheat.Y
 objects();dim(X);dim(Y)
 N<-nrow(X) ; p<-ncol(X)
 y<-Y[,2]
 set.seed(12345)
 tst<-sample(1:N,size=150,replace=FALSE)
 XTRN<-X[-tst,]
 yTRN<-y[-tst]
 XTST<-X[tst,]
 yTST<-y[tst]
```

### 2: Bayesian Regressions


We show here how to fit various Bayesian models using the  [BGLR R-package](https://cran.r-project.org/web/packages/BGLR/index.htm).

Resources: 

  - [BGLR-GitHub](https://github.com/gdlc/BGLR-R)  (include multiple examples)
  - [GENETICS-manuscript](http://www.genetics.org/content/198/2/483)

**2a) Model fitting**

```r
  nIter=1500 # I set this to small value that way it will run quickly, for more serious analyses use longer chains
  burnIn=500 # and longer burnin

 # Gaussian prior ("Bayesian Ridge-Regression")
  LP=list( list(X=XTRN,model='BRR') ) # 2-level list, allows specifying different types of random and fixed effects
  fmBRR=BGLR(y=yTRN,ETA=LP,nIter=nIter,burnIn=burnIn,saveAt='BRR_',verbose=FALSE)
   
  # Scaled-t
   LP[[1]]$model='BayesA'
   fmBA=BGLR(y=yTRN,ETA=LP,nIter=nIter,burnIn=burnIn,saveAt='BA_',verbose=FALSE)
   
  # Double-Exponential
   LP[[1]]$model='BL'
   fmBL=BGLR(y=yTRN,ETA=LP,nIter=nIter,burnIn=burnIn,saveAt='BL_',verbose=FALSE)
   
  # Spike-slab (Gaussian)
   LP[[1]]$model='BayesC'
   fmBC=BGLR(y=yTRN,ETA=LP,nIter=nIter,burnIn=burnIn,saveAt='BC_',verbose=FALSE)
   
  # Spike-slab (Scaled-t)
   LP[[1]]$model='BayesB'
   fmBB=BGLR(y=yTRN,ETA=LP,nIter=nIter,burnIn=burnIn,saveAt='BB_',verbose=FALSE)

```

**2b) Retrieving samples and estimates**

```r
  # Items available to every model
   head(fmBRR$yHat) # predictions (if they were NAs there will be predictions for those points as well)
   fmBRR$mu # intercept (alwasy included by defaul)
   fmBRR$varE # error variance
   fmBRR$fit # DIC and other statistics
   fmBRR$nIter #... other run parameters...
   varE=scan('BRR_varE.dat') # samples from the posterior distribution of the error variance
   plot(varE,type='o',col=4,main='Trace plot of the error variance') 
   
  # Gaussian prior 
   head(fmBRR$ETA[[1]]$b)  # posterior means of effects
   head(fmBRR$ETA[[1]]$SD.b) # posterior SDs
   fmBRR$ETA[[1]]$varB ; fmBRR$ETA[[1]]$SD.varB
  
  # Bayes A
    head(fmBA$ETA[[1]]$b)
    head(fmBA$ETA[[1]]$varB)  
  
  # Bayesian Lasso
    head(fmBL$ETA[[1]]$b)
    head(fmBL$ETA[[1]]$tau) 
    
  # BayesC
    head(fmBC$ETA[[1]]$b)
    head(fmBC$ETA[[1]]$d) # posterior probability of inclussion

  # BayesB
    head(fmBB$ETA[[1]]$b)
    head(fmBB$ETA[[1]]$varB)
    head(fmBB$ETA[[1]]$d) # posterior probability of inclussion

```

**2c) Retrieving samples and estimates**

```r 
   YYHat=cbind( yTST, 
                XTST%*%fmBRR$ETA[[1]]$b ,
                XTST%*%fmBA$ETA[[1]]$b ,
                XTST%*%fmBL$ETA[[1]]$b ,
                XTST%*%fmBC$ETA[[1]]$b ,
                XTST%*%fmBB$ETA[[1]]$b)
   colnames(YYHat)=c('y','BRR','BA','BL','BC','BB')
   cor(YYHat)
```
