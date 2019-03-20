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

## 2: Penalized regressions using glmnet


**2a) Fitting the models**

```r
 library(glmnet)

 # alpha 0 gives Ridge Regression
 fmRR=glmnet(y=yTRN,x=XTRN,alpha=0)
 
 # alpha 1 gives Lassso
 fmL=glmnet(y=yTRN,x=XTRN,alpha=1)
 
 # alpha between 0 and 1 gives elastic net
 fmEN=glmnet(y=yTRN,x=XTRN, alpha=0.5)

 
 COR.RR=rep(NA,100)
 COR.L=rep(NA,100)
 COR.ENet=rep(NA,100)
 
 # evaluating correlation in TST set
 for(i in 1:100){
   COR.RR[i]=cor(yTST,XTST%*%fmRR$beta[,i])
   COR.L[i]=cor(yTST,XTST%*%fmL$beta[,i])
   COR.ENet[i]=cor(yTST,XTST%*%fmEN$beta[,i])
 }
 
 plot(COR.RR,x=log(fmRR$lambda),type='o',col='blue',cex=.5)
 plot(COR.L,x=fmL$df,type='o',  col='blue',cex=.5,xlab='# of active markers')
 plot(COR.ENet,x=fmEN$df,type='o',col='blue',cex=.5,xlab='# of active markers')

 
```


**2b) Extracting estimates**

```r
  dim(fmL$b) #estimated effects over the regularization path
  dim(fmL$a0) # same for intercept
  head(fmL$df) # number of active markers over the regularization path
  fmL$lambda
  fmL$alpha
```

**Problem**: If we chose lambda based on the curves ploted above the accuracy is over estimated because we are choosing a model that is optimal for a particular training set, most likely, averaged over possible training sets, the accuracy will be lower. How do we estimate accuracy when we need to select 'hyper-parameters'.


**Suggested approach**

   - Conduct an CV internal to the training set
   - Select the best model (i.e. the lambda that gives highest accuracy)
   - Refit the model to the entire training set using the chosen lambda
   - Evaluate accuracy in the testing set.
   
   
   
