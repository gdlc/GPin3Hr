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


## 1-Selecting markers using single-marker regression

```r
 pValues<-numeric()
 for(i in 1:p){
	fm<-lsfit(y=yTRN,x=XTRN[,i])
	pValues[i]<-ls.print(fm,print.it=F)$coef[[1]][2,4] # extracts p-value, similar to lm() but a bit faster
 }
 
 plot(-log10(pValues),cex=.5,col=2)
 
####### VARIABLE SELECTION ##############################
 mrk_rank<-order(pValues); corTRN<-numeric(); corTST<-numeric()
 for(i in 1:300){	
	tmpIndex<- mrk_rank[1:i]
	ZTRN=XTRN[,tmpIndex,drop=F]
	ZTST=XTST[,tmpIndex,drop=F]
	
	fm<-lm(yTRN~ZTRN)
	bHat=coef(fm)[-1]
	bHat<-ifelse(is.na(bHat),0,bHat)
	
	yHatTRN=ZTRN%*%bHat
  corTRN[i]<-cor(yTRN,yHatTRN)
  
	yHatTST=ZTST%*%bHat
	corTST[i]<-cor(yTST,yHatTST)
 }
 
 plot(c(0,corTRN),x=0:length(corTRN),type='o',col=2,ylab='Correlation-Training',
       xlab='Number of markers',ylim=c(0,.9))
 
 lines(x=0:length(corTST),y=c(0,corTST),col=4)
 points(x=0:length(corTST),y=c(0,corTST),col=4)
 
 
```


**Suggested problems**:

  - Try forward selection
  - Repeat the above code 100 times each time with a different random split of the data into training and testing.
  
