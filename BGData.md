
## BGData R-package

**[BGData](https://github.com/quantgen/bgdata)** is a suite of R-packages that provide memory mapping for `plink-bed files` and classes and methods for linked arrays that allow you to handle extremely large data sets in R without actually loading the data in memory. The suite also provides basic functionality for genomic analysis (GWAS, G-matrices) in R.

The code below provides basic code demonstrating the use of the package. Fruther examples can be found at the [BGData](https://github.com/quantgen/bgdata) GitHub repository. To illustrate I use data from the UK-Biobank stored in `.bed` files. 

It takes ~6sec to create a memory mapped object for the SNPs in chromosme 1 (n~490,000, p~1.2M SNPs!)


```r
 library(BGData)
 X1<-BEDMatrix('chrom01.bed') 
 # dim(X1) [1]  487395 1208761

```

Once the object is created you can extract data to memory using regular indexing operators.

```r
 X1[1,1]
 head(X1[,1])
 head(colnames(X1))
 head(colnames(X2))
```

We can link data in multiple files in a single array.


```r 
 X2<-BEDMatrix('chrom02.bed') 
 X3<-BEDMatrix('chrom03.bed')  
 X=ColumnLinkedMatrix(X1,X2,X3)
 #ncol(X)~3.6M SNPs
```

Or the whole genome...


```r 
 bedFiles=list.files(pattern='*.bed')
 DATA=list()
 for(i in 1:length(bedFiles)){ DATA[[i]]=BEDMatrix(bedFiles[i])}
 X=X=do.call(ColumnLinkedMatrix,DATA)
 dim(X)
 #[1]   487395 15356050

```




