---
title: "Assignment4"
author: "Tabea"
date: "29 November 2018"
output: 
  html_document: 
    keep_md: yes
---
#2

```r
alleles.diploid <- function(p0,w1,w2,w3,n=100,na.rm=TRUE) {
  p<-rep(NA,n)
  w_aver<-rep(NA,n)
  #starting frequencies
  p[1]<-p0 
  w_aver[1]<-(((p[1]^2)*w1) + (2*p[1]*(1-p[1])*w2) + (((1-p[1])^2)*w3))
  #loop
  for(i in 2:n) {
    w_aver[i-1]<-(((p[i-1]^2)*w1) + (2*p[i-1]*(1-p[i-1])*w2) + (((1-p[i-1])^2)*w3))
    p[i]<-((p[i-1]^2)*(w1/w_aver[i-1])+(p[i-1]*(1-p[i-1])*(w2/w_aver[i-1])))
  }
  #output sentence and plot
  if(any(p>0.9999)) {
    fixationP<-min(which.max(p>0.9999))
    cat("Fixation for A occurs approximately at generation:", fixationP)
  } else if(any(p<0.0001)) {
    fixationQ<-min(which.max(p<0.0001))
    cat("Fixation for a occurs approximately at generation:", fixationQ)
  } else {
    maxAlleleFreq<-max(p)
      cat("Fixation does not occur, max. allele frequency for A is:", print(maxAlleleFreq, digits = 2))
  }
  plot(x = 1:n, y = p,pch=20,ylab="allele frequency", xlab="generations")
}

combi1<-alleles.diploid(p0=0.5,w1=0.9,w2=0.7,w3=1)
```

```
## Fixation for a occurs approximately at generation: 36
```

![](Assignment4_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
combi2<-alleles.diploid(p0=0.4,w1=0.8,w2=0.9,w3=0.6)
```

```
## [1] 0.75
## Fixation does not occur, max. allele frequency for A is: 0.7499841
```

![](Assignment4_files/figure-html/unnamed-chunk-1-2.png)<!-- -->

```r
combi3<-alleles.diploid(p0=0.5,w1=0.9,w2=0.8,w3=0.7)
```

```
## Fixation for A occurs approximately at generation: 79
```

![](Assignment4_files/figure-html/unnamed-chunk-1-3.png)<!-- -->

```r
combi1
```

```
## NULL
```

```r
combi2
```

```
## NULL
```

```r
combi3
```

```
## NULL
```
#3

```r
allele_counts <- function(number_alleles,A,n) {
  allele_freq<-rep(NA,n)
  #1st generation
  allele_freq[1]<-A
  for(i in 2:n) {
   alleles<-sample(c("A","a"),
           size=(number_alleles),
           replace = TRUE,
           prob=c(allele_freq[i-1], (1-allele_freq[i-1])))
   allele_freq[i]<-sum(alleles=="A")/number_alleles
  }
   plot(x = 1:n, y = allele_freq,ylab="allele frequency", xlab="generations",type="l")
   return(allele_freq)
}

allele_counts(20,0.3,50)
```

![](Assignment4_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```
##  [1] 0.30 0.20 0.20 0.05 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
## [15] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
## [29] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
## [43] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
```

```r
allele_counts(20,0.5,100)
```

![](Assignment4_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

```
##   [1] 0.50 0.65 0.65 0.75 0.75 0.70 0.65 0.60 0.55 0.70 0.70 0.65 0.65 0.55
##  [15] 0.75 0.70 0.50 0.60 0.65 0.70 0.50 0.50 0.55 0.55 0.45 0.55 0.40 0.40
##  [29] 0.30 0.20 0.25 0.40 0.40 0.40 0.30 0.30 0.25 0.15 0.10 0.10 0.15 0.05
##  [43] 0.10 0.05 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
##  [57] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
##  [71] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
##  [85] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
##  [99] 0.00 0.00
```

#4

```r
allele_counts_4 <- function(number_alleles,A,n) {
  allele_freq<-rep(NA,n)
  #1st generation
  allele_freq[1]<-A
  for(i in 2:n) {
   alleles<-sample(c("A","a"),
           size=(number_alleles),
           replace = TRUE,
           prob=c(allele_freq[i-1], (1-allele_freq[i-1])))
   allele_freq[i]<-sum(alleles=="A")/number_alleles
  }
  return(allele_freq)
}

alleles_1000<- function(x,number_alleles,A,n) {
  p<-replicate(x,allele_counts_4(number_alleles,A,n))
  A_lost<-sum(p[100,]==0)/x
  return(A_lost)
}

alleles_1000(1000,400,0.5,100)
```

```
## [1] 0.005
```

```r
alleles_1000(1000,400,0.25,100)
```

```
## [1] 0.132
```

```r
alleles_1000(1000,400,0.1,100)
```

```
## [1] 0.411
```

#5

```r
plot(x=1,y=0,type="n",
     xlab="Generation",ylab="Allele frequency",
     xlim = c(0,100), ylim=c(0,1))
allele_counts_5 <- function(number_alleles,A,n,x) {
  allele_freq<-rep(NA,n)
  #1st generation
  allele_freq[1]<-A
  for(i in 2:n) {
   alleles<-sample(c("A","a"),
           size=(number_alleles),
           replace = TRUE,
           prob=c(allele_freq[i-1], (1-allele_freq[i-1])))
   allele_freq[i]<-sum(alleles=="A")/number_alleles
  }
  lines(allele_freq~c(1:x),col=sample(rainbow(6),1))
}
replicate(100,allele_counts_5(400,0.5,100,100))
```

![](Assignment4_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
## 
## [[6]]
## NULL
## 
## [[7]]
## NULL
## 
## [[8]]
## NULL
## 
## [[9]]
## NULL
## 
## [[10]]
## NULL
## 
## [[11]]
## NULL
## 
## [[12]]
## NULL
## 
## [[13]]
## NULL
## 
## [[14]]
## NULL
## 
## [[15]]
## NULL
## 
## [[16]]
## NULL
## 
## [[17]]
## NULL
## 
## [[18]]
## NULL
## 
## [[19]]
## NULL
## 
## [[20]]
## NULL
## 
## [[21]]
## NULL
## 
## [[22]]
## NULL
## 
## [[23]]
## NULL
## 
## [[24]]
## NULL
## 
## [[25]]
## NULL
## 
## [[26]]
## NULL
## 
## [[27]]
## NULL
## 
## [[28]]
## NULL
## 
## [[29]]
## NULL
## 
## [[30]]
## NULL
## 
## [[31]]
## NULL
## 
## [[32]]
## NULL
## 
## [[33]]
## NULL
## 
## [[34]]
## NULL
## 
## [[35]]
## NULL
## 
## [[36]]
## NULL
## 
## [[37]]
## NULL
## 
## [[38]]
## NULL
## 
## [[39]]
## NULL
## 
## [[40]]
## NULL
## 
## [[41]]
## NULL
## 
## [[42]]
## NULL
## 
## [[43]]
## NULL
## 
## [[44]]
## NULL
## 
## [[45]]
## NULL
## 
## [[46]]
## NULL
## 
## [[47]]
## NULL
## 
## [[48]]
## NULL
## 
## [[49]]
## NULL
## 
## [[50]]
## NULL
## 
## [[51]]
## NULL
## 
## [[52]]
## NULL
## 
## [[53]]
## NULL
## 
## [[54]]
## NULL
## 
## [[55]]
## NULL
## 
## [[56]]
## NULL
## 
## [[57]]
## NULL
## 
## [[58]]
## NULL
## 
## [[59]]
## NULL
## 
## [[60]]
## NULL
## 
## [[61]]
## NULL
## 
## [[62]]
## NULL
## 
## [[63]]
## NULL
## 
## [[64]]
## NULL
## 
## [[65]]
## NULL
## 
## [[66]]
## NULL
## 
## [[67]]
## NULL
## 
## [[68]]
## NULL
## 
## [[69]]
## NULL
## 
## [[70]]
## NULL
## 
## [[71]]
## NULL
## 
## [[72]]
## NULL
## 
## [[73]]
## NULL
## 
## [[74]]
## NULL
## 
## [[75]]
## NULL
## 
## [[76]]
## NULL
## 
## [[77]]
## NULL
## 
## [[78]]
## NULL
## 
## [[79]]
## NULL
## 
## [[80]]
## NULL
## 
## [[81]]
## NULL
## 
## [[82]]
## NULL
## 
## [[83]]
## NULL
## 
## [[84]]
## NULL
## 
## [[85]]
## NULL
## 
## [[86]]
## NULL
## 
## [[87]]
## NULL
## 
## [[88]]
## NULL
## 
## [[89]]
## NULL
## 
## [[90]]
## NULL
## 
## [[91]]
## NULL
## 
## [[92]]
## NULL
## 
## [[93]]
## NULL
## 
## [[94]]
## NULL
## 
## [[95]]
## NULL
## 
## [[96]]
## NULL
## 
## [[97]]
## NULL
## 
## [[98]]
## NULL
## 
## [[99]]
## NULL
## 
## [[100]]
## NULL
```

#6

```r
set.seed(720) #confirm to always get the same results

stochastics<-function(intercept,slope,l.o, sd0) {
  x <- seq(from = 1, to = 10, length.out = l.o) 
  a <- intercept
  b <- slope 
  y_deterministic <- a + b*x
  y_simulated <- rnorm(length(x), mean = y_deterministic, sd = sd0)
  mod_sim <- lm(y_simulated ~ x)
  p_val_slope <- summary(mod_sim)$coef[2,4] # extracts the p-value
  p_val_slope
}

stochastics(0.5,0.1,20,2)
```

```
## [1] 0.3620625
```

```r
rep.stoch<-replicate(1000,stochastics(0.5,0.1,20,2))
hist(rep.stoch)
```

![](Assignment4_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
p_less0.05<-sum(rep.stoch<0.05)/length(rep.stoch)
p_less0.05
```

```
## [1] 0.107
```

```r
rep.stoch_2<-replicate(1000,stochastics(0.5,0,20,2))
hist(rep.stoch_2)
```

![](Assignment4_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

```r
p_less0.05<-sum(rep.stoch_2<0.05)/length(rep.stoch_2)
p_less0.05
```

```
## [1] 0.048
```

```r
#The amount of p<0.05 is less

sequence<-seq(10,100,5)
p<-rep(NA,100*length(sequence))
proportion<-sum(p<0.05)/length(p)

for (i in 1:length(sequence)) {
  p<-replicate(100,stochastics(0.5,0.1,sequence[i],1.5))
  proportion[i]<-sum(p<0.05)/length(p)
}

plot(x=sequence,y=proportion)
```

![](Assignment4_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

```r
# The higher the sample  size, the more higher is the probability to have a value below 0.05 
```
