# lnMLE:  Code for Heagerty and Zeger (1999)
Bruce Swihart  
Tuesday, October 21, 2014  

To get started, be sure to have devtools working properly and then run the next three lines.  A lot of warnings will spew forth, but just ignore them.  I included the output for these three lines at the bottom of this readMe in an Appendix.


```r
library(devtools)
install_github("swihart/lnMLE_1.0-2")
library(lnMLE)
```

Then run the example from [Swihart et al 2014](http://onlinelibrary.wiley.com/doi/10.1111/insr.12035/abstract), which concerns relating eye impairment (1-impaired, 0-healthy) to race (White or Black).  Since each individual had two eyes measured, the data contain clustered binary outcomes.  However, since no one changed race, the subject-specific interpretation of the conditional model fixed effects would have a counterfactual context as opposed to a marginal model where the fixed effects coefficients would be comparing the prevalence between the two races.


```r
#
## Eye and Race Example
#
data(eye_race)
attach(eye_race)
marg_model <- logit.normal.mle(meanmodel = value ~ black,
                           logSigma= ~1,
                           id=eye_race$id,
                           model="marginal",
                           data=eye_race,
                           tol=1e-5,
                           maxits=100,
                           r=50)
marg_model
```

```
## 
##   ML Estimation for Logistic-Normal Models
## 
##   model = marginal
## 
##   options:     lambda = 0
##                     r = 50 
## 
## 
##   Mean Parameters: 
## 
##             estimate std. err.      Z  
## (Intercept)  -2.4312   0.05787 -42.0128
## black         0.0822   0.08438   0.9742
## 
## 
##   Variance Components: 
## 
##             estimate std. err.   Z  
## (Intercept)    1.113   0.05315 20.94
## 
## 
##   Maximized logL = -2699.7648
```

```r
cond_model <- logit.normal.mle(meanmodel = value ~ black,
                           logSigma= ~1,
                           id=eye_race$id,
                           model="conditional",
                           data=eye_race,
                           tol=1e-5,
                           maxits=100,
                           r=50)
cond_model
```

```
## 
##   ML Estimation for Logistic-Normal Models
## 
##   model = conditional
## 
##   options:     lambda = 0
##                     r = 50 
## 
## 
##   Mean Parameters: 
## 
##             estimate std. err.      Z  
## (Intercept)  -4.9411    0.2118 -23.3311
## black         0.1463    0.1503   0.9736
## 
## 
##   Variance Components: 
## 
##             estimate std. err.   Z  
## (Intercept)    1.113   0.05315 20.94
## 
## 
##   Maximized logL = -2699.7648
```

```r
compare<-round(cbind(marg_model$beta, cond_model$beta),2)
colnames(compare)<-c("Marginal", "Conditional")
compare
```

```
##             Marginal Conditional
## (Intercept)    -2.43       -4.94
## black           0.08        0.15
```

The conditional logistic regression model race effect differs substantially from the marginalized model.  The marginal effect is attenuated.







Appendix:  The warnings that are generated from the first three lines of code of this readme (just ignore or email me if an easy fix)



```r
# R version 3.1.1 (2014-07-10) -- "Sock it to Me"
# Copyright (C) 2014 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# 
# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.
# 
# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.
# 
# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.
# 
# > library(devtools)
# 
# Attaching package: 'devtools'
# 
# The following objects are masked from 'package:utils':
# 
#     ?, help
# 
# The following object is masked from 'package:base':
# 
#     system.file
# 
# > install_github("swihart/lnMLE_1.0-2")
# Installing github repo lnMLE_1.0-2/master from swihart
# Downloading master.zip from https://github.com/swihart/lnMLE_1.0-2/archive/master.zip
# Installing package from C:\Users\SWIHAR~1\AppData\Local\Temp\RtmpoNMgtK/master.zip
# Installing lnMLE
# "C:/R/R311/bin/x64/R" --vanilla CMD  \
#   INSTALL  \
#   "C:\Users\swihartbj\AppData\Local\Temp\RtmpoNMgtK\devtools13d58775c79b4\lnMLE_1.0-2-master"  \
#   --library="C:/RLIB" --install-tests 
# 
# * installing *source* package 'lnMLE' ...
# ** libs
# 
# *** arch - i386
# gcc -m32 -I"C:/R/R311/include" -DNDEBUG     -I"d:/RCompile/CRANpkg/extralibs64/local/include"     -O3 -Wall  -std=gnu99 -mtune=core2 -c chanmat.c -o chanmat.o
# chanmat.c: In function 'matadd':
# chanmat.c:507:16: warning: unused variable 'z' [-Wunused-variable]
# chanmat.c:507:11: warning: unused variable 'nlen' [-Wunused-variable]
# chanmat.c: In function 'matsub':
# chanmat.c:539:11: warning: unused variable 'nlen' [-Wunused-variable]
# chanmat.c: In function 'matmult':
# chanmat.c:571:17: warning: unused variable 'nlen' [-Wunused-variable]
# chanmat.c: In function 'matxdiagasvec':
# chanmat.c:613:8: warning: unused variable 'tmp' [-Wunused-variable]
# chanmat.c: In function 'matread':
# chanmat.c:793:11: warning: unused variable 'fmt' [-Wunused-variable]
# chanmat.c: In function 'luinv':
# chanmat.c:881:20: warning: variable 'outsub2' set but not used [-Wunused-but-set-variable]
# chanmat.c:881:11: warning: variable 'outsub1' set but not used [-Wunused-but-set-variable]
# chanmat.c:878:33: warning: unused variable 'i' [-Wunused-variable]
# gcc -m32 -I"C:/R/R311/include" -DNDEBUG     -I"d:/RCompile/CRANpkg/extralibs64/local/include"     -O3 -Wall  -std=gnu99 -mtune=core2 -c clinluxxy.c -o clinluxxy.o
# gcc -m32 -I"C:/R/R311/include" -DNDEBUG     -I"d:/RCompile/CRANpkg/extralibs64/local/include"     -O3 -Wall  -std=gnu99 -mtune=core2 -c logit.normal.mle.c -o logit.normal.mle.o
# logit.normal.mle.c: In function 'logit_normal_mle':
# logit.normal.mle.c:56:44: warning: unused variable 'logPi_indep' [-Wunused-variable]
# logit.normal.mle.c:55:33: warning: unused variable 'pm' [-Wunused-variable]
# logit.normal.mle.c:48:17: warning: unused variable 'eta' [-Wunused-variable]
# logit.normal.mle.c:48:12: warning: unused variable 'mu' [-Wunused-variable]
# logit.normal.mle.c: In function 'split':
# logit.normal.mle.c:645:8: warning: unused variable 'j' [-Wunused-variable]
# gcc -m32 -shared -s -static-libgcc -o lnMLE.dll tmp.def chanmat.o clinluxxy.o logit.normal.mle.o -Ld:/RCompile/CRANpkg/extralibs64/local/lib/i386 -Ld:/RCompile/CRANpkg/extralibs64/local/lib -LC:/R/R311/bin/i386 -lR
# installing to C:/RLIB/lnMLE/libs/i386
# 
# *** arch - x64
# gcc -m64 -I"C:/R/R311/include" -DNDEBUG     -I"d:/RCompile/CRANpkg/extralibs64/local/include"     -O2 -Wall  -std=gnu99 -mtune=core2 -c chanmat.c -o chanmat.o
# chanmat.c: In function 'create_matrix':
# chanmat.c:23:20: warning: format '%d' expects argument of type 'int', but argument 3 has type 'long long unsigned int' [-Wformat]
# chanmat.c: In function 'matadd':
# chanmat.c:507:16: warning: unused variable 'z' [-Wunused-variable]
# chanmat.c:507:11: warning: unused variable 'nlen' [-Wunused-variable]
# chanmat.c: In function 'matsub':
# chanmat.c:539:11: warning: unused variable 'nlen' [-Wunused-variable]
# chanmat.c: In function 'matmult':
# chanmat.c:571:17: warning: unused variable 'nlen' [-Wunused-variable]
# chanmat.c: In function 'matxdiagasvec':
# chanmat.c:613:8: warning: unused variable 'tmp' [-Wunused-variable]
# chanmat.c: In function 'matread':
# chanmat.c:793:11: warning: unused variable 'fmt' [-Wunused-variable]
# chanmat.c: In function 'luinv':
# chanmat.c:881:20: warning: variable 'outsub2' set but not used [-Wunused-but-set-variable]
# chanmat.c:881:11: warning: variable 'outsub1' set but not used [-Wunused-but-set-variable]
# chanmat.c:878:33: warning: unused variable 'i' [-Wunused-variable]
# gcc -m64 -I"C:/R/R311/include" -DNDEBUG     -I"d:/RCompile/CRANpkg/extralibs64/local/include"     -O2 -Wall  -std=gnu99 -mtune=core2 -c clinluxxy.c -o clinluxxy.o
# gcc -m64 -I"C:/R/R311/include" -DNDEBUG     -I"d:/RCompile/CRANpkg/extralibs64/local/include"     -O2 -Wall  -std=gnu99 -mtune=core2 -c logit.normal.mle.c -o logit.normal.mle.o
# logit.normal.mle.c: In function 'logit_normal_mle':
# logit.normal.mle.c:56:44: warning: unused variable 'logPi_indep' [-Wunused-variable]
# logit.normal.mle.c:55:33: warning: unused variable 'pm' [-Wunused-variable]
# logit.normal.mle.c:48:17: warning: unused variable 'eta' [-Wunused-variable]
# logit.normal.mle.c:48:12: warning: unused variable 'mu' [-Wunused-variable]
# logit.normal.mle.c: In function 'split':
# logit.normal.mle.c:645:8: warning: unused variable 'j' [-Wunused-variable]
# gcc -m64 -shared -s -static-libgcc -o lnMLE.dll tmp.def chanmat.o clinluxxy.o logit.normal.mle.o -Ld:/RCompile/CRANpkg/extralibs64/local/lib/x64 -Ld:/RCompile/CRANpkg/extralibs64/local/lib -LC:/R/R311/bin/x64 -lR
# installing to C:/RLIB/lnMLE/libs/x64
# ** R
# ** data
# ** inst
# ** preparing package for lazy loading
# ** help
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/logit.normal.mle.Rd:53: unexpected '}'
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/madras.Rd:11: unknown macro '\item'
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/madras.Rd:14: unknown macro '\item'
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/madras.Rd:17: unknown macro '\item'
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/madras.Rd:20: unknown macro '\item'
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/madras.Rd:23: unknown macro '\item'
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/madras.Rd:26: unknown macro '\item'
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/madras.Rd:29: unknown macro '\item'
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/madras.Rd:32: unknown macro '\item'
# Warning: C:/Users/swihartbj/AppData/Local/Temp/RtmpoNMgtK/devtools13d58775c79b4/lnMLE_1.0-2-master/man/madras.Rd:38: unknown macro '\bf'
# *** installing help indices
# ** building package indices
# ** testing if installed package can be loaded
# *** arch - i386
# *** arch - x64
# * DONE (lnMLE)
# > library(lnMLE)
```

