# nlmixr2auto
Automated Population Pharmacokinetic Modelling. This package offers serveral optimisation algorithms designed for automated population pharmacokinetic modeling, serving as a valuable tool for pharmacokinetic model selection.

## R package installation

Installation nlmixr2autoinit and nlmixr2auto:
``` r
library(devtools)
install_github("ucl-pharmacometrics/nlmixr2autoinit")
install_github("ucl-pharmacometrics/nlmixr2auto")
```

Stepwsie example
``` r

library(nlmixr2autoinit)
library(nlmixr2auto)

outs<-sf.operator(dat=pheno_sd,
            search.space = "ivbase",
            filename =  "pheno_sd",
            foldername =   "pheno_sd" )
print(outs)

# Infometrics                               Value          
# ----------------------------------------  ---------------
#   Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Total Number of Subjects                  59             
# Total Number of Observations              155            
# Subjects with First-Dose Interval Data    35             
# Observations in the First-Dose Interval   35             
# Subjects with Multiple-Dose Data          56             
# Observations after Multiple Doses         120            
# ----------------------------------------  ------
#   Estimating half-life....................
# Half-life estimation complete: Estimated t½ = 16.44 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# (hybrid mode: parameters combined across sources)....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 0.0087, estimated Vd: 1.25 
# Run parameter sweeping on nonlinear eliminiation kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters.................... 
# Running Stepwise 1. Structural Model----------------------------------------------------
#   Test number of compartments----------------------------------------------------
#   [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod1.txt
# [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod2.txt
# [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod3.txt
# Analyse elimination type----------------------------------------------------
#   [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod4.txt
# [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod5.txt
# Test IIV on parameters----------------------------------------------------
#   [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod6.txt
# [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod7.txt
# Test Correlation between parameters----------------------------------------------------
#   [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod8.txt
# [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod9.txt
# Explore types of residual errors----------------------------------------------------
#   [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod10.txt
# [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod11.txt
# [Success] Model file created in current working directory:
#   /home/zhonghuihuang/Step_2025-08-11-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod12.txt
# > print(outs)
# 
# 
# === Best Model Code ===
# no.cmpt  eta.km  eta.vc  eta.vp eta.vp2   eta.q  eta.q2      mm   mcorr      rv 
#       1       0       1       0       0       0       0       0       1       1 
# 
# === Best Model Name ===
# iv_1cmpt_etaCLVC_First-order elimination_Eta_correlated_additive 
# 
# === Stepwise Selection History ===
#                   Step                                         Penalty.terms                                                       Model.name          Model.code    Fitness
# 1  No. of compartments                                rse, theta, covariance   iv_1cmpt_etaCL_First-order elimination_No correlation_combined 1,0,0,0,0,0,0,0,0,3 11167.0409
# 2     Elimination type                                rse, theta, covariance   iv_1cmpt_etaCL_First-order elimination_No correlation_combined 1,0,0,0,0,0,0,0,0,3 11167.0409
# 3        IIV (forward)              rse, theta, covariance, shrinkage, omega iv_1cmpt_etaCLVC_First-order elimination_No correlation_combined 1,0,1,0,0,0,0,0,0,3 11076.5477
# 4      ETA correlation rse, theta, covariance, shrinkage, omega, correlation iv_1cmpt_etaCLVC_First-order elimination_Eta_correlated_combined 1,0,1,0,0,0,0,0,1,3  1065.5562
# 5 Residual error types rse, theta, covariance, shrinkage, omega, correlation iv_1cmpt_etaCLVC_First-order elimination_Eta_correlated_additive 1,0,1,0,0,0,0,0,1,1   995.9198
```

Automated Modelling Examples: 
（The following example runs the model with a small quantity of runs and iterations. These numbers can be set according to actual needs or just set simply default values.）
GA example
``` r
ga.operator(dat = d1,
            npopsize = 6,
            max.iter = 6,
            nlocal.search = 3,
            no.cores = getRxThreads(),
            thetalower=c(vp=1,
                         vp2=1),
            control = saemControl(
              seed = 1234,
              print = 5,
              nBurn = 50,
              nEm = 50,
              logLik = T),
            table=tableControl(cwres=T),
            filename =  "pheno_sd",
            foldername =   "pheno_sd" )

# Running GA ---------------------------------------------- Model: 1 Iteration: 1
# Running GA ---------------------------------------------- Model: 2 Iteration: 1
# Running GA ---------------------------------------------- Model: 3 Iteration: 1
# Running GA ---------------------------------------------- Model: 4 Iteration: 1
# Running GA ---------------------------------------------- Model: 5 Iteration: 1
# Running GA ---------------------------------------------- Model: 6 Iteration: 1
# Running GA ---------------------------------------------- Model: 7 Iteration: 2
# Running GA ---------------------------------------------- Model: 8 Iteration: 2
# Running GA ---------------------------------------------- Model: 9 Iteration: 2
# Running GA ---------------------------------------------- Model: 10 Iteration: 2
# Running GA ---------------------------------------------- Model: 11 Iteration: 2
# Running GA ---------------------------------------------- Model: 12 Iteration: 2
# Running GA ---------------------------------------------- Model: 13 Iteration: 3
# Running GA ---------------------------------------------- Model: 14 Iteration: 3
# Running GA ---------------------------------------------- Model: 15 Iteration: 3
# Running GA ---------------------------------------------- Model: 16 Iteration: 3
# Running GA ---------------------------------------------- Model: 17 Iteration: 3
# Running GA ---------------------------------------------- Model: 18 Iteration: 3
# Running GA ---------------------------------------------- Model: 19 Iteration: 3
# Running GA ---------------------------------------------- Model: 20 Iteration: 3
# Running GA ---------------------------------------------- Model: 21 Iteration: 3
# Running GA ---------------------------------------------- Model: 22 Iteration: 3
# Running GA ---------------------------------------------- Model: 23 Iteration: 3
# Running GA ---------------------------------------------- Model: 24 Iteration: 3
# Running GA ---------------------------------------------- Model: 25 Iteration: 3
# Running GA ---------------------------------------------- Model: 26 Iteration: 4
# Running GA ---------------------------------------------- Model: 27 Iteration: 4
# Running GA ---------------------------------------------- Model: 28 Iteration: 4
# Running GA ---------------------------------------------- Model: 29 Iteration: 4
# Running GA ---------------------------------------------- Model: 30 Iteration: 4
# Running GA ---------------------------------------------- Model: 31 Iteration: 4
# Running GA ---------------------------------------------- Model: 32 Iteration: 5
# Running GA ---------------------------------------------- Model: 33 Iteration: 5
# Running GA ---------------------------------------------- Model: 34 Iteration: 5
# Running GA ---------------------------------------------- Model: 35 Iteration: 5
# Running GA ---------------------------------------------- Model: 36 Iteration: 5
# Running GA ---------------------------------------------- Model: 37 Iteration: 5
# Running GA ---------------------------------------------- Model: 38 Iteration: 6
# Running GA ---------------------------------------------- Model: 39 Iteration: 6
# Running GA ---------------------------------------------- Model: 40 Iteration: 6
# Running GA ---------------------------------------------- Model: 41 Iteration: 6
# Running GA ---------------------------------------------- Model: 42 Iteration: 6
# Running GA ---------------------------------------------- Model: 43 Iteration: 6
# Running GA ---------------------------------------------- Model: 44 Iteration: 6
# Running GA ---------------------------------------------- Model: 45 Iteration: 6
# Running GA ---------------------------------------------- Model: 46 Iteration: 6
# Running GA ---------------------------------------------- Model: 47 Iteration: 6
# Running GA ---------------------------------------------- Model: 48 Iteration: 6
# Running GA ---------------------------------------------- Model: 49 Iteration: 6
# Running GA ---------------------------------------------- Model: 50 Iteration: 6
# Running GA ---------------------------------------------- Model: 51 Iteration: 6
# $final.selected.code
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 2        0        1      0      1      0       0     0      0  0     0   0   1
# 
# $final.selected.model
# [1] "1Cmpt,IIV.cl.vc,first-order_elminiation,nocorr,additive"
# 
# $history
# $history$sel.best.code.list
# $history$sel.best.code.list[[1]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        0        1      0      0      0       0     0      0  0     0   0   1
# 
# $history$sel.best.code.list[[2]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 6        0        1      0      0      0       0     0      0  0     0   0   1
# 
# $history$sel.best.code.list[[3]]
# V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12
# 4  0  1  0  1  0  0  0  0  0   0   0   1
# 
# $history$sel.best.code.list[[4]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 4        0        1      0      1      0       0     0      0  0     0   0   1
# 
# $history$sel.best.code.list[[5]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 6        0        1      0      1      0       0     0      0  0     0   0   1
# 
# $history$sel.best.code.list[[6]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 2        0        1      0      1      0       0     0      0  0     0   0   1
# 
# 
# $history$data.pop.list
# $history$data.pop.list[[1]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2    fitness rank
# 1        1        1      0      0      0       1     0      1  0     0   1   1 101273.264    6
# 2        0        1      0      0      0       0     0      0  1     0   1   1  31233.545    2
# 3        1        1      0      1      0       0     0      0  1     0   1   0 101129.345    5
# 4        1        0      0      1      0       0     1      0  1     0   1   0  91137.312    4
# 5        0        1      0      0      0       0     0      0  0     0   0   1   2625.758    1
# 6        1        0      0      1      1       0     1      0  0     1   0   1  61089.019    3
# 
# $history$data.pop.list[[2]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2    fitness rank
# 1        1        1      1      1      1       1     1      1  1     0   1   0 141189.404    6
# 2        0        1      0      0      0       0     0      0  1     0   1   1  31233.545    2
# 3        1        1      1      1      0       0     1      1  1     0   0   1 141153.691    5
# 4        1        0      0      0      0       0     1      0  1     0   1   1  81263.824    3
# 5        1        1      1      0      0       0     0      0  1     1   0   1  82932.734    4
# 6        0        1      0      0      0       0     0      0  0     0   0   1   2625.758    1
# 
# $history$data.pop.list[[3]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2   fitness rank
# 1        0        1      1      0      0       0     0      0  1     1   0   1 12624.591    2
# 2        0        1      0      0      0       0     0      0  1     0   1   1 31233.545    3
# 3        1        0      1      1      0       0     0      0  1     0   0   1 91112.266    5
# 4        0        1      0      0      0       0     0      0  1     0   1   1 31233.545    3
# 5        1        1      0      0      1       0     0      0  1     0   0   1 92896.057    6
# 6        0        1      0      0      0       0     0      0  0     0   0   1  2625.758    1
# 
# $history$data.pop.list[[4]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2   fitness rank
# 1        0        1      0      0      0       0     0      0  0     0   0   1  2625.758    3
# 2        0        1      1      0      0       0     0      0  1     1   0   1 12624.591    4
# 3        1        0      1      0      0       0     0      0  1     0   0   1 52548.431    6
# 4        0        1      0      1      0       0     0      0  0     0   0   1  1033.996    1
# 5        1        0      0      1      0       0     0      0  0     0   0   1 51084.718    5
# 6        0        1      0      1      0       0     0      0  0     0   0   1  1033.996    1
# 
# $history$data.pop.list[[5]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2   fitness rank
# 1        1        1      0      0      0       0     0      0  0     0   0   1 72587.502    6
# 2        1        0      0      0      0       0     0      0  0     0   0   1 32467.136    4
# 3        0        1      0      1      0       0     0      0  1     0   0   1 21075.829    3
# 4        0        1      0      0      0       0     0      0  0     0   0   1  2625.758    2
# 5        1        0      0      1      0       0     0      0  0     1   0   1 51067.656    5
# 6        0        1      0      1      0       0     0      0  0     0   0   1  1033.996    1
# 
# $history$data.pop.list[[6]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2   fitness rank
# 1        0        1      0      0      0       0     0      0  0     0   0   1  2625.758    4
# 2        0        1      0      1      0       0     0      0  0     0   0   1  1033.996    1
# 3        0        1      0      0      0       0     0      0  1     0   0   1 22423.416    6
# 4        0        1      0      1      0       0     0      0  0     1   0   1  1033.476    1
# 5        0        1      0      1      0       0     0      0  1     0   0   1 21075.829    5
# 6        0        1      0      1      0       0     0      0  0     0   0   1  1033.996    1
# 
# 
# $history$sel.population.list
# $history$sel.population.list[[1]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        1        1      0      1      0       0     0      0  1     0   1   0
# [2,]        0        1      0      0      0       0     0      0  1     0   1   1
# [3,]        1        1      0      1      0       0     0      0  1     0   1   0
# [4,]        1        0      0      1      0       0     1      0  1     0   1   0
# [5,]        0        1      0      0      0       0     0      0  0     0   0   1
# [6,]        0        1      0      0      0       0     0      0  1     0   1   1
# 
# $history$sel.population.list[[2]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 5        1        1      1      0      0       0     0      0  1     1   0   1
# 2        0        1      0      0      0       0     0      0  1     0   1   1
# 4        1        0      0      0      0       0     1      0  1     0   1   1
# 2        0        1      0      0      0       0     0      0  1     0   1   1
# 5        1        1      1      0      0       0     0      0  1     1   0   1
# 6        0        1      0      0      0       0     0      0  0     0   0   1
# 
# $history$sel.population.list[[3]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 1        0        1      1      0      0       0     0      0  1     1   0   1
# 6        0        1      0      0      0       0     0      0  0     0   0   1
# 3        1        0      1      1      0       0     0      0  1     0   0   1
# 4        0        1      0      0      0       0     0      0  1     0   1   1
# 3        1        0      1      1      0       0     0      0  1     0   0   1
# 6        0        1      0      0      0       0     0      0  0     0   0   1
# 
# $history$sel.population.list[[4]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 1        0        1      0      0      0       0     0      0  0     0   0   1
# 1        0        1      0      0      0       0     0      0  0     0   0   1
# 4        0        1      0      1      0       0     0      0  0     0   0   1
# 4        0        1      0      1      0       0     0      0  0     0   0   1
# 5        1        0      0      1      0       0     0      0  0     0   0   1
# 4        0        1      0      1      0       0     0      0  0     0   0   1
# 
# $history$sel.population.list[[5]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 5        1        0      0      1      0       0     0      0  0     1   0   1
# 2        1        0      0      0      0       0     0      0  0     0   0   1
# 4        0        1      0      0      0       0     0      0  0     0   0   1
# 4        0        1      0      0      0       0     0      0  0     0   0   1
# 3        0        1      0      1      0       0     0      0  1     0   0   1
# 6        0        1      0      1      0       0     0      0  0     0   0   1
# 
# $history$sel.population.list[[6]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 1        0        1      0      0      0       0     0      0  0     0   0   1
# 2        0        1      0      1      0       0     0      0  0     0   0   1
# 6        0        1      0      1      0       0     0      0  0     0   0   1
# 4        0        1      0      1      0       0     0      0  0     1   0   1
# 4        0        1      0      1      0       0     0      0  0     1   0   1
# 2        0        1      0      1      0       0     0      0  0     0   0   1
# 
# 
# $history$ls.population.list
# $history$ls.population.list[[1]]
# NULL
# 
# $history$ls.population.list[[2]]
# NULL
# 
# $history$ls.population.list[[3]]
# V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12   fitness rank
# ls.population.s     1  1  0  0  0  0  0  0  0   0   0   1 72587.502    7
# ls.population.s.1   1  0  0  0  0  0  0  0  0   0   0   1 32467.136    6
# ls.population.s.2   0  1  0  0  0  0  0  0  0   0   0   1  2625.758    3
# ls.population.s.3   0  1  0  1  0  0  0  0  0   0   0   1  1033.996    1
# ls.population.s.8   0  1  0  0  0  0  0  0  1   0   0   1 22423.416    5
# ls.population.s.10  0  1  0  0  0  0  0  0  0   0   1   1 11180.126    4
# ls.population.s.11  0  1  0  0  0  0  0  0  0   0   1   0  1883.545    2
# 
# $history$ls.population.list[[4]]
# NULL
# 
# $history$ls.population.list[[5]]
# NULL
# 
# $history$ls.population.list[[6]]
# V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12   fitness rank
# ls.population.s     1  1  0  1  0  0  0  0  0   0   0   1 81062.130    8
# ls.population.s.1   1  0  0  1  0  0  0  0  0   0   0   1 51084.718    7
# ls.population.s.2   0  1  0  1  0  0  0  0  0   0   0   1  1033.996    1
# ls.population.s.3   0  1  0  0  0  0  0  0  0   0   0   1  2625.758    5
# ls.population.s.8   0  1  0  1  0  0  0  0  1   0   0   1 21075.829    6
# ls.population.s.9   0  1  0  1  0  0  0  0  0   1   0   1  1033.476    1
# ls.population.s.10  0  1  0  1  0  0  0  0  0   0   1   1  1176.734    4
# ls.population.s.11  0  1  0  1  0  0  0  0  0   0   1   0  1034.986    1
# 
# 
# $history$children.cross.list
# $history$children.cross.list[[1]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        1        1      0      1      0       0     0      0  1     0   1   0
# [2,]        0        1      0      0      0       0     0      0  1     0   1   1
# [3,]        1        1      0      1      0       0     0      0  1     0   1   0
# [4,]        1        0      0      1      0       0     1      0  1     0   1   0
# [5,]        0        1      0      0      0       0     0      0  0     0   0   1
# [6,]        0        1      0      0      0       0     0      0  1     0   1   1
# 
# $history$children.cross.list[[2]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        1        1      1      0      0       0     0      0  1     1   0   1
# [2,]        0        1      0      0      0       0     0      0  1     0   1   1
# [3,]        1        0      0      0      0       0     0      0  1     0   1   1
# [4,]        0        1      0      0      0       0     1      0  1     0   1   1
# [5,]        1        1      0      0      0       0     0      0  1     1   0   1
# [6,]        0        1      1      0      0       0     0      0  0     0   0   1
# 
# $history$children.cross.list[[3]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        0        1      1      0      0       0     0      0  0     0   0   1
# [2,]        0        1      0      0      0       0     0      0  1     1   0   1
# [3,]        1        0      1      0      0       0     0      0  1     0   0   1
# [4,]        0        1      0      1      0       0     0      0  1     0   1   1
# [5,]        1        0      1      1      0       0     0      0  0     0   0   1
# [6,]        0        1      0      0      0       0     0      0  1     0   0   1
# 
# $history$children.cross.list[[4]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 0        1      0      0      0       0     0      0  0     0   0   1
# 0        1      0      0      0       0     0      0  0     0   0   1
# 4        0        1      0      1      0       0     0      0  0     0   0   1
# 4        0        1      0      1      0       0     0      0  0     0   0   1
# 1        0      0      1      0       0     0      0  0     0   0   1
# 0        1      0      1      0       0     0      0  0     0   0   1
# 
# $history$children.cross.list[[5]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 1        0      0      0      0       0     0      0  0     0   0   1
# 1        0      0      1      0       0     0      0  0     1   0   1
# 4        0        1      0      0      0       0     0      0  0     0   0   1
# 4        0        1      0      0      0       0     0      0  0     0   0   1
# 0        1      0      1      0       0     0      0  1     0   0   1
# 0        1      0      1      0       0     0      0  0     0   0   1
# 
# $history$children.cross.list[[6]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        0        1      0      1      0       0     0      0  0     0   0   1
# [2,]        0        1      0      0      0       0     0      0  0     0   0   1
# [3,]        0        1      0      1      0       0     0      0  0     0   0   1
# [4,]        0        1      0      1      0       0     0      0  0     1   0   1
# [5,]        0        1      0      1      0       0     0      0  0     0   0   1
# [6,]        0        1      0      1      0       0     0      0  0     1   0   1
# 
# 
# $history$children.mutation.list
# $history$children.mutation.list[[1]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        1        1      1      1      1       1     1      1  1     0   1   0
# [2,]        0        1      0      0      0       0     0      0  1     0   1   1
# [3,]        1        1      1      1      0       0     1      1  1     0   0   1
# [4,]        1        0      0      0      0       0     1      0  1     0   1   1
# [5,]        1        1      1      0      0       0     0      0  1     1   0   1
# [6,]        0        1      0      0      0       0     0      0  1     0   1   0
# 
# $history$children.mutation.list[[2]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        0        1      1      0      0       0     0      0  1     1   0   1
# [2,]        0        0      0      0      0       0     0      0  1     0   1   1
# [3,]        1        0      1      1      0       1     0      0  1     0   0   0
# [4,]        0        1      0      0      0       0     0      1  1     0   1   1
# [5,]        1        1      0      0      1       0     0      0  1     0   0   1
# [6,]        1        1      0      0      1       0     0      0  0     0   0   1
# 
# $history$children.mutation.list[[3]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        0        1      0      0      1       1     0      0  0     0   0   1
# [2,]        0        1      1      0      0       0     0      1  1     1   0   1
# [3,]        1        0      1      0      0       0     0      0  1     0   0   0
# [4,]        0        1      0      1      0       0     0      0  0     0   0   1
# [5,]        1        0      0      1      0       0     0      1  0     0   0   1
# [6,]        0        1      0      0      0       0     0      1  1     0   0   1
# 
# $history$children.mutation.list[[4]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 1        1      0      0      0       0     0      0  0     0   0   1
# 1        0      1      0      0       0     0      0  0     1   0   1
# 4        0        1      0      1      0       0     0      0  1     0   0   1
# 4        0        1      0      0      0       0     0      0  0     1   0   1
# 1        0      0      1      0       1     0      0  0     1   0   1
# 1        1      0      1      0       0     0      0  0     0   1   1
# 
# $history$children.mutation.list[[5]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 0        0      0      0      0       0     0      0  0     0   0   1
# 0        1      1      1      0       0     0      0  0     0   0   1
# 4        0        1      0      0      0       0     0      0  1     0   0   1
# 4        0        1      0      1      0       0     0      0  0     1   0   1
# 0        1      0      1      0       1     0      0  1     1   0   1
# 1        1      0      1      0       0     0      0  0     0   1   1
# 
# $history$children.mutation.list[[6]]
# cmpt.iv1 cmpt.iv2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# [1,]        0        1      0      1      0       0     0      0  0     1   0   1
# [2,]        0        1      1      0      0       1     1      0  0     0   1   0
# [3,]        0        1      1      1      0       0     0      0  0     0   0   1
# [4,]        0        1      0      1      0       0     0      0  0     0   0   1
# [5,]        0        0      0      1      0       0     0      1  1     1   0   1
# [6,]        0        1      0      1      0       1     0      0  0     1   1   1
```
ACO example
``` r
aco.operator(dat = d1,
             no.ants<-6,
             max.iter<-6,
             search.space = 1,
             no.cores = getRxThreads(),
             thetalower=c(vp=1,
                          vp2=1),
             control = saemControl(
               seed = 1234,
               print = 5,
               nBurn = 50,
               nEm = 50,
               logLik = T),
             table=tableControl(cwres=T),
             filename =  "pheno_sd",
             foldername =   "pheno_sd" )
# [1] "Output directory for ACO analysis is created"
# Running ACO ---------------------------------------------- Model: 1 Iteration: 1
# Running ACO ---------------------------------------------- Model: 2 Iteration: 1
# Running ACO ---------------------------------------------- Model: 3 Iteration: 1
# Running ACO ---------------------------------------------- Model: 4 Iteration: 1
# Running ACO ---------------------------------------------- Model: 5 Iteration: 1
# Running ACO ---------------------------------------------- Model: 6 Iteration: 1
# Running ACO ---------------------------------------------- Model: 7 Iteration: 2
# Running ACO ---------------------------------------------- Model: 8 Iteration: 2
# Running ACO ---------------------------------------------- Model: 9 Iteration: 2
# Running ACO ---------------------------------------------- Model: 10 Iteration: 2
# Running ACO ---------------------------------------------- Model: 11 Iteration: 2
# Running ACO ---------------------------------------------- Model: 12 Iteration: 2
# Running ACO ---------------------------------------------- Model: 13 Iteration: 3
# Running ACO ---------------------------------------------- Model: 14 Iteration: 3
# Running ACO ---------------------------------------------- Model: 15 Iteration: 3
# Running ACO ---------------------------------------------- Model: 16 Iteration: 3
# Running ACO ---------------------------------------------- Model: 17 Iteration: 3
# Running ACO ---------------------------------------------- Model: 18 Iteration: 3
# Running ACO ---------------------------------------------- Model: 19 Iteration: 4
# Running ACO ---------------------------------------------- Model: 20 Iteration: 4
# Running ACO ---------------------------------------------- Model: 21 Iteration: 4
# Running ACO ---------------------------------------------- Model: 22 Iteration: 4
# Running ACO ---------------------------------------------- Model: 23 Iteration: 4
# Running ACO ---------------------------------------------- Model: 24 Iteration: 4
# Running ACO ---------------------------------------------- Model: 25 Iteration: 5
# Running ACO ---------------------------------------------- Model: 26 Iteration: 5
# Running ACO ---------------------------------------------- Model: 27 Iteration: 5
# Running ACO ---------------------------------------------- Model: 28 Iteration: 5
# Running ACO ---------------------------------------------- Model: 29 Iteration: 5
# Running ACO ---------------------------------------------- Model: 30 Iteration: 5
# Running ACO ---------------------------------------------- Model: 31 Iteration: 6
# Running ACO ---------------------------------------------- Model: 32 Iteration: 6
# Running ACO ---------------------------------------------- Model: 33 Iteration: 6
# Running ACO ---------------------------------------------- Model: 34 Iteration: 6
# Running ACO ---------------------------------------------- Model: 35 Iteration: 6
# Running ACO ---------------------------------------------- Model: 36 Iteration: 6
# $bestmodelcode
# cmpt.iv eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 14       1      0      1      0       0     0      0  0     1  1
# 
# $best_model_name
# [1] "1Cmpt,IIV.cl.vc,first-order_elminiation,full_omega_matrix,additive"
# 
# $history.list
# $history.list$node.list.all
# travel node.no local.node.no  node.names node.group      phi  delta_phi         p
# 1        0       1             1       1Cmpt          1 1.000000 0.00000000 0.3330000
# 2        0       2             2       2Cmpt          1 1.000000 0.00000000 0.3330000
# 3        0       3             3       3Cmpt          1 1.000000 0.00000000 0.3330000
# 4        0       4             0  eta.vp2.no          2 1.000000 0.00000000 0.5000000
# 5        0       5             1 eta.vp2.yes          2 1.000000 0.00000000 0.5000000
# 6        0       6             0   eta.q2.no          3 1.000000 0.00000000 0.5000000
# 7        0       7             1  eta.q2.yes          3 1.000000 0.00000000 0.5000000
# 8        0       8             0   eta.vp.no          4 1.000000 0.00000000 0.5000000
# 9        0       9             1  eta.vp.yes          4 1.000000 0.00000000 0.5000000
# 10       0      10             0    eta.q.no          5 1.000000 0.00000000 0.5000000
# 11       0      11             1   eta.q.yes          5 1.000000 0.00000000 0.5000000
# 12       0      12             0   eta.vc.no          6 1.000000 0.00000000 0.5000000
# 13       0      13             1  eta.vc.yes          6 1.000000 0.00000000 0.5000000
# 14       0      14             0       mm.no          7 1.000000 0.00000000 0.5000000
# 15       0      15             1      mm.yes          7 1.000000 0.00000000 0.5000000
# 16       0      16             0   eta.km.no          8 1.000000 0.00000000 0.5000000
# 17       0      17             1  eta.km.yes          8 1.000000 0.00000000 0.5000000
# 18       0      18             0    mcorr.no          9 1.000000 0.00000000 0.5000000
# 19       0      19             1   mcorr.yes          9 1.000000 0.00000000 0.5000000
# 20       0      20             1         add         10 1.000000 0.00000000 0.3330000
# 21       0      21             2        prop         10 1.000000 0.00000000 0.3330000
# 22       0      22             3        comb         10 1.000000 0.00000000 0.3330000
# 23       1       1             1       1Cmpt          1 2.300000 1.50000000 0.4742268
# 24       1       2             2       2Cmpt          1 1.383333 0.58333333 0.2852234
# 25       1       3             3       3Cmpt          1 1.166667 0.36666667 0.2500000
# 26       1       4             0  eta.vp2.no          2 1.166667 0.36666667 0.5384615
# 27       1       5             1 eta.vp2.yes          2 1.000000 0.00000000 0.4615385
# 28       1       6             0   eta.q2.no          3 1.166667 0.36666667 0.5384615
# 29       1       7             1  eta.q2.yes          3 1.000000 0.00000000 0.4615385
# 30       1       8             0   eta.vp.no          4 1.750000 0.95000000 0.6363636
# 31       1       9             1  eta.vp.yes          4 1.000000 0.00000000 0.3636364
# 32       1      10             0    eta.q.no          5 1.750000 0.95000000 0.6363636
# 33       1      11             1   eta.q.yes          5 1.000000 0.00000000 0.3636364
# 34       1      12             0   eta.vc.no          6 3.250000 2.45000000 0.7647059
# 35       1      13             1  eta.vc.yes          6 1.000000 0.00000000 0.2352941
# 36       1      14             0       mm.no          7 3.250000 2.45000000 0.7647059
# 37       1      15             1      mm.yes          7 1.000000 0.00000000 0.2352941
# 38       1      16             0   eta.km.no          8 1.000000 0.00000000 0.5000000
# 39       1      17             1  eta.km.yes          8 1.000000 0.00000000 0.5000000
# 40       1      18             0    mcorr.no          9 3.250000 2.45000000 0.7647059
# 41       1      19             1   mcorr.yes          9 1.000000 0.00000000 0.2352941
# 42       1      20             1         add         10 1.716667 0.91666667 0.3399340
# 43       1      21             2        prop         10 2.333333 1.53333333 0.4620462
# 44       1      22             3        comb         10 1.000000 0.00000000 0.2500000
# 45       2       1             1       1Cmpt          1 4.506667 2.66666667 0.6671601
# 46       2       2             2       2Cmpt          1 1.231667 0.12500000 0.2500000
# 47       2       3             3       3Cmpt          1 1.016667 0.08333333 0.2500000
# 48       2       4             0  eta.vp2.no          2 1.683333 0.75000000 0.6273292
# 49       2       5             1 eta.vp2.yes          2 1.000000 0.00000000 0.3726708
# 50       2       6             0   eta.q2.no          3 1.683333 0.75000000 0.6273292
# 51       2       7             1  eta.q2.yes          3 1.000000 0.00000000 0.3726708
# 52       2       8             0   eta.vp.no          4 2.191667 0.79166667 0.6866841
# 53       2       9             1  eta.vp.yes          4 1.000000 0.08333333 0.3133159
# 54       2      10             0    eta.q.no          5 2.275000 0.87500000 0.6946565
# 55       2      11             1   eta.q.yes          5 1.000000 0.00000000 0.3053435
# 56       2      12             0   eta.vc.no          6 3.391667 0.79166667 0.5405046
# 57       2      13             1  eta.vc.yes          6 2.883333 2.08333333 0.4594954
# 58       2      14             0       mm.no          7 5.475000 2.87500000 0.8000000
# 59       2      15             1      mm.yes          7 1.000000 0.00000000 0.2000000
# 60       2      16             0   eta.km.no          8 1.466667 0.66666667 0.5945946
# 61       2      17             1  eta.km.yes          8 1.000000 0.00000000 0.4054054
# 62       2      18             0    mcorr.no          9 5.475000 2.87500000 0.8000000
# 63       2      19             1   mcorr.yes          9 1.000000 0.00000000 0.2000000
# 64       2      20             1         add         10 3.498333 2.12500000 0.4916842
# 65       2      21             2        prop         10 2.616667 0.75000000 0.3677676
# 66       2      22             3        comb         10 1.000000 0.00000000 0.2500000
# 67       3       1             1       1Cmpt          1 5.796242 2.19090909 0.7385402
# 68       3       2             2       2Cmpt          1 1.052000 0.06666667 0.2500000
# 69       3       3             3       3Cmpt          1 1.000000 0.00000000 0.2500000
# 70       3       4             0  eta.vp2.no          2 2.346667 1.00000000 0.7011952
# 71       3       5             1 eta.vp2.yes          2 1.000000 0.00000000 0.2988048
# 72       3       6             0   eta.q2.no          3 2.346667 1.00000000 0.7011952
# 73       3       7             1  eta.q2.yes          3 1.000000 0.00000000 0.2988048
# 74       3       8             0   eta.vp.no          4 2.820000 1.06666667 0.7382199
# 75       3       9             1  eta.vp.yes          4 1.000000 0.00000000 0.2617801
# 76       3      10             0    eta.q.no          5 2.886667 1.06666667 0.7427101
# 77       3      11             1   eta.q.yes          5 1.000000 0.00000000 0.2572899
# 78       3      12             0   eta.vc.no          6 2.904242 0.19090909 0.3990673
# 79       3      13             1  eta.vc.yes          6 4.373333 2.06666667 0.6009327
# 80       3      14             0       mm.no          7 6.546667 2.16666667 0.8000000
# 81       3      15             1      mm.yes          7 1.000000 0.09090909 0.2000000
# 82       3      16             0   eta.km.no          8 2.264242 1.09090909 0.6936502
# 83       3      17             1  eta.km.yes          8 1.000000 0.00000000 0.3063498
# 84       3      18             0    mcorr.no          9 5.637576 1.25757576 0.7579857
# 85       3      19             1   mcorr.yes          9 1.800000 1.00000000 0.2420143
# 86       3      20             1         add         10 4.798667 2.00000000 0.6029486
# 87       3      21             2        prop         10 2.160000 0.06666667 0.2714022
# 88       3      22             3        comb         10 1.000000 0.19090909 0.2500000
# 89       4       1             1       1Cmpt          1 6.761994 2.12500000 0.7500000
# 90       4       2             2       2Cmpt          1 1.000000 0.05263158 0.2500000
# 91       4       3             3       3Cmpt          1 1.000000 0.09166667 0.2500000
# 92       4       4             0  eta.vp2.no          2 3.927333 2.05000000 0.7970505
# 93       4       5             1 eta.vp2.yes          2 1.000000 0.04166667 0.2029495
# 94       4       6             0   eta.q2.no          3 3.969000 2.09166667 0.7987523
# 95       4       7             1  eta.q2.yes          3 1.000000 0.00000000 0.2012477
# 96       4       8             0   eta.vp.no          4 4.400298 2.14429825 0.8000000
# 97       4       9             1  eta.vp.yes          4 1.000000 0.00000000 0.2000000
# 98       4      10             0    eta.q.no          5 4.401000 2.09166667 0.8000000
# 99       4      11             1   eta.q.yes          5 1.000000 0.05263158 0.2000000
# 100      4      12             0   eta.vc.no          6 2.501026 0.17763158 0.3090983
# 101      4      13             1  eta.vc.yes          6 5.590333 2.09166667 0.6909017
# 102      4      14             0       mm.no          7 7.506632 2.26929825 0.8000000
# 103      4      15             1      mm.yes          7 1.000000 0.00000000 0.2000000
# 104      4      16             0   eta.km.no          8 3.811394 2.00000000 0.7921600
# 105      4      17             1  eta.km.yes          8 1.000000 0.00000000 0.2078400
# 106      4      18             0    mcorr.no          9 4.779359 0.26929825 0.5814759
# 107      4      19             1   mcorr.yes          9 3.440000 2.00000000 0.4185241
# 108      4      20             1         add         10 5.941565 2.10263158 0.6724094
# 109      4      21             2        prop         10 1.894667 0.16666667 0.2500000
# 110      4      22             3        comb         10 1.000000 0.00000000 0.2500000
# 111      5       1             1       1Cmpt          1 7.476262 2.06666667 0.7500000
# 112      5       2             2       2Cmpt          1 1.000000 0.04545455 0.2500000
# 113      5       3             3       3Cmpt          1 1.000000 0.07500000 0.2500000
# 114      5       4             0  eta.vp2.no          2 5.183533 2.04166667 0.8000000
# 115      5       5             1 eta.vp2.yes          2 1.000000 0.03333333 0.2000000
# 116      5       6             0   eta.q2.no          3 5.250200 2.07500000 0.8000000
# 117      5       7             1  eta.q2.yes          3 1.000000 0.00000000 0.2000000
# 118      5       8             0   eta.vp.no          4 5.640693 2.12045455 0.8000000
# 119      5       9             1  eta.vp.yes          4 1.000000 0.00000000 0.2000000
# 120      5      10             0    eta.q.no          5 5.595800 2.07500000 0.8000000
# 121      5      11             1   eta.q.yes          5 1.000000 0.04545455 0.2000000
# 122      5      12             0   eta.vc.no          6 2.112942 0.11212121 0.2439828
# 123      5      13             1  eta.vc.yes          6 6.547267 2.07500000 0.7560172
# 124      5      14             0       mm.no          7 8.192426 2.18712121 0.8000000
# 125      5      15             1      mm.yes          7 1.000000 0.00000000 0.2000000
# [ reached 'max' / getOption("max.print") -- omitted 29 rows ]
# 
# $history.list$ant.travel.list
# $history.list$ant.travel.list[[1]]
# X1 X2 X3 X4 X5 X6
# 1   1  1  2  2  3  3
# 2  -1 -1 -1 -1  0  0
# 3  -1 -1 -1 -1  0  0
# 4  -1 -1  0  0  0  0
# 5  -1 -1  0  0  0  0
# 6   0  0  0  0  0  0
# 7   0  0  0  0  0  0
# 8  -1 -1 -1 -1 -1 -1
# 9   0  0  0  0  0  0
# 10  1  2  1  2  1  2
# 
# $history.list$ant.travel.list[[2]]
# X1 X2 X3 X4 X5 X6
# 1   3  2  1  1  2  1
# 2   0 -1 -1 -1 -1 -1
# 3   0 -1 -1 -1 -1 -1
# 4   1  0 -1 -1  0 -1
# 5   0  0 -1 -1  0 -1
# 6   1  0  1  1  1  0
# 7   0  0  0  0  0  1
# 8  -1 -1 -1 -1 -1  0
# 9   0  0  0  0  0  1
# 10  2  1  1  1  3  2
# 
# $history.list$ant.travel.list[[3]]
# X1 X2 X3 X4 X5 X6
# 1   1  1  1  2  2  1
# 2  -1 -1 -1 -1 -1 -1
# 3  -1 -1 -1 -1 -1 -1
# 4  -1 -1 -1  0  0 -1
# 5  -1 -1 -1  0  0 -1
# 6   0  0  1  1  0  0
# 7   1  0  0  0  0  0
# 8   0 -1 -1 -1 -1 -1
# 9   0  0  1  0  0  1
# 10  3  3  1  2  1  2
# 
# $history.list$ant.travel.list[[4]]
# X1 X2 X3 X4 X5 X6
# 1   3  3  1  2  1  2
# 2   0  1 -1 -1 -1 -1
# 3   0  0 -1 -1 -1 -1
# 4   0  0 -1  0 -1  1
# 5   0  0 -1  1 -1  0
# 6   1  1  0  0  0  1
# 7   0  0  0  0  0  0
# 8  -1 -1 -1 -1 -1 -1
# 9   0  0  0  0  0  0
# 10  1  2  2  1  3  3
# 
# $history.list$ant.travel.list[[5]]
# X1 X2 X3 X4 X5 X6
# 1   3  3  1  2  1  2
# 2   0  1 -1 -1 -1 -1
# 3   0  0 -1 -1 -1 -1
# 4   0  0 -1  0 -1  1
# 5   0  0 -1  1 -1  0
# 6   1  1  0  0  0  1
# 7   0  0  0  0  0  0
# 8  -1 -1 -1 -1 -1 -1
# 9   1  0  1  0  0  0
# 10  1  3  3  1  2  2
# 
# $history.list$ant.travel.list[[6]]
# X1 X2 X3 X4 X5 X6
# 1   1  1  1  1  2  1
# 2  -1 -1 -1 -1 -1 -1
# 3  -1 -1 -1 -1 -1 -1
# 4  -1 -1 -1 -1  0 -1
# 5  -1 -1 -1 -1  1 -1
# 6   1  0  1  0  0  1
# 7   0  0  0  0  0  0
# 8  -1 -1 -1 -1 -1 -1
# 9   0  1  1  0  1  1
# 10  2  1  1  3  1  1
``` 


Tabu example
``` r
tabu.operator(dat = d1,
              no.cores = getRxThreads(),
              tabu.duration=2,
              max.round=6,
              search.space = 1,
              thetalower=c(vp=1,
                           vp2=1),
              control = saemControl(
                seed = 1234,
                print = 5,
                nBurn = 50,
                nEm = 50,
                logLik = T),
              table=tableControl(cwres=T),
              filename =  "pheno_sd",
              foldername =   "pheno_sd" )
# 
# [1] "Output directory for Tabu analysis is created"
# Running Tabu ---------------------------------------------- Model: 1 Iteration: 0
# Running Tabu ---------------------------------------------- Model: 2 Iteration: 1
# Running Tabu ---------------------------------------------- Model: 3 Iteration: 1
# Running Tabu ---------------------------------------------- Model: 4 Iteration: 1
# Running Tabu ---------------------------------------------- Model: 5 Iteration: 1
# Running Tabu ---------------------------------------------- Model: 6 Iteration: 1
# Running Tabu ---------------------------------------------- Model: 7 Iteration: 1
# Running Tabu ---------------------------------------------- Model: 8 Iteration: 1
# Running Tabu ---------------------------------------------- Model: 9 Iteration: 1
# Running Tabu ---------------------------------------------- Model: 10 Iteration: 1
# Running Tabu ---------------------------------------------- Model: 11 Iteration: 2
# Running Tabu ---------------------------------------------- Model: 12 Iteration: 2
# Running Tabu ---------------------------------------------- Model: 13 Iteration: 2
# Running Tabu ---------------------------------------------- Model: 14 Iteration: 2
# Running Tabu ---------------------------------------------- Model: 15 Iteration: 2
# Running Tabu ---------------------------------------------- Model: 16 Iteration: 3
# Running Tabu ---------------------------------------------- Model: 17 Iteration: 3
# Running Tabu ---------------------------------------------- Model: 18 Iteration: 4
# Running Tabu ---------------------------------------------- Model: 19 Iteration: 4
# Running Tabu ---------------------------------------------- Model: 20 Iteration: 5
# Running Tabu ---------------------------------------------- Model: 21 Iteration: 5
# Running Tabu ---------------------------------------------- Model: 22 Iteration: 5
# Running Tabu ---------------------------------------------- Model: 23 Iteration: 5
# Running Tabu ---------------------------------------------- Model: 24 Iteration: 5
# Running Tabu ---------------------------------------------- Model: 25 Iteration: 6
# Running Tabu ---------------------------------------------- Model: 26 Iteration: 6
# Running Tabu ---------------------------------------------- Model: 27 Iteration: 6
# Running Tabu ---------------------------------------------- Model: 28 Iteration: 6
# Running Tabu ---------------------------------------------- Model: 29 Iteration: 6
# $best_model_name
# [1] "1Cmpt,IIV.cl.vc,first-order_elminiation,full_omega_matrix,additive"
# 
# $best_model_code
# [1] 1 0 1 0 0 0 0 0 1 1
# 
# $history
# $history$starting.points.history.list
# $history$starting.points.history.list[[1]]
# cmpt.iv  eta.km  eta.vc  eta.vp eta.vp2   eta.q  eta.q2      mm   mcorr      rv 
# 2       0       1       0       0       0       0       0       0       1 
# 
# $history$starting.points.history.list[[2]]
# cmpt.iv  eta.km  eta.vc  eta.vp eta.vp2   eta.q  eta.q2      mm   mcorr      rv 
# 1       0       1       0       0       0       0       0       0       1 
# 
# $history$starting.points.history.list[[3]]
# cmpt.iv  eta.km  eta.vc  eta.vp eta.vp2   eta.q  eta.q2      mm   mcorr      rv 
# 1       0       1       0       0       0       0       0       1       1 
# 
# $history$starting.points.history.list[[4]]
# cmpt.iv  eta.km  eta.vc  eta.vp eta.vp2   eta.q  eta.q2      mm   mcorr      rv 
# 1       0       1       0       0       0       0       0       1       2 
# 
# $history$starting.points.history.list[[5]]
# cmpt.iv  eta.km  eta.vc  eta.vp eta.vp2   eta.q  eta.q2      mm   mcorr      rv 
# 2       0       1       0       0       0       0       0       1       2 
# 
# $history$starting.points.history.list[[6]]
# cmpt.iv  eta.km  eta.vc  eta.vp eta.vp2   eta.q  eta.q2      mm   mcorr      rv 
# 2       0       0       0       0       0       0       0       0       2 
# 
# 
# $history$tabu.elements.history.list
# $history$tabu.elements.history.list[[1]]
# NULL
# 
# $history$tabu.elements.history.list[[2]]
# tabu.num elements elements.value tabu.iteration.left
# 1        3  cmpt.iv              1                   2
# 
# $history$tabu.elements.history.list[[3]]
# tabu.num elements elements.value tabu.iteration.left
# 1        3  cmpt.iv              1                   1
# 2        4    mcorr              1                   2
# 
# $history$tabu.elements.history.list[[4]]
# tabu.num elements elements.value tabu.iteration.left
# 2        4    mcorr              1                   1
# 3        5       rv              2                   2
# 
# $history$tabu.elements.history.list[[5]]
# tabu.num elements elements.value tabu.iteration.left
# 3        5       rv              2                   1
# 1        6  cmpt.iv              2                   2
# 
# $history$tabu.elements.history.list[[6]]
# tabu.num elements elements.value tabu.iteration.left
# 1         6  cmpt.iv              2                   1
# 11        7   eta.vc              0                   2
# 2         7    mcorr              0                   2
# 
# 
# $history$tried.neighbors
# $history$tried.neighbors[[1]]
# cmpt.iv eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 1       1      0      1      0       0     0      0  0     0  1
# 2       3      0      1      0       0     0      0  0     0  1
# 3       2      0      0      0       0     0      0  0     0  1
# 4       2      0      1      1       0     0      0  0     0  1
# 5       2      0      1      0       0     1      0  0     0  1
# 6       2      0      1      0       0     0      0  1     0  1
# 7       2      0      1      0       0     0      0  0     1  1
# 8       2      0      1      0       0     0      0  0     0  2
# 9       2      0      1      0       0     0      0  0     0  3
# 
# $history$tried.neighbors[[2]]
# cmpt.iv eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 1       1      0      0      0       0     0      0  0     0  1
# 2       1      0      1      0       0     0      0  1     0  1
# 3       1      0      1      0       0     0      0  0     1  1
# 4       1      0      1      0       0     0      0  0     0  2
# 5       1      0      1      0       0     0      0  0     0  3
# 
# $history$tried.neighbors[[3]]
# cmpt.iv eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 1       1      0      1      0       0     0      0  0     1  2
# 2       1      0      1      0       0     0      0  0     1  3
# 
# $history$tried.neighbors[[4]]
# cmpt.iv eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 1       2      0      1      0       0     0      0  0     1  2
# 2       3      0      1      0       0     0      0  0     1  2
# 
# $history$tried.neighbors[[5]]
# cmpt.iv eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 1       2      0      0      0       0     0      0  0     0  2
# 2       2      0      1      1       0     0      0  0     1  2
# 3       2      0      1      0       0     1      0  0     1  2
# 4       2      0      1      0       0     0      0  1     0  2
# 5       2      0      1      0       0     0      0  0     0  2
# 
# $history$tried.neighbors[[6]]
# cmpt.iv eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 1       2      0      0      1       0     0      0  0     0  2
# 2       2      0      0      0       0     1      0  0     0  2
# 3       2      0      0      0       0     0      0  1     0  2
# 4       2      0      0      0       0     0      0  0     0  1
# 5       2      0      0      0       0     0      0  0     0  3

```


