# nlmixr2auto
Automated Population Pharmacokinetic Modelling. This package offers serveral optimisation algorithms designed for automated population pharmacokinetic modeling, serving as a valuable tool for pharmacokinetic model selection.

## Installation of nlmixr2autoinit
```r
install.packages("nlmixr2autoinit", dependencies = TRUE)
```

## Installation of nlmixr2auto
``` r
library(devtools)
install_github("ucl-pharmacometrics/nlmixr2auto")
```

## Examples
Stepwsie example
``` r
library(nlmixr2)
library(nlmixr2autoinit)
library(nlmixr2auto)

outs <- sf.operator(
  dat = pheno_sd,                    # Dataset used for model fitting
  search.space = "ivbase",           # Structural search space for IV PK models
  filename = "pheno_sd",             # Prefix for output files
  foldername = "pheno_sd",           # Folder where results will be stored
  saem.control = saemControl(        # SAEM estimation control settings
    seed = 1234,                     # Random seed
    nBurn = 200,                     # SAEM burn-in iterations 
    nEm   = 300,                     # SAEM EM-phase iterations 
    rxControl = rxControl(cores = 4),# CPU cores for ODE solving
    logLik    = TRUE                 # Compute log-likelihood
  ),
  table.control = tableControl(
    cwres = TRUE                     # Compute conditional weighted residuals (CWRES)
  ),
  max_wall_time = 2 * 60 * 60        # Maximum allowed wall-clock time (seconds) per model; here: 2 hours
)

print(outs)

# Infometrics                               Value          
# ----------------------------------------  ---------------
#   Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Number of Subjects                        59             
# Number of Observations                    155            
# Subjects with First-Dose Interval Data    35             
# Observations in the First-Dose Interval   35             
# Subjects with Multiple-Dose Data          56             
# Observations after Multiple Doses         120            
# ----------------------------------------  ------
#   Estimating half-life....................
# Half-life estimation complete: Estimated t1/2 = 16.44 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 0.0087, estimated Vd: 1.25 
# Run parameter sweeping on nonlinear elimination kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters....................
# Running Step 1: Number of compartments ----------------------------------------------------
#   [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod1.txt
# SAEM control (core) = niter=200|300; nBurn=200; nEm=300; seed=1234; print=1
# [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod2.txt
# [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod3.txt
# Running Step 2: Elimination type ----------------------------------------------------
#   [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod4.txt
# [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod5.txt
# Running Step 3: IIV on Km ----------------------------------------------------
#   Step 3: IIV on Km skipped.
# Running Step 5: IIV (forward selection) ----------------------------------------------------
#   [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod6.txt
# [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod7.txt
# Running Step 6: ETA correlation ----------------------------------------------------
#   [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod8.txt
# [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod9.txt
# Running Step 7: Residual error model ----------------------------------------------------
#   [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod10.txt
# [Success] Model file created:
#   /home/zhonghuihuang/pheno_sd/mod11.txt
# [Success] Model file created:
  # /home/zhonghuihuang/pheno_sd/mod12.txt
# > print(outs)
# === Best Model Code ===
#   no.cmpt eta.vmax eta.km eta.cl eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 1       1        0      0      1      1      0       0     0      0  0     1  3
# 
# === Best Model Name ===
#   bolus_1cmpt_etaCLVC_FOelim_correlated_combined 
# 
# === Stepwise Selection History ===
#   Step                                                Penalty.terms                                       Model.name          Model.code  Fitness      AIC      BIC      OFV
# 1  No. of compartments                                       rse, theta, covariance   bolus_1cmpt_etaCL_FOelim_uncorrelated_combined 1,0,0,0,0,0,0,0,0,3 11182.26 1167.041 1182.258 872.1700
# 2     Elimination type                                       rse, theta, covariance   bolus_1cmpt_etaCL_FOelim_uncorrelated_combined 1,0,0,0,0,0,0,0,0,3 11182.26 1167.041 1182.258 872.1700
# 5        IIV (forward)                     rse, theta, covariance, shrinkage, omega bolus_1cmpt_etaCLVC_FOelim_uncorrelated_combined 1,0,1,0,0,0,0,0,0,3 11094.81 1076.548 1094.808 779.6768
# 6      Eta correlation        rse, theta, covariance, shrinkage, omega, correlation   bolus_1cmpt_etaCLVC_FOelim_correlated_combined 1,0,1,0,0,0,0,0,1,3  1086.86 1055.556 1076.860 756.6853
# 7 Residual error types rse, theta, covariance, shrinkage, omega, correlation, sigma   bolus_1cmpt_etaCLVC_FOelim_correlated_combined 1,0,1,0,0,0,0,0,1,3  1086.86 1055.556 1076.860 756.6853

GA example
``` r
outs<-ga.operator(dat=pheno_sd,
                  search.space = "ivbase",
                  filename =  "pheno_sd",
                  foldername =   "pheno_sd" )
# Infometrics                               Value          
# ----------------------------------------  ---------------
# Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Total Number of Subjects                  59             
# Total Number of Observations              155            
# Subjects with First-Dose Interval Data    35             
# Observations in the First-Dose Interval   35             
# Subjects with Multiple-Dose Data          56             
# Observations after Multiple Doses         120            
# ----------------------------------------  ------
# Estimating half-life....................
# Half-life estimation complete: Estimated t½ = 16.44 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# (hybrid mode: parameters combined across sources)....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 0.0087, estimated Vd: 1.25 
# Run parameter sweeping on nonlinear eliminiation kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters.................... 
# [Success] Model file created in current working directory:                      
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod1.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod2.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod3.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod4.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod5.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod6.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod7.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod8.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod9.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod10.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod11.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod12.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod13.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod14.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod15.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod16.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod17.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod18.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod19.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod20.txt
#  GA Search [==>---------------------]  13% (iteration 2/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod21.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod22.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod23.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod24.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod25.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod26.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod27.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod28.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod29.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod30.txt
#  GA Search [====>-------------------]  20% (iteration 3/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod31.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod32.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod33.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod34.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod35.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod36.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod37.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod38.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod39.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod40.txt
#  GA Search [=====>------------------]  27% (iteration 4/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod41.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod42.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod43.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod44.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod45.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod46.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod47.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod48.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod49.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod50.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod51.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod52.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod53.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod54.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod55.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod56.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod57.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod58.txt
#  GA Search [=======>----------------]  33% (iteration 5/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod59.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod60.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod61.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod62.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod63.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod64.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod65.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod66.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod67.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod68.txt
#  GA Search [=========>--------------]  40% (iteration 6/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod69.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod70.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod71.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod72.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod73.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod74.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod75.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod76.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod77.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod78.txt
#  GA Search [==========>-------------]  47% (iteration 7/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod79.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod80.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod81.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod82.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod83.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod84.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod85.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod86.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod87.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod88.txt
#  GA Search [============>-----------]  53% (iteration 8/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod89.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod90.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod91.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod92.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod93.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod94.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod95.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod96.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod97.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod98.txt
#  GA Search [=============>----------]  60% (iteration 9/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod99.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod100.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod101.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod102.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod103.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod104.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod105.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod106.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod107.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod108.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod109.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod110.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod111.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod112.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod113.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod114.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod115.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod116.txt
#  GA Search [==============>--------]  67% (iteration 10/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod117.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod118.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod119.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod120.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod121.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod122.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod123.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod124.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod125.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod126.txt
#  GA Search [================>------]  73% (iteration 11/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod127.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod128.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod129.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod130.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod131.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod132.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod133.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod134.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod135.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod136.txt
#  GA Search [=================>-----]  80% (iteration 12/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod137.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod138.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod139.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod140.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod141.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod142.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod143.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod144.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod145.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod146.txt
#  GA Search [===================>---]  87% (iteration 13/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod147.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod148.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod149.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod150.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod151.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod152.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod153.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod154.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod155.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod156.txt
#  GA Search [====================>--]  93% (iteration 14/15)
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod157.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod158.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod159.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod160.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod161.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod162.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod163.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod164.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod165.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod166.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod167.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod168.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod169.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod170.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod171.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod172.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod173.txt
# [Success] Model file created in current working directory:
# /home/zhonghuihuang/Desktop/test/Step_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/GA_2025-09-04-pheno_sd_359e255ca82e284a7b41aefd444b39a1_temp/mod174.txt
#  GA Search [=======================] 100% (iteration 15/15)
```
GA results
```r
# > print(outs)
# 
# === Final Selected Model Code ===
#   no.cmpt1 no.cmpt2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 1        0        1      0      1      0       0     0      0  0     0   0   1
# 
# === Final Selected Model Name ===
#   iv_1cmpt_etaCLVC_First-order elimination_No correlation_additive 
```

ACO example
``` r
outs<-aco.operator(dat=pheno_sd,
                  search.space = "ivbase",
                  filename =  "pheno_sd",
                  foldername =   "pheno_sd" )
print(outs)
``` 

Tabu example
``` r
outs<-tabu.operator(dat=pheno_sd,
                   search.space = "ivbase",
                   filename =  "pheno_sd",
                   foldername =   "pheno_sd" )
print(outs)

```


