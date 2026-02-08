

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
```

<details>
<summary>📊 Click to view Stepwise output</summary>

```r
# Infometrics                               Value          
# ----------------------------------------  ---------------
# Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Number of Subjects                        59             
# Number of Observations                    155            
# Subjects with First-Dose Interval Data    35             
# Observations in the First-Dose Interval   35             
# Subjects with Multiple-Dose Data          56             
# Observations after Multiple Doses         120            
# ----------------------------------------  ------
# Estimating half-life....................
# Half-life estimation complete: Estimated t1/2 = 16.44 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 0.0087, estimated Vd: 1.25 
# Run parameter sweeping on nonlinear elimination kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters....................
# Running Step 1: Number of compartments ----------------------------------------------------
# Loaded cached model results from 'pheno_sd_sf/pheno_sd_sf.csv'.
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod1.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod2.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod3.txt
# Running Step 2: Elimination type ----------------------------------------------------
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod4.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod5.txt
# Running Step 3: IIV on Km ----------------------------------------------------
# Step 3: IIV on Km skipped.
# Running Step 5: IIV (forward selection) ----------------------------------------------------
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod6.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod7.txt
# Running Step 6: ETA correlation ----------------------------------------------------
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod8.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod9.txt
# Running Step 7: Residual error model ----------------------------------------------------
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod10.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod11.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_sf/mod12.txt
# > print(pheno_sd_sf)
# 
# === Best Model Code ===
#   no.cmpt eta.vmax eta.km eta.cl eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 1       1        0      0      1      1      0       0     0      0  0     1  1
# 
# === Best Model Name ===
# bolus_1cmpt_etaCLVC_FOelim_correlated_add 
# 
# === Stepwise Selection History ===
#                   Step                                                Penalty.terms                                       Model.name          Model.code   Fitness       AIC      BIC      OFV
# 1  No. of compartments                                       rse, theta, covariance   bolus_1cmpt_etaCL_FOelim_uncorrelated_combined 1,0,0,0,0,0,0,0,0,3 21185.358 1170.1409 1185.358 875.2699
# 2     Elimination type                                       rse, theta, covariance   bolus_1cmpt_etaCL_FOelim_uncorrelated_combined 1,0,0,0,0,0,0,0,0,3 21185.358 1170.1409 1185.358 875.2699
# 5        IIV (forward)                     rse, theta, covariance, shrinkage, omega bolus_1cmpt_etaCLVC_FOelim_uncorrelated_combined 1,0,1,0,0,0,0,0,0,3 11104.215 1085.9549 1104.215 789.0840
# 6      Eta correlation        rse, theta, covariance, shrinkage, omega, correlation   bolus_1cmpt_etaCLVC_FOelim_correlated_combined 1,0,1,0,0,0,0,0,1,3  1097.977 1066.6730 1087.977 767.8021
# 7 Residual error types rse, theta, covariance, shrinkage, omega, correlation, sigma        bolus_1cmpt_etaCLVC_FOelim_correlated_add 1,0,1,0,0,0,0,0,1,1  1014.616  986.3554 1004.616 689.4845
```

</details>


GA example
``` r
pheno_sd_ga <- ga.operator(
  seed = 20,
  dat = pheno_sd,                    # Dataset used for model fitting
  search.space = "ivbase",           # Structural search space for IV PK models
  filename = "pheno_sd_ga",          # Prefix for output files
  foldername = "pheno_sd_ga",        # Folder where results will be stored
  saem.control = saemControl(        # SAEM estimation control settings
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

print(pheno_sd_ga)

# Infometrics                               Value          
# ----------------------------------------  ---------------
# Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Number of Subjects                        59             
# Number of Observations                    155            
# Subjects with First-Dose Interval Data    35             
# Observations in the First-Dose Interval   35             
# Subjects with Multiple-Dose Data          56             
# Observations after Multiple Doses         120            
# ----------------------------------------  ------
# Estimating half-life....................
# Half-life estimation complete: Estimated t1/2 = 16.44 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 0.0087, estimated Vd: 1.25 
# Run parameter sweeping on nonlinear elimination kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters....................
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod1.txt
# SAEM control (core) = niter=200|300; nBurn=200; nEm=300; seed=99; print=1
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod2.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod3.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod4.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod5.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod6.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod7.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod8.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod9.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod10.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod11.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod12.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod13.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod14.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod15.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod16.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod17.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod18.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod19.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod20.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod21.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod22.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod23.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod24.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod25.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod26.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod27.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod28.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod29.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod30.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod31.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod32.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod33.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod34.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod35.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod36.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod37.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod38.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod39.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod40.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod41.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod42.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod43.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod44.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod45.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod46.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod47.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod48.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod49.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod50.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod51.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod52.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod53.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod54.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod55.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod56.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod57.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod58.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod59.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod60.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod61.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod62.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod63.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod64.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod65.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod66.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod67.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod68.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod69.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod70.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod71.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod72.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod73.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod74.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod75.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod76.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod77.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod78.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod79.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod80.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod81.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod82.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod83.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod84.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod85.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod86.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod87.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod88.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod89.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod90.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod91.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod92.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod93.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod94.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod95.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod96.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod97.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod98.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod99.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod100.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod101.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod102.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod103.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod104.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod105.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod106.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod107.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod108.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod109.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod110.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod111.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod112.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod113.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod114.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod115.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod116.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod117.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod118.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod119.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod120.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod121.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod122.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod123.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod124.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod125.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod126.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod127.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod128.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod129.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod130.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod131.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod132.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod133.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod134.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod135.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod136.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod137.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod138.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod139.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod140.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod141.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod142.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod143.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod144.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod145.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod146.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod147.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod148.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod149.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod150.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod151.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod152.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod153.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod154.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod155.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod156.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod157.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod158.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod159.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod160.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod161.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod162.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod163.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod164.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod165.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod166.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod167.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod168.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod169.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod170.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod171.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod172.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod173.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod174.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod175.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod176.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod177.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod178.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod179.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod180.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod181.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod182.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod183.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod184.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod185.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod186.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod187.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod188.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod189.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod190.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod191.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod192.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod193.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod194.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod195.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod196.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod197.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod198.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod199.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod200.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod201.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod202.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod203.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod204.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod205.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod206.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod207.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod208.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod209.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod210.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod211.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod212.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod213.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod214.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod215.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod216.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod217.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod218.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod219.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod220.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod221.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod222.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod223.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod224.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod225.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod226.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod227.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod228.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod229.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod230.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod231.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod232.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod233.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod234.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod235.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod236.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod237.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod238.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod239.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod240.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod241.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod242.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod243.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod244.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod245.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod246.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod247.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod248.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod249.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod250.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod251.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod252.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod253.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod254.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod255.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod256.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod257.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod258.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod259.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod260.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod261.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod262.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod263.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod264.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod265.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod266.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod267.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod268.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod269.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod270.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod271.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod272.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod273.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod274.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod275.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod276.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod277.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod278.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod279.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod280.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod281.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod282.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod283.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod284.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod285.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod286.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod287.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod288.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod289.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod290.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod291.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod292.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod293.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod294.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod295.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod296.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod297.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod298.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod299.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod300.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod301.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod302.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod303.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod304.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod305.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod306.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod307.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod308.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod309.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod310.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod311.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod312.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod313.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod314.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod315.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod316.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod317.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod318.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod319.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod320.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod321.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod322.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod323.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod324.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod325.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod326.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod327.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod328.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod329.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod330.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod331.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod332.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod333.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod334.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod335.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod336.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod337.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod338.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod339.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod340.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod341.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod342.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod343.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod344.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod345.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod346.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod347.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod348.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod349.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod350.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod351.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod352.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod353.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod354.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod355.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod356.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod357.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod358.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod359.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod360.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod361.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod362.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod363.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod364.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod365.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod366.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod367.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod368.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod369.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod370.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod371.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod372.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod373.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod374.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod375.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod376.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod377.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod378.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod379.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod380.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod381.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod382.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod383.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod384.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod385.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod386.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod387.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod388.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod389.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod390.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod391.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod392.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod393.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod394.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod395.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod396.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod397.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod398.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod399.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod400.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod401.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod402.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod403.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod404.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod405.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod406.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod407.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod408.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod409.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod410.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod411.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod412.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod413.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod414.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod415.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod416.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod417.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod418.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod419.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod420.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod421.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod422.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod423.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod424.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod425.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod426.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod427.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod428.txt
# [Success] Model file created:                                                                                                                                                                
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod429.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod430.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod431.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod432.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod433.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod434.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod435.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod436.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod437.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod438.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod439.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod440.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod441.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod442.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod443.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod444.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod445.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod446.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod447.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_ga/mod448.txt
# >                                                                                                                                                                                            
# > print(pheno_sd_ga)
# 
# === Final Selected Model Code ===
#   no.cmpt1 no.cmpt2 eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv1 rv2
# 1        0        1      0      1      0       0     0      0  0     1   0   1
# 
# === Final Selected Model Name ===
# bolus_1cmpt_etaCLVC_FOelim_correlated_add 
```

ACO example
``` r
pheno_sd_aco <-  aco.operator(
  seed = 20,                         # Random seed for aco
  dat = pheno_sd,                    # Dataset used for model fitting
  search.space = "ivbase",           # Structural search space for IV PK models
  filename = "pheno_sd_aco",         # Prefix for output files
  foldername = "pheno_sd_aco",       # Folder where results will be stored
  saem.control = saemControl(        # SAEM estimation control settings
    nBurn = 200,                     # SAEM burn-in iterations
    nEm   = 300,                     # SAEM EM-phase iterations
    rxControl = rxControl(cores = 4),# Number of cores for ODE solving.
    logLik    = TRUE                 # Compute log-likelihood
  ),
  table.control = tableControl(
    cwres = TRUE                     # Compute conditional weighted residuals (CWRES)
  ),
  max_wall_time = 2 * 60 * 60        # Maximum allowed wall-clock time (seconds) per model; here: 2 hours
)
print(pheno_sd_aco) 

# Infometrics                               Value          
# ----------------------------------------  ---------------
# Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Number of Subjects                        59             
# Number of Observations                    155            
# Subjects with First-Dose Interval Data    35             
# Observations in the First-Dose Interval   35             
# Subjects with Multiple-Dose Data          56             
# Observations after Multiple Doses         120            
# ----------------------------------------  ------
# Estimating half-life....................
# Half-life estimation complete: Estimated t1/2 = 16.44 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 0.0087, estimated Vd: 1.25 
# Run parameter sweeping on nonlinear elimination kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters....................
# Loaded cached model results from 'pheno_sd_aco/pheno_sd_aco.csv'.                                                                                                                             
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod1.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod2.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod3.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod4.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod5.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod6.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod7.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod8.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod9.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod10.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod11.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod12.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod13.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod14.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod15.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod16.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod17.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod18.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod19.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod20.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod21.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod22.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod23.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod24.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod25.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod26.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod27.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod28.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod29.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod30.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod31.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod32.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod33.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod34.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod35.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod36.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod37.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod38.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod39.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod40.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod41.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod42.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod43.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod44.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod45.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod46.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod47.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod48.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod49.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod50.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod51.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod52.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod53.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod54.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod55.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod56.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod57.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod58.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod59.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod60.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod61.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod62.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod63.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod64.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod65.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod66.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod67.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod68.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod69.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod70.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod71.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod72.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod73.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod74.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod75.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod76.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod77.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod78.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod79.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod80.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod81.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod82.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod83.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod84.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod85.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod86.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod87.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod88.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod89.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod90.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod91.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod92.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod93.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod94.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod95.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod96.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod97.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod98.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod99.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod100.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod101.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod102.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod103.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod104.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod105.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod106.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod107.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod108.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod109.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod110.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod111.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod112.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod113.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod114.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod115.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod116.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod117.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod118.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod119.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod120.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod121.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod122.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod123.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod124.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod125.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod126.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod127.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod128.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod129.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod130.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod131.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod132.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod133.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod134.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod135.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod136.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod137.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod138.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod139.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod140.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod141.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod142.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod143.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod144.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod145.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod146.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod147.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod148.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod149.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod150.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod151.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod152.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod153.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod154.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod155.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod156.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod157.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod158.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod159.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod160.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod161.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod162.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod163.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod164.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod165.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod166.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod167.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod168.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod169.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod170.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod171.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod172.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod173.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod174.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod175.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod176.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod177.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod178.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod179.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod180.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod181.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod182.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod183.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod184.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod185.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod186.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod187.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod188.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod189.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod190.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod191.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod192.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod193.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod194.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod195.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod196.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod197.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod198.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod199.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod200.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod201.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod202.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod203.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod204.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod205.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod206.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod207.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod208.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod209.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod210.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod211.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod212.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod213.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod214.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod215.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod216.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod217.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod218.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod219.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod220.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod221.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod222.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod223.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod224.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod225.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod226.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod227.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod228.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod229.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod230.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod231.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod232.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod233.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod234.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod235.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod236.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod237.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod238.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod239.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod240.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod241.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod242.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod243.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod244.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod245.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod246.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod247.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod248.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod249.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod250.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod251.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod252.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod253.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod254.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod255.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod256.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod257.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod258.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod259.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod260.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod261.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod262.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod263.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod264.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod265.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod266.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod267.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod268.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod269.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod270.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod271.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod272.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod273.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod274.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod275.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod276.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod277.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod278.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod279.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod280.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod281.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod282.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod283.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod284.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod285.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod286.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod287.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod288.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod289.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod290.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod291.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod292.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod293.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod294.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod295.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod296.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod297.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod298.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod299.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_aco/mod300.txt
# > print(pheno_sd_aco)                                                                                                                                                                         
# 
# === Final Selected Model Code (ACO) ===
#    no.cmpt eta.km eta.vc eta.vp eta.vp2 eta.q eta.q2 mm mcorr rv
# 67       1      0      1      0       0     0      0  0     1  1
# 
# === Final Selected Model Name (ACO) ===
# bolus_1cmpt_etaCLVC_FOelim_correlated_add 
``` 

Tabu example
``` r
pheno_sd_tabu <- tabu.operator(
  dat = pheno_sd,                    # Dataset used for model fitting
  search.space = "ivbase",           # Structural search space for IV PK models
  filename = "pheno_sd_tabu",        # Prefix for output files
  foldername = "pheno_sd_tabu",      # Folder where results will be stored
  saem.control = saemControl(        # SAEM estimation control settings
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

print(pheno_sd_tabu)

# Infometrics                               Value          
# ----------------------------------------  ---------------
# Dose Route                                bolus          
# Dose Type                                 combined_doses 
# Number of Subjects                        59             
# Number of Observations                    155            
# Subjects with First-Dose Interval Data    35             
# Observations in the First-Dose Interval   35             
# Subjects with Multiple-Dose Data          56             
# Observations after Multiple Doses         120            
# ----------------------------------------  ------
# Estimating half-life....................
# Half-life estimation complete: Estimated t1/2 = 16.44 h
# Evaluating the predictive performance of calculated one-compartment model parameters....................
# Base PK parameter analysis finished. Estimated ka: NA, estimated CL: 0.0087, estimated Vd: 1.25 
# Run parameter sweeping on nonlinear elimination kinetics PK parameters....................
# Run parameter sweeping on multi-compartmental PK parameters....................
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod1.txt
# SAEM control (core) = niter=200|300; nBurn=200; nEm=300; seed=99; print=1
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod2.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod3.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod4.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod5.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod6.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod7.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod8.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod9.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod10.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod11.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod12.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod13.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod14.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod15.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod16.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod17.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod18.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod19.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod20.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod21.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod22.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod23.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod24.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod25.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod26.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod27.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod28.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod29.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod30.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod31.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod32.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod33.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod34.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod35.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod36.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod37.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod38.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod39.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod40.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod41.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod42.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod43.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod44.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod45.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod46.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod47.txt
# Iteration 6: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod48.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod49.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod50.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod51.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod52.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod53.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod54.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod55.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod56.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod57.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod58.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod59.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod60.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod61.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod62.txt
# Iteration 8: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod63.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod64.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod65.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod66.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod67.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod68.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod69.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod70.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod71.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod72.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod73.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod74.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod75.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod76.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod77.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod78.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod79.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod80.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod81.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod82.txt
# Iteration 10: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod83.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod84.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod85.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod86.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod87.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod88.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod89.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod90.txt
# Iteration 11: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod91.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod92.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod93.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod94.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod95.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod96.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod97.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod98.txt
# Iteration 12: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod99.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod100.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod101.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod102.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod103.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod104.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod105.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod106.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod107.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod108.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod109.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod110.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod111.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod112.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod113.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod114.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod115.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod116.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod117.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod118.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod119.txt
# Iteration 14: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod120.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod121.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod122.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod123.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod124.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod125.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod126.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod127.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod128.txt
# Iteration 15: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod129.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod130.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod131.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod132.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod133.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod134.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod135.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod136.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod137.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod138.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod139.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod140.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod141.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod142.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod143.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod144.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod145.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod146.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod147.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod148.txt
# Iteration 17: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod149.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod150.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod151.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod152.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod153.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod154.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod155.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod156.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod157.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod158.txt
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod159.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod160.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod161.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod162.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod163.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod164.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod165.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod166.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod167.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod168.txt
# Iteration 19: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# [Success] Model file created:                                                                                                                                                                 
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod169.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod170.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod171.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod172.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod173.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod174.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod175.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod176.txt
# [Success] Model file created:
# /home/zhonghuihuang/Desktop/nlmixr2autotest/pheno_sd_tabu/mod177.txt
# Iteration 20: candidate already used as a starting point. Applying 2-bit perturbation to avoid cycling.
# >                                                                                                                                                                                             
# > print(pheno_sd_tabu)
# 
# === Final Selected Model Code (Tabu Search) ===
# no.cmpt  eta.km  eta.vc  eta.vp eta.vp2   eta.q  eta.q2      mm   mcorr      rv 
#       1       0       1       0       0       0       0       0       1       1 
# 
# === Final Selected Model Name (Tabu Search) ===
# bolus_1cmpt_etaCLVC_FOelim_correlated_add 
```


