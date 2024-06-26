#' Run stepwise algorithm for automated model selection
#'
#' Perform a stepwise algorithm to select the best model for pharmacokinetic modelling
#'
#' @param dat A data frame containing the pharmacokinetic data to be used for model fitting.
#' @param search.space An integer specifying the search space type (1 for standard).
#' @param no.cores The number of CPU cores to use for parallel processing.
#' @param thetalower A named vector specifying lower bounds for model parameters.
#' @param crse Constraint value for the relative standard error (RSE).
#' @param cshrink Constraint value for the shrinkage.
#' @param cadd Constraint value for the additive error.
#' @param cprop Constraint value for the proportional error.
#' @param ccorr Constraint value for the correlation between parameters.
#' @param penalty.type The type of penalty to apply during the fitness evaluation:
#' \itemize{
#'   \item \code{penalty.type = 1}: Considers constraints on RSE.
#'   \item \code{penalty.type = 2}: Considers constraints on RSE and shrinkage.
#'   \item \code{penalty.type = 3}: Considers constraints on RSE, shrinkage, and sigma values (variance of residual error).
#' }
#' @param penalty.value The value of the penalty to apply during the fitness evaluation.
#' @param foldername The name of the folder for output storage. Defaults to "test".
#' @param filename The name of the output file. Defaults to "test".
#' @param ... Additional arguments passed to the model fitting function.
#'
#' @return A list include the best model selection and stepwise history running.
#'
#' @examples
#' \dontrun{
#' library(nlmixr2autoinit)
#' dat<-pheno_sd
#' filename =  "pheno_sd"
#' foldername =   "pheno_sd"
#' inits.out<-getppkinit(dat = dat,runnpd = 0)
#' autosetinit(dat = dat,
#' inits= inits.out$Recommended_initial_estimates)
#'sf.operator(dat,
#' no.cores = 4,
#' thetalower=c(vp=1,
#'  vp2=1),
#'  control = saemControl(
#'   seed = 1234,
#'  print = 5,
#'  nBurn = 20,
#'  nEm = 30,
#'  logLik = T,
#'  rxControl = rxControl(cores = 4,
#'  maxsteps =70000)),
#'  table=tableControl(cwres=T),
#'  filename =  "Bolus_1CPT",
#'  foldername =   "Bolus_1CPT"  )
#'}

#' @export

sf.operator <- function(dat,
                        search.space,
                        no.cores,
                        thetalower,
                        crse,
                        cshrink,
                        cadd,
                        cprop,
                        ccorr,
                        penalty.type,
                        penalty.value,
                        foldername,
                        filename,
                        ...) {
  # If exists, then create a new one
  # Output directory
  if (missing(foldername)) {
    foldername <- "test"
  }
  if (missing(filename)) {
    filename <- "test"
  }
  
  foldername <<- foldername
  filename <<- filename
  
  current.date <- Sys.Date()
  
  outputdir <-
    paste0("Step_",
           current.date,
           "-",
           foldername,
           "_",
           digest::digest(dat),
           "_temp")
  
  # for linux ubuntu
  if (dir.exists(outputdir)) {
    print("Warning: current directory for stepwsie analysis already exists")
    unlink(outputdir, recursive = T, force = T)
    print("Warning: a new one was created and replace the previous one")
    dir.create(outputdir, showWarnings = F, recursive = T)
  }
  
  if (!dir.exists(outputdir)) {
    print("Output directory for stepwise analysis is created")
    dir.create(outputdir, showWarnings = F, recursive = T)
  }
  setwd(paste0(getwd(), "/", outputdir))
  storage.path <- getwd()
  
  if (missing(search.space)) {
    search.space <- 1
  }
  
  ################################## Default Modelling setting###########################
  # Default of constraint setting
  # RSE<=20%
  # Shrinkage <=30%
  # Theta lower boundary, default=0
  # SD of additive >=1 (approximately ) 1/1000-
  # SD of proportional >= 5%
  # R square for omega correlation coefficient < 0.1
  thetalower0 <- c(
    lbka = 0,
    lbvc = 0,
    lbcl = 0,
    lbvp = 0,
    lbvp2 = 0,
    lbq = 0,
    lbq2 = 0,
    lbtlag = 0,
    lbvmax = 0,
    lbkm = 0
  )
  
  if (missing(thetalower)) {
    thetalower = thetalower0
  }
  
  if (is.na(thetalower["ka"]) == F) {
    thetalower0["lbka"] = thetalower["ka"]
  }
  if (is.na(thetalower["cl"]) == F) {
    thetalower0["lbcl"] = thetalower["cl"]
  }
  if (is.na(thetalower["vc"]) == F) {
    thetalower0["lbvc"] = thetalower["vc"]
  }
  if (is.na(thetalower["vp"]) == F) {
    thetalower0["lbvp"] = thetalower["vp"]
  }
  if (is.na(thetalower["vp2"]) == F) {
    thetalower0["lbvp2"] = thetalower["vp2"]
  }
  if (is.na(thetalower["q"]) == F) {
    thetalower0["lbq"] = thetalower["q"]
  }
  if (is.na(thetalower["q2"]) == F) {
    thetalower0["lbq2"] = thetalower["q2"]
  }
  if (is.na(thetalower["tlag"]) == F) {
    thetalower0["lbtlag"] = thetalower["tlag"]
  }
  if (is.na(thetalower["vmax"]) == F) {
    thetalower0["lbvmax"] = thetalower["vmax"]
  }
  if (is.na(thetalower["km"]) == F) {
    thetalower0["lbkm"] = thetalower["km"]
  }
  
  lbka = thetalower0["lbka"]
  lbcl = thetalower0["lbcl"]
  lbvc = thetalower0["lbvc"]
  lbvp = thetalower0["lbvp"]
  lbvp2 = thetalower0["lbvp2"]
  lbq = thetalower0["lbq"]
  lbq2 = thetalower0["lbq2"]
  lbtlag = thetalower0["lbtlag"]
  lbvmax = thetalower0["lbvmax"]
  lbkm = thetalower0["lbkm"]
  
  
  if (missing(crse)) {
    crse <- 20
  }
  if (missing(cshrink)) {
    cshrink <- 30
  }
  
  
  if (missing(cadd) & search.space == 1) {
    dat.obs <- dat[dat$EVID == 0, ]
    pop.cmax <- aggregate(dat.obs, by = list(dat.obs$ID), FUN = max)
    dat.cmax <<- median(pop.cmax$DV)
    cadd <- round(dat.cmax / 1000, 0)
  }
  
  if (missing(cprop)) {
    cprop <- 0.05
  }
  
  if (missing(ccorr)) {
    ccorr <- 0 # default 0, no constraint on the correlation
  }
  
  if (missing(no.cores)) {
    no.cores <-
      getRxThreads() # default 0, no constraint on the correlation
  }
  
  #Initial Blank setting
  mod.record.best.all <- NULL
  mod.record.all <- NULL
  model.tried. <- 0
  mod.record.best <- NULL
  mod.record <- NULL
  
  #################################Step1. No. of compartment##########################
  message(blue(
    paste0(
      "Running Stepwise 1. Structural Model----------------------------------------------------"
    )
  ))
  
  message(blue(
    paste0(
      "Test number of compartments----------------------------------------------------"
    )
  ))
  
  # set as the variable in the global environment.
  modi <<- 1
  r <<- 1
  dat <<- dat
  search.space <<- search.space
  crse <<- crse
  cshrink <<- cshrink
  lbcl <<- lbcl
  lbvc <<- lbvc
  lbvp <<- lbvp
  lbq <<- lbq
  lbvp2 <<- lbvp2
  lbq2 <<- lbq2
  cadd <<- cadd
  cprop <<- cprop
  ccorr <<- ccorr
  
  penalty.type <<- 1
  penalty.value <<- 10000
  
  
  # Start from no. of compartment, only random effects on CL, combined residual error model
  mod.string1 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 3)
  mod.string2 <- c(2, 0, 0, 0, 0, 0, 0, 0, 0, 3)
  mod.string3 <- c(3, 0, 0, 0, 0, 0, 0, 0, 0, 3)
  
  mod1.fit <- sf.mod.run(modi = modi, string = mod.string1, ...)
  mod2.fit <- sf.mod.run(modi = modi, string = mod.string2, ...)
  mod3.fit <- sf.mod.run(modi = modi, string = mod.string3, ...)
  
  mod.step1 <-
    data.frame(
      modname = c(
        read.code2(search.space = search.space, sel.best.code = mod.string1),
        read.code2(search.space = search.space, sel.best.code =
                     mod.string2),
        read.code2(search.space = search.space, sel.best.code =
                     mod.string3)
      ),
      cmpt = c(mod.string1[1], mod.string2[1], mod.string3[1]),
      fitness = c(mod1.fit,
                  mod2.fit,
                  mod3.fit),
      mod = c(
        toString(mod.string1),
        toString(mod.string2),
        toString(mod.string3)
      )
    )
  
  mod.step1$step = "no. of compartments"
  mod.record <- rbind(mod.record, mod.step1)
  # Select the minimum fitness values
  mod.step1 <-
    mod.step1[mod.step1$fitness == min(mod.step1$fitness), ]
  # Record the best solution for current steps
  
  mod.record.best <- rbind(mod.record.best, mod.step1)
  
  model.tried. <- model.tried. + 3
  
  
  r <<- 1.2
  ##################### Step 1.2. Analyse the nonlinear elimination###############
  
  message(blue(
    paste0(
      "Analyse elimination type----------------------------------------------------"
    )
  ))
  
  
  mod.step1.string <- as.numeric(strsplit(mod.step1$mod, ",")[[1]])
  mod.step1.string.fit <-
    sf.mod.run(modi = modi, string = mod.step1.string, ...)
  mod.string.r <- mod.step1.string
  mod.string.r[8] <- 1
  
  
  mod.string.r.fit <-
    sf.mod.run(modi = modi, string = mod.string.r, ...)
  mod.step1.2 <-
    data.frame(
      modname = c(
        read.code2(search.space = search.space, sel.best.code = mod.step1.string),
        read.code2(search.space = search.space, sel.best.code = mod.string.r)
      ),
      cmpt = c(mod.step1.string[1], mod.string.r[1]),
      fitness = c(mod.step1.string.fit,
                  mod.string.r.fit),
      mod = c(toString(mod.step1.string),
              toString(mod.string.r))
    )
  
  mod.step1.2$step = "elimination type"
  mod.record <- rbind(mod.record, mod.step1.2)
  mod.step1.2 <-
    mod.step1.2[mod.step1.2$fitness == min(mod.step1.2$fitness), ]
  mod.record.best <- rbind(mod.record.best, mod.step1.2)
  
  model.tried. <- model.tried. + 1
  
  
  r <<- 1.3
  ##################### Step 1.3. Introduce Km #################################
  
  type <<- 2
  penalty <<- 10000
  
  message(blue(
    paste0(
      "Test IIV on Km----------------------------------------------------"
    )
  ))
  
  
  mod.step1.2.string <-
    as.numeric(strsplit(mod.step1.2$mod, ",")[[1]])
  mod.step1.2.string.fit <-
    sf.mod.run(modi = modi, string = mod.step1.2.string, ...)
  mod.string.r <- mod.step1.2.string
  
  if (mod.string.r[8] == 1) {
    message(blue(
      paste0(
        "Running Stepwise 2. Analyse the elimination type---------------------------------------- ",
      )
    ))
    
    
    mod.string.r[2] <- 1
    
    
    mod.string.r.fit <-
      sf.mod.run(modi = modi, string = mod.string.r, ...)
    mod.step1.3 <-
      data.frame(
        modname = c(
          read.code2(search.space = search.space, mod.step1.2.string),
          read.code2(search.space = search.space, mod.string.r)
        ),
        cmpt = c(mod.step1.2.string[1], mod.string.r[1]),
        fitness = c(mod.step1.2.string.fit,
                    mod.string.r.fit),
        mod = c(toString(mod.step1.2.string),
                toString(mod.string.r))
      )
    
    mod.step1.3$step = "IIV on Km"
    mod.record <- rbind(mod.record, mod.step1.3)
    mod.step1.3 <-
      mod.step1.3[mod.step1.3$fitness == min(mod.step1.3$fitness), ]
    mod.record.best <- rbind(mod.record.best, mod.step1.3)
    model.tried. <- model.tried. + 1
    
  }
  
  
  
  # Temporarily use model.step1 name
  mod.step1 <- mod.step1.2
  if (exists("mod.step1.3")) {
    mod.step1 <- mod.step1.3
  }
  
  ############ Step 3. Introduce random effects by stepwise method###############
  
  message(blue(
    paste0(
      "Introduce IIV on parameters----------------------------------------------------"
    )
  ))
  
  r <<- 2
  mod.step1.string <- as.numeric(strsplit(mod.step1$mod, ",")[[1]])
  # Re-run the best solution by the new fitness function
  mod.step1.string.fit <-
    sf.mod.run(modi = modi, string = mod.step1.string, ...)
  # 1Cmpt
  
  if (mod.step1$cmpt == 1) {
    mod.string.r <- mod.step1.string
    mod.string.r[3] <- 1
    
    mod.string.r.fit <-
      sf.mod.run(modi = modi, string = mod.string.r, ...)
    mod.step2 <-
      data.frame(
        modname = c(
          read.code2(1, mod.step1.string),
          read.code2(1, mod.string.r)
        ),
        cmpt = c(mod.step1.string[1], mod.string.r[1]),
        fitness = c(mod.step1.string.fit, mod.string.r.fit),
        mod = c(toString(mod.step1.string), toString(mod.string.r))
      )
    
    mod.step2$step = "IIV on parameters"
    mod.record <- rbind(mod.record, mod.step2)
    mod.step2 <-
      mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
    mod.record.best <- rbind(mod.record.best, mod.step2)
    
    model.tried. <- model.tried. + 1
  }
  
  
  if (mod.step1$cmpt > 1) {
    mod.test <- mod.step1
    mod.test.string.fit <-
      sf.mod.run(modi = modi, string = mod.step1.string, ...)
  }
  
  # 2Cmpt
  if (mod.step1$cmpt == 2) {
    mod.string.r1 <- mod.step1.string
    mod.string.r2 <- mod.step1.string
    mod.string.r3 <- mod.step1.string
    
    mod.string.r1[3] <- 1
    mod.string.r2[4] <- 1
    mod.string.r3[6] <- 1
    
    mod.string.r1.fit <-
      sf.mod.run(modi = modi, string = mod.string.r1, ...)
    mod.string.r2.fit <-
      sf.mod.run(modi = modi, string = mod.string.r2, ...)
    mod.string.r3.fit <-
      sf.mod.run(modi = modi, string = mod.string.r3, ...)
    
    
    mod.step2 <-
      data.frame(
        modname = c(
          read.code2(search.space = search.space, mod.step1.string),
          read.code2(search.space = search.space, mod.string.r1),
          read.code2(search.space = search.space, mod.string.r2),
          read.code2(search.space = search.space, mod.string.r3)
        ),
        cmpt = c(
          mod.step1.string[1],
          mod.string.r1[1],
          mod.string.r2[1],
          mod.string.r3[1]
        ),
        fitness = c(
          mod.step1.string.fit,
          mod.string.r1.fit,
          mod.string.r2.fit,
          mod.string.r3.fit
        ),
        mod = c(
          toString(mod.test$mod),
          toString(mod.string.r1),
          toString(mod.string.r2),
          toString(mod.string.r3)
        )
      )
    
    mod.step2$step = "IIV on parameters"
    mod.record <- rbind(mod.record, mod.step2)
    
    mod.step2 <-
      mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
    mod.record.best <- rbind(mod.record.best, mod.step2)
    model.tried. <- model.tried. + 3
    
    if (mod.step2$mod != mod.test$mod) {
      mod.test <- mod.step2
      mod.step2.string <-
        as.numeric(strsplit(mod.step2$mod, ",")[[1]])
      
      if (mod.step2.string [3] == 1) {
        mod.string.r1 <- mod.step2.string
        mod.string.r2 <- mod.step2.string
        
        mod.string.r1[4] <- 1
        mod.string.r2[6] <- 1
        
        mod.string.r1.fit <-
          sf.mod.run(modi = modi, string = mod.string.r1, ...)
        mod.string.r2.fit <-
          sf.mod.run(modi = modi, string = mod.string.r2, ...)
        
        mod.step2 <-
          data.frame(
            modname = c(
              read.code2(search.space = search.space, mod.step2.string),
              read.code2(search.space = search.space, mod.string.r1),
              read.code2(search.space = search.space, mod.string.r2)
            ),
            cmpt = c(mod.step2.string[1], mod.string.r1[1], mod.string.r2[1]),
            fitness = c(
              mod.test$fitness,
              mod.string.r1.fit,
              mod.string.r2.fit
            ),
            mod = c(
              toString(mod.test$mod),
              toString(mod.string.r1),
              toString(mod.string.r2)
            )
          )
        
        mod.step2$step = "IIV on parameters"
        mod.record <- rbind(mod.record, mod.step2)
        
        mod.step2 <-
          mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
        mod.record.best <- rbind(mod.record.best, mod.step2)
        model.tried. <- model.tried. + 2
      }
      
      
      if (mod.step2.string [4] == 1) {
        mod.string.r1 <- mod.step2.string
        mod.string.r2 <- mod.step2.string
        
        mod.string.r1[3] <- 1
        mod.string.r2[6] <- 1
        
        mod.string.r1.fit <-
          sf.mod.run(modi = modi, string = mod.string.r1, ...)
        mod.string.r2.fit <-
          sf.mod.run(modi = modi, string = mod.string.r2, ...)
        
        
        mod.step2 <-
          data.frame(
            modname = c(
              read.code2(search.space = search.space, mod.step2.string),
              read.code2(search.space = search.space, mod.string.r1),
              read.code2(search.space = search.space, mod.string.r2)
            ),
            cmpt = c(mod.step2.string, mod.string.r1[1], mod.string.r2[1]),
            fitness = c(
              mod.test$fitness,
              mod.string.r1.fit,
              mod.string.r2.fit
            ),
            mod = c(
              toString(mod.test$mod),
              toString(mod.string.r1),
              toString(mod.string.r2)
            )
          )
        
        mod.step2$step = "IIV on parameters"
        mod.record <- rbind(mod.record, mod.step2)
        
        mod.step2 <-
          mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
        mod.record.best <- rbind(mod.record.best, mod.step2)
        model.tried. <- model.tried. + 2
      }
      
      
      if (mod.step2.string [6] == 1) {
        mod.string.r1 <- mod.step2.string
        mod.string.r2 <- mod.step2.string
        
        mod.string.r1[3] <- 1
        mod.string.r2[4] <- 1
        
        mod.string.r1.fit <-
          sf.mod.run(modi = modi, string = mod.string.r1, ...)
        mod.string.r2.fit <-
          sf.mod.run(modi = modi, string = mod.string.r2, ...)
        
        mod.step2 <-
          data.frame(
            modname = c(
              read.code2(search.space = search.space, mod.step2.string),
              read.code2(search.space = search.space, mod.string.r1),
              read.code2(search.space = search.space, mod.string.r2)
            ),
            cmpt = c(mod.step2.string[1], mod.string.r1[1], mod.string.r2[1]),
            fitness = c(
              mod.test$fitness,
              mod.string.r1.fit,
              mod.string.r2.fit
            ),
            mod = c(
              toString(mod.test$mod),
              toString(mod.string.r1),
              toString(mod.string.r2)
            )
          )
        
        mod.step2$step = "IIV on parameters"
        mod.record <- rbind(mod.record, mod.step2)
        
        mod.step2 <-
          mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
        mod.record.best <- rbind(mod.record.best, mod.step2)
        model.tried. <- model.tried. + 2
      }
      
      
      if (mod.step2$mod != mod.test$mod) {
        mod.test <- mod.step2
        
        mod.step2.string <-
          as.numeric(strsplit(mod.step2$mod, ",")[[1]])
        
        if (mod.step2.string [3] != 1) {
          mod.string.r1 <- mod.step2.string
          
          mod.string.r1[3] <- 1
          
          mod.string.r1.fit <-
            sf.mod.run(modi = modi, string = mod.string.r1, ...)
          
          mod.step2 <-
            data.frame(
              modname = c(
                read.code2(search.space = search.space, mod.step2.string),
                read.code2(search.space = search.space, mod.string.r1)
              ),
              cmpt = c(mod.step2.string[1], mod.string.r1[1]),
              fitness = c(mod.test$fitness, mod.string.r1.fit),
              mod = c(toString(mod.test$mod), toString(mod.string.r1))
            )
          
          mod.step2$step = "IIV on parameters"
          mod.record <- rbind(mod.record, mod.step2)
          
          mod.step2 <-
            mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
          mod.record.best <- rbind(mod.record.best, mod.step2)
          model.tried. <- model.tried. + 1
        }
        
        
        if (mod.step2.string [4] != 1) {
          mod.string.r1 <- mod.step2.string
          
          mod.string.r1[4] <- 1
          
          mod.string.r1.fit <-
            sf.mod.run(modi = modi, string = mod.string.r1, ...)
          
          mod.step2 <-
            data.frame(
              modname = c(
                read.code2(search.space = search.space, mod.step2.string),
                read.code2(search.space = search.space, mod.string.r1)
              ),
              cmpt = c(mod.step2.string[1], mod.string.r1[1]),
              fitness = c(mod.step2$fitness, mod.string.r1.fit),
              mod = c(toString(mod.step2$mod), toString(mod.string.r1))
            )
          
          mod.step2$step = "IIV on parameters"
          mod.record <- rbind(mod.record, mod.step2)
          
          mod.step2 <-
            mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
          mod.record.best <- rbind(mod.record.best, mod.step2)
          model.tried. <- model.tried. + 1
        }
        
        
        if (mod.step2.string [6] != 1) {
          mod.string.r1 <- mod.step2.string
          
          mod.string.r1[6] <- 1
          
          mod.string.r1.fit <-
            sf.mod.run(modi = modi, string = mod.string.r1, ...)
          
          mod.step2 <-
            data.frame(
              modname = c(
                read.code2(search.space = search.space, mod.step2.string),
                read.code2(search.space = search.space, mod.string.r1)
              ),
              cmpt = c(mod.step2.string[1], mod.string.r1[1]),
              fitness = c(mod.test$fitness, mod.string.r1.fit),
              mod = c(toString(mod.test$mod), toString(mod.string.r1))
            )
          
          mod.step2$step = "IIV on parameters"
          mod.record <- rbind(mod.record, mod.step2)
          
          mod.step2 <-
            mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
          mod.record.best <- rbind(mod.record.best, mod.step2)
          model.tried. <- model.tried. + 1
        }
        
        
      }
      
      
    }
    
    
  }
  
  
  # 3Cmpt
  
  if (mod.step1$cmpt == 3) {
    mod.string.r1 <- mod.step1.string
    mod.string.r2 <- mod.step1.string
    mod.string.r3 <- mod.step1.string
    mod.string.r4 <- mod.step1.string
    mod.string.r5 <- mod.step1.string
    
    
    mod.string.r1[3] <- 1
    mod.string.r2[4] <- 1
    mod.string.r3[6] <- 1
    mod.string.r4[5] <- 1
    mod.string.r5[7] <- 1
    
    mod.string.r1.fit <-
      sf.mod.run(modi = modi, string = mod.string.r1, ...)
    mod.string.r2.fit <-
      sf.mod.run(modi = modi, string = mod.string.r2, ...)
    mod.string.r3.fit <-
      sf.mod.run(modi = modi, string = mod.string.r3, ...)
    mod.string.r4.fit <-
      sf.mod.run(modi = modi, string = mod.string.r4, ...)
    mod.string.r5.fit <-
      sf.mod.run(modi = modi, string = mod.string.r5, ...)
    
    mod.step2 <-
      data.frame(
        modname = c(
          read.code2(search.space = search.space, mod.step1.string),
          read.code2(search.space = search.space, mod.string.r1),
          read.code2(search.space = search.space, mod.string.r2),
          read.code2(search.space = search.space, mod.string.r3),
          read.code2(search.space = search.space, mod.string.r4),
          read.code2(search.space = search.space, mod.string.r5)
        ),
        
        cmpt = c(
          mod.step1.string[1],
          mod.string.r1[1],
          mod.string.r2[1],
          mod.string.r3[1],
          mod.string.r4[1],
          mod.string.r5[1]
        ),
        
        fitness = c(
          mod.step1.string.fit,
          mod.string.r1.fit,
          mod.string.r2.fit,
          mod.string.r3.fit,
          mod.string.r4.fit,
          mod.string.r5.fit
        ),
        
        mod = c(
          toString(mod.test$mod),
          toString(mod.string.r1),
          toString(mod.string.r2),
          toString(mod.string.r3),
          toString(mod.string.r4),
          toString(mod.string.r5)
        )
      )
    
    mod.step2$step = "IIV on parameters"
    mod.record <- rbind(mod.record, mod.step2)
    
    mod.step2 <-
      mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
    mod.record.best <- rbind(mod.record.best, mod.step2)
    model.tried. <- model.tried. + 5
    
    if (mod.step2$mod != mod.test$mod) {
      mod.test <- mod.step2
      
      mod.step2.string <-
        as.numeric(strsplit(mod.step2$mod, ",")[[1]])
      
      if (mod.step2.string [3] == 1) {
        mod.string.r1 <- mod.step2.string
        mod.string.r2 <- mod.step2.string
        mod.string.r3 <- mod.step2.string
        mod.string.r4 <- mod.step2.string
        
        mod.string.r1[4] <- 1
        mod.string.r2[6] <- 1
        mod.string.r3[5] <- 1
        mod.string.r4[7] <- 1
        
        mod.string.r1.fit <-
          sf.mod.run(modi = modi, string = mod.string.r1, ...)
        mod.string.r2.fit <-
          sf.mod.run(modi = modi, string = mod.string.r2, ...)
        mod.string.r3.fit <-
          sf.mod.run(modi = modi, string = mod.string.r3, ...)
        mod.string.r4.fit <-
          sf.mod.run(modi = modi, string = mod.string.r4, ...)
        
        mod.step2 <-
          data.frame(
            modname = c(
              read.code2(search.space = search.space, mod.step2.string),
              read.code2(search.space = search.space, mod.string.r1),
              read.code2(search.space = search.space, mod.string.r2),
              read.code2(search.space = search.space, mod.string.r3),
              read.code2(search.space = search.space, mod.string.r4)
            ),
            cmpt = c(
              mod.step2.string[1],
              mod.string.r1[1],
              mod.string.r2[1],
              mod.string.r3[1],
              mod.string.r4[1]
            ),
            
            fitness = c(
              mod.test$fitness,
              mod.string.r1.fit,
              mod.string.r2.fit,
              mod.string.r3.fit,
              mod.string.r4.fit
            ),
            
            mod = c(
              toString(mod.test$mod),
              toString(mod.string.r1),
              toString(mod.string.r2),
              toString(mod.string.r3),
              toString(mod.string.r4)
            )
          )
        
        mod.step2$step = "IIV on parameters"
        mod.record <- rbind(mod.record, mod.step2)
        
        mod.step2 <-
          mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
        mod.record.best <- rbind(mod.record.best, mod.step2)
        model.tried. <- model.tried. + 4
      }
      
      
      if (mod.step2.string [4] == 1) {
        mod.string.r1 <- mod.step2.string
        mod.string.r2 <- mod.step2.string
        mod.string.r3 <- mod.step2.string
        mod.string.r4 <- mod.step2.string
        
        mod.string.r1[3] <- 1
        mod.string.r2[6] <- 1
        mod.string.r3[5] <- 1
        mod.string.r4[7] <- 1
        
        mod.string.r1.fit <-
          sf.mod.run(modi = modi, string = mod.string.r1, ...)
        mod.string.r2.fit <-
          sf.mod.run(modi = modi, string = mod.string.r2, ...)
        mod.string.r3.fit <-
          sf.mod.run(modi = modi, string = mod.string.r3, ...)
        mod.string.r4.fit <-
          sf.mod.run(modi = modi, string = mod.string.r4, ...)
        
        mod.step2 <- data.frame(
          modname = c(
            read.code2(search.space = search.space, mod.step2.string),
            read.code2(search.space = search.space, mod.string.r1),
            read.code2(search.space = search.space, mod.string.r2),
            read.code2(search.space = search.space, mod.string.r3),
            read.code2(search.space = search.space, mod.string.r4)
          ),
          cmpt = c(
            mod.step2.string[1],
            mod.string.r1[1],
            mod.string.r2[1],
            mod.string.r3[1],
            mod.string.r4[1]
          ),
          
          fitness = c(
            mod.test$fitness,
            mod.string.r1.fit,
            mod.string.r2.fit,
            mod.string.r3.fit,
            mod.string.r4.fit
          ),
          
          mod = c(
            toString(mod.test$mod),
            toString(mod.string.r1),
            toString(mod.string.r2),
            toString(mod.string.r3),
            toString(mod.string.r4)
          )
        )
        
        mod.step2$step = "IIV on parameters"
        mod.record <- rbind(mod.record, mod.step2)
        
        mod.step2 <-
          mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
        mod.record.best <- rbind(mod.record.best, mod.step2)
        model.tried. <- model.tried. + 4
      }
      
      
      if (mod.step2.string [6] == 1) {
        mod.string.r1 <- mod.step2.string
        mod.string.r2 <- mod.step2.string
        mod.string.r3 <- mod.step2.string
        mod.string.r4 <- mod.step2.string
        
        mod.string.r1[3] <- 1
        mod.string.r2[4] <- 1
        mod.string.r3[5] <- 1
        mod.string.r4[7] <- 1
        
        mod.string.r1.fit <-
          sf.mod.run(modi = modi, string = mod.string.r1, ...)
        mod.string.r2.fit <-
          sf.mod.run(modi = modi, string = mod.string.r2, ...)
        mod.string.r3.fit <-
          sf.mod.run(modi = modi, string = mod.string.r3, ...)
        mod.string.r4.fit <-
          sf.mod.run(modi = modi, string = mod.string.r4, ...)
        
        mod.step2 <-
          data.frame(
            modname = c(
              read.code2(search.space = search.space, mod.step2.string),
              read.code2(search.space = search.space, mod.string.r1),
              read.code2(search.space = search.space, mod.string.r2),
              read.code2(search.space = search.space, mod.string.r3),
              read.code2(search.space = search.space, mod.string.r4)
            ),
            cmpt = c(
              mod.step2.string[1],
              mod.string.r1[1],
              mod.string.r2[1],
              mod.string.r3[1],
              mod.string.r4[1]
            ),
            
            fitness = c(
              mod.test$fitness,
              mod.string.r1.fit,
              mod.string.r2.fit,
              mod.string.r3.fit,
              mod.string.r4.fit
            ),
            
            mod = c(
              toString(mod.test$mod),
              toString(mod.string.r1),
              toString(mod.string.r2),
              toString(mod.string.r3),
              toString(mod.string.r4)
            )
          )
        
        mod.step2$step = "IIV on parameters"
        mod.record <- rbind(mod.record, mod.step2)
        
        mod.step2 <-
          mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
        mod.record.best <- rbind(mod.record.best, mod.step2)
        model.tried. <- model.tried. + 4
      }
      
      
      
      if (mod.step2.string [7] == 1) {
        mod.string.r1 <- mod.step2.string
        mod.string.r2 <- mod.step2.string
        mod.string.r3 <- mod.step2.string
        mod.string.r4 <- mod.step2.string
        
        mod.string.r1[3] <- 1
        mod.string.r2[4] <- 1
        mod.string.r3[6] <- 1
        mod.string.r4[5] <- 1
        
        mod.string.r1.fit <-
          sf.mod.run(modi = modi, string = mod.string.r1, ...)
        mod.string.r2.fit <-
          sf.mod.run(modi = modi, string = mod.string.r2, ...)
        mod.string.r3.fit <-
          sf.mod.run(modi = modi, string = mod.string.r3, ...)
        mod.string.r4.fit <-
          sf.mod.run(modi = modi, string = mod.string.r4, ...)
        
        mod.step2 <-
          data.frame(
            modname = c(
              read.code2(search.space = search.space, mod.step2.string),
              read.code2(search.space = search.space, mod.string.r1),
              read.code2(search.space = search.space, mod.string.r2),
              read.code2(search.space = search.space, mod.string.r3),
              read.code2(search.space = search.space, mod.string.r4)
            ),
            cmpt = c(
              mod.step2.string[1],
              mod.string.r1[1],
              mod.string.r2[1],
              mod.string.r3[1],
              mod.string.r4[1]
            ),
            
            fitness = c(
              mod.test$fitness,
              mod.string.r1.fit,
              mod.string.r2.fit,
              mod.string.r3.fit,
              mod.string.r4.fit
            ),
            
            mod = c(
              toString(mod.test$mod),
              toString(mod.string.r1),
              toString(mod.string.r2),
              toString(mod.string.r3),
              toString(mod.string.r4)
            )
          )
        
        mod.step2$step = "IIV on parameters"
        mod.record <- rbind(mod.record, mod.step2)
        
        mod.step2 <-
          mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
        mod.record.best <- rbind(mod.record.best, mod.step2)
        model.tried. <- model.tried. + 4
      }
      
      
      if (mod.step2$mod != mod.test$mod) {
        mod.test <- mod.step2
        
        mod.step2.string <-
          as.numeric(strsplit(mod.step2$mod, ",")[[1]])
        
        mod.string.r1 <- mod.step2.string
        mod.string.r2 <- mod.step2.string
        mod.string.r3 <- mod.step2.string
        
        count1 <- NULL
        
        for (check.cycle in 3:7) {
          if (mod.string.r1[check.cycle] == 1) {
            count1 <- c(count1, check.cycle)
          }
          
        }
        
        all.string <- c(3, 4, 5, 6, 7)
        all_string <-
          all.string[all.string != count1[1] &
                       all.string != count1[2]]
        
        mod.string.r1[all_string[1]] <- 1
        mod.string.r2[all_string[2]] <- 1
        mod.string.r3[all_string[3]] <- 1
        
        mod.string.r1.fit <-
          sf.mod.run(modi = modi, string = mod.string.r1, ...)
        mod.string.r2.fit <-
          sf.mod.run(modi = modi, string = mod.string.r2, ...)
        mod.string.r3.fit <-
          sf.mod.run(modi = modi, string = mod.string.r3, ...)
        
        mod.step2 <-
          data.frame(
            modname = c(
              read.code2(search.space = search.space, mod.step2.string),
              read.code2(search.space = search.space, mod.string.r1),
              read.code2(search.space = search.space, mod.string.r2),
              read.code2(search.space = search.space, mod.string.r3)
            ),
            
            cmpt = c(
              mod.step2.string[1],
              mod.string.r1[1],
              mod.string.r2[1],
              mod.string.r3[1]
            ),
            
            fitness = c(
              mod.test$fitness,
              mod.string.r1.fit,
              mod.string.r2.fit,
              mod.string.r3.fit
            ),
            
            mod = c(
              toString(mod.test$mod),
              toString(mod.string.r1),
              toString(mod.string.r2),
              toString(mod.string.r3)
            )
          )
        
        mod.step2$step = "IIV on parameters"
        mod.record <- rbind(mod.record, mod.step2)
        
        mod.step2 <-
          mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
        mod.record.best <- rbind(mod.record.best, mod.step2)
        model.tried. <- model.tried. + 3
        
        if (mod.step2$mod != mod.test$mod) {
          mod.step2.string <- as.numeric(strsplit(mod.step2$mod, ",")[[1]])
          
          mod.string.r1 <- mod.step2.string
          mod.string.r2 <- mod.step2.string
          
          count1 <- NULL
          
          for (check.cycle in 3:7) {
            if (mod.string.r1[check.cycle] == 1) {
              count1 <- c(count1, check.cycle)
            }
            
          }
          
          # all string
          all.string <- c(3, 4, 5, 6, 7)
          all_string <-
            all.string[all.string != count1[1] &
                         all.string != count1[2] &
                         all.string != count1[3]]
          
          mod.string.r1[all_string[1]] <- 1
          mod.string.r2[all_string[2]] <- 1
          
          
          mod.string.r1.fit <-
            sf.mod.run(modi = modi, string = mod.string.r1, ...)
          mod.string.r2.fit <-
            sf.mod.run(modi = modi, string = mod.string.r2, ...)
          
          
          mod.step2 <-
            data.frame(
              modname = c(
                read.code2(search.space = search.space, mod.step2.string),
                read.code2(search.space = search.space, mod.string.r1),
                read.code2(search.space = search.space, mod.string.r2)
              ),
              
              cmpt = c(mod.step2.string[1],
                       mod.string.r1[1],
                       mod.string.r2[1]),
              
              fitness = c(
                mod.test$fitness,
                mod.string.r1.fit,
                mod.string.r2.fit
              ),
              
              mod = c(
                toString(mod.test$mod),
                toString(mod.string.r1),
                toString(mod.string.r2)
              )
            )
          
          mod.step2$step = "IIV on parameters"
          mod.record <- rbind(mod.record, mod.step2)
          
          mod.step2 <-
            mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
          mod.record.best <- rbind(mod.record.best, mod.step2)
          model.tried. <- model.tried. + 2
          
          if (mod.step2$mod != mod.test$mod) {
            mod.step2.string <- as.numeric(strsplit(mod.step2$mod, ",")[[1]])
            
            if (mod.step2.string [3] != 1) {
              mod.string.r1 <- mod.step2.string
              
              mod.string.r1[3] <- 1
              mod.string.r1.fit <-
                sf.mod.run(modi = modi, string = mod.string.r1, ...)
              
              
              mod.step2 <-
                data.frame(
                  modname = c(
                    read.code2(search.space = search.space, mod.step2.string),
                    read.code2(search.space = search.space, mod.string.r1)
                  ),
                  
                  cmpt = c(mod.step2.string[1],
                           mod.string.r1[1]),
                  
                  fitness = c(mod.test$fitness, mod.string.r1.fit),
                  mod = c(
                    toString(mod.test$mod),
                    toString(mod.string.r1)
                  )
                )
              
              mod.step2$step = "IIV on parameters"
              mod.record <- rbind(mod.record, mod.step2)
              
              mod.step2 <-
                mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
              mod.record.best <- rbind(mod.record.best, mod.step2)
              
            }
            
            
            if (mod.step2.string [4] != 1) {
              mod.string.r1 <- mod.step2.string
              
              mod.string.r1[4] <- 1
              
              mod.string.r1.fit <-
                sf.mod.run(modi = modi, string = mod.string.r1, ...)
              
              mod.step2 <-
                data.frame(
                  modname = c(
                    read.code2(search.space = search.space, mod.step2.string),
                    read.code2(search.space = search.space, mod.string.r1)
                  ),
                  
                  cmpt = c(mod.step2.string[1],
                           mod.string.r1[1]),
                  
                  fitness = c(mod.test$fitness, mod.string.r1.fit),
                  mod = c(
                    toString(mod.test$mod),
                    toString(mod.string.r1)
                  )
                )
              
              mod.step2$step = "IIV on parameters"
              mod.record <- rbind(mod.record, mod.step2)
              
              mod.step2 <-
                mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
              mod.record.best <- rbind(mod.record.best, mod.step2)
              model.tried. <- model.tried. + 1
            }
            
            if (mod.step2.string [5] != 1) {
              mod.string.r1 <- mod.step2.string
              
              mod.string.r1[5] <- 1
              mod.string.r1.fit <-
                sf.mod.run(modi = modi, string = mod.string.r1, ...)
              
              
              mod.step2 <-
                data.frame(
                  modname = c(
                    read.code2(search.space = search.space, mod.step2.string),
                    read.code2(search.space = search.space, mod.string.r1)
                  ),
                  
                  cmpt = c(mod.step2.string[1],
                           mod.string.r1[1]),
                  
                  fitness = c(mod.test$fitness, mod.string.r1.fit),
                  mod = c(
                    toString(mod.test$mod),
                    toString(mod.string.r1)
                  )
                )
              
              mod.step2$step = "IIV on parameters"
              mod.record <- rbind(mod.record, mod.step2)
              
              mod.step2 <-
                mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
              mod.record.best <- rbind(mod.record.best, mod.step2)
              
            }
            
            if (mod.step2.string [6] != 1) {
              mod.string.r1 <- mod.step2.string
              
              mod.string.r1[6] <- 1
              
              mod.string.r1.fit <-
                sf.mod.run(modi = modi, string = mod.string.r1, ...)
              
              mod.step2 <-
                data.frame(
                  modname = c(
                    read.code2(search.space = search.space, mod.step2.string),
                    read.code2(search.space = search.space, mod.string.r1)
                  ),
                  
                  cmpt = c(mod.step2.string[1],
                           mod.string.r1[1]),
                  
                  fitness = c(mod.test$fitness, mod.string.r1.fit),
                  mod = c(
                    toString(mod.test$mod),
                    toString(mod.string.r1)
                  )
                )
              
              mod.step2$step = "IIV on parameters"
              mod.record <- rbind(mod.record, mod.step2)
              
              mod.step2 <-
                mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
              mod.record.best <- rbind(mod.record.best, mod.step2)
              model.tried. <- model.tried. + 1
            }
            
            
            if (mod.step2.string [7] != 1) {
              mod.string.r1 <- mod.step2.string
              
              mod.string.r1[7] <- 1
              
              mod.string.r1.fit <-
                sf.mod.run(modi = modi, string = mod.string.r1, ...)
              
              mod.step2 <-
                data.frame(
                  modname = c(
                    read.code2(search.space = search.space, mod.step2.string),
                    read.code2(search.space = search.space, mod.string.r1)
                  ),
                  
                  cmpt = c(mod.step2.string[1],
                           mod.string.r1[1]),
                  
                  fitness = c(mod.test$fitness, mod.string.r1.fit),
                  mod = c(
                    toString(mod.test$mod),
                    toString(mod.string.r1)
                  )
                )
              
              mod.step2$step = "IIV on parameters"
              mod.record <- rbind(mod.record, mod.step2)
              
              mod.step2 <-
                mod.step2[mod.step2$fitness == min(mod.step2$fitness), ]
              mod.record.best <- rbind(mod.record.best, mod.step2)
              model.tried. <- model.tried. + 1
            }
            
          }
          
          
        }
        
        
        
      }
      
      
    }
    
  }
  
  
  ######################## Step 3. Explore correlation. ##########################
  message(blue(
    paste0(
      "Test Correlation between parameters----------------------------------------------------"
    )
  ))
  
  r <<- 3
  mod.step2.string <- as.numeric(strsplit(mod.step2$mod, ",")[[1]])
  type <<- 3
  penalty <<- 10000
  
  
  if (mod.step2.string[9] == 0) {
    mod.step3.string <- mod.step2.string
    mod.step3.string.fit <-
      sf.mod.run(modi = modi, string = mod.step3.string, ...)
    mod.step2$fitness <- mod.step3.string.fit
    mod.step3.string[9] = 1
    
    mod.step3.string.fit <<-
      sf.mod.run(modi = modi, string = mod.step3.string, ...)
    
    if (mod.step3.string.fit < mod.step2$fitness) {
      mod.step3 <-
        data.frame(
          modname = read.code2(search.space = search.space, mod.step3.string),
          cmpt = c(mod.step3.string[1]),
          
          fitness = mod.step3.string.fit,
          mod = toString(mod.step3.string)
        )
      
      mod.step3$step = "Correlation between parameters"
      mod.record <- rbind(mod.record, mod.step3)
      
      mod.record.best <- rbind(mod.record.best, mod.step3)
      model.tried. <- model.tried. + 1
    }
    
    else{
      mod.step3 <- mod.step2
      
      mod.step3$step = "Correlation between parameters"
      mod.record <- rbind(mod.record, mod.step3)
      
      mod.record.best <- rbind(mod.record.best, mod.step3)
      model.tried. <- model.tried. + 1
    }
    
  }
  
  
  # Step 4. Explore RV model
  #################################################################################
  
  message(blue(
    paste0(
      "Explore types of residual errors----------------------------------------------------"
    )
  ))
  
  r <<- 4
  type <<- 3
  penalty <<- 10000
  
  mod.step3.string <- as.numeric(strsplit(mod.step3$mod, ",")[[1]])
  
  
  if (mod.step3.string[10] == 3) {
    mod.step4.string1 <- mod.step3.string
    mod.step4.string2 <- mod.step3.string
    
    mod.step4.string1[10] <- 1
    mod.step4.string2[10] <- 2
    
    mod.step4.string1.fit <-
      sf.mod.run(modi = modi, string = mod.step4.string1, ...)
    mod.step4.string2.fit <-
      sf.mod.run(modi = modi, string = mod.step4.string2, ...)
    
    
    mod.step4 <-
      data.frame(
        modname = c(
          read.code2(search.space = search.space, mod.step3.string),
          read.code2(search.space = search.space, mod.step4.string1),
          read.code2(search.space = search.space, mod.step4.string2)
        ),
        
        cmpt = c(
          mod.step3.string[1],
          mod.step4.string1[1],
          mod.step4.string2[1]
        ),
        
        fitness = c(
          mod.step3$fitness,
          mod.step4.string1.fit,
          mod.step4.string2.fit
        ),
        mod = c(
          toString(mod.step3.string),
          toString(mod.step4.string1),
          toString(mod.step4.string2)
        )
      )
    
  }
  
  mod.step4$step = "Residual error types"
  mod.record <- rbind(mod.record, mod.step4)
  
  mod.step4 <-
    mod.step4[mod.step4$fitness == min(mod.step4$fitness), ]
  
  mod.record.best <- rbind(mod.record.best, mod.step4)
  model.tried. <- model.tried. + 2
  
  mod.record$penalty.type = NA
  
  mod.record[mod.record$step == "no. of compartments", ]$penalty.type <-
    "RSE"
  
  mod.record[mod.record$step == "elimination type", ]$penalty.type <-
    "RSE, Shrinkage"
  
  if (nrow(mod.record[mod.record$step == "IIV on Km", ]) > 0) {
    mod.record[mod.record$step == "IIV on Km", ]$penalty.type <-
      "RSE, Shrinkage, RV model"
  }
  
  mod.record[mod.record$step == "IIV on parameters", ]$penalty.type <-
    "RSE, Shrinkage, RV model"
  mod.record[mod.record$step == "Correlation between parameters", ]$penalty.type <-
    "RSE, Shrinkage, RV model"
  mod.record[mod.record$step == "Residual error types", ]$penalty.type <-
    "RSE, Shrinkage, RV model"
  
  mod.record <- mod.record[, c(5, 1, 4, 3, 6)]
  colnames(mod.record) <-
    c("Step", "Model name", "Model code", "Fitness", "Penalty terms")
  
  mod.record.best$penalty.type = NA
  mod.record.best[mod.record.best$step == "no. of compartments", ]$penalty.type <-
    "RSE"
  mod.record.best[mod.record.best$step == "elimination type", ]$penalty.type <-
    "RSE"
  if (nrow(mod.record.best[mod.record.best$step == "IIV on Km", ]) > 0) {
    mod.record.best[mod.record.best$step == "IIV on Km", ]$penalty.type <-
      "RSE, Shrinkage"
  }
  
  mod.record.best[mod.record.best$step == "IIV on parameters", ]$penalty.type <-
    "RSE, Shrinkage"
  mod.record.best[mod.record.best$step == "Correlation between parameters", ]$penalty.type <-
    "RSE, Shrinkage, RV model"
  mod.record.best[mod.record.best$step == "Residual error types", ]$penalty.type <-
    "RSE, Shrinkage, RV model"
  
  mod.record.best <- mod.record.best[, c(5, 1, 4, 3, 6)]
  colnames(mod.record.best) <-
    c("Step", "Model name", "Model code", "Fitness", "Penalty terms")
  
  rownames(mod.record)<-c(seq(1,nrow(mod.record),1))
  rownames(mod.record.best)<-c(seq(1,nrow(mod.record.best),1))
  
  history.list <- list(local.best.models = mod.record.best,
                       allmodels = mod.record)
  
  finalstep = mod.record.best[mod.record.best$Step=="Residual error types", ]
  
  bestmodel = finalstep[finalstep$Fitness == min(finalstep$Fitness), ]
  
  return(list(bestmodel = bestmodel,
              history = history.list))
  
  
}
