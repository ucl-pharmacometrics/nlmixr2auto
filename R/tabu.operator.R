#' Run Tabu search for automated pharmacokinetic model selection
#'
#' Perform Tabu search to select the best model for pharmacokinetic modelling
#'
#' @param dat A data frame containing the pharmacokinetic data to be used for model fitting.
#' @param tabu.duration Integer specifying the number of iterations a target remains tabu list.
#' @param max.round Integer specifying the maximum number of rounds for the search.
#' @param start.point Named vector specifying a model as starting points for the search.
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
#' @return A list include the best solution and Tabu running history
#'
#' @examples
#' \dontrun{
#'
#' d1<-pheno_sd
#' inits.out<-getppkinit(dat = d1,runnpd = 0)
#' autosetinit(dat = d1,
#'             inits= inits.out$Recommended_initial_estimates)
#' tabu.operator(dat = d1,
#'               tabu.duration=2,
#'               max.round=4,
#'               search.space = 1,
#'               no.cores = 4,
#'               thetalower=c(vp=1,
#'                            vp2=1),
#'               control = saemControl(
#'                 seed = 1234,
#'                 print = 5,
#'                 nBurn = 20,
#'                 nEm = 30,
#'                 logLik = T,
#'                 rxControl = rxControl(cores = 4,
#'                                       maxsteps =70000)),
#'                table=tableControl(cwres=T),
#'                filename =  "pheno_sd",
#'                foldername =   "pheno_sd" )
#'
#'
#' }
#' @export

tabu.operator <- function(dat,
                          tabu.duration,
                          max.round,
                          start.point,
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
  #####################Create temporary storage folder for output###############
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
    paste0("TB_",
           current.date,
           "-",
           foldername,
           "_",
           digest::digest(dat),
           "_temp")
  
  # for linux ubuntu
  if (dir.exists(outputdir)) {
    print("Warning: current directory for Tabu analysis already exists")
    unlink(outputdir, recursive = T, force = T)
    print("Warning: a new one was created and replace the previous one")
    dir.create(outputdir, showWarnings = F, recursive = T)
  }
  
  if (!dir.exists(outputdir)) {
    print("Output directory for Tabu analysis is created")
    dir.create(outputdir, showWarnings = F, recursive = T)
  }
  setwd(paste0(getwd(), "/", outputdir))
  storage.path <- getwd()
  
  ################################ Default Tabu setting###########################
  
  if (missing(tabu.duration)) {
    tabu.duration <- 2
  }
  
  if (missing(max.round)) {
    max.round <- 30
  }
  
  
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
    dat.obs <- dat[dat$EVID == 0,]
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
  
  if (missing(penalty.type)) {
    penalty.type <- 3
  }
  
  if (missing(penalty.value)) {
    penalty.value <- 10000
  }
  
  
  
  r <<- 0
  modi <<- 1
  
  # Initial solution as the starting points
  if (missing(start.point)) {
    cmpt.iv = 2
    eta.km = 0
    eta.vc = 1
    eta.vp2 = 0
    eta.q2 = 0
    eta.vp = 0
    eta.q = 0
    mm = 0
    mcorr = 0
    rv  = 1
  }
  
  
  string <- c(cmpt.iv,
              eta.km,
              eta.vc,
              eta.vp,
              eta.vp2,
              eta.q,
              eta.q2,
              mm,
              mcorr,
              rv)
  
  fitness.values <- try(tabu.mod.run(
    modi = modi,
    r = r,
    dat = dat,
    search.space = search.space,
    string = string,
    crse = crse,
    cshrink = cshrink,
    lbcl = lbcl,
    lbvc = lbvc,
    lbvp = lbvp,
    lbq = lbq,
    lbvp2 = lbvp2,
    lbq2 = lbq2,
    cadd = cadd,
    cprop = cprop,
    ccorr = ccorr,
    penalty.type = penalty.type,
    penalty.value = penalty.value,
    ...
  ),
  
  silent = T)
  
  
  
  #  Start with an empty tabu list
  r <<- 1
  tabu.elements <- NULL
  tabu.elements.all <- NULL
  tabu.num <- 0
  local.best.all <<- NULL
  tabu.elements.history.list <- NULL
  starting.points.history.list <- NULL
  neighbors_df.list <- NULL
  
  # current local best model
  localbest <-
    Store.all[Store.all$fitness == min(Store.all$fitness),][1,]
  
  
  for (step.r in 1:max.round) {
    r <<- step.r
    tabu.num <- step.r
    # tabu.list<-Store.all
    localbest0 <- localbest
    
    if (r > 1) {
      local.store <- Store.all[Store.all$round.num == (r - 1),]
      localbest <-
        local.store[local.store$fitness == min(local.store$fitness),][1,]
    }
    
    compare. <-
      data.frame(t(rbind(localbest0[, (ncol(Store.all) - 9):ncol(Store.all)], localbest[, (ncol(Store.all) -
                                                                                             9):ncol(Store.all)])))
    colnames(compare.) <- c("local0", "local1")
    compare.$check <- compare.$local0 == compare.$local1
    
    if (nrow(compare.[compare.$check == F,]) > 0) {
      tabu.num <- tabu.num + 1
      tabu.elements <- data.frame(
        tabu.num = tabu.num,
        elements = rownames(compare.[compare.$check ==
                                       FALSE,]),
        elements.value = compare.[compare.$check == F,]$local1,
        tabu.iteration.left = tabu.duration
      )
      
    }
    
    if (is.null(tabu.elements.all) == F) {
      tabu.elements.all$tabu.iteration.left <-
        tabu.elements.all$tabu.iteration.left - 1
    }
    tabu.elements.all <- rbind(tabu.elements.all, tabu.elements)
    tabu.elements.all <-
      tabu.elements.all[tabu.elements.all$tabu.iteration.left > 0,]
    local.best.all <- rbind(local.best.all, localbest)
    
    tabu.elements.history.list [[r]] <- tabu.elements.all
  
    
    
    # current local best, also the starting points for the next generation
    current_string <- c(
      cmpt.iv = localbest[1, (ncol(localbest) - 9)],
      eta.km = localbest[1, (ncol(localbest) - 8)],
      eta.vc = localbest[1, (ncol(localbest) - 7)],
      eta.vp = localbest[1, (ncol(localbest) - 6)],
      eta.vp2 = localbest[1, (ncol(localbest) - 5)],
      eta.q = localbest[1, (ncol(localbest) - 4)],
      eta.q2 = localbest[1, (ncol(localbest) - 3)],
      mm = localbest[1, (ncol(localbest) - 2)],
      mcorr = localbest[1, (ncol(localbest) - 1)],
      rv = localbest[1, ncol(localbest)]
    )
    
    starting.points.history.list [[r]] <- current_string
    
    # Generate the all of neighbor models
    neighbors_df <- generate_neighbors_df(current_string)
    
    colnames(neighbors_df) <- c(
      "cmpt.iv",
      "eta.km",
      "eta.vc",
      "eta.vp",
      "eta.vp2",
      "eta.q",
      "eta.q2",
      "mm",
      "mcorr",
      "rv"
    )
    
    # Remove the invalid models
    neighbors_df$eta.vp[neighbors_df$cmpt.iv == 1] <- 0
    neighbors_df$eta.vp2[neighbors_df$cmpt.iv == 1] <- 0
    neighbors_df$eta.q[neighbors_df$cmpt.iv == 1] <- 0
    neighbors_df$eta.q2[neighbors_df$cmpt.iv == 1] <- 0
    
    neighbors_df$eta.vp2[neighbors_df$cmpt.iv == 2] <- 0
    neighbors_df$eta.q2[neighbors_df$cmpt.iv == 2] <- 0
    neighbors_df$eta.km[neighbors_df$mm == 0] <- 0
    
    neighbors_df$mcorr[neighbors_df$mm == 0 &
                         (
                           neighbors_df$eta.vc + neighbors_df$eta.vp + neighbors_df$eta.vp2 + neighbors_df$eta.q + neighbors_df$eta.q2
                         ) == 0] <- 0
    neighbors_df$mcorr[neighbors_df$mm == 1 &
                         neighbors_df$eta.km == 0 &
                         (
                           neighbors_df$eta.vc + neighbors_df$eta.vp + neighbors_df$eta.vp2 + neighbors_df$eta.q + neighbors_df$eta.q2
                         ) < 2] <- 0
    
    neighbors_df$flag <- 0
    
    # Exclude model in the tabu list
    for (cycle in 1:nrow(neighbors_df)) {
      string <- neighbors_df[cycle,]
      for (cycle2 in 1:nrow(local.best.all)) {
        if (toString(string) == toString(local.best.all[cycle2, (ncol(Store.all) - 9):ncol(Store.all)])) {
          neighbors_df[cycle,]$flag <- 1
        }

      }
    }
    
    #remove itselfs
    
    for (cycle in 1:nrow(neighbors_df)) {
     string <- neighbors_df[cycle,]
     if (current_string[1]==string[1] &
         current_string[2]==string[2] &
         current_string[3]==string[3] &
         current_string[4]==string[4] &
         current_string[5]==string[5] &
         current_string[6]==string[6] &
         current_string[7]==string[7] &
         current_string[8]==string[8] &
         current_string[9]==string[9] &
         current_string[10]==string[10] )
       neighbors_df[cycle,]$flag <- 1
    }
 
    
    neighbors_df <- neighbors_df[neighbors_df$flag == 0,]
    neighbors_df <- neighbors_df[, 1:10]
    
    
    # exclude the models/elements in the tabu list
    if (is.null(tabu.elements.all) == F) {
      # Remove the model, which includes the elements in the tabu.list
      for (tabu.elements.num in 1:nrow(tabu.elements.all)) {
        tabu.remove <- tabu.elements.all[tabu.elements.num,]$elements
        neighbors_df <-
          neighbors_df[neighbors_df[which(names(neighbors_df) == tabu.remove)] == tabu.elements.all[tabu.elements.num,]$elements.value,]
      }
    }
    
    # Remove the repeat models
    neighbors_df <- distinct(neighbors_df)
    neighbors_df.list [[r]] <- neighbors_df
    
    # Start local exhaustive search
    for (local.cycle in 1:nrow(neighbors_df)) {
      string <- neighbors_df[local.cycle,]
      string <- strsplit(toString(string), " ")[[1]]
      string <- as.numeric(gsub(",", "", string))
  
      fitness.values.s <- try( tabu.mod.run(
        modi = modi,
        r = r,
        dat = dat,
        search.space = search.space,
        string = string,
        crse = crse,
        cshrink = cshrink,
        lbcl = lbcl,
        lbvc = lbvc,
        lbvp = lbvp,
        lbq = lbq,
        lbvp2 = lbvp2,
        lbq2 = lbq2,
        cadd = cadd,
        cprop = cprop,
        ccorr = ccorr,
        penalty.type = penalty.type,
        penalty.value = penalty.value,
        ...
      ),silent = T)
      
      
    }
    
  }
  
  history.list = list(
    starting.points.history.list=starting.points.history.list,
    tabu.elements.history.list = tabu.elements.history.list,
    tried.neighbors = neighbors_df.list
  )
  
  localbestf <-
    Store.all[Store.all$fitness == min(Store.all$fitness),][1,]
  
  best_model_code <-
    toString(localbestf[, (ncol(Store.all) - 9):ncol(Store.all)])
  best_model_code <- strsplit(toString( best_model_code), " ")[[1]]
  best_model_code <- as.numeric(gsub(",", "",  best_model_code))
  best_model_name <- read.code2(1, best_model_code)
  
  return(
    list(
      best_model_name = best_model_name,
      best_model_code = best_model_code,
      history = history.list
    )
  )
}
