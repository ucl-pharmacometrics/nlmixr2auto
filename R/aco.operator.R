#' Run ant colony optimisation for automated model selection
#'
#' Perform an ant colony optimisation (ACO) algorithm to select the best model for pharmacokinetic modelling
#'
#' @param dat A data frame containing the pharmacokinetic data to be used for model fitting.
#' @param no.ants An integer specifying the number of ants to use. Default is 20.
#' @param max.iter An integer specifying the maximum number of iterations. Default is 15.
#' @param rho A numeric value for the pheromone evaporation rate. Default is 0.2.
#' @param initial.phi A numeric value for the initial pheromone level. Default is 1.
#' @param lower.limit.phi A numeric value for the lower limit of the pheromone level. Default is 1.
#' @param upper.limit.phi A numeric value for the upper limit of the pheromone level. Default is Inf.
#' @param alpha.value A numeric value for the pheromone influence parameter. Default is 1.
#' @param elitism.percentage A numeric value specifying the percentage of elitism. Default is 1/3
#' @param sig.diff The significance difference threshold for ranking.
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
#' @return A list include the best solution and ACO running history
#'
#' @examples
#' \dontrun{
#' d1<-pheno_sd
#' inits.out<-getppkinit(dat = d1,runnpd = 0)
#' autosetinit(dat = d1,
#'            inits= inits.out$Recommended_initial_estimates)
#' result <- aco.operator(
#'     dat = d1,
#'     no.ants = 6,
#'     max.iter = 6,
#'     rho = 0.2,
#'     sig.diff = 1,
#'     search.space = 1,
#'     no.cores = 4,
#'     thetalower = c(vp = 1, vp2 = 1),
#'     control = saemControl(
#'     seed = 1234,
#'     print = 5,
#'     nBurn = 20,
#'     nEm = 30,
#'     logLik = TRUE,
#'     rxControl = rxControl(cores = 4, maxsteps = 70000)
#'    ),
#'    table = tableControl(cwres = TRUE),
#'    filename = "pheno_sd",
#'    foldername = "pheno_sd"
#'  )
#' result
#' }
#'
#' @export



aco.operator <- function(dat,
                         no.ants,
                         max.iter,
                         rho,
                         initial.phi,
                         lower.limit.phi,
                         upper.limit.phi,
                         alpha.value,
                         elitism.percentage,
                         sig.diff,
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
    paste0("ACO_",
           current.date,
           "-",
           foldername,
           "_",
           digest::digest(dat),
           "_temp")
  
  # for linux ubuntu
  if (dir.exists(outputdir)) {
    print("Warning: current directory for ACO analysis already exists")
    unlink(outputdir, recursive = T, force = T)
    print("Warning: a new one was created and replace the previous one")
    dir.create(outputdir, showWarnings = F, recursive = T)
  }
  
  if (!dir.exists(outputdir)) {
    print("Output directory for ACO analysis is created")
    dir.create(outputdir, showWarnings = F, recursive = T)
  }
  setwd(paste0(getwd(), "/", outputdir))
  storage.path <- getwd()
  
  ################################ Default ACO setting###########################
  if (missing(no.ants)) {
    no.ants = 20
  }
  if (missing(max.iter)) {
    max.iter = 15
  }
  
  if (missing(rho)) {
    rho = 0.2
  }
  
  if (missing(initial.phi)) {
    initial.phi = 1
  }
  
  if (missing(lower.limit.phi)) {
    lower.limit.phi = 1
  }
  
  if (missing(upper.limit.phi)) {
    upper.limit.phi = Inf
  }
  
  if (missing(alpha.value)) {
    alpha.value = 1
  }
  
  if (missing(sig.diff)) {
    sig.diff = 1
  }
  
  if (missing(search.space)) {
    search.space <- 1
  }
  
  if (missing(elitism.percentage)) {
    elitism.percentage <- 0.33
  }
  
  if (search.space == 1) {
    no.nodes = 22
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
  
  if (missing(penalty.type)) {
    penalty.type <- 3
  }
  
  if (missing(penalty.value)) {
    penalty.value <- 10000
  }
  
  
  cycle.all.list = NULL
  
  ############################### Initial conditions#############################
  node.list.all <- data.frame(
    travel = 0,
    node.no = seq(1, no.nodes, 1),
    local.node.no = c(
      seq(1, 3, 1),
      seq(0, 1, 1),
      seq(0, 1, 1),
      seq(0, 1, 1),
      seq(0, 1, 1),
      seq(0, 1, 1),
      seq(0, 1, 1),
      seq(0, 1, 1),
      seq(0, 1, 1),
      seq(1, 3, 1)
    ),
    node.names = c(
      "1Cmpt",
      "2Cmpt",
      "3Cmpt",
      "eta.vp2.no",
      "eta.vp2.yes",
      "eta.q2.no",
      "eta.q2.yes",
      "eta.vp.no",
      "eta.vp.yes",
      "eta.q.no",
      "eta.q.yes",
      "eta.vc.no",
      "eta.vc.yes",
      "mm.no",
      "mm.yes",
      "eta.km.no",
      "eta.km.yes",
      "mcorr.no",
      "mcorr.yes",
      "add",
      "prop",
      "comb"
    ),
    node.group = c(rep(1, 3),
                   sort(rep(seq(
                     2, 9, 1
                   ), 2)),
                   rep(10, 3)),
    phi = rep(initial.phi, no.nodes),
    delta_phi = rep(0, no.nodes),
    p = c(rep(round(1 / 3, 3), 3),
          rep(0.5, 16),
          rep(round(1 / 3, 3), 3))
  )
  
  node.list.0 <<- node.list.all
  pre.name <- colnames(node.list.0)
  
  # generate the initial ants for the first travel.
  initial.all <- create.ant(
    search.space = search.space,
    no.ants = no.ants,
    initialize = T,
    node.list = node.list.0
  )
  
  cycle.all.list <- NULL
  cycle.all.list [[1]] <- initial.all
  
  
  modi <<- 1 # Set i as the global variable
  setRxThreads(no.cores)
  #########################Run ACO################################################
  # First travel
  bestmodelcode <- NULL
  aco.iter <- 1
  
  #single individual
  for (si in 1:ncol(initial.all)) {
    run.dat <- initial.all
    cmpt.iv = run.dat[1, si]
    eta.vp2 = run.dat[2, si]
    eta.q2 = run.dat[3, si]
    eta.vp = run.dat[4, si]
    eta.q = run.dat[5, si]
    eta.vc = run.dat[6, si]
    mm = run.dat[7, si]
    eta.km = run.dat[8, si]
    mcorr = run.dat[9, si]
    rv  = run.dat[10, si]
    
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
    
    # change -1 to 0
    string <- pmax(string, 0)
    
    results <- try(aco.mod.run(
      modi = modi,
      r = aco.iter,
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
  }
  
  
  # Part4: Calculate phi for each ants
  
  node.list.1 <- phi.calculate(
    r = 1,
    search.space = search.space,
    input.ants.dat = initial.all,
    all.lib.dat = Store.all,
    # all fitness
    node.list.0 = node.list.0,
    node.list.all = node.list.all,
    alpha.value = alpha.value,
    rho = rho,
    sig.diff = sig.diff,
    lower.limit.phi = lower.limit.phi,
    upper.limit.phi = upper.limit.phi
  )
  
  node.list.s <-
    p.calculation(search.space = 1, node.list.cal = node.list.1)
  node.list.all <- rbind(node.list.all, node.list.s)
  
  
  for (aco.iter in 2:max.iter) {
    r <<- aco.iter
    bestmodel <- Store.all[Store.all$fitness == min(Store.all$fitness), ]
    bestmodelcode <- bestmodel[1, (ncol(Store.all) - 9):ncol(Store.all)]
    
    cycle.all <- create.ant(
      search.space = search.space,
      no.ants = no.ants,
      initialize = F,
      node.list = node.list.s
    )
    
    cycle.all.list [[r]] <- cycle.all
    
    # elitism
    #==============================================================================#
    # keep a fixed number of percentage of elitism
    if (search.space == 1) {
      bestmodelcode.r <- c(
        bestmodelcode[1],
        bestmodelcode[5],
        bestmodelcode[7],
        bestmodelcode[4],
        bestmodelcode[6],
        bestmodelcode[3],
        bestmodelcode[8],
        bestmodelcode[2],
        bestmodelcode[9],
        bestmodelcode[10]
      )
      
      no.elitism <- round(no.ants * elitism.percentage, 0)
      # at least one
      no.elitism <- pmax(no.elitism, 1)
      for (loop.no.elitism in 1:no.elitism) {
        cycle.all[1:10, (no.ants - loop.no.elitism + 1)] <-
          as.numeric(bestmodelcode.r)
      }
    }
    
    
    #==============================================================================#
    for (si in 1:ncol(cycle.all)) {
      run.dat <- cycle.all
      cmpt.iv = run.dat[1, si]
      eta.vp2 = run.dat[2, si]
      eta.q2 = run.dat[3, si]
      eta.vp = run.dat[4, si]
      eta.q = run.dat[5, si]
      eta.vc = run.dat[6, si]
      mm = run.dat[7, si]
      eta.km = run.dat[8, si]
      mcorr = run.dat[9, si]
      rv  = run.dat[10, si]
      
      
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
      
      # change -1 to 0
      string <- pmax(string, 0)
      
      results <- try(aco.mod.run(
        modi = modi,
        r = aco.iter,
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
    }
    
    node.list.s <- phi.calculate(
      r = r,
      search.space = search.space,
      input.ants.dat = cycle.all,
      all.lib.dat = Store.all,
      # all fitness
      node.list.0 = node.list.0,
      node.list.all = node.list.all,
      alpha.value = alpha.value,
      rho = rho,
      sig.diff = sig.diff,
      lower.limit.phi = lower.limit.phi,
      upper.limit.phi = upper.limit.phi
    )
    
    node.list.s <-
      p.calculation(search.space = 1, node.list.cal = node.list.s)
    node.list.all <- rbind(node.list.all, node.list.s)
    
  }
  
  best_model_name <- read.code2(1, bestmodelcode)
  
  history.list <- list(node.list.all = node.list.all,
                       ant.travel.list = cycle.all.list)
  
  return(
    list(
      bestmodelcode = bestmodelcode,
      best_model_name = best_model_name,
      history.list = history.list
    )
  )
  
}
