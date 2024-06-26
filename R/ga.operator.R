#' Run genetic algorithm (GA) for automated model selection
#'
#' Perform a genetic algorithm to select the best model for pharmacokinetic modelling
#'
#' @param dat A data frame containing the pharmacokinetic data to be used for model fitting.
#' @param npopsize An integer specifying the size of the population.
#' @param max.iter An integer specifying the maximum number of iterations for the genetic algorithm.
#' @param prob.crossover The probability of crossover in the genetic algorithm.
#' @param prob.mutation The probability of mutation in the genetic algorithm.
#' @param nlocal.search The frequency of local exhaustive search iterations.
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
#' @return A matrix representing the best solution found by the genetic algorithm.
#' @import rxode2
#' @import nlmixr2
#' @import crayon
#' @import stringr
#' @examples
#' \dontrun{
#' library(nlmixr2autoinit)
#' library(nlmixr2autoGA)
#' d1<-pheno_sd
#' inits.out<-getppkinit(dat = d1,runnpd = 0)
#' autosetinit(dat = d1,
#'            inits= inits.out$Recommended_initial_estimates)
#' best.model<-ga.operator(dat = d1,
#'                         npopsize = 8,
#'                         max.iter = 6,
#'                         prob.crossover = 0.8,
#'                         prob.mutation = 0.2,
#'                         nlocal.search = 3,
#'                         sig.diff = 1,
#'                         search.space = 1,
#'                         no.cores = 4,
#'                         thetalower=c(vp=1,
#'                                      vp2=1),
#'
#'                         control = saemControl(
#'                           seed = 1234,
#'                           print = 5,
#'                           nBurn = 20,
#'                           nEm = 30,
#'                           logLik = T,
#'                          rxControl = rxControl(cores = 4,
#'                                                 maxsteps =70000)),
#'                         table=tableControl(cwres=T),
#'                         filename = "test",
#'                         foldername = "test"
#' )
#' best.model
#' }
#'
#' @export
ga.operator <- function(dat,
                        npopsize,
                        max.iter,
                        prob.crossover,
                        prob.mutation,
                        nlocal.search,
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
    paste0("GA_",
           current.date,
           "-",
           foldername,
           "_",
           digest::digest(dat),
           "_temp")

  # for linux ubuntu
  if (dir.exists(outputdir)) {
    print("Warning: current directory for GA analysis already exists")
    unlink(outputdir, recursive = T, force = T)
    print("Warning: a new one was created and replace the previous one")
    dir.create(outputdir, showWarnings = F, recursive = T)
  }

  if (!dir.exists(outputdir)) {
    print("Output directory for GA analysis is created")
    dir.create(outputdir, showWarnings = F, recursive = T)
  }
  setwd(paste0(getwd(), "/", outputdir))
  storage.path <- getwd()
  ################################ Default setting#################################
  if (missing(npopsize)) {
    npopsize = 10
  }
  if (missing(max.iter)) {
    max.iter = 15
  }

  if (missing(prob.crossover)) {
    prob.crossover = 0.8
  }

  if (missing(prob.mutation)) {
    prob.mutation = 0.2
  }

  if (missing(sig.diff)) {
    sig.diff = 1
  }

  if (missing(search.space)) {
    search.space = 1
  }

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

  sel.best.code.list = NULL
  data.pop.list = NULL
  sel.population.list = NULL
  ls.population.list = NULL
  children.cross.list = NULL
  children.mutation.list = NULL

  ##########Search space#########################################################
  # Search space (BASE). intravenous case, ,Michaelis–Menten elimination included.
  if (search.space == 1) {
    bit.names <- c(
      "cmpt.iv1",
      "cmpt.iv2",
      "eta.km",
      "eta.vc",
      "eta.vp",
      "eta.vp2",
      "eta.q",
      "eta.q2",
      "mm",
      "mcorr",
      "rv1",
      "rv2"
    )
    nbits <- length(bit.names)
  }


  ################################ GA running###################################
  ######Population generation.====================================================
  modi <<- 1 # Set i as the global variable
  setRxThreads(no.cores)
  for (ga.iter in 1:max.iter) {
    if (ga.iter == 1) {
      population = create.pop(npopsize, nbits)
    }
    # Use the children chromosomes as the new population in the next iterations
    else{
      population <- children.all
    }
    # Revise the column names to the model element names.
    colnames(population) <- c(bit.names)

    # Decode- interpret invalid model code into valid model code.
    for (decode.num in 1:nrow(population)) {
      population[decode.num,] <-
        de.code(string = population[decode.num,], search.space = search.space)
    }

    ######Fitness calculation.===================================================

    data.pop <- as.data.frame(population)
    data.pop$fitness <- NA

    for (k in 1:nrow(data.pop)) {
      results <- try(ga.mod.run(
        modi = modi,
        r = ga.iter,
        dat = dat,
        search.space = search.space,
        string = population[k, 1:nbits],
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

      # results<-try(ga.mod.run(  modi=modi,
      #                           r = ga.iter,
      #                           dat=dat,
      #                           search.space=search.space,
      #                           string=population[k,1:nbits],
      #                           crse=crse,
      #                           cshrink=cshrink,
      #                           lbcl=lbcl,
      #                           lbvc=lbvc,
      #                           lbvp=lbvp,
      #                           lbq=lbq,
      #                           lbvp2=lbvp2,
      #                           lbq2=lbq2,
      #                           cadd=cadd,
      #                           cprop=cprop,
      #                           ccorr=ccorr,
      #                           penalty.type=penalty.type,
      #                           penalty.value=penalty.value,
      #                           control = saemControl(
      #                             seed = 1234,
      #                             print = 5,
      #                             nBurn = 20,
      #                             nEm = 30,
      #                             logLik = T,
      #                             rxControl = rxControl(cores = 22,
      #                                                   maxsteps =70000))),
      #
      #              silent = T)
      data.pop[k,]$fitness <- results
    }

    # Sort out the fitness for each individuals
    # A customized ranking function is used where very small differences between numbers
    # (less than 1, as default setting) do not change their rank.
    # data.pop$fitness<-as.numeric(data.pop$fitness)
    data.pop$rank <- rank_new(data.pop$fitness, sig.diff)
    data.pop.list[[ga.iter]] <- data.pop
    # write.csv(data.pop, file=paste0(storage.path,'/pop.iter',ga.iter,'.fitness..csv'),row.names = F)

    ######Calculate global best and local best======================================
    # Local best, the best individual of the current iteration
    # Global best, the best individual including models run by local exhaustive search
    # Local best is used as reference model for the local exhaustive search.
    # Current elitism strategy of GA is to keep (global best) in each generation.
    data.pop$local.num <- seq(1, nrow(data.pop), 1)
    rownames(data.pop) <- seq(1, nrow(data.pop), 1)

    # If there is no local exhaustive search, then the local best is the global best
    # as well given of the elitism strategy (keep best always)
    sel.best <- data.pop[data.pop$rank == min(data.pop$rank),][1,]
    sel.best <- as.numeric (sel.best$local.num)
    sel.best.code = population[sel.best, , drop = FALSE]

    ######Local exhaustive search===================================================
    # Generate models based on the 1-bit rule, and the local best individual as reference
    # Binary code was used as reference
    if (ga.iter %% nlocal.search == 0) {
      ls.population <- runlocal(
        ga.iter,
        dat,
        search.space,
        sel.best.code,
        nbits,
        sig.diff,
        crse,
        cshrink,
        lbcl,
        lbvc,
        lbvp,
        lbq,
        lbvp2,
        lbq2,
        cadd,
        cprop,
        ccorr,
        penalty.type,
        penalty.value,
        ...
      )

      ls.population.list[[ga.iter]] <- ls.population
      # ls.population<-runlocal( ga.iter,
      #                         dat,
      #                         sel.best.code,
      #                         nbits,
      #                         sig.diff,
      #                         crse,
      #                         cshrink,
      #                         lbcl,
      #                         lbvc,
      #                         lbvp,
      #                         lbq,
      #                         lbvp2,
      #                         lbq2,
      #                         cadd,
      #                         cprop,
      #                         ccorr,
      #                         penalty.type,
      #                         penalty.value,
      #                         control = saemControl(
      #                           seed = 1234,
      #                           print = 5,
      #                           nBurn = 20,
      #                           nEm = 30,
      #                           logLik = T,
      #                           rxControl = rxControl(cores = 22,
      #                                                 maxsteps =70000)))
      if (min(data.pop$fitness) > min(ls.population$fitness)) {
        rownames(ls.population) <- seq(1, nrow(ls.population), 1)
        sel.best <-
          as.numeric (rownames(ls.population[ls.population$rank == min(ls.population$rank),][1,]))
        ls.population2 <- data.matrix(ls.population[, 1:nbits])
        sel.best.code = ls.population2[sel.best, , drop = FALSE]
      }
    }

    sel.best.code.list[[ga.iter]] <- sel.best.code

    ###### Selection================================================================
    # Tournament method used
    sel.population <-
      ga.sel.tournament(data.pop, population, npopsize, nbits)
    sel.population.list[[ga.iter]] <- sel.population

    ###### Crossover================================================================
    children.all <-
      ga.crossover(sel.population, prob.crossover, npopsize, nbits, bit.names)
    # iter0 means before the crossover.
    # write.csv(children.all,file=paste0('children.cross.iter0',ga.iter,'.csv'),row.names = F)
    children.cross.list[[ga.iter]] <- children.all

    #decode
    for (decode.num in 1:nrow(children.all)) {
      children.all[decode.num,] <-
        de.code(string = children.all[decode.num,], search.space = search.space)
    }



    ###### Mutation==================================================================
    children.all <-
      ga.mutation(children.all, prob.mutation, npopsize, nbits)
    children.mutation.list[[ga.iter]] <- children.all

    # covert to the valid model structure.
    for (decode.num in 1:nrow(children.all)) {
      children.all[decode.num,] <-
        de.code(string = children.all[decode.num,], search.space = search.space)
    }


    ######Elitism strategy===========================================================
    # Current elitism strategy of GA is to keep the best individual in that local iterations.
    children.all <-
      rbind(children.all[1:(nrow(children.all) - 1),], sel.best.code)
    rownames(children.all) <- seq(1, nrow(children.all), 1)


  }

  # analysis sel.best.code
  best_model_name <- read.code(1, sel.best.code)

  history_list <- list(
    sel.best.code.list = sel.best.code.list,
    data.pop.list = data.pop.list,
    sel.population.list = sel.population.list,
    ls.population.list = ls.population.list,
    children.cross.list = children.cross.list,
    children.mutation.list = children.mutation.list
  )

  return(
    list(
      final.selected.code = sel.best.code,
      final.selected.model = best_model_name,
      history = history_list
    )
  )
}
