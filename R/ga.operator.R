
gaControl <- function(
    npopsize       = 10,
    max.iter       = 15,
    prob.crossover = 0.8,
    prob.mutation  = 0.2,
    sig.diff       = 1,
    nlocal.search  = 5
) {
  list(
    npopsize       = npopsize,
    max.iter       = max.iter,
    prob.crossover = prob.crossover,
    prob.mutation  = prob.mutation,
    sig.diff       = sig.diff,
    nlocal.search  = nlocal.search
  )
}


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
                        param_table=NULL,
                        search.space="ivbase",
                        no.cores = rxode2::getRxThreads(),
                        foldername = "test",
                        filename = "test",
                        ga.control=gaControl(),
                        penalty.control=penaltyControl(),
                        precomputed_results_file=NULL,
                        seed.no=1234,
                        ...) {
  #####################Create temporary storage folder for output###############
  current.date <- Sys.Date()
  set.seed(seed.no)

  outputdir <-
    paste0("GA_",
           current.date,
           "-",
           foldername,
           "_",
           digest::digest(dat),
           "_temp")

  if (!dir.exists(outputdir)) {
    dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
  } else {
    message(
      sprintf(
        "Output directory '%s' already exists. Using existing directory.",
        outputdir
      )
    )
  }

  npopsize       <- ga.control$npopsize
  max.iter       <- ga.control$max.iter
  prob.crossover <- ga.control$prob.crossover
  prob.mutation  <- ga.control$prob.mutation
  sig.diff       <- ga.control$sig.diff
  nlocal.search  <- ga.control$nlocal.search

  setwd(paste0(getwd(), "/", outputdir))
  storage.path <- getwd()
  #################################Initial estimate###############################

  param_table <- auto_param_table(
    dat = dat,
    param_table = param_table,
    nlmixr2autoinits = T,
    foldername = foldername
  )
  # Global variable for model and round number run index

  ##########Search Space Definiation ###########################################
  if (search.space == "ivbase") {
    bit.names <- c(
      "no.cmpt1",   # IV compartment code 1
      "no.cmpt2",   # IV compartment code 2
      "eta.km",     # Eta on Michaelis–Menten constant Km
      "eta.vc",     # Eta on central volume of distribution
      "eta.vp",     # Eta on peripheral volume 1
      "eta.vp2",    # Eta on peripheral volume 2
      "eta.q",      # Eta on intercompartmental clearance 1
      "eta.q2",     # Eta on intercompartmental clearance 2
      "mm",         # Michaelis–Menten elimination included
      "mcorr",      # Random effect correlation included
      "rv1",        # Residual error code 1
      "rv2"         # Residual error code 2
    )

  } else if (search.space == "oralbase") {
    bit.names <- c(
      "no.cmpt1",   # Oral compartment 1
      "no.cmpt2",   # Oral compartment 2
      "eta.km",     # Eta on Michaelis–Menten constant Km
      "eta.vc",     # Eta on central volume of distribution
      "eta.vp",     # Eta on peripheral volume 1
      "eta.vp2",    # Eta on peripheral volume 2
      "eta.q",      # Eta on intercompartmental clearance 1
      "eta.q2",     # Eta on intercompartmental clearance 2
      "eta.ka",     # Eta on absorption rate constant ka
      "mm",         # Michaelis–Menten elimination included
      "mcorr",      # Random effect correlation included
      "rv1",        # Residual error code 1
      "rv2"         # Residual error code 2
    )

  } else {
    stop("Unknown search.space type: must be 'ivbase' or 'oralbase'")
  }

  nbits <- length(bit.names)  # Total number of bits in the chromosome

  # --- GA Main Loop ---
  # Runs genetic algorithm over multiple generations:
  # 1) Init/update population
  # 2) Decode to valid models
  # 3) Evaluate fitness (mod.run + penalty.control)
  # 4) Local search every nlocal.search gens
  # 5) Selection → crossover → mutation
  # 6) Elitism: keep best model
  # 7) Save iteration results in history

  history <- vector("list", max.iter)  # Store results for each iteration

  for (ga.iter in 1:max.iter) {
    # 1. Initialize or update population
    if (ga.iter == 1) {
      population <-
        create.pop(npopsize, nbits)   # Initial random population
    } else {
      population <-
        children.all                  # Offspring from previous generation
    }

    colnames(population) <- bit.names

    population <- t(vapply(seq_len(nrow(population)),
                           function(i)
                             validateModels(
                               string = population[i,],
                               search.space = search.space,
                               code.source = "GA"
                             ),
                           numeric(nbits)))

    colnames(population) <- bit.names
    # 3. Model running and fitness evaluation
    data.pop <- as.data.frame(population)

    data.pop$fitness <- vapply(seq_len(nrow(data.pop)),
                               function(k) {
                                 result <- try(mod.run(
                                   r                = ga.iter,
                                   dat              = dat,
                                   search.space     = search.space,
                                   string           = parseCode(population[k, 1:nbits]),
                                   param_table      = param_table,
                                   penalty.control  = penalty.control,
                                   precomputed_results_file = precomputed_results_file,
                                   filename         = filename
                                 ),
                                 silent = TRUE)
                                 if (is.numeric(result) && length(result) == 1)
                                   result
                                 else
                                   NA_real_
                               },
                               numeric(1))

    # Rank individuals
    data.pop$rank <- rank_new(data.pop$fitness, sig.diff)

    # 4. Chromosome selection
    data.pop$local.num <- seq_len(nrow(data.pop))
    sel.best <- data.pop[data.pop$rank == min(data.pop$rank),][1,]
    sel.best <- as.numeric(sel.best$local.num)
    sel.best.code <- population[sel.best, , drop = FALSE]

    # 5. Local exhaustive search
    # ls.population <- NULL
    if (ga.iter %% nlocal.search == 0) {
      ls.population <- runlocal(
        ga.iter                   = ga.iter,
        dat                       = dat,
        search.space              = search.space,
        sel.best.code             = sel.best.code,
        sig.diff                  = sig.diff,
        param_table               = param_table,
        penalty.control           = penalty.control,
        precomputed_results_file  = precomputed_results_file,
        filename                  = filename,
        bit.names                 = bit.names,
        ...
      )
      if (min(data.pop$fitness, na.rm = TRUE) > min(ls.population$fitness, na.rm = TRUE)) {
        sel.best <-
          as.numeric(rownames(ls.population[ls.population$rank == min(ls.population$rank), ][1, ]))
        ls.population2 <- data.matrix(ls.population[, 1:nbits])
        sel.best.code <- ls.population2[sel.best, , drop = FALSE]
      }
    }

    # 6. Selection (tournament)
    sel.population <-
      ga.sel.tournament(data.pop = data.pop, population=population, npopsize=npopsize, nbits=nbits)

    # 7. Crossover
    children.cross <-
      ga.crossover(sel.population=sel.population,
                   prob.crossover=prob.crossover,
                   npopsize=npopsize,
                   nbits=nbits,
                   bit.names=bit.names)

    children.cross <- t(apply(children.cross, 1, function(x) {
      validateModels(
        string       = x,
        search.space = search.space,
        code.source  = "GA"
      )
    }))


    # 8. Mutation
    children.mutation <-
      ga.mutation(children.cross=children.cross,
                  prob.mutation=prob.mutation)

    children.mutation <- t(apply(children.mutation, 1, function(x) {
      validateModels(
        string       = x,
        search.space = search.space,
        code.source  = "GA"
      )
    }))


    # 9. Elitism: Keep best model in next generation
    children.all <- rbind(children.mutation[1:(nrow(children.mutation) - 1),],
                          sel.best.code)
    rownames(children.all) <- seq_len(nrow(children.all))


    # 10. Save results to history
    history[[ga.iter]] <- list(
      iteration              = ga.iter,
      population        = population,
      data.pop          = data.pop,
      sel.best.code     = sel.best.code,
      sel.population    = sel.population,
      ls.population     = ls.population,
      children.cross    = children.cross,
      children.mutation = children.all
    )
  }

  # ----------------------------
  # Final output
  # ----------------------------
  names(sel.best.code) <- bit.names
  best_model_name <- CodetoMod(sel.best.code = parseCode(sel.best.code,search.space=search.space),search.space=search.space)


  out <- new.env(parent = emptyenv())
  class(out) <- "gaOperatorResult"
  out[["Final Selected Code"]] <- sel.best.code
  out[["Final Selected Model Name"]] <- best_model_name
  out[["Model Run History"]] <- as.data.frame(Store.all, stringsAsFactors = FALSE)
  out[["Selection History"]] <- history

  on.exit({
  rm(modi, r, Store.all, precomputed_cache_loaded, envir = .GlobalEnv)
  }, add = TRUE)

  return(out)
}


#' Print method for gaOperatorResult objects
#'
#' Custom print method for results returned by the GA operator.
#' Displays only:
#'   - Final selected model code
#'   - Final selected model name
#'
#' @param x An object containing GA operator output (class gaOperatorResult).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#' @export
print.gaOperatorResult <- function(x, ...) {
  # Print final selected model code
  cat(crayon::green$bold("\n=== Final Selected Model Code ===\n"))
  print(x$`Final Selected Code`)

  # Print final selected model name
  cat(crayon::green$bold("\n=== Final Selected Model Name ===\n"))
  cat(x$`Final Selected Model Name`, "\n")

  invisible(x)
}













