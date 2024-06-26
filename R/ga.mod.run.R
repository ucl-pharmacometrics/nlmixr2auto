#' Run and evaluate a model (GA)
#'
#' Conduct a model run and evaluate its fitness within the context of a genetic algorithm.
#' Converts binary strings to model parameters, fits the model using `nlmixr2`,
#' and evaluates the model's performance based on `ga.fitness`.
#'
#' @param modi The current model number.
#' @param r An integer representing the current round number.
#' @param dat A data frame containing the data to be used for model fitting.
#' @param search.space An integer indicating the search space.
#' @param string A binary string representing the genetic algorithm's chromosome.
#' @param crse Constraint value for the relative standard error (RSE).
#' @param cshrink Constraint value for the shrinkage.
#' @param lbcl Lower bound for the clearance parameter.
#' @param lbvc Lower bound for the central volume of distribution parameter.
#' @param lbvp Lower bound for the peripheral volume of distribution parameter.
#' @param lbq Lower bound for the intercompartmental clearance parameter.
#' @param lbvp2 Lower bound for the second peripheral volume of distribution parameter.
#' @param lbq2 Lower bound for the second intercompartmental clearance parameter.
#' @param cadd Constraint value for the additive error.
#' @param cprop Constraint value for the proportional error.
#' @param ccorr Constraint value for the correlation between parameters.
#' @param penalty.type The type of penalty to apply during the fitness evaluation:
#' \itemize{
#'   \item \code{penalty.type = 1}: Considers constraints on RSE.
#'   \item \code{penalty.type = 2}: Considers constraints on RSE and shrinkage.
#'   \item \code{penalty.type = 3}: Considers constraints on RSE, shrinkage, and sigma values (variance of residual error).
#' }
#' All types also consider whether the model covariance step was successful.
#' @param penalty.value The value of the penalty to apply during the fitness evaluation.
#' @param ... Additional arguments to be passed to the `nlmixr2` function.
#'
#' @return A numeric value representing the fitness of the model.
#' @import rxode2
#' @import nlmixr2
#' @import crayon
#' @import stringr
#'
#' @examples
#' \dontrun{
#' dat <- pheno_sd
#' inits.out<-getppkinit(dat = dat)
#' autosetinit(dat = dat,
#'            inits= inits.out$Recommended_initial_estimates)
#' string <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
#' fitness <- ga.mod.run(modi=1, r = 1, dat = dat, search.space = 1, string = string,
#'                       crse = 20, cshrink = 30, lbcl = 1, lbvc = 1, lbvp = 1,
#'                       lbq = 0.01, lbvp2 = 1, lbq2 = 0.01, cadd = 0.1, cprop = 0.01,
#'                       ccorr = 0.2, penalty.type = 3, penalty.value = 10000,control = saemControl(logLik = T))
#' }
#'
#' @export


ga.mod.run <- function(modi,
                       r,
                       dat,
                       search.space,
                       string,
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
                       ...) {
  message(red(
    paste0(
      "Running GA ---------------------------------------------- ",
      "Model: ",
      modi,
      " Iteration: ",
      r
    )
  ))

  if (file.exists(filename) == F) {
    Store.all <- NULL
  }

  if (search.space == 1) {
    string <- ga.mod.create  (
      modi = modi,
      cmpt.iv1 = string[1],
      cmpt.iv2 = string[2],
      eta.km  = string[3],
      eta.vc   = string[4],
      eta.vp   = string[5],
      eta.vp2  = string[6],
      eta.q    = string[7],
      eta.q2   = string[8],
      mm       = string[9],
      mcorr    = string[10],
      rv1     = string[11],
      rv2     = string[12]
    )
    source(file = paste0('mod', modi, '.txt'))

    ##############################Avoid running repeat model########################
    norow <- 1
    Store.flag <- 0
    if (modi > 1) {
      for (norow in 1:nrow(Store.all)) {
        # 9 =nbins-1
        if (toString(string) == toString(Store.all[norow, (ncol(Store.all) -
                                                           9):ncol(Store.all)])) {
          Store2 <- Store.all[norow, ]
          Store2$model.num = modi
          Store2$current.time = Sys.time()
          Store2$round.num <- r
          Store.flag <- 1
          break
        }
      }
    }
    ##############################Avoid running repeat model########################
    # for the first one model or the model which was not run before
    if (modi == 1 || Store.flag == 0) {
      loadError <<- FALSE
      # fit.s <- tryCatch(nlmixr2(f,dat,est="saem",
      #                           control = saemControl(seed = seed.no,
      #                                                 print = 5,
      #                                                 nBurn = nBurn,
      #                                                 nEm = nEm,
      #                                                 logLik = T,
      #                                                 rxControl=rxControl(cores=no.cores,
      #                                                                     maxsteps =max.steps)),
      #                           table = tableControl(cwres = T)),
      #                   error=function(e) loadError <<- TRUE)

      fit.s <-
        suppressMessages(suppressWarnings(tryCatch(
          nlmixr2(f, dat, est = "saem", ...),
          error = function(e)
            loadError <<- TRUE
        )))
      ##############################Avoid running repeat model########################
      # Get the model parameters and performance
      Store. <- get.mod.lst(fit.s, modi)
      Store2 <- ga.fitness(
        search.space = search.space,
        fit.results = Store.,
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
        penalty.value = penalty.value
      )

      Store2$cmpt.iv = string[1]
      Store2$eta.km = string[2]
      Store2$eta.vc = string[3]
      Store2$eta.vp = string[4]
      Store2$eta.vp2 = string[5]
      Store2$eta.q = string[6]
      Store2$eta.q2 = string[7]
      Store2$mm = string[8]
      Store2$mcorr = string[9]
      Store2$rv = string[10]
      Store2$round.num <- r
      Store2 <- Store2[, c(ncol(Store2), 1:(ncol(Store2) - 1))]
    }
  } # close search.space=1



  # Fitness evaluation by rules.
  FitnessValue <- Store2$fitness
  Store.all <<- rbind(Store.all, Store2)

  write.csv(Store.all, file = filename, row.names = F)

  modi <<- modi + 1

  return(c(FitnessValue))
}
