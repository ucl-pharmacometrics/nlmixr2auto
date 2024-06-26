#' Run and evaluate a model (stepwise)
#'
#' Conduct a model run and evaluate its fitness within the context of a genetic algorithm.
#' Converts binary strings to model parameters, fits the model using `nlmixr2`,
#' and evaluates the model's performance based on `fitness`.
#'
#' @param modi The current model number.
#' @param string A binary string representing the genetic algorithm's chromosome.
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
#' string <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1)
#' fitness <- sf.mod.run(modi=1, r = 1, dat = dat, search.space = 1, string = string,
#'                       crse = 20, cshrink = 30, lbcl = 1, lbvc = 1, lbvp = 1,
#'                       lbq = 0.01, lbvp2 = 1, lbq2 = 0.01, cadd = 0.1, cprop = 0.01,
#'                       ccorr = 0.2, penalty.type = 3, penalty.value = 10000,control = saemControl(logLik = T))
#' }
#'
#' @export


sf.mod.run <- function(modi,
                       string,
                        ...) {
  message(red(
    paste0(
      "Model: ",
      modi,
      "-------------------------------------------------------------------------------------------------------------"
    )
  ))
  
  if (file.exists(filename) == F) {
    # create a new one for the start
    Store.all <<- NULL
  }
  
  if (search.space == 1) {
    string <- de.code2(string, search.space = search.space)
    ex.mod.create.iv.mm(
      modi = modi,
      cmpt.iv = string[1],
      eta.km = string[2],
      eta.vc = string[3],
      eta.vp = string[4],
      eta.vp2 = string[5],
      eta.q = string[6],
      eta.q2 = string[7],
      mm = string[8],
      mcorr = string[9],
      rv = string[10]
    )
    
    source(file = paste0('mod', modi, '.txt'))
    
    
    ##############################Avoid running repeat model########################
    norow <- 1
    Store.flag <- 0
    if (modi > 1) {
      for (norow in 1:nrow(Store.all)) {
        if (toString(string) == toString(Store.all[norow, (ncol(Store.all) - 9):ncol(Store.all)])) {
          Store2 <- Store.all[norow,]
          Store2$model.num = modi
          Store2$current.time = Sys.time()
          Store2$round.num <- r
          Store.flag <- 1
          break
        }
      }
    }
    
    if (modi == 1 || Store.flag == 0) {
      loadError <<- FALSE
      fit.s <-
        suppressMessages(suppressWarnings(tryCatch(
          nlmixr2(f, dat, est = "saem", ...),
          error = function(e)
            loadError <<- TRUE
        )))
      
      
   # tryCatch(
   #      nlmixr2(f, dat, est = "saem", ...),
   #      error = function(e)
   #        loadError <<- TRUE
   # )
      ##############################Avoid running repeat model########################
      # Get the model parameters and performance
      Store. <- get.mod.lst(fit.s, modi)
      Store2 <- fitness(
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
  
  # Fitness evaluation by pre-defined rules.
  FitnessValue <- Store2$fitness
  Store.all <<- rbind(Store.all, Store2)
  
  write.csv(Store.all, file = filename, row.names = F)
  
  modi <<- modi + 1
  
  return(c(FitnessValue))
}
