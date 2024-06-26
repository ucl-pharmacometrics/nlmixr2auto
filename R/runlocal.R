#' Perform local exhaustive search for genetic algorithm
#'
#' Generates a population of chromosomes for local exhaustive search based on the best individual
#' from the previous iteration, and valuates the fitness of each new chromosome and ranks them.
#'
#' @param modi Current model iteration number.
#' @param ga.iter Current genetic algorithm iteration number.
#' @param dat Data frame containing the input data for model fitting.
#' @param search.space An integer indicating the search space.
#' @param sel.best.code The best chromosome from the previous iteration.
#' @param nbits Number of bits (parameters) in the chromosome.
#' @param sig.diff Significant difference threshold for ranking.
#' @param crse Constraint value for relative standard error.
#' @param cshrink Constraint value for shrinkage.
#' @param lbcl Lower bound for clearance (CL).
#' @param lbvc Lower bound for volume of central compartment (Vc).
#' @param lbvp Lower bound for volume of peripheral compartment (Vp).
#' @param lbq Lower bound for inter-compartmental clearance (Q).
#' @param lbvp2 Lower bound for volume of second peripheral compartment (Vp2).
#' @param lbq2 Lower bound for second inter-compartmental clearance (Q2).
#' @param cadd Constraint value for additive error.
#' @param cprop Constraint value for proportional error.
#' @param ccorr Constraint value for correlation.
#' @param penalty.type Type of penalty applied during fitness evaluation.
#' @param penalty.value Value of the penalty applied during fitness evaluation.
#' @param ... Additional arguments passed to `ga.mod.run`.
#'
#' @return A data frame containing the population of chromosomes generated during the local exhaustive search and their fitness values.
#' @import rxode2
#' @import nlmixr2
#' @import stringr
#' @examples
#' \dontrun{
#'
#' dat<-pheno_sd
#' ga.iter<-1
#' search.space<-1
#' sel.best.code <- data.frame(0,
#'                            0,
#'                             0,
#'                            0,
#'                            0,
#'                            0,
#'                            0,
#'                            0,
#'                            0,
#'                            0,
#'                            1,
#'                            1)
#' sig.diff<-1
#' nbits<-12
#' crse <- 20
#' cshrink <- 20
#' lbcl <- 0.1
#' lbvc <- 0.1
#' lbvp <- 0.1
#' lbq <- 0.1
#' lbvp2 <- 0.1
#' lbq2 <- 0.1
#' cadd <- 0.1
#' cprop <- 0.1
#' ccorr <- 0.1
#' penalty.type <- 3
#' penalty.value <- 10000
#'
#'  ls.population<-runlocal( ga.iter,
#'                        dat,
#'                        search.space,
#'                        sel.best.code,
#'                         nbits,
#'                        sig.diff,
#'                        crse,
#'                        cshrink,
#'                        lbcl,
#'                        lbvc,
#'                        lbvp,
#'                        lbq,
#'                        lbvp2,
#'                        lbq2,
#'                        cadd,
#'                        cprop,
#'                         ccorr,
#'                        penalty.type,
#'                        penalty.value,
#'                         control = saemControl(
#'                          seed = 1234,
#'                          print = 5,
#'                          nBurn = 20,
#'                          nEm = 30,
#'                          logLik = T,
#'                          rxControl = rxControl(cores = 2,
#'                                                maxsteps =70000)))
#'
#' }
#'
#' @export
runlocal <- function(ga.iter,
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
                     ...) {
  if (search.space == 1) {
    # Generate population chromonsomes for local exhaustive search
    ls.population <- NULL
    for (ls.cycle in 1:ncol(sel.best.code)) {
      if (sel.best.code[ls.cycle] == 0) {
        ls.population.s = sel.best.code
        ls.population.s[ls.cycle] = 1
      }
      if (sel.best.code[ls.cycle] == 1) {
        ls.population.s = sel.best.code
        ls.population.s[ls.cycle] = 0
      }
      # Three scenarios
      if (ls.cycle == 2) {
        if (sel.best.code[ls.cycle - 1] == 0) {
          ls.population.s = sel.best.code
          ls.population.s[ls.cycle] = 0
          ls.population.s[ls.cycle - 1] = 1
        }
      }
      # Three scenarios
      if (ls.cycle == 12) {
        if (sel.best.code[ls.cycle - 1] == 0) {
          ls.population.s = sel.best.code
          ls.population.s[ls.cycle] = 0
          ls.population.s[ls.cycle - 1] = 1
        }
      }
      ls.population.s <- de.code(ls.population.s, search.space)
      ls.population = rbind(ls.population, ls.population.s)
    }

    ls.population <- as.data.frame(ls.population)
    ls.population <- unique(ls.population)
    ls.population$fitness <- NA

    for (k in 1:nrow(ls.population)) {
      # format revision
      ls.population <- as.matrix(ls.population)
      # results<-try(ga.mod.run(modi=modi,
      #                         r = ga.iter,
      #                         dat=dat,
      #                         search.space=search.space,
      #                         string = as.numeric(ls.population[k,1:nbits]),
      #                         crse=crse,
      #                         cshrink=cshrink,
      #                         lbcl=lbcl,
      #                         lbvc=lbvc,
      #                         lbvp=lbvp,
      #                         lbq=lbq,
      #                         lbvp2=lbvp2,
      #                         lbq2=lbq2,
      #                         cadd=cadd,
      #                         cprop=cprop,
      #                         ccorr=ccorr,
      #                         penalty.type=penalty.type,
      #                         penalty.value=penalty.value,...),
      #              silent = T)

      results <- try(ga.mod.run(
        modi = modi,
        r = ga.iter,
        dat = dat,
        search.space = search.space,
        string = as.numeric(ls.population[k, 1:nbits]),
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
        control = saemControl(
          seed = 1234,
          print = 5,
          nBurn = 20,
          nEm = 30,
          logLik = T,
          rxControl = rxControl(cores = 22,
                                maxsteps = 70000)
        )
      ),
      silent = T)
      ls.population[k, (nbits + 1)] <- results
    }
    # print(ls.population)
    ls.population <- as.data.frame(ls.population)
    ls.population$rank <- rank_new(ls.population$fitness, sig.diff)

  } # close search space 1

  return(ls.population)


}
