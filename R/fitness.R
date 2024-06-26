#' Evaluate fitness of a population pharmacokinetic model
#'
#' Evaluate the fitness of model performance based on various constraints and penalty criteria.
#'
#' @param search.space An integer indicating the search space.
#' @param fit.results A data frame containing the results of the model fitting, obtained from `get.mod.lst`.
#' @param string A binary string representing the model.
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
#'
#' @return A data frame with updated fitness values and constraints flags.
#'
#' @examples
#' \dontrun{
#' # Define the model
#' pheno <- function() {
#'   ini({
#'     tcl <- log(0.008) # typical value of clearance
#'     tv <-  log(0.6)   # typical value of volume
#'     eta.cl + eta.v ~ c(1,
#'                        0.01, 1)
#'     add.err <- 0.1    # residual variability
#'   })
#'   model({
#'     cl <- exp(tcl + eta.cl) # individual value of clearance
#'     v <- exp(tv + eta.v)    # individual value of volume
#'     ke <- cl / v            # elimination rate constant
#'     d/dt(A1) = - ke * A1    # model differential equation
#'     cp = A1 / v             # concentration in plasma
#'     cp ~ add(add.err)       # define error model
#'   })
#' }
#' fit <- nlmixr(pheno, pheno_sd, "saem",control=list(print=0),
#' table=list(cwres=TRUE, npde=TRUE))
#'
#' # Fit the model using nlmixr2
#' fit.results <- get.mod.lst(fit,1)
#' string <- c(0, 0, 0, 1, 0, 0, 0, 0, 1, 1)
#' search.space<-1
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
#' fitness_results <- aco.fitness(search.space,fit.results, string, crse, cshrink, lbcl, lbvc, lbvp, lbq, lbvp2, lbq2, cadd, cprop, ccorr, penalty.type, penalty.value)
#' print(fitness_results)
#' }
#' @export
fitness <- function(search.space,
                       fit.results,
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
                       penalty.value) {
  if (search.space == 1) {
    lastcol0 <- ncol(fit.results)
    fit.results[, (lastcol0 + 1):(lastcol0 + 10)] <- 0
    
    colnames(fit.results)[(lastcol0 + 1):(lastcol0 + 10)] <-
      paste0("flag.", colnames(fit.results)[15:24])
    
    try(fit.results[fit.results$thetacl < lbcl &
                      is.na(fit.results$thetacl) == F,]$flag.thetacl <-
          1, silent = T)
    try(fit.results[fit.results$thetavc < lbvc &
                      is.na(fit.results$thetavc) == F,]$flag.thetavc <-
          1, silent = T)
    try(fit.results[fit.results$thetavp < lbvp &
                      is.na(fit.results$thetavp) == F,]$flag.thetavp <-
          1, silent = T)
    try(fit.results[fit.results$thetavp2 < lbvp2 &
                      is.na(fit.results$thetavp2) == F,]$flag.thetavp2 <-
          1, silent = T)
    try(fit.results[fit.results$thetaq < lbq &
                      is.na(fit.results$thetaq) == F,]$flag.thetaq <-
          1, silent = T)
    try(fit.results[fit.results$thetaq2 < lbq2 &
                      is.na(fit.results$thetaq2) == F,]$flag.thetaq2 <-
          1, silent = T)
    try(fit.results[fit.results$thetavmax < lbvmax &
                      is.na(fit.results$thetavmax) == F,]$flag.thetavmax <-
          1, silent = T)
    try(fit.results[fit.results$thetakm < lbkm &
                      is.na(fit.results$thetakm) == F,]$flag.thetakm <-
          1, silent = T)
    
    
    lastcol <- ncol(fit.results)
    fit.results[, (lastcol + 1):(lastcol + 10)] <- 0
    
    colnames(fit.results)[(lastcol + 1):(lastcol + 10)] <-
      paste0("flag.", colnames(fit.results)[25:34])
    
    try(fit.results[fit.results$rsecl > crse &
                      is.na(fit.results$rsecl) == F,]$flag.rsecl <-
          1, silent = T)
    try(fit.results[fit.results$rsevc > crse &
                      is.na(fit.results$rsevc) == F,]$flag.rsevc <-
          1, silent = T)
    try(fit.results[fit.results$rsevp > crse &
                      is.na(fit.results$rsevp) == F,]$flag.rsevp <-
          1, silent = T)
    try(fit.results[fit.results$rsevp2 > crse &
                      is.na(fit.results$rsevp2) == F,]$flag.rsevp2 <-
          1, silent = T)
    try(fit.results[fit.results$rseq > crse &
                      is.na(fit.results$rseq) == F,]$flag.rseq <-
          1, silent = T)
    try(fit.results[fit.results$rseq2 > crse &
                      is.na(fit.results$rseq2) == F,]$flag.rseq2 <-
          1, silent = T)
    try(fit.results[fit.results$rsevmax > crse &
                      is.na(fit.results$rsevmax) == F,]$flag.rsevmax <-
          1, silent = T)
    try(fit.results[fit.results$rsekm > crse &
                      is.na(fit.results$rsekm) == F,]$flag.rsekm <-
          1, silent = T)
    
    lastcol <- ncol(fit.results)
    fit.results[, (lastcol + 1):(lastcol + 10)] <- 0
    
    colnames(fit.results)[(lastcol + 1):(lastcol + 10)] <-
      paste0("flag.", colnames(fit.results)[45:54])
    
    try(fit.results[fit.results$shrinkcl > cshrink &
                      is.na(fit.results$shrinkcl) == F,]$flag.shrinkcl <-
          1, silent = T)
    try(fit.results[fit.results$shrinkvc > cshrink &
                      is.na(fit.results$shrinkvc) == F,]$flag.shrinkvc <-
          1, silent = T)
    try(fit.results[fit.results$shrinkvp > cshrink &
                      is.na(fit.results$shrinkvp) == F,]$flag.shrinkvp <-
          1, silent = T)
    try(fit.results[fit.results$shrinkvp2 > cshrink &
                      is.na(fit.results$shrinkvp2) == F,]$flag.shrinkvp2 <-
          1, silent = T)
    try(fit.results[fit.results$shrinkq > cshrink &
                      is.na(fit.results$shrinkq) == F,]$flag.shrinkq <-
          1, silent = T)
    try(fit.results[fit.results$shrinkq2 > cshrink &
                      is.na(fit.results$shrinkq2) == F,]$flag.shrinkq2 <-
          1, silent = T)
    try(fit.results[fit.results$shrinkvmax > cshrink &
                      is.na(fit.results$shrinkvmax) == F,]$flag.shrinkvmax <-
          1, silent = T)
    try(fit.results[fit.results$shrinkkm > cshrink &
                      is.na(fit.results$shrinkkm) == F,]$flag.shrinkkm <-
          1, silent = T)
    
    
    lastcol <- ncol(fit.results)
    fit.results[, (lastcol + 1):(lastcol + 16)] <- 0
    
    colnames(fit.results)[(lastcol + 1):(lastcol + 16)] <-
      paste0("flag.", colnames(fit.results)[101:116])
    
    try(fit.results[abs(fit.results$cor.eta.vc.cl) < ccorr &
                      is.na(fit.results$cor.eta.vc.cl) == F &
                      abs(fit.results$cor.eta.vc.cl) > 0 ,]$flag.cor.eta.vc.cl <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.vp.cl) < ccorr &
                      is.na(fit.results$cor.eta.vp.cl) == F &
                      abs(fit.results$cor.eta.vp.cl) > 0,]$flag.cor.eta.vp.cl <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.q.cl) < ccorr &
                      is.na(fit.results$cor.eta.q.cl) == F &
                      abs(fit.results$cor.eta.q.cl) > 0,]$flag.cor.eta.q.cl <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.vp2.cl) < ccorr &
                      is.na(fit.results$cor.eta.vp2.cl) == F &
                      abs(fit.results$cor.eta.vp2.cl) > 0,]$flag.cor.eta.vp2.cl <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.q2.cl) < ccorr &
                      is.na(fit.results$cor.eta.q2.cl) == F &
                      abs(fit.results$cor.eta.q2.cl) > 0,]$flag.cor.eta.q2.cl <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.vp.vc) < ccorr &
                      is.na(fit.results$cor.eta.vp.vc) == F &
                      abs(fit.results$cor.eta.vp.vc) > 0,]$flag.cor.eta.vp.vc <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.q.vc) < ccorr &
                      is.na(fit.results$cor.eta.q.vc) == F &
                      abs(fit.results$cor.eta.q.vc) > 0,]$flag.cor.eta.q.vc <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.vp2.vc) < ccorr &
                      is.na(fit.results$cor.eta.vp2.vc) == F &
                      abs(fit.results$cor.eta.vp2.vc) > 0,]$flag.cor.eta.vp2.vc <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.q2.vc) < ccorr &
                      is.na(fit.results$cor.eta.q2.vc) == F &
                      abs(fit.results$cor.eta.q2.vc) > 0,]$flag.cor.eta.q2.vc <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.q.vp) < ccorr &
                      is.na(fit.results$cor.eta.q.vp) == F &
                      abs(fit.results$cor.eta.q.vp) > 0,]$flag.cor.eta.q.vp <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.vp2.vp) < ccorr &
                      is.na(fit.results$cor.eta.vp2.vp) == F &
                      abs(fit.results$cor.eta.vp2.vp) > 0,]$flag.cor.eta.vp2.vp <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.q2.vp) < ccorr &
                      is.na(fit.results$cor.eta.q2.vp) == F &
                      abs(fit.results$cor.eta.q2.vp) > 0,]$flag.cor.eta.q2.vp <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.vp2.q) < ccorr &
                      is.na(fit.results$cor.eta.vp2.q) == F &
                      abs(fit.results$cor.eta.vp2.q) > 0,]$flag.cor.eta.vp2.q <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.q2.q) < ccorr &
                      is.na(fit.results$cor.eta.q2.q) == F  &
                      abs(fit.results$cor.eta.q2.q) > 0,]$flag.cor.eta.q2.q <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.q2.vp2) < ccorr &
                      is.na(fit.results$cor.eta.q2.vp2) == F &
                      abs(fit.results$cor.eta.q2.vp2) > 0,]$flag.cor.eta.q2.vp2 <-
          1, silent = T)
    try(fit.results[abs(fit.results$cor.eta.vmax.km) < ccorr &
                      is.na(fit.results$cor.eta.vmax.km) == F &
                      abs(fit.results$cor.eta.vmax.km) > 0,]$flag.cor.eta.vmax.km <-
          1, silent = T)
    
    
    lastcol <- ncol(fit.results)
    fit.results[, (lastcol + 1):(lastcol + 2)] <- 0
    
    colnames(fit.results)[(lastcol + 1):(lastcol + 2)] <-
      paste0("flag.", colnames(fit.results)[117:118])
    
    # string[no.] is based on the search space
    if (string[10] != 4) {
      try(fit.results[fit.results$add < cadd &
                        is.na(fit.results$add) == F,]$flag.add <-
            1, silent = T)
    }
    
    try(fit.results[fit.results$prop < cprop &
                      is.na(fit.results$prop) == F,]$flag.prop <-
          1, silent = T)
    
    fit.results$flag.covariance <- Inf
    try(fit.results[fit.results$model.covMethod %in% c("r", "s"), ]$flag.covariance <-
          penalty.value, silent = T)
    # For focei
    try(fit.results[fit.results$model.covMethod %in% c("r,s"), ]$flag.covariance <-
          0, silent = T)
    # For focei
    try(fit.results[fit.results$model.covMethod %in% c("linFim"), ]$flag.covariance <-
          0, silent = T)
    # For saem algorithm
    
    fit.results$no.theta.constraints <-
      rowSums(fit.results[, (lastcol0 + 1):(lastcol0 + 10)])
    fit.results$no.rse.constraints <-
      rowSums(fit.results[, (lastcol0 + 11):(lastcol0 + 20)])
    fit.results$no.shrink.constraints <-
      rowSums(fit.results[, (lastcol0 + 21):(lastcol0 + 30)])
    fit.results$no.corr.constraints <-
      rowSums(fit.results[, (lastcol0 + 31):(lastcol0 + 46)])
    fit.results$no.add.prop.constraints <-
      rowSums(fit.results[, (lastcol0 + 47):(lastcol0 + 48)])
    # Consideration of all constraints
    fit.results$no.constraints <-
      rowSums(fit.results[, (lastcol0 + 1):(ncol(fit.results) - 5)])
    
    # Static method: only calculate the number of violated constraints and regardless of its amount of constraint violation
    if (penalty.type == 1) {
      fit.results$fitness = fit.results$AIC + penalty.value * (
        fit.results$no.theta.constraints + fit.results$no.rse.constraints + fit.results$flag.covariance
      )
    }
    if (penalty.type == 2) {
      fit.results$fitness = fit.results$AIC + penalty.value * (
        fit.results$no.rse.constraints + fit.results$no.shrink.constraints + fit.results$no.corr.constraints + fit.results$flag.covariance
      )
    }
    if (penalty.type == 3) {
      fit.results$fitness = fit.results$AIC + penalty.value * fit.results$no.constraints
    }
    
  }
  return(fit.results)
}
