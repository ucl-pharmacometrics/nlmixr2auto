#' Extract model parameters and statistics
#'
#' Extract model parameters and statistics from a fitted `nlmixr2` model object.
#'
#' @param modi The current model number.
#' @param fit.s A fitted `nlmixr2` model object.
#' @return A data frame containing the extracted model parameters and statistics.
#'
#' @import stringr
#' @examples
#' \dontrun{
#'
#' pheno <- function() {
#'   ini({
#'     tcl <- log(0.008) # typical value of clearance
#'     tv <-  log(0.6)   # typical value of volume
#'     eta.cl + eta.v ~ c(1,
#'                        0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
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
#'
#' # Fit the model using nlmixr2
#' fit <- nlmixr2(one.compartment, pheno_sd, est="saem", saemControl(print=0))
#'
#' # Extract model results
#' model_results <- get.mod.lst(fit,1)
#' print(model_results)
#' }
#'
#' @export

get.mod.lst <- function(fit.s,
                        modi) {
  current.time <- Sys.time()
  AIC <- Inf
  BIC <- Inf
  OBJFV <- Inf
  ll <- NA
  npar <- NA

  model.covMethod <- NA
  model.time.setup <- NA
  model.time.covariance <- NA
  model.time.saem <- NA
  model.time.table <- NA
  model.time.compress <- NA
  model.time.other <- NA

  thetaka <- NA
  thetacl <- NA
  thetavc <- NA
  thetavp <- NA
  thetavp2 <- NA
  thetaq <- NA
  thetaq2 <- NA
  thetatlag <- NA
  thetavmax <- NA
  thetakm <- NA

  rseka <- NA
  rsecl <- NA
  rsevc <- NA
  rsevp <- NA
  rsevp2 <- NA
  rseq <- NA
  rseq2 <- NA
  rsetlag <- NA
  rsevmax <- NA
  rsekm <- NA

  bsvka <- NA
  bsvcl <- NA
  bsvvc <- NA
  bsvvp <- NA
  bsvvp2 <- NA
  bsvq <- NA
  bsvq2 <- NA
  bsvtlag <- NA
  bsvvmax <- NA
  bsvkm <- NA

  shrinkka <- NA
  shrinkcl <- NA
  shrinkvc <- NA
  shrinkvp <- NA
  shrinkvp2 <- NA
  shrinkq <- NA
  shrinkq2 <- NA
  shrinktlag <- NA
  shrinkvmax <- NA
  shrinkkm <- NA

  CIlowerka <- NA
  CIlowercl <- NA
  CIlowervc <- NA
  CIlowervp <- NA
  CIlowervp2 <- NA
  CIlowerq <- NA
  CIlowerq2 <- NA
  CIlowertlag <- NA
  CIlowervmax <- NA
  CIlowerkm <- NA

  CIupperka <- NA
  CIuppercl <- NA
  CIuppervc <- NA
  CIuppervp <- NA
  CIuppervp2 <- NA
  CIupperq <- NA
  CIupperq2 <- NA
  CIuppertlag <- NA
  CIuppervmax <- NA
  CIupperkm <- NA

  omegaka <- NA
  omegacl <- NA
  omegavc <- NA
  omegavp <- NA
  omegavp2 <- NA
  omegaq <- NA
  omegaq2 <- NA
  omegatlag <- NA
  omegavmax <- NA
  omegakm <- NA

  omega.vc.cl <- NA
  omega.vp.cl <- NA
  omega.q.cl <- NA
  omega.vp2.cl <- NA
  omega.q2.cl <- NA
  omega.vp.vc <- NA
  omega.q.vc <- NA

  omega.vp2.vc <- NA
  omega.q2.vc <- NA
  omega.q.vp <- NA
  omega.vp2.vp <- NA
  omega.q2.vp <- NA
  omega.vp2.q <- NA
  omega.q2.q <- NA
  omega.q2.vp2 <- NA
  omega.vmax.km <- NA

  cor.eta.vc.cl <- NA
  cor.eta.vp.cl <- NA
  cor.eta.q.cl <- NA
  cor.eta.vp2.cl <- NA
  cor.eta.q2.cl <- NA
  cor.eta.vp.vc <- NA
  cor.eta.q.vc <- NA

  cor.eta.vp2.vc <- NA
  cor.eta.q2.vc <- NA
  cor.eta.q.vp <- NA
  cor.eta.vp2.vp <- NA
  cor.eta.q2.vp <- NA
  cor.eta.vp2.q <- NA
  cor.eta.q2.q <- NA
  cor.eta.q2.vp2 <- NA
  cor.eta.vmax.km <- NA

  add <- NA
  prop <- NA

  if (length(fit.s) > 1) {
    parlst <- as.data.frame(fit.s$parFixedDf)
    parlst$param.names <- rownames(parlst)

    omegalst <- as.data.frame(fit.s$omega)
    omegalst$omega.names <- rownames(omegalst)

    if (is.null(fit.s$objDf[2]$AIC) == F) {
      AIC <- fit.s$objDf[2][1,]
    }

    if (is.null(fit.s$objDf[3]$BIC) == F) {
      BIC <- fit.s$objDf[3][1,]
    }

    if (is.null(fit.s$objDf[1]$OBJF) == F) {
      OBJFV <- fit.s$objDf[1][1,]
    }

    if (is.null(fit.s$objDf[4]$`Log-likelihood`) == F) {
      ll <- fit.s$objDf[4][1,]
      npar <- (AIC - 2 * (-ll)) / 2
    }


    # read message of saem
    if (is.null(fit.s$covMethod) == F) {
      model.covMethod <- fit.s$covMethod
    }

    if (is.null(fit.s$time$setup) == F) {
      model.time.setup <- fit.s$time$setup
    }

    if (is.null(fit.s$time$covariance) == F) {
      model.time.covariance <- fit.s$time$covariance
    }

    if (is.null(fit.s$time$saem) == F) {
      model.time.saem <- fit.s$time$saem
    }

    if (is.null(fit.s$time$table) == F) {
      model.time.table <- fit.s$time$table
    }

    if (is.null(fit.s$time$compress) == F) {
      model.time.compress <- fit.s$time$compress
    }

    if (is.null(fit.s$time$other) == F) {
      model.time.other <- fit.s$time$other
    }

    # sometimes there was no 95CI
    if (nrow(parlst[parlst$param.names == "lcl", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "lcl", ]$`Back-transformed`) == F) {
        thetacl <- parlst[parlst$param.names == "lcl", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "lcl", ]$`%RSE`) == F) {
        rsecl <- parlst[parlst$param.names == "lcl", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "lcl", ]$`BSV(CV%)`) == F) {
        bsvcl <- parlst[parlst$param.names == "lcl", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "lcl", ]$`Shrink(SD)%`) == F) {
        shrinkcl <- parlst[parlst$param.names == "lcl", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "lcl", ]$`CI Lower`) == F) {
        CIlowercl <- parlst[parlst$param.names == "lcl", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "lcl", ]$`CI Upper`) == F) {
        CIuppercl <- parlst[parlst$param.names == "lcl", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "lvmax", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "lvmax", ]$`Back-transformed`) ==
          F) {
        thetavmax <-
          parlst[parlst$param.names == "lvmax", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "lvmax", ]$`%RSE`) == F) {
        rsevmax <- parlst[parlst$param.names == "lvmax", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "lvmax", ]$`BSV(CV%)`) == F) {
        bsvvmax <- parlst[parlst$param.names == "lvmax", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "lvmax", ]$`Shrink(SD)%`) ==
          F) {
        shrinkvmax <- parlst[parlst$param.names == "lvmax", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "lvmax", ]$`CI Lower`) == F) {
        CIlowervmax <- parlst[parlst$param.names == "lvmax", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "lvmax", ]$`CI Upper`) == F) {
        CIuppervmax <- parlst[parlst$param.names == "lvmax", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "lkm", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "lkm", ]$`Back-transformed`) == F) {
        thetakm <- parlst[parlst$param.names == "lkm", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "lkm", ]$`%RSE`) == F) {
        rsekm <- parlst[parlst$param.names == "lkm", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "lkm", ]$`BSV(CV%)`) == F) {
        bsvkm <- parlst[parlst$param.names == "lkm", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "lkm", ]$`Shrink(SD)%`) == F) {
        shrinkkm <- parlst[parlst$param.names == "lkm", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "lkm", ]$`CI Lower`) == F) {
        CIlowerkm <- parlst[parlst$param.names == "lkm", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "lkm", ]$`CI Upper`) == F) {
        CIupperkm <- parlst[parlst$param.names == "lkm", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "lvc", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "lvc", ]$`Back-transformed`) == F) {
        thetavc <- parlst[parlst$param.names == "lvc", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "lvc", ]$`%RSE`) == F) {
        rsevc <- parlst[parlst$param.names == "lvc", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "lvc", ]$`BSV(CV%)`) == F) {
        bsvvc <- parlst[parlst$param.names == "lvc", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "lvc", ]$`Shrink(SD)%`) == F) {
        shrinkvc <- parlst[parlst$param.names == "lvc", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "lvc", ]$`CI Lower`) == F) {
        CIlowervc <- parlst[parlst$param.names == "lvc", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "lvc", ]$`CI Upper`) == F) {
        CIuppervc <- parlst[parlst$param.names == "lvc", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "lvp", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "lvp", ]$`Back-transformed`) == F) {
        thetavp <- parlst[parlst$param.names == "lvp", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "lvp", ]$`%RSE`) == F) {
        rsevp <- parlst[parlst$param.names == "lvp", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "lvp", ]$`BSV(CV%)`) == F) {
        bsvvp <- parlst[parlst$param.names == "lvp", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "lvp", ]$`Shrink(SD)%`) == F) {
        shrinkvp <- parlst[parlst$param.names == "lvp", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "lvp", ]$`CI Lower`) == F) {
        CIlowervp <- parlst[parlst$param.names == "lvp", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "lvp", ]$`CI Upper`) == F) {
        CIuppervp <- parlst[parlst$param.names == "lvp", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "lvp2", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "lvp2", ]$`Back-transformed`) ==
          F) {
        thetavp2 <- parlst[parlst$param.names == "lvp2", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "lvp2", ]$`%RSE`) == F) {
        rsevp2 <- parlst[parlst$param.names == "lvp2", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "lvp2", ]$`BSV(CV%)`) == F) {
        bsvvp2 <- parlst[parlst$param.names == "lvp2", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "lvp2", ]$`Shrink(SD)%`) == F) {
        shrinkvp2 <- parlst[parlst$param.names == "lvp2", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "lvp2", ]$`CI Lower`) == F) {
        CIlowervp2 <- parlst[parlst$param.names == "lvp2", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "lvp2", ]$`CI Upper`) == F) {
        CIuppervp2 <- parlst[parlst$param.names == "lvp2", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "lq", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "lq", ]$`Back-transformed`) == F) {
        thetaq <- parlst[parlst$param.names == "lq", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "lq", ]$`%RSE`) == F) {
        rseq <- parlst[parlst$param.names == "lq", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "lq", ]$`BSV(CV%)`) == F) {
        bsvq <- parlst[parlst$param.names == "lq", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "lq", ]$`Shrink(SD)%`) == F) {
        shrinkq <- parlst[parlst$param.names == "lq", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "lq", ]$`CI Lower`) == F) {
        CIlowerq <- parlst[parlst$param.names == "lq", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "lq", ]$`CI Upper`) == F) {
        CIupperq <- parlst[parlst$param.names == "lq", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "lq2", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "lq2", ]$`Back-transformed`) == F) {
        thetaq2 <- parlst[parlst$param.names == "lq2", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "lq2", ]$`%RSE`) == F) {
        rseq2 <- parlst[parlst$param.names == "lq2", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "lq2", ]$`BSV(CV%)`) == F) {
        bsvq2 <- parlst[parlst$param.names == "lq2", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "lq2", ]$`Shrink(SD)%`) == F) {
        shrinkq2 <- parlst[parlst$param.names == "lq2", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "lq2", ]$`CI Lower`) == F) {
        CIlowerq2 <- parlst[parlst$param.names == "lq2", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "lq2", ]$`CI Upper`) == F) {
        CIupperq2 <- parlst[parlst$param.names == "lq2", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "ltlag", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "ltlag", ]$`Back-transformed`) ==
          F) {
        thetatlag <-
          parlst[parlst$param.names == "ltlag", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "ltlag", ]$`%RSE`) == F) {
        rsetlag <- parlst[parlst$param.names == "ltlag", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "ltlag", ]$`BSV(CV%)`) == F) {
        bsvtlag <- parlst[parlst$param.names == "ltlag", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "ltlag", ]$`Shrink(SD)%`) ==
          F) {
        shrinktlag <- parlst[parlst$param.names == "ltlag", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "ltlag", ]$`CI Lower`) == F) {
        CIlowertlag <- parlst[parlst$param.names == "ltlag", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "ltlag", ]$`CI Upper`) == F) {
        CIuppertlag <- parlst[parlst$param.names == "ltlag", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "lka", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "lka", ]$`Back-transformed`) == F) {
        thetaka <- parlst[parlst$param.names == "lka", ]$`Back-transformed`
      }
      if (is.null(parlst[parlst$param.names == "lka", ]$`%RSE`) == F) {
        rseka <- parlst[parlst$param.names == "lka", ]$`%RSE`
      }
      if (is.null(parlst[parlst$param.names == "lka", ]$`BSV(CV%)`) == F) {
        bsvka <- parlst[parlst$param.names == "lka", ]$`BSV(CV%)`
      }

      if (is.null(parlst[parlst$param.names == "lka", ]$`Shrink(SD)%`) == F) {
        shrinkka <- parlst[parlst$param.names == "lka", ]$`Shrink(SD)%`
      }
      if (is.null(parlst[parlst$param.names == "lka", ]$`CI Lower`) == F) {
        CIlowerka <- parlst[parlst$param.names == "lka", ]$`CI Lower`
      }

      if (is.null(parlst[parlst$param.names == "lka", ]$`CI Upper`) == F) {
        CIupperka <- parlst[parlst$param.names == "lka", ]$`CI Upper`
      }
    }

    if (nrow(parlst[parlst$param.names == "add.err", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "add.err", ]$`Back-transformed`) ==
          F) {
        add <- parlst[parlst$param.names == "add.err", ]$`Back-transformed`
      }
    }

    if (nrow(parlst[parlst$param.names == "prop.err", ]) > 0) {
      if (is.null(parlst[parlst$param.names == "prop.err", ]$`Back-transformed`) ==
          F) {
        prop <- parlst[parlst$param.names == "prop.err", ]$`Back-transformed`
      }
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.ka", ]) > 0) {
      omegaka <- omegalst[omegalst$omega.names == "eta.ka", ]$eta.ka
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.cl", ]) > 0) {
      omegacl <- omegalst[omegalst$omega.names == "eta.cl", ]$eta.cl
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.vc", ]) > 0) {
      omegavc <- omegalst[omegalst$omega.names == "eta.vc", ]$eta.vc
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.vp", ]) > 0) {
      omegavp <- omegalst[omegalst$omega.names == "eta.vp", ]$eta.vp
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.vp2", ]) > 0) {
      omegavp2 <- omegalst[omegalst$omega.names == "eta.vp2", ]$eta.vp2
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.q", ]) > 0) {
      omegaq <- omegalst[omegalst$omega.names == "eta.q", ]$eta.q
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.q2", ]) > 0) {
      omegaq2 <- omegalst[omegalst$omega.names == "eta.q2", ]$eta.q2
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.tlag", ]) > 0) {
      omegatlag <- omegalst[omegalst$omega.names == "eta.tlag", ]$eta.tlag
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.vmax", ]) > 0) {
      omegavmax <- omegalst[omegalst$omega.names == "eta.vmax", ]$eta.vmax
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.km", ]) > 0) {
      omegavmax <- omegalst[omegalst$omega.names == "eta.km", ]$eta.km
    }

    if (nrow(omegalst[omegalst$omega.names == "eta.tlag", ]) > 0) {
      omegatlag <- omegalst[omegalst$omega.names == "eta.tlag", ]$eta.tlag
    }


    if (nrow(omegalst[omegalst$omega.names %in% c("eta.cl", "eta.vc"), ]) ==
        2) {
      omega.vc.cl <-
        omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.vc
      cor.eta.vc.cl <-
        omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.vc / (sqrt(omegacl) *
                                                                      sqrt(omegavc))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.cl", "eta.vp"), ]) ==
        2) {
      omega.vp.cl <-
        omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.vp
      cor.eta.vp.cl <-
        omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.vp / (sqrt(omegacl) *
                                                                      sqrt(omegavp))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.cl", "eta.q"), ]) ==
        2) {
      omega.q.cl <- omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.q
      cor.eta.q.cl <-
        omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.q / (sqrt(omegacl) *
                                                                     sqrt(omegaq))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.cl", "eta.vp2"), ]) ==
        2) {
      omega.vp2.cl <-
        omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.vp2
      cor.eta.vp2.cl <-
        omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.vp2 / (sqrt(omegacl) *
                                                                       sqrt(omegavp2))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.cl", "eta.q2"), ]) ==
        2) {
      omega.q2.cl <-
        omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.q2
      cor.eta.q2.cl <-
        omegalst[omegalst$omega.names %in% c("eta.cl"), ]$eta.q2 / (sqrt(omegacl) *
                                                                      sqrt(omegaq2))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.vp", "eta.vc"), ]) ==
        2) {
      omega.vp.vc <-
        omegalst[omegalst$omega.names %in% c("eta.vc"), ]$eta.vp
      cor.eta.vp.vc <-
        omegalst[omegalst$omega.names %in% c("eta.vc"), ]$eta.vp / (sqrt(omegavc) *
                                                                      sqrt(omegavp))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.q", "eta.vc"), ]) ==
        2) {
      omega.q.vc <- omegalst[omegalst$omega.names %in% c("eta.vc"), ]$eta.q
      cor.eta.q.vc <-
        omegalst[omegalst$omega.names %in% c("eta.vc"), ]$eta.q / (sqrt(omegavc) *
                                                                     sqrt(omegaq))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.vp2", "eta.vc"), ]) ==
        2) {
      omega.vp2.vc <-
        omegalst[omegalst$omega.names %in% c("eta.vc"), ]$eta.vp2
      cor.eta.vp2.vc <-
        omegalst[omegalst$omega.names %in% c("eta.vc"), ]$eta.vp2 / (sqrt(omegavc) *
                                                                       sqrt(omegavp2))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.q2", "eta.vc"), ]) ==
        2) {
      omega.q2.vc <-
        omegalst[omegalst$omega.names %in% c("eta.vc"), ]$eta.q2
      cor.eta.q2.vc <-
        omegalst[omegalst$omega.names %in% c("eta.vc"), ]$eta.q2 / (sqrt(omegavc) *
                                                                      sqrt(omegaq2))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.q", "eta.vp"), ]) ==
        2) {
      omega.q.vp <- omegalst[omegalst$omega.names %in% c("eta.vp"), ]$eta.q
      cor.eta.q.vp <-
        omegalst[omegalst$omega.names %in% c("eta.vp"), ]$eta.q / (sqrt(omegavp) *
                                                                     sqrt(omegaq))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.vp2", "eta.vp"), ]) ==
        2) {
      omega.vp2.vp <-
        omegalst[omegalst$omega.names %in% c("eta.vp"), ]$eta.vp2
      cor.eta.vp2.vp <-
        omegalst[omegalst$omega.names %in% c("eta.vp"), ]$eta.vp2 / (sqrt(omegavp) *
                                                                       sqrt(omegavp2))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.q2", "eta.vp"), ]) ==
        2) {
      omega.q2.vp <-
        omegalst[omegalst$omega.names %in% c("eta.vp"), ]$eta.q2
      cor.eta.q2.vp <-
        omegalst[omegalst$omega.names %in% c("eta.vp"), ]$eta.q2 / (sqrt(omegavp) *
                                                                      sqrt(omegaq2))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.vp2", "eta.q"), ]) ==
        2) {
      omega.vp2.q <-
        omegalst[omegalst$omega.names %in% c("eta.q"), ]$eta.vp2
      cor.eta.vp2.q <-
        omegalst[omegalst$omega.names %in% c("eta.q"), ]$eta.vp2 / (sqrt(omegaq) *
                                                                      sqrt(omegavp2))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.q2", "eta.q"), ]) ==
        2) {
      omega.q2.q <- omegalst[omegalst$omega.names %in% c("eta.q"), ]$eta.q2
      cor.eta.q2.q <-
        omegalst[omegalst$omega.names %in% c("eta.q"), ]$eta.q2 / (sqrt(omegaq) *
                                                                     sqrt(omegaq2))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.q2", "eta.vp2"), ]) ==
        2) {
      omega.q2.vp2 <-
        omegalst[omegalst$omega.names %in% c("eta.vp2"), ]$eta.q2
      cor.eta.q2.vp2 <-
        omegalst[omegalst$omega.names %in% c("eta.vp2"), ]$eta.q2 / (sqrt(omegavp2) *
                                                                       sqrt(omegaq2))
    }

    if (nrow(omegalst[omegalst$omega.names %in% c("eta.vmax", "eta.km"), ]) ==
        2) {
      omega.vmax.km <-
        omegalst[omegalst$omega.names %in% c("eta.vmax"), ]$eta.km
      cor.eta.vmax.km <-
        omegalst[omegalst$omega.names %in% c("eta.vmax"), ]$eta.km / (sqrt(omegavmax) *
                                                                        sqrt(omegakm))
    }
  }

  Store. <- data.frame(
    model.num = modi,
    current.time = current.time,
    AIC = AIC,
    BIC = BIC,
    OBJFV = OBJFV,
    ll = ll,
    npar = npar,
    model.covMethod = model.covMethod,
    model.time.setup = model.time.setup,
    model.time.covariance = model.time.covariance,
    model.time.saem = model.time.saem,
    model.time.table = model.time.table,
    model.time.compress = model.time.compress,
    model.time.other = model.time.other,
    thetaka = thetaka,
    thetacl = thetacl,
    thetavc = thetavc,
    thetavp = thetavp,
    thetavp2 = thetavp2,
    thetaq = thetaq,
    thetaq2 = thetaq2,
    thetatlag = thetatlag,
    thetavmax = thetavmax,
    thetakm = thetakm,
    rseka = rseka,
    rsecl = rsecl,
    rsevc = rsevc,
    rsevp = rsevp,
    rsevp2 = rsevp2,
    rseq = rseq,
    rseq2 = rseq2,
    rsetlag = rsetlag,
    rsevmax = rsevmax,
    rsekm = rsekm,
    bsvka = bsvka,
    bsvcl = bsvcl,
    bsvvc = bsvvc,
    bsvvp = bsvvp,
    bsvvp2 = bsvvp2,
    bsvq = bsvq,
    bsvq2 = bsvq2,
    bsvtlag = bsvtlag,
    bsvvmax = bsvvmax,
    bsvkm = bsvkm,
    shrinkka = shrinkka,
    shrinkcl = shrinkcl,
    shrinkvc = shrinkvc,
    shrinkvp = shrinkvp,
    shrinkvp2 = shrinkvp2,
    shrinkq = shrinkq,
    shrinkq2 = shrinkq2,
    shrinktlag = shrinktlag,
    shrinkvmax = shrinkvmax,
    shrinkkm = shrinkkm,
    CIlowerka = CIlowerka,
    CIlowercl = CIlowercl,
    CIlowervc = CIlowervc,
    CIlowervp = CIlowervp,
    CIlowervp2 = CIlowervp2,
    CIlowerq = CIlowerq,
    CIlowerq2 = CIlowerq2,
    CIlowertlag = CIlowertlag,
    CIlowervmax = CIlowervmax,
    CIlowerkm = CIlowerkm,
    CIupperka = CIupperka,
    CIuppercl = CIuppercl,
    CIuppervc = CIuppervc,
    CIuppervp = CIuppervp,
    CIuppervp2 = CIuppervp2,
    CIupperq = CIupperq,
    CIupperq2 = CIupperq2,
    CIuppertlag = CIuppertlag,
    CIuppervmax = CIuppervmax,
    CIupperkm = CIupperkm,
    omegaka = omegaka,
    omegacl = omegacl,
    omegavc = omegavc,
    omegavp = omegavp,
    omegavp2 = omegavp2,
    omegaq = omegaq,
    omegaq2 = omegaq2,
    omegatlag = omegatlag,
    omegavmax = omegavmax,
    omegakm = omegakm,
    omega.vc.cl = omega.vc.cl,
    omega.vp.cl = omega.vp.cl,
    omega.q.cl = omega.q.cl,
    omega.vp2.cl = omega.vp2.cl,
    omega.q2.cl = omega.q2.cl,
    omega.vp.vc = omega.vp.vc,
    omega.q.vc = omega.q.vc,
    omega.vp2.vc = omega.vp2.vc,
    omega.q2.vc = omega.q2.vc,
    omega.q.vp = omega.q.vp,
    omega.vp2.vp = omega.vp2.vp,
    omega.q2.vp = omega.q2.vp,
    omega.vp2.q = omega.vp2.q,
    omega.q2.q = omega.q2.q,
    omega.q2.vp2 = omega.q2.vp2,
    omega.vmax.km = omega.vmax.km,
    cor.eta.vc.cl = cor.eta.vc.cl,
    cor.eta.vp.cl = cor.eta.vp.cl,
    cor.eta.q.cl = cor.eta.q.cl,
    cor.eta.vp2.cl = cor.eta.vp2.cl,
    cor.eta.q2.cl = cor.eta.q2.cl,
    cor.eta.vp.vc = cor.eta.vp.vc,
    cor.eta.q.vc = cor.eta.q.vc,
    cor.eta.vp2.vc = cor.eta.vp2.vc,
    cor.eta.q2.vc = cor.eta.q2.vc,
    cor.eta.q.vp = cor.eta.q.vp,
    cor.eta.vp2.vp = cor.eta.vp2.vp,
    cor.eta.q2.vp = cor.eta.q2.vp,
    cor.eta.vp2.q = cor.eta.vp2.q,
    cor.eta.q2.q = cor.eta.q2.q,
    cor.eta.q2.vp2 = cor.eta.q2.vp2,
    cor.eta.vmax.km = cor.eta.vmax.km,
    add = add,
    prop = prop
  )
  return(Store.)

}
