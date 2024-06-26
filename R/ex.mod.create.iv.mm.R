#' Generate model files for population pharmacokinetic modelling
#'
#' Generate a model file based on the specified parameters for use in automated modelling,
#' and the model's initial estimates, residual errors, and differential equations based on the
#' number of compartments and the presence of random effects, correlations, and Michaelis-Menten elimination.
#'
#' @param modi The current model number.
#' @param cmpt.iv The number of compartments, represented by a categorical variable:
#' \itemize{
#'   \item \code{cmpt.iv = 1}: One compartment model
#'   \item \code{cmpt.iv = 2}: Two compartment model
#'   \item \code{cmpt.iv = 3}: Three compartment model
#' }
#' @param eta.km A binary value indicating the presence (1) or absence (0) of random effects for the Michaelis-Menten constant (Km).
#' @param eta.vc A binary value indicating the presence (1) or absence (0) of random effects for the central volume of distribution (Vc).
#' @param eta.vp A binary value indicating the presence (1) or absence (0) of random effects for the peripheral volume of distribution (Vp).
#' @param eta.vp2 A binary value indicating the presence (1) or absence (0) of random effects for the second peripheral volume of distribution (Vp2).
#' @param eta.q A binary value indicating the presence (1) or absence (0) of random effects for the intercompartmental clearance (Q).
#' @param eta.q2 A binary value indicating the presence (1) or absence (0) of random effects for the second intercompartmental clearance (Q2).
#' @param mm A parameter representing the Michaelis-Menten elimination.
#' @param mcorr A parameter representing the existence of correlation between parameters.
#' @param rv A categorical variable representing the type of residual error:
#' \itemize{
#'   \item \code{rv = 1}: Additive error
#'   \item \code{rv = 2}: Proportional error
#'   \item \code{rv = 3}: Combined error
#' }
#'
#' @return None. The function creates a model file.
#' @import stringr
#'
#' @examples
#' \dontrun{
#' ex.mod.create.iv.mm(modi = 1, cmpt.iv = 2, eta.km = 1, eta.vc = 1, eta.vp = 1, eta.vp2 = 0,
#'                     eta.q = 1, eta.q2 = 0, mm = 1, mcorr = 1, rv = 3)
#' }
#'
#' @export

ex.mod.create.iv.mm <- function(modi,
                                cmpt.iv,
                                eta.km,
                                eta.vc,
                                eta.vp,
                                eta.vp2,
                                eta.q,
                                eta.q2,
                                mm,
                                mcorr,
                                rv)
{
  netas <-
    eta.vc + eta.vp + eta.vp2 + eta.q + eta.q2 # netas. etas on parameters
  # Note: IIV on cl is always considered
  clini = paste0('lcl  = round(log(', input_cl, '),2)')
  modcl = c("cl = exp(lcl+eta.cl)")
  etaclini = paste0('eta.cl ~ ', input_eta.cl)

  vmaxini = ''
  kmini = ''
  etavmaxini = ''
  etakmini = ''
  modvmax = ''
  modkm = ''
  vcini = paste0('lvc = round(log(', input_vc, '),2)')
  modvc = c("vc = exp(lvc")
  etavcini = ''

  vpini       = ''
  vp2ini      = ''
  qini        = ''
  q2ini       = ''

  # initial estimate for eta
  etavpini    = ''
  etavp2ini   = ''
  etaqini     = ''
  etaq2ini    = ''

  # other parameter mod initial setting
  modvp     = ''
  modvp2    = ''
  modq      = ''
  modq2     = ''
  correta   = ''     # for v,q
  correta2  = ''     # for vmax,km


  ##############################Base setting for one-compartment model#############

  if (cmpt.iv == 1) {
    odeline2 <- "d/dt(centr) <- - cl / vc * centr"
    odeline3 <- ""
    odeline4 <- ""
    odelinef <- "cp <- centr / vc"

    if (eta.vc == 1) {
      modvc = paste0(modvc, '+eta.vc')
      etavcini = paste0('eta.vc ~ ', input_eta.vc)
    }

    modvc = paste0(modvc, ')')

    #########Correlation############################################################
    if (mcorr == 1 & netas > 0 & mm == 0) {
      etaclini = ""
      etavcini = ""

      etacov = sqrt(input_eta.cl) * sqrt(input_eta.vc) * input_etarsquare # R = 0.1 weak correlation assumed
      correta <-
        paste0("eta.cl + eta.vc ~ c(",
               input_eta.cl,
               ",",
               etacov,
               ",",
               input_eta.vc,
               ")")
    }

  }

  #######################Base setting for two-compartment model####################

  # correlation parameter candidate
  corr.candidate <- data.frame(corr.param = "cl")

  if (cmpt.iv == 2) {
    odeline2 <-
      "d/dt(centr) <- - cl / vc * centr - q/vc* centr + q/vp *peri"
    odeline3 <- "d/dt(peri) <- q/vc* centr - q/vp *peri"
    odeline4 <- ""
    odelinef <- "cp <- centr / vc"

    # Ratio of vc to vp was set 0.5
    vcini  = paste0('lvc = round(log(', input_vc2cpt, '),1)')
    vpini = paste0('lvp = round(log(', input_vp2cpt, '),1)')
    qini  = paste0('lq = round(log(', input_q2cpt, '),2)')

    modvp = c("vp = exp(lvp")
    modq = c("q  = exp(lq")

    # Check which parameters were introduced IIV
    eta.param.all <- data.frame(eta.vc, eta.vp, eta.q)

    eta.ans1 <- str_detect(eta.param.all$eta.vc, "1")
    eta.ans2 <- str_detect(eta.param.all$eta.vp, "1")
    eta.ans3 <- str_detect(eta.param.all$eta.q, "1")


    if (eta.ans1 == T) {
      modvc = paste0(modvc, '+eta.vc')
      corr.candidate.s <- data.frame(corr.param = "eta.vc")
      corr.candidate <- rbind(corr.candidate, corr.candidate.s)
      etavcini = paste0('eta.vc ~ ', input_eta.vc)
    }

    if (eta.ans2 == T) {
      modvp = paste0(modvp, '+eta.vp')
      corr.candidate.s <- data.frame(corr.param = "eta.vp")
      corr.candidate <- rbind(corr.candidate, corr.candidate.s)
      etavpini = paste0('eta.vp ~ ', input_eta.vp)
    }

    if (eta.ans3 == T) {
      modq = paste0(modq, '+eta.q')
      corr.candidate.s <- data.frame(corr.param = "eta.q")
      corr.candidate <- rbind(corr.candidate, corr.candidate.s)
      etaqini = paste0('eta.q ~ ', input_eta.q)
    }

    modvc = paste0(modvc, ')')
    modvp = paste0(modvp, ')')
    modq = paste0(modq, ')')

    #########Correlation###########################################################
    # cl exists.
    if (mcorr == 1 & netas > 0 & mm == 0) {
      etaclini = ""
      etavcini = ""
      etavpini = ""
      etaqini = ""

      etacov = sqrt(input_eta) * sqrt(input_eta) * input_etarsquare  # weak correlation assumed
      if (nrow(corr.candidate) == 2) {
        varname <- corr.candidate$corr.param[2]
        correta <-
          paste0("eta.cl",
                 "+",
                 varname,
                 " ~ c(",
                 input_eta.cl,
                 ",",
                 etacov,
                 ",",
                 input_eta,
                 ")")
      }

      if (nrow(corr.candidate) == 3) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]

        correta <- paste0(
          "eta.cl",
          "+",
          varname1,
          "+",
          varname2,
          " ~ c(",
          input_eta,
          ",",
          etacov,
          ",",
          input_eta,
          ",",
          etacov,
          ",",
          etacov,
          ",",
          input_eta,
          ")"
        )
      }
      if (nrow(corr.candidate) == 4) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]
        varname3 <- corr.candidate$corr.param[4]

        correta <-
          paste0(
            "eta.cl",
            " + ",
            varname1,
            " + ",
            varname2,
            " + ",
            varname3,
            " ~ c(",
            input_eta,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ")"
          )
      }
    }
    # cl does not exists.
    if (mcorr == 1 & netas > 1 & mm == 1) {
      etaclini = ""
      etavcini = ""
      etavpini = ""
      etaqini = ""

      etacov = sqrt(input_eta) * sqrt(input_eta) * input_etarsquare # weak correlation assumed

      if (nrow(corr.candidate) == 2) {
        varname <- corr.candidate$corr.param[2]
        correta <-
          paste0("eta.cl",
                 "+",
                 varname,
                 " ~ c(",
                 input_eta.cl,
                 ",",
                 etacov,
                 ",",
                 input_eta,
                 ")")
      }

      if (nrow(corr.candidate) == 3) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]

        correta <-
          paste0(varname1,
                 "+",
                 varname2,
                 " ~ c(",
                 input_eta.cl,
                 ",",
                 etacov,
                 ",",
                 input_eta,
                 ")")
      }
      if (nrow(corr.candidate) == 4) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]
        varname3 <- corr.candidate$corr.param[4]

        correta <- paste0(
          varname1,
          " + ",
          varname2,
          " + ",
          varname3,
          " ~ c(",
          input_eta,
          ",",
          etacov,
          ",",
          input_eta,
          ",",
          etacov,
          ",",
          etacov,
          ",",
          input_eta,
          ")"
        )
      }
    }
  }

  ################Base setting for three-compartment model#######################
  if (cmpt.iv == 3) {
    odeline2 <-
      "d/dt(centr) <- - cl / vc * centr - q/vc* centr + q/vp *peri  - q2/vc* centr + q2/vp2 *peri2 "
    odeline3 <- "d/dt(peri) <- q/vc* centr - q/vp *peri"
    odeline4 <- "d/dt(peri2) <- q2/vc* centr - q2/vp2 *peri2"
    odelinef <- "cp <- centr / vc"

    vcini  = paste0('lvc = round(log(', input_vc3cpt, '),1)')
    vpini =  paste0('lvp = round(log(', input_vp3cpt, '),1)')
    qini  =  paste0('lq = round(log(', input_q3cpt, '),2)')
    vp2ini = paste0('lvp2 = round(log(', input_vp23cpt, '),1)')
    q2ini  = paste0('lq2 = round(log(', input_q23cpt, '),2)')

    modvp = c("vp = exp(lvp")
    modq = c("q  = exp(lq")
    modvp2 = c("vp2 = exp(lvp2")
    modq2 = c("q2  = exp(lq2")

    # Check which parameters were introduced IIV

    eta.param.all <- data.frame(eta.vc, eta.vp, eta.q, eta.vp2, eta.q2)

    eta.ans1 <- str_detect(eta.param.all$eta.vc, "1")
    eta.ans2 <- str_detect(eta.param.all$eta.vp, "1")
    eta.ans3 <- str_detect(eta.param.all$eta.q, "1")
    eta.ans4 <- str_detect(eta.param.all$eta.vp2, "1")
    eta.ans5 <- str_detect(eta.param.all$eta.q2, "1")

    if (eta.ans1 == T) {
      modvc = paste0(modvc, '+eta.vc')
      corr.candidate.s <- data.frame(corr.param = "eta.vc")
      corr.candidate <- rbind(corr.candidate, corr.candidate.s)
      etavcini = paste0('eta.vc ~ ', input_eta.vc)
    }

    if (eta.ans2 == T) {
      modvp = paste0(modvp, '+eta.vp')
      corr.candidate.s <- data.frame(corr.param = "eta.vp")
      corr.candidate <- rbind(corr.candidate, corr.candidate.s)
      etavpini = paste0('eta.vp ~ ', input_eta.vp)
    }

    if (eta.ans3 == T) {
      modq = paste0(modq, '+eta.q')
      corr.candidate.s <- data.frame(corr.param = "eta.q")
      corr.candidate <- rbind(corr.candidate, corr.candidate.s)
      etaqini = paste0('eta.q ~ ', input_eta.q)
    }

    if (eta.ans4 == T) {
      modvp2 = paste0(modvp2, '+eta.vp2')
      corr.candidate.s <- data.frame(corr.param = "eta.vp2")
      corr.candidate <- rbind(corr.candidate, corr.candidate.s)
      etavp2ini = paste0('eta.vp2 ~ ', input_eta.vp2)
    }

    if (eta.ans5 == T) {
      modq2 = paste0(modq2, '+eta.q2')
      corr.candidate.s <- data.frame(corr.param = "eta.q2")
      corr.candidate <- rbind(corr.candidate, corr.candidate.s)
      etaq2ini = paste0('eta.q2 ~ ', input_eta.q2)
    }

    modvc = paste0(modvc, ')')
    modvp = paste0(modvp, ')')
    modq = paste0(modq, ')')
    modvp2 = paste0(modvp2, ')')
    modq2 = paste0(modq2, ')')

    #########Correlation###########################################################
    # cl exists.
    if (mcorr == 1 & netas > 0 & mm == 0) {
      etaclini = ""
      etavcini = ""
      etavpini = ""
      etaqini = ""
      etavp2ini = ""
      etaq2ini = ""

      etacov = sqrt(input_eta) * sqrt(input_eta) * input_etarsquare # weak correlation assumed

      if (nrow(corr.candidate) == 2) {
        varname <- corr.candidate$corr.param[2]
        correta <-
          paste0("eta.cl",
                 "+",
                 varname,
                 "~ c(",
                 input_eta.cl,
                 ",",
                 etacov,
                 ",",
                 input_eta,
                 ")")
      }
      if (nrow(corr.candidate) == 3) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]


        correta <- paste0(
          "eta.cl",
          "+",
          varname1,
          "+",
          varname2,
          " ~ c(",
          input_eta,
          ",",
          etacov,
          ",",
          input_eta,
          ",",
          etacov,
          ",",
          etacov,
          ",",
          input_eta,
          ")"
        )

      }
      if (nrow(corr.candidate) == 4) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]
        varname3 <- corr.candidate$corr.param[4]

        correta <-
          paste0(
            "eta.cl",
            " +",
            varname1,
            " +",
            varname2,
            " +",
            varname3,
            " ~ c(",
            input_eta,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta ,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ")"
          )
      }
      if (nrow(corr.candidate) == 5) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]
        varname3 <- corr.candidate$corr.param[4]
        varname4 <- corr.candidate$corr.param[5]

        correta <-
          paste0(
            "eta.cl",
            " +",
            varname1,
            " +",
            varname2,
            " +",
            varname3,
            " +",
            varname4,
            " ~ c(",
            input_eta,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ")"
          )
      }
      if (nrow(corr.candidate) == 6) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]
        varname3 <- corr.candidate$corr.param[4]
        varname4 <- corr.candidate$corr.param[5]
        varname5 <- corr.candidate$corr.param[6]

        correta <-
          paste0(
            "eta.cl",
            " +",
            varname1,
            " +",
            varname2,
            " +",
            varname3,
            " +",
            varname4,
            " +",
            varname5,
            " ~ c(",
            input_eta,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ")"
          )
      }

    }
    if (mcorr == 1 & netas > 1 & mm == 1) {
      etaclini = ""
      etavcini = ""
      etavpini = ""
      etaqini = ""
      etavp2ini = ""
      etaq2ini = ""
      etacov = sqrt(input_eta) * sqrt(input_eta) * input_etarsquare # 0.1-0.2 weak correlation assumed


      if (nrow(corr.candidate) == 3) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]
        correta <- paste0(varname1,
                          "+",
                          varname2,
                          "~ c(",
                          input_eta.cl,
                          ",",
                          etacov,
                          ",",
                          input_eta,
                          ")")

      }
      if (nrow(corr.candidate) == 4) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]
        varname3 <- corr.candidate$corr.param[4]
        correta <- paste0(
          varname1,
          " +",
          varname2,
          " +",
          varname3,
          " ~ c(",
          input_eta,
          ",",
          etacov,
          ",",
          input_eta,
          ",",
          etacov,
          ",",
          etacov,
          ",",
          input_eta,
          ")"
        )
      }
      if (nrow(corr.candidate) == 5) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]
        varname3 <- corr.candidate$corr.param[4]
        varname4 <- corr.candidate$corr.param[5]
        correta <-
          paste0(
            varname1,
            " +",
            varname2,
            " +",
            varname3,
            " +",
            varname4,
            " ~ c(",
            input_eta,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta ,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ")"
          )
      }
      if (nrow(corr.candidate) == 6) {
        varname1 <- corr.candidate$corr.param[2]
        varname2 <- corr.candidate$corr.param[3]
        varname3 <- corr.candidate$corr.param[4]
        varname4 <- corr.candidate$corr.param[5]
        varname5 <- corr.candidate$corr.param[6]
        correta <-
          paste0(
            varname1,
            " +",
            varname2,
            " +",
            varname3,
            " +",
            varname4,
            " +",
            varname5,
            " ~ c(",
            input_eta,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            etacov,
            ",",
            input_eta,
            ")"
          )
      }

    }
  }

  ########################### Statistical Model part#############################
  modresline = ""
  iniresline1 = ''
  iniresline2 = ''

  if (rv == 1) {
    iniresline2 = paste0('add.err  <- ', input_add.err)
    modresline = paste0(modresline, "cp~add(add.err)")
  }

  if (rv == 2) {
    iniresline1 = paste0('prop.err  <- ', input_prop.err)
    modresline = paste0(modresline, "cp~prop(prop.err)")
  }

  if (rv == 3) {
    iniresline1 = paste0('prop.err  <- ', input_prop.err)
    modresline = paste0(modresline, "cp~prop(prop.err)")

    iniresline2 = paste0('add.err  <- ', input_add.err)
    modresline = paste0(modresline, "+add(add.err)")
  }

  if (mm == 1) {
    clini = ""
    etaclini = ""
    modcl = ""
    vmaxini = paste0('lvmax  = round(log(', input_vmax, '),2)')
    # Note: IIV on vmax is always considered
    modvmax = c("vmax = exp(lvmax+eta.vmax)")
    etavmaxini = paste0('eta.vmax ~ ', input_eta.vmax)

    kmini = paste0('lkm  = round(log(', input_km, '),2)')
    modkm = c("km = exp(lkm")

    if (eta.km == 1) {
      modkm = paste0(modkm, '+eta.km')
      etakmini = paste0('eta.km ~ ', input_eta.km)
    }
    modkm = paste0(modkm, ')')

    if (cmpt.iv == 1) {
      odeline2 <- "d/dt(centr) <- - (vmax/(km+ centr/vc)) / vc * centr"
    }
    if (cmpt.iv == 2) {
      odeline2 <-
        "d/dt(centr) <- -  (vmax/(km+ centr/vc)) / vc * centr - q/vc* centr + q/vp *peri"
    }
    if (cmpt.iv == 3) {
      odeline2 <-
        "d/dt(centr) <- -  (vmax/(km+ centr/vc)) / vc * centr - q/vc* centr + q/vp *peri - q2/vc* centr + q2/vp2 *peri2"
    }

    if (mcorr == 1 & eta.km == 1) {
      etavmaxini <- ""
      etakmini <- ""
      r = sqrt(input_eta) * sqrt(input_eta) * input_etarsquare # weak correlation assumed
      correta2 <-
        paste0("eta.vmax + eta.km ~ c(",
               input_eta.vmax,
               ",",
               r,
               ",",
               input_eta.km,
               ")")
    }
  }
  ##########################Write model file part################################

  fileConn <- file(paste0("mod", modi, '.txt'))

  writeLines(
    c(
      "f <- function(){",
      "ini({",

      clini,
      #     lcl      = log(input)
      vmaxini,
      kmini,
      vcini,
      #     lvc      = log(input)
      vpini,
      #     lvp      = log(input)
      vp2ini,
      #     lvp2     = log(input)
      qini,
      #     lq       = log(input)
      q2ini,
      #     lq2      = log(input)

      iniresline1,
      #    "add.err <- c(0, 0.1, 1)",
      iniresline2,
      #    "prop.err <- c(0, 0.1, 1)",


      etaclini,
      #     "eta.cl      ~ 0.01",
      etavmaxini,
      #     "eta.vmax    ~ 0.01",
      etakmini,
      #     "eta.km      ~ 0.01",
      etavcini,
      #     "eta.vc      ~ 0.01",
      etavpini,
      #     "eta.vp      ~ 0.01",
      etavp2ini,
      #     "eta.vp2     ~ 0.01",
      etaqini,
      #     "eta.q       ~ 0.01",
      etaq2ini,
      #     "eta.q2      ~ 0.01",

      correta,
      #  eta.cl  + eta.vc  ~c (0.01,-0.04,0.01)
      correta2,
      # eta.vmax + eta.km ~c (0.01,-0.04,0.01)
      "})",

      "model({",

      modcl,
      #    cl    =   exp(lcl+eta.cl)",
      modvmax,
      modkm,
      modvc,
      #    vc    =   exp(lv2+eta.vc)",
      modvp,
      #    vp    =   exp(lv3+eta.vp)"
      modvp2,
      #    vp2   =   exp(lq+eta.vp2)"
      modq,
      #    q     =   exp(lq+eta.q)"
      modq2,
      #    q2    =   exp(lq+eta.q2)"


      odeline2,
      odeline3,
      odeline4,
      odelinef,

      modresline,
      "})",
      "}"



    ),

    fileConn
  )
  close(fileConn)

}
