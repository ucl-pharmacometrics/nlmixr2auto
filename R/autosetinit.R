#' Automatically set initial values for model parameters
#'
#' Sets initial values for various model parameters based on provided data and initial estimates.
#' It updates global variables with the specified or default initial values.
#'
#' @param dat A data frame containing the data to be used for setting initial values.
#' @param inits A data frame containing the initial estimates of the parameters, which is obtained from `nlmixr2autoinit`
#' @param input_vmax Optional. Initial value for Vmax. If not provided, it will be set based on `inits`.
#' @param input_km Optional. Initial value for Km. If not provided, it will be set based on `inits`.
#' @param input_cl Optional. Initial value for clearance (CL). If not provided, it will be set based on `inits`.
#' @param input_vc Optional. Initial value for volume of distribution (Vc). If not provided, it will be set based on `inits`.
#' @param input_vc2cpt Optional. Initial value for Vc in a two-compartment model. If not provided, it will be set based on `inits`.
#' @param input_vp2cpt Optional. Initial value for Vp in a two-compartment model. If not provided, it will be set based on `inits`.
#' @param input_q2cpt Optional. Initial value for Q in a two-compartment model. If not provided, it will be set based on `inits`.
#' @param input_vc3cpt Optional. Initial value for Vc in a three-compartment model. If not provided, it will be set based on `inits`.
#' @param input_vp3cpt Optional. Initial value for Vp in a three-compartment model. If not provided, it will be set based on `inits`.
#' @param input_vp23cpt Optional. Initial value for Vp2 in a three-compartment model. If not provided, it will be set based on `inits`.
#' @param input_q3cpt Optional. Initial value for Q in a three-compartment model. If not provided, it will be set based on `inits`.
#' @param input_q23cpt Optional. Initial value for Q2 in a three-compartment model. If not provided, it will be set based on `inits`.
#' @param input_prop.err Optional. Initial value for proportional error. Default is 0.1.
#' @param input_add.err Optional. Initial value for additive error. If not provided, it will be set based on observed data.
#' @param input_add.err Optional. Initial value for log-normal distributed additive error. If not provided, it will be set based on observed data.
#' @param input_eta.cl Optional. Initial value for ETA on clearance. Default is 0.1.
#' @param input_eta.vc Optional. Initial value for ETA on Vc. Default is 0.1.
#' @param input_eta.vp Optional. Initial value for ETA on Vp. Default is 0.1.
#' @param input_eta.vp2 Optional. Initial value for ETA on Vp2. Default is 0.1.
#' @param input_eta.q Optional. Initial value for ETA on Q. Default is 0.1.
#' @param input_eta.q2 Optional. Initial value for ETA on Q2. Default is 0.1.
#' @param input_eta.vmax Optional. Initial value for ETA on Vmax. Default is 0.1.
#' @param input_eta.km Optional. Initial value for ETA on Km. Default is 0.1.
#' @param input_etarsquare Optional. Initial value for the square of the correlation coefficient for ETA. Default is 0.2.
#'
#' @return None. This function updates global variables.
#'
#' @examples
#' \dontrun{
#'
#' library(nlmixr2autoinit)
#' dat<-pheno_sd
#' inits.out<-getppkinit(dat = dat,runnpd = 0)
#' autosetinit(dat = dat,
#'            inits= inits.out$Recommended_initial_estimates)
#' }
#'
#' @export

autosetinit <-  function(dat,
                         inits,
                         input_vmax = NULL,
                         input_km = NULL,
                         input_cl = NULL,
                         input_vc = NULL,
                         input_vc2cpt = NULL,
                         input_vp2cpt = NULL,
                         input_q2cpt = NULL,
                         input_vc3cpt = NULL,
                         input_vp3cpt = NULL,
                         input_vp23cpt = NULL,
                         input_q3cpt = NULL,
                         input_q23cpt = NULL,
                         input_prop.err = NULL,
                         input_add.err = NULL,
                         input_logn.sd = NULL,
                         input_eta.cl = NULL,
                         input_eta.vc = NULL,
                         input_eta.vp = NULL,
                         input_eta.vp2 = NULL,
                         input_eta.q = NULL,
                         input_eta.q2 = NULL,
                         input_eta.vmax = NULL,
                         input_eta.km = NULL,
                         input_etarsquare = NULL) {
  # For allometric scalling
  # mean weight of population
  mean_wt <- 70 # default values
  if (!is.null(dat$WT)) {
    mean_wt <- round(mean(dat[!duplicated(dat$ID),]$WT, na.rm = T), 1)
  }

  if (!missing(input_vmax)) {
    input_vmax <- input_vmax
    input_vmax_70 <- input_vmax * (70 / mean_wt) ^ 0.75
  } else {
    input_vmax <- signif(inits[inits$Parameters == "Vmax", ]$Values, 3)
    input_vmax_70 <-
      signif(inits[inits$Parameters == "Vmax", ]$Values * (70 / mean_wt) ^ 0.75, 3)
  }

  if (!missing(input_km)) {
    input_km <- input_km
  } else {
    input_km <- signif (inits[inits$Parameters == "Km", ]$Values, 3)
  }

  if (!missing(input_cl)) {
    input_cl <- input_cl
    input_cl_70 <- input_vmax * (70 / mean_wt) ^ 0.75
  } else {
    input_cl <- signif(inits[inits$Parameters == "CL", ]$Values, 3)
    input_cl_70 <-
      signif(inits[inits$Parameters == "CL", ]$Values * (70 / mean_wt) ^ 0.75, 3)
  }

  if (!missing(input_vc)) {
    input_vc <- input_vc
    input_vc_70 <- input_vc * (70 / mean_wt)
  } else {
    input_vc <- signif(inits[inits$Parameters == "Vd", ]$Values, 3)
    input_vc_70 <-
      signif(inits[inits$Parameters == "Vd", ]$Values * (70 / mean_wt), 3)
  }

  if (!missing(input_vc2cpt)) {
    input_vc2cpt <- input_vc2cpt
    input_vc2cpt_70 <- input_vc2cpt * (70 / mean_wt)
  } else {
    input_vc2cpt <-
      signif(inits[inits$Parameters == "Vc(2CMPT)", ]$Values, 1)
    input_vc2cpt_70 <-
      signif(inits[inits$Parameters == "Vc(2CMPT)", ]$Values * (70 / mean_wt), 3)
  }

  if (!missing(input_vp2cpt)) {
    input_vp2cpt <- input_vp2cpt
    input_vp2cpt_70 <- input_vp2cpt * (70 / mean_wt)
  } else {
    input_vp2cpt <-
      signif(inits[inits$Parameters == "Vp(2CMPT)", ]$Values, 3)
    input_vp2cpt_70 <-
      signif(inits[inits$Parameters == "Vp(2CMPT)", ]$Values * (70 / mean_wt), 3)
  }

  if (!missing(input_q2cpt)) {
    input_q2cpt <- input_q2cpt
    input_q2cpt_70 <- input_q2cpt * (70 / mean_wt) ^ 0.75
  } else {
    input_q2cpt <- signif(inits[inits$Parameters == "CL", ]$Values, 3)
    input_q2cpt_70 <-
      signif(inits[inits$Parameters == "CL", ]$Values * (70 / mean_wt) ^ 0.75, 3)
  }

  if (!missing(input_vc3cpt)) {
    input_vc3cpt <- input_vc3cpt
    input_vc3cpt_70 <- input_vc3cpt * (70 / mean_wt)
  } else {
    input_vc3cpt <-
      signif(inits[inits$Parameters == "Vc(3CMPT)", ]$Values, 3)
    input_vc3cpt_70 <-
      signif(inits[inits$Parameters == "Vc(3CMPT)", ]$Values * (70 / mean_wt), 3)
  }

  if (!missing(input_vp3cpt)) {
    input_vp3cpt <- input_vp3cpt
    input_vp3cpt_70 <- input_vp3cpt * (70 / mean_wt)
  } else {
    input_vp3cpt <-
      signif(inits[inits$Parameters == "Vp(3CMPT)", ]$Values, 3)
    input_vp3cpt_70 <-
      signif(inits[inits$Parameters == "Vp(3CMPT)", ]$Values * (70 / mean_wt), 3)
  }

  if (!missing(input_vp23cpt)) {
    input_vp23cpt <- input_vp23cpt
    input_vp23cpt_70 <- input_vp3cpt * (70 / mean_wt)
  } else {
    input_vp23cpt <-
      signif(inits[inits$Parameters == "Vp2(3CMPT)", ]$Values, 3)
    input_vp23cpt_70 <-
      signif(inits[inits$Parameters == "Vp2(3CMPT)", ]$Values * (70 / mean_wt), 3)
  }

  if (!missing(input_q3cpt)) {
    input_q3cpt <- input_q3cpt
    input_q3cpt_70 <- input_q3cpt * (70 / mean_wt) ^ 0.75
  } else {
    input_q3cpt <-  signif(inits[inits$Parameters == "CL", ]$Values, 3)
    input_q3cpt_70 <-
      signif(inits[inits$Parameters == "CL", ]$Values * (70 / mean_wt) ^ 0.75, 3)
  }

  if (!missing(input_q23cpt)) {
    input_q23cpt <- input_q23cpt
    input_q23cpt_70 <- input_q23cpt * (70 / mean_wt) ^ 0.75
  } else {
    input_q23cpt <-  signif(inits[inits$Parameters == "CL", ]$Values, 3)
    input_q23cpt_70 <-
      signif(inits[inits$Parameters == "CL", ]$Values * (70 / mean_wt) ^ 0.75, 3)
  }

  if (!missing(input_prop.err)) {
    input_prop.err <- input_prop.err
  } else {
    input_prop.err <- 0.1
  }

  if (!missing(input_add.err)) {
    input_add.err <- input_add.err
    input_logn.sd <- input_logn.sd
  } else {
    dat.obs <- dat[dat$EVID == 0, ]
    pop.cmax <- aggregate(dat.obs, by = list(dat.obs$ID), FUN = max)
    dat.cmax <- median(pop.cmax$DV)
    input_add.err <- round(dat.cmax / 1000, 1) # cmax/1000

    input_logn.sd <- round(log(dat.cmax / 1000), 3)

    if (input_add.err < 0.001) {
      input_add.err <- 0.001
    }

    if (input_logn.sd < 0.001) {
      input_logn.sd <- 0.001
    }
  }

  if (!missing(input_eta.cl)) {
    input_eta.cl <- input_eta.cl
  } else {
    input_eta.cl <- 0.1
  }

  if (!missing(input_eta.vc)) {
    input_eta.vc <- input_eta.vc
  } else {
    input_eta.vc <- 0.1
  }

  if (!missing(input_eta.vp)) {
    input_eta.vp <- input_eta.vp
  } else {
    input_eta.vp <- 0.1
  }

  if (!missing(input_eta.vp2)) {
    input_eta.vp2 <- input_eta.vp2
  } else {
    input_eta.vp2 <- 0.1
  }

  if (!missing(input_eta.q)) {
    input_eta.q <- input_eta.q
  } else {
    input_eta.q <- 0.1
  }

  if (!missing(input_eta.q2)) {
    input_eta.q2 <- input_eta.q2
  } else {
    input_eta.q2 <- 0.1
  }

  if (!missing(input_eta.vmax)) {
    input_eta.vmax <- input_eta.vmax
  } else {
    input_eta.vmax <- 0.1
  }

  if (!missing(input_eta.km)) {
    input_eta.km <- input_eta.km
  } else {
    input_eta.km <- 0.1
  }


  if (!missing(input_etarsquare)) {
    input_etarsquare <- input_etarsquare
  } else {
    input_etarsquare <- 0.2
  }


  input_eta <<- 0.1

  # Function to check and replace 0 with 0.001
  replace_zero <- function(x) {
    if (x == 0) {
      return(0.001)
    } else {
      return(x)
    }
  }

  # Apply the function to each variable
  input_vmax <- replace_zero(input_vmax)
  input_km <- replace_zero(input_km)
  input_cl <- replace_zero(input_cl)
  input_vc <- replace_zero(input_vc)
  input_vc2cpt <- replace_zero(input_vc2cpt)
  input_vp2cpt <- replace_zero(input_vp2cpt)
  input_q2cpt <- replace_zero(input_q2cpt)
  input_vc3cpt <- replace_zero(input_vc3cpt)
  input_vp3cpt <- replace_zero(input_vp3cpt)
  input_vp23cpt <- replace_zero(input_vp23cpt)
  input_q3cpt <- replace_zero(input_q3cpt)
  input_q23cpt <- replace_zero(input_q23cpt)
  input_vmax_70 <- replace_zero(input_vmax_70)
  input_cl_70 <- replace_zero(input_cl_70)
  input_vc_70 <- replace_zero(input_vc_70)
  input_vc2cpt_70 <- replace_zero(input_vc2cpt_70)
  input_vp2cpt_70 <- replace_zero(input_vp2cpt_70)
  input_q2cpt_70 <- replace_zero(input_q2cpt_70)
  input_vc3cpt_70 <- replace_zero(input_vc3cpt_70)
  input_vp3cpt_70 <- replace_zero(input_vp3cpt_70)
  input_vp23cpt_70 <- replace_zero(input_vp23cpt_70)
  input_q3cpt_70 <- replace_zero(input_q3cpt_70)
  input_q23cpt_70 <- replace_zero(input_q23cpt_70)

  input_vmax <<- input_vmax
  input_km <<- input_km
  input_cl <<- input_cl
  input_vc <<- input_vc
  input_vc2cpt <<-   input_vc2cpt
  input_vp2cpt <<-  input_vp2cpt
  input_q2cpt <<-   input_q2cpt
  input_vc3cpt <<-   input_vc3cpt
  input_vp3cpt <<-   input_vp3cpt
  input_vp23cpt <<- input_vp23cpt
  input_q3cpt <<- input_q3cpt
  input_q23cpt <<-   input_q23cpt
  input_vmax_70 <<- input_vmax_70
  input_cl_70 <<- input_cl_70
  input_vc_70 <<- input_vc_70
  input_vc2cpt_70 <<-   input_vc2cpt_70
  input_vp2cpt_70 <<-  input_vp2cpt_70
  input_q2cpt_70 <<-   input_q2cpt_70
  input_vc3cpt_70 <<-   input_vc3cpt_70
  input_vp3cpt_70 <<-   input_vp3cpt_70
  input_vp23cpt_70 <<- input_vp23cpt_70
  input_q3cpt_70 <<- input_q3cpt_70
  input_q23cpt_70 <<-   input_q23cpt_70
  input_prop.err <<- input_prop.err
  input_add.err <<-   input_add.err
  input_logn.sd <<-    input_logn.sd
  input_eta.cl <<- input_eta.cl
  input_eta.vc <<- input_eta.vc
  input_eta.vp <<-   input_eta.vp
  input_eta.vp2 <<-   input_eta.vp2
  input_eta.q <<- input_eta.q
  input_eta.q2 <<-   input_eta.q2
  input_eta.vmax <<-  input_eta.vmax
  input_eta.km <<-   input_eta.km
  input_etarsquare <<-   input_etarsquare

}
