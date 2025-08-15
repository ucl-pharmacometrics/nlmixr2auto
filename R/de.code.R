#' Revise invalid model codes to valid models (GA)
#'
#' Revise invalid model codes to ensure strings generated from GA are valid based on the specified rules.
#'
#' @param string A binary string representing the model configuration.
#' @param search.space An integer indicating the search space.
#'
#' @return A numeric vector representing the revised model string.
#'
#' @examples
#' \dontrun{
#' string <- c(0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1)
#' de.code(string,1)
#' }
#'
#' @export

de.code <- function(string,
                    search.space) {
  if (search.space == 1) {
    string <- as.numeric(as.vector(as.matrix(string)))
    # de.code, revise the invalid model into valid models
    no.cmpt1 = string[1]
    no.cmpt2 = string[2]
    eta.km   = string[3]
    eta.vc   = string[4]
    eta.vp   = string[5]
    eta.vp2  = string[6]
    eta.q    = string[7]
    eta.q2   = string[8]
    mm       = string[9]
    mcorr    = string[10]
    rv1      = string[11]
    rv2      = string[12]

    if (cmpt.iv1 == 0 & cmpt.iv2 == 0) {
      cmpt.iv2 = 1
    }

    if (cmpt.iv1 == 0) {
      eta.vp = 0
      eta.vp2 = 0
      eta.q = 0
      eta.q2 = 0
    }

    if (cmpt.iv2 == 0) {
      eta.vp2 = 0
      eta.q2 = 0
    }

    if (mm == 0) {
      eta.km = 0
    }

    if ((eta.vc + eta.vp + eta.vp2 + eta.q + eta.q2) == 0 & mm == 0) {
      mcorr = 0
    }

    if ((eta.vc + eta.vp + eta.vp2 + eta.q + eta.q2) <2 &
        mm == 1 & eta.km == 0) {
      mcorr = 0
    }

    if (rv1 == 0 & rv2 == 0) {
      rv2 = 1
    }
    string <- c(
      cmpt.iv1,
      cmpt.iv2,
      eta.km,
      eta.vc,
      eta.vp,
      eta.vp2,
      eta.q,
      eta.q2,
      mm,
      mcorr,
      rv1,
      rv2
    )
  }
  return(string)
}



#' Revise invalid model codes to valid models  (Type2)
#'
#' Revise invalid model codes to ensure strings generated from ACO are valid based on the specified rules.
#'
#' @param string A binary string representing the model configuration.
#' @param search.space An integer indicating the search space.
#'
#' @return A numeric vector representing the revised model string.
#'
#' @examples
#' \dontrun{
#' string <- c(1, 0, 0, 0, 0, 0, 0, 1, 1, 1)
#' aco.de.code(string,1)
#' }
#'
#' @export


de.code2 <- function(string,
                    search.space) {
  if (search.space == 1) {

    # de.code, revise the invalid model into valid models
    cmpt.iv = string[1]
    eta.km   = string[2]
    eta.vc   = string[3]
    eta.vp   = string[4]
    eta.vp2  = string[5]
    eta.q    = string[6]
    eta.q2   = string[7]
    mm       = string[8]
    mcorr    = string[9]
    rv    = string[10]

    string<-c(cmpt.iv,
              eta.km,
              eta.vc,
              eta.vp,
              eta.vp2,
              eta.q,
              eta.q2,
              mm,
              mcorr,
              rv)

    if ((eta.vc + eta.vp + eta.vp2 + eta.q + eta.q2) == 0 & mm == 0) {
      mcorr = 0
    }

    if (( eta.vc +eta.vp + eta.vp2 + eta.q + eta.q2) < 2 &
        mm == 1 & eta.km == 0) {
      mcorr = 0
    }


    string <- c(
      cmpt.iv,
      eta.km,
      eta.vc,
      eta.vp,
      eta.vp2,
      eta.q,
      eta.q2,
      mm,
      mcorr,
      rv
    )
  }

  return(string)
}
