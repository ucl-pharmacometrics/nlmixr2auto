
maptoValid <-
  function(string,
           search.space = 1,
           code.source = NULL) {
    string <- as.numeric(as.vector(as.matrix(string)))

    # --- Auto-detect code.source if not provided ---
    if (is.null(code.source)) {
      if (length(string) == 12) {
        code.source <- "GA"
      } else if (length(string) == 10) {
        code.source <- "ACO"
      } else {
        stop(
          "Cannot infer code.source from string length. Please specify code.source explicitly."
        )
      }
    }

    # --- GA Code Mapping ---
    if (toupper(code.source) == "GA") {
      cmpt.iv1 = string[1]
      cmpt.iv2 = string[2]
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

      # GA Specific Validity Corrections
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

      if ((eta.vc + eta.vp + eta.vp2 + eta.q + eta.q2) == 0 &
          mm == 0) {
        mcorr = 0
      }

      if ((eta.vc + eta.vp + eta.vp2 + eta.q + eta.q2) < 2 &
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

    # --- ACO Code Mapping ---
    if (toupper(code.source) == "ACO") {
      cmpt.iv = string[1]
      eta.km   = string[2]
      eta.vc   = string[3]
      eta.vp   = string[4]
      eta.vp2  = string[5]
      eta.q    = string[6]
      eta.q2   = string[7]
      mm       = string[8]
      mcorr    = string[9]
      rv       = string[10]

      # ACO Specific Validity Corrections
      if ((eta.vc + eta.vp + eta.vp2 + eta.q + eta.q2) == 0 &
          mm == 0) {
        mcorr = 0
      }

      if ((eta.vc + eta.vp + eta.vp2 + eta.q + eta.q2) < 2 &
          mm == 1 & eta.km == 0) {
        mcorr = 0
      }

      string <- c(cmpt.iv,
                  eta.km,
                  eta.vc,
                  eta.vp,
                  eta.vp2,
                  eta.q,
                  eta.q2,
                  mm,
                  mcorr,
                  rv)
    }

    return(string)
  }
