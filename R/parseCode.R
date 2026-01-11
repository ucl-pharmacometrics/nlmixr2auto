#' 2-bit code helper
#'
#' Internal utility used by decodeBinary() and encodeBinary() to convert between
#' a 2-bit representation and categorical values via lookup tables.
#'
#' @param mode Character, either "decode" or "encode".
#' @param value For decode: the first bit (0/1). For encode: the categorical
#'   value to encode.
#' @param bit2 For decode: the second bit (0/1). Ignored for encode.
#' @param decode_map Numeric vector of length 4 used for decoding, in the order
#'   corresponding to 00, 01, 10, 11. Use NA for illegal codes.
#' @param encode_map Named integer/numeric vector mapping categorical values
#'   (as names) to integer codes 0..3 (corresponding to 00..11).
#' @param param_name Character, parameter name used in error messages.
#'
#' @details
#' The helper supports two modes:
#' \itemize{
#'   \item decode: converts (bit1, bit2) to a categorical value using a length-4
#'     lookup table decode_map corresponding to 00, 01, 10, 11.
#'   \item encode: converts a categorical value to a 2-bit code using a named
#'     lookup table encode_map mapping values to integer codes 0..3
#'     (corresponding to 00..11).
#' }
#'
#' @return
#' For decode: a single numeric categorical value. For encode: an integer vector
#' length 2 containing bits c(bit1, bit2).
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Decode example (00/01 map to level 1).
#' .twoBitCode("decode", 0, 0, decode_map = c(1, 1, 2, 3))  # 1
#' .twoBitCode("decode", 0, 1, decode_map = c(1, 1, 2, 3))  # 1
#' .twoBitCode("decode", 1, 0, decode_map = c(1, 1, 2, 3))  # 2
#' .twoBitCode("decode", 1, 1, decode_map = c(1, 1, 2, 3))  # 3
#'
#' # Encode example (level 1 emits 01).
#' encode_map <- stats::setNames(c(1, 2, 3), c(1, 2, 3))
#' .twoBitCode("encode", 1, encode_map = encode_map)  # c(0, 1)
#' .twoBitCode("encode", 2, encode_map = encode_map)  # c(1, 0)
#' .twoBitCode("encode", 3, encode_map = encode_map)  # c(1, 1)
#'
#' # Decode 4-level example (00..11 map to 1..4).
#' .twoBitCode("decode", 0, 0, decode_map = c(1, 2, 3, 4))  # 1
#'
#' @seealso
#' \code{\link{decodeBinary}}, \code{\link{encodeBinary}}
#'
#' @export

.twoBitCode <- function(mode = c("decode", "encode"),
                        value,
                        bit2 = NULL,
                        decode_map = NULL,
                        encode_map = NULL,
                        param_name = "param") {
  mode <- match.arg(mode)

  if (mode == "decode") {
    lookup_index <- as.integer(value) * 2L + as.integer(bit2) + 1L
    decoded_value <- decode_map[lookup_index]
    if (is.na(decoded_value)) {
      stop(sprintf("Invalid %s encoding: (%d, %d)", param_name, value, bit2),
           call. = FALSE)
    }
    return(decoded_value)
  }

  code_value <- encode_map[as.character(value)]
  if (is.na(code_value)) {
    stop(sprintf("Invalid %s value: %s", param_name, value), call. = FALSE)
  }

  code_value <- as.integer(code_value)
  bit1 <- code_value %/% 2L
  bit2 <- code_value %% 2L
  return(c(bit1, bit2))
}


#' Decode binary encoding to categorical encoding
#'
#' Converts a binary-encoded GA chromosome (0/1 vector) into a categorical
#' parameter vector.
#'
#' @param binary_string Numeric vector of bits (0/1). Length must match the
#'   expected layout for the selected search.space.
#' @param search.space Character string specifying which search space to use.
#'   Options are "ivbase", "oralbase", or "custom". Default is "ivbase".
#' @param custom_config Optional named list defining a custom parameter structure.
#'   If provided, the parameter names are taken from the names of this list.
#'   If NULL, a default parameter structure is used based on the selected
#'   search space.
#'
#' @details
#' Supported search spaces:
#' \itemize{
#'   \item "ivbase": binary string has 12 bits and decodes to 10 values.
#'   \item "oralbase": binary string has 13 bits and decodes to 11 values.
#'   \item "custom": binary string has 29 bits and decodes to 24 values.
#' }
#'
#' Legacy layout ("ivbase" and "oralbase"):
#' \itemize{
#'   \item The first two bits encode the number of compartments (no.cmpt) as
#'     00 or 01 for 1, 10 for 2, and 11 for 3.
#'   \item The middle parameters are copied directly and preserve values such as
#'     0, 1, and -1.
#'   \item The last two bits encode the residual error model (rv) as 00 or 01
#'     for 1, 10 for 2, and 11 for 3.
#' }
#'
#' Custom layout ("custom"):
#' \itemize{
#'   \item Multi-level parameters are stored as 2-bit fields (for example,
#'     no.cmpt, absorption type, absorption delay, residual error model, and
#'     allometric scaling).
#'   \item Binary flags are stored as single bits, including absolute
#'     bioavailability (abs.bio), the Michaelis-Menten indicator (mm), and the
#'     correlation indicator (mcorr).
#'   \item Inter-individual variability indicators (eta.*) are stored as 16
#'     single-bit flags in a fixed order.
#' }
#'
#' For "custom", the categorical output order is:
#' \preformatted{
#' no.cmpt, abs.type, abs.delay, abs.bio,
#' eta.vmax, eta.km, eta.cl, eta.vc, eta.vp, eta.vp2, eta.q, eta.q2,
#' eta.ka, eta.tlag, eta.D2, eta.F1, eta.Fr, eta.mtt, eta.n, eta.bio,
#' mm, mcorr, rv, allometric_scaling
#' }
#'
#' @return Numeric vector of categorical parameter values.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' binary_iv <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1)
#' decodeBinary(binary_iv, "ivbase")
#'
#' binary_oral <- c(0, 1, rep(0, 9), 1, 1)
#' decodeBinary(binary_oral, "oralbase")
#'
#' binary_custom <- c(
#'   0, 0, 0, 1, 1, 1, 1, rep(0, 16), 0, 1, 1, 0, 0, 1
#' )
#' decodeBinary(binary_custom, "custom")
#'
#' @seealso
#' \code{\link{.twoBitCode}}, \code{\link{encodeBinary}}
#'
#' @export

decodeBinary <- function(binary_string,
                         search.space = "ivbase",
                         custom_config) {
  if (missing(custom_config)){
    custom_config <-NULL
  }
  if (search.space %in% c("ivbase", "oralbase")) {
    legacy_decode_map <- c(1, 1, 2, 3) # 00/01->1, 10->2, 11->3
    is_oral <- (search.space == "oralbase")
    expected_len <- if (is_oral)
      13L
    else
      12L

    if (length(binary_string) != expected_len) {
      stop(
        sprintf(
          "Length mismatch for %s: got %d bits, expected %d.",
          search.space,
          length(binary_string),
          expected_len
        ),
        call. = FALSE
      )
    }

    no.cmpt <- .twoBitCode(
      "decode",
      binary_string[1],
      binary_string[2],
      decode_map = legacy_decode_map,
      param_name = "no.cmpt"
    )

    if (is_oral) {
      middle_params <- binary_string[3:11]
      rv_bits <- binary_string[12:13]
    } else {
      middle_params <- binary_string[3:10]
      rv_bits <- binary_string[11:12]
    }

    rv <- .twoBitCode("decode",
                      rv_bits[1],
                      rv_bits[2],
                      decode_map = legacy_decode_map,
                      param_name = "rv")

    return(as.numeric(c(no.cmpt, middle_params, rv)))
  }

  if (search.space == "custom") {
    default_params <- c(
      "no.cmpt",
      "abs.type",
      "abs.delay",
      "abs.bio",
      "eta.vmax",
      "eta.km",
      "eta.cl",
      "eta.vc",
      "eta.vp",
      "eta.vp2",
      "eta.q",
      "eta.q2",
      "eta.ka",
      "eta.tlag",
      "eta.D2",
      "eta.F1",
      "eta.Fr",
      "eta.mtt",
      "eta.n",
      "eta.bio",
      "mm",
      "mcorr",
      "rv",
      "allometric_scaling"
    )

    params <- if (!is.null(custom_config) &&
                  !is.null(custom_config$params) &&
                  length(custom_config$params) > 0) {
      custom_config$params
    } else {
      default_params
    }

    # Which params use 2 bits in custom encoding
    twobit_params <- c("no.cmpt",
                       "abs.type",
                       "abs.delay",
                       "rv",
                       "allometric_scaling")

    # expected_bits = sum(2 for twobit, 1 otherwise) without nested functions
    is_twobit <- !is.na(match(params, twobit_params))
    expected_bits <- sum(1L + as.integer(is_twobit))

    if (length(binary_string) != expected_bits) {
      stop(
        sprintf(
          "Length mismatch for custom: got %d bits, expected %d.",
          length(binary_string),
          expected_bits
        ),
        call. = FALSE
      )
    }

    out <- numeric(length(params))
    i <- 1L

    for (k in seq_along(params)) {
      p <- params[k]

      if (!is.na(match(p, twobit_params))) {
        b1 <- binary_string[i]
        b2 <- binary_string[i + 1L]
        i <- i + 2L

        out[k] <- switch(
          p,
          "no.cmpt" = .twoBitCode(
            "decode",
            b1,
            b2,
            decode_map = c(1, 1, 2, 3),
            param_name = "no.cmpt"
          ),
          "abs.type" = .twoBitCode(
            "decode",
            b1,
            b2,
            decode_map = c(1, 2, 3, 4),
            param_name = "abs.type"
          ),
          "abs.delay" = .twoBitCode(
            "decode",
            b1,
            b2,
            decode_map = c(0, 1, 2, 0),
            param_name = "abs.delay"
          ),
          "rv" = .twoBitCode(
            "decode",
            b1,
            b2,
            decode_map = c(1, 2, 3, 4),
            param_name = "rv"
          ),
          "allometric_scaling" = .twoBitCode(
            "decode",
            b1,
            b2,
            decode_map = c(0, 1, 2, 3),
            param_name = "allometric_scaling"
          ),
          stop(
            sprintf("Unknown 2-bit parameter in custom space: %s", p),
            call. = FALSE
          )
        )

      } else {
        out[k] <- as.numeric(binary_string[i])
        i <- i + 1L
      }
    }

    return(as.numeric(out))
  }

  stop('search.space must be "ivbase", "oralbase", or "custom".',
       call. = FALSE)
}


#' Encode categorical encoding to binary encoding
#'
#' Converts a categorical parameter vector into a binary-encoded GA chromosome.
#'
#' @param categorical_string Numeric vector of categorical parameter values.
#'   Length must match the expected layout for the selected \code{search.space}.
#' @param search.space Character string specifying which search space to use.
#'   Options are "ivbase", "oralbase", or "custom". Default is "ivbase".
#' @param custom_config Optional named list defining a custom parameter structure.
#'   If provided, the parameter names are taken from the names of this list.
#'   If NULL, a default parameter structure is used based on the selected
#'   search space.
#'
#' @details
#' This function converts a vector of categorical parameter values into a 0/1
#' bit string (a GA chromosome). The required input layout and the encoding rules
#' depend on the selected search space.
#'
#' \strong{Built-in search spaces: ivbase / oralbase}
#'
#' \itemize{
#'   \item Expected input length:
#'     \itemize{
#'       \item ivbase: 10 categorical values
#'       \item oralbase: 11 categorical values
#'     }
#'   \item Structure: the first value is no.cmpt and the last value is rv.
#'         All middle values are copied as-is (they are expected to be 0/1 flags).
#'         Specifically, ivbase copies indices 2..9 and oralbase copies indices 2..10.
#'   \item no.cmpt and rv use the legacy 3-level 2-bit encoding:
#'     \itemize{
#'       \item 1 -> 01
#'       \item 2 -> 10
#'       \item 3 -> 11
#'     }
#'     The code 00 is not used in this mapping.
#'   \item Output length:
#'     \itemize{
#'       \item ivbase: 12 bits (2 + 8 + 2)
#'       \item oralbase: 13 bits (2 + 9 + 2)
#'     }
#' }
#'
#' \strong{Custom search space: custom}
#'
#' For search.space = "custom", the function reads categorical values in the
#' order given by params. If custom_config$params is provided and non-empty,
#' that vector defines the parameter names and order; otherwise, the default
#' 24-parameter layout is used:
#' \preformatted{
#' no.cmpt, abs.type, abs.delay, abs.bio,
#' eta.vmax, eta.km, eta.cl, eta.vc, eta.vp, eta.vp2, eta.q, eta.q2,
#' eta.ka, eta.tlag, eta.D2, eta.F1, eta.Fr, eta.mtt, eta.n, eta.bio,
#' mm, mcorr, rv, allometric_scaling
#' }
#'
#' \itemize{
#'   \item Parameters encoded with 2 bits: no.cmpt, abs.type, abs.delay, rv,
#'         and allometric_scaling.
#'   \item All other parameters must be single-bit flags (0 or 1) and are appended
#'         directly to the chromosome.
#' }
#'
#' \strong{2-bit encoding rules}
#' \itemize{
#'   \item no.cmpt (1..3): 1->01, 2->10, 3->11 (00 unused)
#'   \item abs.type (1..4): 1->00, 2->01, 3->10, 4->11
#'   \item abs.delay (0..2): 0->00, 1->01, 2->10 (11 unused)
#'   \item rv (1..4): 1->00, 2->01, 3->10, 4->11
#'   \item allometric_scaling (0..3): 0->00, 1->01, 2->10, 3->11
#' }

#' @return Numeric vector of bits (0/1).
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # ivbase: 10 categorical -> 12 bits
#' cat_iv <- c(1, rep(0, 8), 3)
#' encodeBinary(cat_iv, "ivbase")
#'
#' # Custom: 24 categorical -> 29 bits
#' cat_custom <- c(
#'   1, 2, 0, 1,
#'   rep(0, 16),
#'   0, 1, 3, 1
#' )
#' encodeBinary(cat_custom, "custom")
#'
#' @seealso
#' \code{\link{.twoBitCode}}, \code{\link{decodeBinary}}
#'
#' @export
#'
encodeBinary <- function(categorical_string,
                         search.space = "ivbase",
                         custom_config) {

  if (missing(custom_config)){
    custom_config <-NULL
  }
  if (search.space %in% c("ivbase", "oralbase")) {
    legacy_encode_map <-
      stats::setNames(c(1, 2, 3), c(1, 2, 3)) # 1->01,2->10,3->11
    is_oral <- (search.space == "oralbase")
    expected_len <- if (is_oral)
      11L
    else
      10L

    if (length(categorical_string) != expected_len) {
      stop(
        sprintf(
          "Length mismatch for %s: got %d values, expected %d.",
          search.space,
          length(categorical_string),
          expected_len
        ),
        call. = FALSE
      )
    }

    no.cmpt <- categorical_string[1]
    rv <- categorical_string[length(categorical_string)]

    no_bits <- .twoBitCode("encode",
                           no.cmpt,
                           encode_map = legacy_encode_map,
                           param_name = "no.cmpt")
    rv_bits <- .twoBitCode("encode",
                           rv,
                           encode_map = legacy_encode_map,
                           param_name = "rv")

    middle_params <-
      if (is_oral)
        categorical_string[2:10]
    else
      categorical_string[2:9]
    return(as.numeric(c(no_bits, middle_params, rv_bits)))
  }

  if (search.space == "custom") {
    default_params <- c(
      "no.cmpt",
      "abs.type",
      "abs.delay",
      "abs.bio",
      "eta.vmax",
      "eta.km",
      "eta.cl",
      "eta.vc",
      "eta.vp",
      "eta.vp2",
      "eta.q",
      "eta.q2",
      "eta.ka",
      "eta.tlag",
      "eta.D2",
      "eta.F1",
      "eta.Fr",
      "eta.mtt",
      "eta.n",
      "eta.bio",
      "mm",
      "mcorr",
      "rv",
      "allometric_scaling"
    )

    params <- if (!is.null(custom_config) &&
                  !is.null(custom_config$params) &&
                  length(custom_config$params) > 0) {
      custom_config$params
    } else {
      default_params
    }
    if (length(categorical_string) != length(params)) {
      stop(
        sprintf(
          "Length mismatch for custom: got %d values, expected %d.",
          length(categorical_string),
          length(params)
        ),
        call. = FALSE
      )
    }
    twobit_params <- c("no.cmpt",
                       "abs.type",
                       "abs.delay",
                       "rv",
                       "allometric_scaling")

    bits <- integer(0)
    for (k in seq_along(params)) {
      p <- params[k]
      x <- categorical_string[k]

      if (!is.na(match(p, twobit_params))) {
        bits_k <- switch(
          p,
          "no.cmpt" = .twoBitCode(
            "encode",
            x,
            encode_map = stats::setNames(c(1, 2, 3), c(1, 2, 3)),
            param_name = "no.cmpt"
          ),
          "abs.type" = .twoBitCode(
            "encode",
            x,
            encode_map = stats::setNames(0:3, 1:4),
            param_name = "abs.type"
          ),
          "abs.delay" = .twoBitCode(
            "encode",
            x,
            encode_map = stats::setNames(c(0, 1, 2), c(0, 1, 2)),
            param_name = "abs.delay"
          ),
          "rv" = .twoBitCode(
            "encode",
            x,
            encode_map = stats::setNames(0:3, 1:4),
            param_name = "rv"
          ),
          "allometric_scaling" = .twoBitCode(
            "encode",
            x,
            encode_map = stats::setNames(0:3, 0:3),
            param_name = "allometric_scaling"
          ),
          stop(
            sprintf("Unknown 2-bit parameter in custom space: %s", p),
            call. = FALSE
          )
        )
        bits <- c(bits, as.integer(bits_k))
      } else {
        if (!(x %in% c(0, 1))) {
          stop(sprintf("%s must be 0 or 1 in custom encoding.", p),
               call. = FALSE)
        }
        bits <- c(bits, as.integer(x))
      }
    }
    return(as.numeric(bits))
  }
  stop('search.space must be "ivbase", "oralbase", or "custom".',
       call. = FALSE)
}
