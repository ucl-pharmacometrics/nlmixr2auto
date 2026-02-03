#' Generate initial parameter table for pharmacometric model estimation
#'
#' Creates a structured parameter table containing initial estimates with constraints for
#' base parameters, inter-individual variability (ETA), residual errors (SIGMA), and
#' correlation terms (OMEGA) to initialize nonlinear mixed-effects model fitting.
#'
#' @details
#' This table includes:
#' \itemize{
#'   \item Base PK parameters (absorption, clearance, volumes, etc.) in log-scale
#'   \item Michaelis-Menten kinetics parameters (vmax, km)
#'   \item Absorption parameters including zero-order, mixed-order, and transit compartment models
#'   \item Residual variability components (additive and proportional error)
#'   \item Inter-individual variability (ETA) terms with variance parameters
#'   \item Correlation parameters between ETA terms in two blocks:
#'   \itemize{
#'     \item Block 1: vmax and km parameters
#'     \item Block 2: clearance, volumes, and inter-compartmental clearance
#'   }
#' }
#'
#' Parameters are organized with:
#' \itemize{
#'   \item Name: Parameter name following standard nomenclature
#'   \item init: Initial estimate for model fitting
#'   \item lb/ub: Lower/upper bounds for parameter estimation
#'   \item fixed: Flag indicating fixed parameters (1) vs estimated (0)
#'   \item Description: Plain-text explanation of parameter meaning
#' }
#'
#' @return A data.frame with 29 columns containing parameter specifications. The structure includes:
#' \describe{
#'   \item{Name}{Parameter name (e.g., "lcl", "eta.vc", "cor.eta_cl_vc")}
#'   \item{init}{Numeric initial value for parameter estimation}
#'   \item{lb}{Lower bound constraint (use -Inf for unconstrained)}
#'   \item{ub}{Upper bound constraint (use Inf for unconstrained)}
#'   \item{fixed}{Integer flag indicating whether parameter should be fixed (1) or estimated (0)}
#'   \item{Description}{Text description of parameter's biological/pharmacometric meaning}
#' }
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Generate default parameter table
#' initialize_param_table()
#'
#' @export
#'
initialize_param_table <- function() {
  base_params <- list(
    # Base parameters
    list(
      Name = "lka",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of the absorption rate constant"
    ),
    list(
      Name = "lcl",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of clearance"
    ),
    # list(
    #   Name = "lvc",
    #   init = 1,
    #   lb = -Inf,
    #   ub = Inf,
    #   fixed = 0,
    #   Description = "Log of central volume of distribution"
    # ),
    # list(
    #   Name = "lvp",
    #   init = 1,
    #   lb = -Inf,
    #   ub = Inf,
    #   fixed = 0,
    #   Description = "Log of peripheral volume of distribution for the first peripheral compartment"
    # ),

    # Modified central volume of distribution
    list(
      Name = "lvc1cmpt",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of central volume of distribution for 1-compartment model"
    ),
    list(
      Name = "lvc2cmpt",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of central volume of distribution for 2-compartment model"
    ),
    list(
      Name = "lvc3cmpt",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of central volume of distribution for 3-compartment model"
    ),

    # Modified peripheral volume of distribution
    list(
      Name = "lvp2cmpt",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of peripheral volume for 2-compartment model"
    ),
    list(
      Name = "lvp3cmpt",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of peripheral volume for 3-compartment model"
    ),
    list(
      Name = "lq2cmpt",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of intercompartmental clearance between central and first peripheral compartments for 2-compartment model"
    ),
    list(
      Name = "lq3cmpt",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of intercompartmental clearance between central and first peripheral compartments for 3-compartment model"
    ),
    list(
      Name = "lvp2",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of volume for the second peripheral compartment"
    ),
    list(
      Name = "lq2",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of intercompartmental clearance between central and second peripheral compartments"
    ),

    # Michaelis-Menten parameters
    list(
      Name = "lvmax",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of the maximum elimination rate"
    ),
    list(
      Name = "lkm",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of the concentration at which the rate is half of Vmax"
    ),

    # Absorption parameters
    list(
      Name = "lD2",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of zero-order absorption duration"
    ),
    list(
      Name = "lF1",
      init = 0.01,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of absolute bioavailability, typically used in IV and oral mixed dosing"
    ),
    list(
      Name = "lFr",
      init = 0.01,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of the fraction absorbed through first-order kinetics in a mixed absorption model"
    ),

    # absorption delay parameter
    list(
      Name = "ltlag",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of lag time before absorption begins"
    ),
    list(
      Name = "lmtt",
      init = 1,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of mean transit time"
    ),
    list(
      Name = "ln",
      init = 0.01,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of the number of transit compartments"
    ),
    list(
      Name = "lbio",
      init = 0.01,
      lb = -Inf,
      ub = Inf,
      fixed = 0,
      Description = "Log of the bioavailability of transit compartment model"
    ),

    # Residual variability
    list(
      Name = "sigma_add",
      init = 1,
      lb = 0,
      ub = Inf,
      fixed = 0,
      Description = "Additive residual error component, constant across all concentrations."
    ),
    list(
      Name = "sigma_prop",
      init = 0.1,
      lb = 0,
      ub = 1,
      fixed = 0,
      Description = "Proportional residual error component, scales with drug concentration."
    )
  )

  # Define eta parameters
  eta_names <-
    c(
      "ka",
      "cl",
      "vc",
      "vp",
      "q",
      "vp2",
      "q2",
      "vmax",
      "km",
      "tlag",
      "D2",
      "mtt",
      "n",
      "bio",
      "F1",
      "Fr"
    )
  eta_params <- lapply(eta_names, function(p) {
    list(
      Name = paste0("eta.", p),
      init = 0.1,
      lb = 0,
      ub = 1,
      fixed = 0,
      Description = paste0(
        "Inter-individual variability in",
        p,
        ", modeled with variance as the initial estimate"
      )
    )
  })

  # Combine base and eta parameters
  base_params <- append(base_params, eta_params)

  # Define omega blocks
  omega_block1 <- c("eta.vmax", "eta.km")
  omega_block2 <-
    c("eta.cl", "eta.vc", "eta.vp", "eta.vp2", "eta.q", "eta.q2")

  # Add correlation parameters for block 1
  combination_block1 <- utils::combn(omega_block1, 2, simplify = FALSE)
  for (combo in combination_block1) {
    combined_name <-
      paste0("cor.eta_",
             gsub("eta.", "", combo[1]),
             "_",
             gsub("eta.", "", combo[2]))
    base_params <- append(base_params, list(
      list(
        Name = combined_name,
        init = 0.1,
        lb = 0,
        ub = 1,
        fixed = 0,
        Description = paste0("Correlation between ", combo[1], " and ", combo[2])
      )
    ))
  }

  # Add correlation parameters for block 2
  combination_block2 <- utils::combn(omega_block2, 2, simplify = FALSE)
  for (combo in combination_block2) {
    combined_name <-
      paste0("cor.eta_",
             gsub("eta.", "", combo[1]),
             "_",
             gsub("eta.", "", combo[2]))
    base_params <- append(base_params, list(
      list(
        Name = combined_name,
        init = 0.1,
        lb = 0,
        ub = 1,
        fixed = 0,
        Description = paste0("Correlation between ", combo[1], " and ", combo[2])
      )
    ))
  }

  # Convert to data.frame
  param_table <- do.call(rbind, lapply(base_params, as.data.frame))
  return(param_table)
}


#' Automatically generate a parameter table with initial estimates
#'
#' Constructs a parameter table for nlmixr2 model fitting. It supports:
#' \itemize{
#'   \item Direct use of a user-provided parameter table.
#'   \item Automatic initialization of parameters from data using
#'         \code{getPPKinits()}.
#'   \item Fallback to a default parameter table created by
#'         \code{initialize_param_table()}.
#' }
#'
#' @details
#' When `nlmixr2autoinits = TRUE`, this function estimates initial values
#' from data, applies a name mapping to internal model parameters,
#' performs log transformations where appropriate, and replaces
#' problematic log values (e.g. log(0) or `NA`) with `log(0.01)` for
#' numerical stability.
#'
#' @param dat A data frame containing observed data (required if
#'   `nlmixr2autoinits = TRUE`).
#' @param param_table Optional. A user-provided parameter table (if
#'   provided, all other logic is skipped).
#' @param nlmixr2autoinits Logical. Whether to automatically estimate
#'   initial values using \code{getPPKinits()}. Default is `TRUE`.
#' @param foldername Character string specifying the folder name for storing
#'   `nlmixr2autoinits` outputs. If \code{NULL} (default), \code{tempdir()}
#'   is used for temporary storage. If specified, a cache directory
#'   is created in the current working directory.
#' @param filename Character string specifying the base name for model output
#'   files generated during evaluation.
#' @param out.inits Logical flag indicating whether the results returned
#'   by the automated initialization procedure should be saved to an RDS
#'   file. When TRUE, the output of the initialization step is written to
#'   disk for reproducibility or debugging purposes.
#' @param ... Additional arguments passed to \code{getPPKinits()}.
#'
#' @return A `data.frame` representing the parameter table with initial
#'   estimates, ready for use in `nlmixr2()`.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' \donttest{
#' auto_param_table(dat = pheno_sd)
#' }
#' @seealso \code{\link[nlmixr2autoinit:getPPKinits]{getPPKinits}},
#' \code{\link{initialize_param_table}}
#'
#' @export

auto_param_table <- function(dat = NULL,
                             param_table = NULL,
                             nlmixr2autoinits = TRUE,
                             foldername= NULL,
                             filename = "test",
                             out.inits = TRUE,
                             ...) {
  # Case 1: User has provided a param_table — use it as-is
  if (!is.null(param_table)) {
    return(param_table)
  }

  # Case 2: Automatically generate initial estimates using getPPKinits()
  if (nlmixr2autoinits) {
    if (is.null(dat)) {
      stop("`dat` must be provided when `nlmixr2autoinits = TRUE`")
    }
    # Run automated initialization procedure
    getinits. <- nlmixr2autoinit::getPPKinits(dat, ...)

    if (is.null(foldername) || !nzchar(foldername)) {
      # foldername <-
      #   paste0("modRunCache_", filename, "_", digest::digest(dat))
      foldername <- tempdir()
    }
    if (isTRUE(out.inits)) {
      rds_file <- file.path(foldername,paste0(filename, ".inits.RDS"))
      saveRDS(getinits., rds_file)
    }
    # Extract the recommended initial estimates
    inits.out <- getinits.$Recommended_initial_estimates
    inits.out$Values <-
      suppressWarnings(as.numeric(inits.out$Values))

    # Remove rows with missing values
    inits.out <- inits.out[!is.na(inits.out$Values),]

    # Map external parameter names to internal model variable names
    param_map <- c(
      "Ka" = "lka",
      "CL" = "lcl",
      "Vd" = "lvc1cmpt",
      "Vmax" = "lvmax",
      "Km" = "lkm",
      "Vc(2CMPT)" = "lvc2cmpt",
      "Vp(2CMPT)" = "lvp2cmpt",
      "Q(2CMPT)" = "lq2cmpt",
      "Vc(3CMPT)" = "lvc3cmpt",
      "Vp(3CMPT)" = "lvp3cmpt",
      "Vp2(3CMPT)" = "lvp2",
      "Q(3CMPT)" = "lq3cmpt",
      "Q2(3CMPT)" = "lq2",
      "Sigma additive" = "sigma_add",
      "Sigma proportional" = "sigma_prop"
    )

    # Keep only recognized parameters
    valid_rows <- inits.out$Parameters %in% names(param_map)
    filtered_vals <- inits.out[valid_rows,]

    # Rename parameters using the internal naming convention
    vals <- filtered_vals$Values
    names(vals) <- param_map[filtered_vals$Parameters]

    # Identify parameters that require log transformation
    log_transform_names <-
      setdiff(names(vals), c("sigma_add", "sigma_prop"))
    log_vals <- vals

    # Apply log transformation (base e) to required parameters
    log_vals[log_transform_names] <- log(vals[log_transform_names])

    # Replace log(0), -Inf, or NA with log(0.01) to ensure numerical stability
    log_vals[log_transform_names][is.na(log_vals[log_transform_names]) |
                                    log_vals[log_transform_names] == -Inf] <-
      log(0.01)

    # Generate a default parameter table
    param_table <- initialize_param_table()

    # Replace default initial values with estimated/log-transformed values
    # param_table <- param_table %>%
    #   dplyr::mutate(init = dplyr::if_else(Name %in% names(log_vals), log_vals[Name], init))
    param_table <- dplyr::mutate(
      param_table,
      init = dplyr::if_else(Name %in% names(log_vals), log_vals[Name], init)
    )
    return(param_table)
  }

  # Case 3: Fallback — use default parameter table if no inits provided and auto is off
  return(initialize_param_table())
}


#' Initialize model parameters from parameter table
#'
#' Generates parameter initialization code based on a parameter table, handling
#' both fixed and estimated parameters.
#'
#' @param param_name Character, name of the parameter to initialize (without
#'   "l" prefix)
#' @param param_table Dataframe containing parameter specifications, must
#'   include:
#'   - `Name`: Character parameter names with "l" prefix (e.g., "lka"
#'     corresponds to param_name="ka")
#'   - `init`: Numeric initial values
#'   - `fixed`: Integer flag (0/1) indicating fixed status (1 = fixed)
#'
#' @return Character vector containing generated initialization code line.
#'   Format:
#'   - Fixed parameters: `<param_name> <- fix(initial_value)`
#'   - Estimated parameters: `l<param_name> <- initial_value`
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Create sample parameter table
#' param_table <- initialize_param_table()
#'
#' # Generate initialization code
#' initialize_param("ka", param_table)  # Returns "ka <- fix(0.500)"
#' initialize_param("cl", param_table)  # Returns "lcl <- 1.200"
#' @export

initialize_param <- function(param_name, param_table) {
  # Retrieve the corresponding row from the parameter table
  param_row <-
    param_table[param_table$Name == paste0("l", param_name),]

  # Check if the parameter exists in the table
  if (nrow(param_row) == 0) {
    stop(paste("Parameter", paste0("l", param_name), "not found in parameter table."))
  }
  # Check for duplicate entries in the table
  if (nrow(param_row) > 1) {
    stop(paste("Multiple entries found for", paste0("l", param_name)))
  }

  # Extract the initial value and fixed status from the row
  init_value <- param_row$init
  fixed_flag <- param_row$fixed

  # Round the initial value to 3 decimal places
  rounded_value <- round(init_value, 3)

  # Determine output variable name
  output_name <- if (param_name %in% c("vc1cmpt", "vc2cmpt", "vc3cmpt")) {
    "vc"
  } else if (param_name %in% c("vp2cmpt", "vp3cmpt")) {
    "vp"
  } else if (param_name %in% c("q2cmpt", "q3cmpt")) {
    "q"
  } else {
    param_name
  }

  # Generate the initialization code based on the fixed flag
  if (fixed_flag == 1) {
    # For fixed parameters: use the parameter name without 'l' prefix and wrap in fix()
    code_line <-
      sprintf("%s <- fix(%s)",
              output_name,
              format(rounded_value, nsmall = 1))
  } else {
    # For estimated parameters: use the parameter name with 'l' prefix and assign directly
    code_line <-
      sprintf("l%s <- %s", output_name, format(rounded_value, nsmall = 1))
  }

  return(code_line)
}

#' Add inter-individual variability to a parameter
#'
#' Defines a model string for a parameter, optionally adding inter-individual
#' variability.
#'
#' @param param_name Character. The name of the parameter.
#' @param eta_flag Integer. If 1, inter-individual variability is added;
#'   otherwise, it is not.
#' @param param_table Data frame. A table containing parameter details with
#'   columns `Name`, `init`, and optionally bounds like `lb` and `ub`.
#' @param param.type Integer. Transformation type: 1=Exponential, 2=Logistic.
#'   Defaults to 1.
#'
#' @return A list containing:
#'   \item{mod}{Character. The model string for the parameter.}
#'   \item{eta_init}{Character. The initialization string for the variability
#'   parameter (if applicable).}
#'
#' @author Zhonghui Huang
#'
#' @examples
#' param_table <- initialize_param_table()
#' add_variability("cl", 1, param_table)
#'
#' @export
#'
add_variability <- function(param_name, eta_flag, param_table,param.type=1) {

  eta_init <- ""
  if (param.type == 1) {
    # Exponential Transformation
    mod <- paste0(param_name, " = exp(l", param_name)
  }
  if (param.type == 2) {
    # Logistic Transformation
    mod <- paste0(param_name, " = 1 / (1 + exp(-(l", param_name)
  }

  # If eta_flag is enabled, add inter-individual variability
  if (eta_flag == 1) {
    mod <- paste0(mod, "+eta.", param_name)
    # Retrieve `eta.<param_name>` row from the parameter table
    eta_row <-
      param_table[param_table$Name == paste0("eta.", param_name), ]
    if (nrow(eta_row) == 0) {
      stop(paste0(
        "Parameter 'eta.",
        param_name,
        "' not found in parameter table."
      ))
    }

    # Format the eta initialization statement
    eta_init <- paste0("eta.", param_name, " ~ ", eta_row$init)
  }

  # Close the parameter model definition
  if (param.type == 1) {
    mod <- paste0(mod, ")")
  }
  if (param.type == 2) {
    mod <- paste0(mod, ")))")
  }

  return(list(mod = mod, eta_init = eta_init))
}

#' Generate omega block Code for nlmixr2 model
#'
#' Generates the code for the omega block matrix in nlmixr2 syntax,
#' supporting both independent variance terms and correlated covariance structures.
#'
#' @param param_list A character vector of parameter names requiring inter-individual
#'   variability (IIV) terms.
#' @param mcorr Integer flag indicating covariance structure:
#'   - `0`: Generate independent variance terms only
#'   - `1`: Generate full block covariance structure
#' @param eta_table A data frame containing eta initialization values and correlation
#'   coefficients. Must contain columns:
#'   - `Name`: Parameter names (format "eta.X" for variances, "cor.eta_X_Y" for correlations)
#'   - `init`: Initialization values for variance/covariance terms
#'
#' @return A character string containing nlmixr2 omega matrix specification code.
#'   - When `mcorr = 0`: Returns individual variance terms in formula syntax
#'   - When `mcorr = 1`: Returns covariance block structure in matrix syntax
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Example eta table structure
#' eta_table <- initialize_param_table()
#'
#' # Generate independent terms
#' omega_block(c("eta.cl", "eta.vc"), mcorr = 0, eta_table)
#'
#' # Generate covariance block
#' omega_block(c("eta.cl", "eta.vc"), mcorr = 1, eta_table)
#'
#' @export
#'
omega_block <- function(param_list, mcorr, eta_table) {

  # Handling mcorr=0 case: output only independent variance terms
  if (mcorr == 0) {
    eta_values <-
      stats::setNames(eta_table$init[grep("^eta\\.", eta_table$Name)],
               eta_table$Name[grep("^eta\\.", eta_table$Name)])

    # Generate individual lines for each parameter in formula syntax
    single_lines <- sapply(param_list, function(p) {
      sprintf("%s ~ %s", p, eta_values[[p]])
    })

    # Combine all lines with newline separators
    return(paste(single_lines, collapse = "\n"))
  }


  # If the parameter list has one or no parameters left, return an empty string
  if (length(param_list) <= 1)
    return(stop("mcorr was set to 1 but only one IIV on parameter was found"))

  # Create lookup for eta values
  eta_values <-
    stats::setNames(eta_table$init[grep("^eta\\.", eta_table$Name)],
             eta_table$Name[grep("^eta\\.", eta_table$Name)])

  # Create lookup for cor.eta values
  cor_eta_values <-
    stats::setNames(eta_table$init[grep("^cor\\.eta_", eta_table$Name)],
             eta_table$Name[grep("^cor\\.eta_", eta_table$Name)])

  # Generate a lower triangular matrix representation
  correlation_lines <- c()
  for (i in seq_along(param_list)) {
    line_values <- c()  # Store values for the current line
    for (j in seq_len(i)) {
      if (i == j) {
        # Diagonal elements: variance
        line_values <- c(line_values, eta_values[[param_list[i]]])
      } else {
        # Off-diagonal elements: covariance
        cov_name <-
          paste0(
            "cor.eta_",
            gsub("eta\\.", "", param_list[j]),
            "_",
            gsub("eta\\.", "", param_list[i])
          )
        if (!cov_name %in% names(cor_eta_values)) {
          # Try reversed order
          cov_name_rev <- paste0(
            "cor.eta_",
            gsub("eta\\.", "", param_list[i]),
            "_",
            gsub("eta\\.", "", param_list[j])
          )
          if (cov_name_rev %in% names(cor_eta_values)) {
            cov_name <- cov_name_rev
          } else {
            stop(paste0("Missing correlation value for: ", cov_name))
          }
        }
        cov <-
          sqrt(eta_values[[param_list[j]]]) * sqrt(eta_values[[param_list[i]]]) * cor_eta_values[[cov_name]]
        line_values <- c(line_values, cov)
      }
    }
    # Combine the line into a string and add it to correlation_lines
    correlation_lines <-
      c(correlation_lines, paste(line_values, collapse = ", "))
  }

  # Generate the nlmixr2 code in the desired format
  correlation_code <- paste0(
    paste(param_list, collapse = " + "),
    " ~ c(\n  ",
    paste(correlation_lines, collapse = ",\n  "),
    "\n)"
  )

  return(correlation_code)
}


#' Build ODE model lines for pharmacokinetic modeling
#'
#' Constructs a system of ordinary differential equations (ODEs) for pharmacokinetic modeling
#' with various configurations including different absorption models, compartmental structures,
#' and elimination kinetics.
#'
#' @param mm Michaelis-Menten elimination flag. 0 = linear elimination (default),
#'           1 = Michaelis-Menten elimination.
#' @param no.cmpt Number of compartments. Supported values: 1 (central only),
#'                2 (central + peripheral), or 3 (central + 2 peripheral).
#' @param route Administration route. Options: "bolus" (IV), "oral",
#'              or "mixed_iv_oral".
#' @param abs.bio Bioavailability estimation flag. 0 = no bioavailability estimation (default),
#'                1 = include bioavailability parameter (F1).
#' @param abs.type Absorption type for oral route:
#'                 \itemize{
#'                   \item 1 = First-order absorption (default)
#'                   \item 2 = Zero-order absorption
#'                   \item 3 = Sequential zero-first order absorption
#'                   \item 4 = Dual zero-first order absorption
#'                 }
#' @param abs.delay Absorption delay type:
#'                  \itemize{
#'                    \item 0 = No delay (default)
#'                    \item 1 = Lag time (tlag)
#'                    \item 2 = Transit compartment model
#'                  }
#'
#' @details
#'  Parameter Constraints:
#' The function includes error checking for incompatible parameter combinations:
#' \itemize{
#'   \item abs.bio=1 cannot be used with abs.type=4 or abs.delay=3
#'   \item Dual absorption (abs.type=4) not supported for mixed routes
#' }
#'
#' @return A character vector containing ODE equations for the specified configuration.
#'         Includes differential equations for drug compartments, absorption models,
#'         and derived parameters.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Two-compartment model with first-order absorption
#' build_odeline(no.cmpt = 2, route = "oral")
#'
#' # One-compartment IV model with Michaelis-Menten elimination
#' build_odeline(mm = 1, route = "bolus")
#' @export

build_odeline <- function(mm = 0,
                          no.cmpt = 1,
                          route = "bolus",
                          abs.bio = 0,
                          abs.type = 1,
                          abs.delay = 0) {

  # Initialize the central ODE line based on mm (Michaelis-Menten)
  if (mm == 0) {
    central_odeline <- "- cl / vc * central"

  } else {
    central_odeline <- "- (vmax / (km + central / vc)) / vc * central"
  }

  # Add terms for compartments based on no.cmpt
  if (no.cmpt >= 2) {
    central_odeline <-
      paste0(central_odeline, " - q / vc * central + q / vp * peri")
    peri_odeline <- "d/dt(peri) <- q / vc * central - q / vp * peri"
  } else {
    peri_odeline <- NULL
  }

  if (no.cmpt == 3) {
    central_odeline <-
      paste0(central_odeline, " - q2 / vc * central + q2 / vp2 * peri2")
    peri2_odeline <-
      "d/dt(peri2) <- q2 / vc * central - q2 / vp2 * peri2"
  } else {
    peri2_odeline <- NULL
  }

  depot_odeline <- NULL
  transit_line <- NULL
  durline <- NULL
  F1line <- NULL
  Frline <- NULL

  # Add dosing route logic for central compartment
  if (route == "oral" || route == "mixed_iv_oral") {
    if (abs.bio == 1 & abs.type == 4) {
      stop(
        "Error: can not do parallel zero-order and first order absorption model and bioavailalibity estimation simultaneously"
      )
    }

    if (abs.bio == 1 & abs.delay == 3) {
      stop(
        "Error: can not do transit compartment absorption delay and bioavailalibity estimation simultaneously"
      )
    }

    # first-order absorption
    if (abs.type == 1) {
      depot_odeline <- "d/dt(depot) <- - ka * depot"
      central_odeline <- paste0(central_odeline, " + ka * depot")

      if (abs.delay == 1) {
        alagline <-  "alag(depot) <- tlag"
      } else if (abs.delay == 2) {
        transit_line <- "ktr <- (n+1)/mtt"
        depot_odeline <-
          "d/dt(depot) <-  transit(n,mtt,bio) - ka*depot"
      }
    }

    # zero-order absorption
    if (abs.type == 2) {
      durline <-  "dur(central) <- D2"
      if (abs.delay == 1) {
        alagline <-  "alag(central) <- tlag"
      } else if (abs.delay == 2) {
        transit_line <- "ktr <- (n+1)/mtt"
        depot_odeline <- paste0("transit(n,mtt,bio)", central_odeline)
      }
    }

    # sequential zero and first order
    if (abs.type == 3) {
      depot_odeline <- "d/dt(depot) <- - ka * depot"
      durline <-  "dur(depot) <- D2"
      central_odeline <- paste0(central_odeline, " + ka * depot")
      if (abs.delay == 1) {
        alagline <-  "alag(depot) <- tlag"
      } else if (abs.delay == 2) {
        transit_line <- "ktr <- (n+1)/mtt"
        depot_odeline <-
          "d/dt(depot) <-  transit(n,mtt,bio) - ka*depot"
      }
    }

    # dual zero-order and first-order
    if (abs.type == 4) {
      durline <-  "dur(central) <- D2"
      Frline <- c("f(depot) <- Fr", "f(central) <- 1 - Fr")
      depot_odeline <- "d/dt(depot) <- - ka * depot"
      central_odeline <- paste0(central_odeline, " + ka * depot")

      if (abs.delay == 1) {
        alagline <-  "alag(central) <- tlag"
      } else if (abs.delay == 2) {
        transit_line <- "ktr <- (n+1)/mtt"
        central_odeline <-
          paste0("transit(n,mtt,bio)", central_odeline)
      }
    }

    # Currently dont support dual absorption for mixed route
    if (abs.bio == 1) {
      if (abs.type %in% c(1, 3)) {
        F1line <- "f(depot) <-F1"
      }
      if (abs.type == 2) {
        F1line <- "f(central) <-F1"
      }
    }
  }

  # Combine all ODE lines
  central_odeline <- paste0("d/dt(central) <- ", central_odeline)

  dvline <- "cp <- central/vc"
  odelines <-
    c(
      transit_line,
      depot_odeline,
      central_odeline,
      peri_odeline,
      peri2_odeline,
      durline,
      F1line,
      Frline,
      dvline
    )

  # Return all ODE lines as a character vector
  return(odelines)
}

#' Add a covariate effect to a parameter model
#'
#' Automates the creation of covariate effects in pharmacometric models by
#' generating appropriate beta coefficients and modifying model expressions.
#' Supports both standard allometric scaling rules and custom covariate effects.
#'
#' @param param_name Character. Target parameter name (e.g., "cl", "vc").
#' @param covariate_var Character. Covariate variable name (e.g., "WT", "BMI").
#' @param param_model Character. Current parameter model expression (e.g., "cl = exp(tcl)").
#' @param beta_value Numeric. Optional fixed beta value. If NULL, uses built-in rules.
#' @param existing_betas Character vector. Existing beta definitions to append to.
#' @param use_fix Logical. Use `fix()` for beta values? Default TRUE.
#'
#' @return List with two elements:
#' \itemize{
#'   \item betas - Updated character vector of beta definitions
#'   \item mod - Modified model expression with covariate term
#' }
#'
#' @details
#' Automatic beta selection rules:
#' \itemize{
#'   \item Standard covariates ("wt"/"ffm"/"bmi"/"bsa"):
#'     \itemize{
#'       \item 0.75 for clearance parameters (cl/q/q2)
#'       \item 1.0 for volume parameters (vc/vp/vp2)
#'     }
#'   \item Other covariates: Default beta = -0.1 with message
#' }
#'
#' @author Zhonghui Huang
#'
#' @examples
#' # Add weight effect to clearance
#'  add_covariate( "cl", "WT", "cl = exp(tcl)")
#'
#' # Custom beta value for BMI effect
#' add_covariate(
#'   "vc", "BMI", "vc = exp(tvc)",
#'   beta_value = -0.2, use_fix = FALSE
#' )
#' @export
add_covariate <- function(param_name,
                          covariate_var,
                          param_model,
                          beta_value = NULL,
                          existing_betas = c(),
                          use_fix = TRUE) {
  # Determine beta value based on covariate type and parameter category
  if (is.null(beta_value)) {
    predefined_covs <-
      c("wt", "ffm", "bmi", "bsa")  # Supported covariates
    cov_lower <- tolower(covariate_var)

    if (cov_lower %in% predefined_covs) {
      # Apply standard PK allometric scaling rules
      if (param_name %in% c("vmax","cl", "q", "q2")) {
        beta_value <- 0.75  # Clearance parameters
      } else if (param_name %in% c("vc", "vp", "vp2")) {
        beta_value <- 1     # Volume parameters
      } else {
        message(
          paste(
            "No default beta for parameter type",
            param_name,
            "with standard covariate",
            covariate_var,
            "- using 1.0"
          )
        )
        beta_value <- 1
      }
    } else {
      # Default negative effect for non-standard covariates
      beta_value <- -0.1
      message("Using default beta (-0.1) for [",
              covariate_var, "]")
    }
  }

  # Generate standardized covariate names
  l_cov <- toupper(paste0("log", covariate_var))
  beta_name <-
    paste("beta", tolower(covariate_var), param_name, sep = ".")

  # Create beta definition line
  if (use_fix) {
    new_beta <- paste0(beta_name, " <- fix(", beta_value, ")")
  } else {
    new_beta <- paste0(beta_name, " <- ", beta_value)
  }

  # Update beta list and model
  updated_betas <- c(existing_betas, new_beta)
  cov_term <- paste0(" + ", beta_name, "*", l_cov)

  # Pattern matching for parameter model structure
  pattern <- "^\\s*(\\w+)\\s*=\\s*exp\\(\\s*([^)]+)\\s*\\)\\s*$"
  if (grepl(pattern, param_model)) {
    new_model <-
      sub(pattern, paste0("\\1 = exp(\\2", cov_term, ")"), param_model)
    return(list(betas = updated_betas, mod = new_model))
  }

  list(betas = existing_betas, mod = param_model)
}


#' Generate a Pharmacokinetic (PK) Model for nlmixr2
#'
#' Constructs a PK model based on specified parameters, absorption characteristics,
#' variability components, and residual error models. The model is generated as a text file compatible
#' with nlmixr syntax. The function handles various absorption types, multi-compartment models,
#' Michaelis-Menten kinetics, and different residual variability structures.
#'
#' @param modi Model identification number (default: 1). Used for generating unique model filenames.
#' @param route Administration route. Valid options: "bolus", "oral", "mixed_iv_oral" (default: "bolus").
#' @param no.cmpt Number of compartments in the model (1, 2, or 3) (default: 1).
#' @param abs.bio Bioavailability flag (0 = no bioavailability, 1 = with bioavailability) (default: 0).
#' @param abs.type Absorption type (1 = first-order, 2 = zero-order,
#' 3 = sequential first-order and zero-order absorption,
#' 4 = dual first-order and zero-order absorption) (default: 1).
#' @param abs.delay Absorption delay type (0 = none, 1 = lag time, 2 = transit compartments) (default: 0).
#' @param eta.ka Variability flag for absorption rate (ka) (0 = no variability, 1 = include variability).
#' @param eta.cl Variability flag for clearance (CL) (0 = no variability, 1 = include variability).
#' @param eta.vc Variability flag for central volume (Vc) (0 = no variability, 1 = include variability).
#' @param eta.vp Variability flag for peripheral volume (Vp) in multi-compartment models.
#' @param eta.vp2 Variability flag for second peripheral volume (Vp2) in 3-compartment models.
#' @param eta.q Variability flag for intercompartmental clearance (Q) in multi-compartment models.
#' @param eta.q2 Variability flag for second intercompartmental clearance (Q2) in 3-compartment models.
#' @param mm Michaelis-Menten kinetics flag (0 = linear kinetics, 1 = Michaelis-Menten kinetics).
#' @param eta.vmax Variability flag for Vmax when using Michaelis-Menten kinetics.
#' @param eta.km Variability flag for Km when using Michaelis-Menten kinetics.
#' @param eta.tlag Variability flag for lag time (tlag) when abs.delay=1.
#' @param eta.n Variability flag for number of transit compartments when abs.delay=2.
#' @param eta.mtt Variability flag for mean transit time when abs.delay=2.
#' @param eta.bio Variability flag for bioavailability when abs.delay=2.
#' @param eta.D2 Variability flag for zero-order duration (D2) when abs.type=2 or 3.
#' @param eta.F1 Variability flag for bioavailability fraction (F1) when abs.bio=1.
#' @param eta.Fr Variability flag for absorption fraction (Fr) when abs.type=4.
#' @param mcorr Correlation flag for omega blocks (0 = no correlation, 1 = include correlations).
#' @param rv Residual variability type (1 = additive, 2 = proportional, 3 = combined, 4 = log-normal).
#' @param allometric_scaling Allometric scaling type (0 = none, 1 = weight, 2 = BMI, 3 = FFM).
#' @param param_table Data frame containing parameter initial values and variability components.
#' Should contain columns: Name (parameter name), init (initial value),
#' eta (TRUE/FALSE for variability inclusion), cov (covariate relationships).
#' @param return.func Logical, whether to return a compiled function (default `FALSE` returns model code as text).
#' @param out.dir Directory where model files and results are written. Defaults to
#'   the current working directory when not provided.
#' @param verbose Logical; if `TRUE`, progress messages are printed.
#'
#' @return Generates a text file ('modX.txt' where X = modi) containing the nlmixr-compatible model code.
#' The file is written to the current working directory. No explicit return value.
#' If `return.func = TRUE`, returns a compiled model function object.
#'
#' @author Zhonghui Huang
#'
#' @examples
#' withr::with_dir(tempdir(), {
#' #' # Create a 1-compartment oral model with first-order absorption
#'  ppkmodGen( no.cmpt = 1, abs.type = 1,return.func = TRUE,param_table = initialize_param_table())
#' })
#' @export

ppkmodGen<- function(modi=1,
                     route="bolus",
                     no.cmpt=1,
                     abs.bio=0,
                     abs.type=1,
                     abs.delay=0,
                     eta.ka=0,
                     eta.cl=0,
                     eta.vc=0,
                     eta.vp=0,
                     eta.vp2=0,
                     eta.q=0,
                     eta.q2=0,
                     mm=0,
                     eta.vmax=0,
                     eta.km=0,
                     eta.tlag=0,
                     eta.n=0,
                     eta.mtt=0,
                     eta.bio=0,
                     eta.D2=0,
                     eta.F1=0,
                     eta.Fr=0,
                     mcorr=0,
                     rv=1,
                     allometric_scaling=0,
                     param_table=NULL,
                     return.func=FALSE,
                     out.dir = NULL,
                     verbose=TRUE) {

 # Initialize tokens
  ka_init <- ""
  ka_var <- list(mod = "", eta_init = "")

  cl_init <- ""
  cl_var <- list(mod = "", eta_init = "")

  vc_init <- ""
  vc_var <- list(mod = "", eta_init = "")

  vp_init <- ""
  vp_var <- list(mod = "", eta_init = "")

  q_init <- ""
  q_var <- list(mod = "", eta_init = "")

  vp2_init <- ""
  vp2_var <- list(mod = "", eta_init = "")

  q2_init <- ""
  q2_var <- list(mod = "", eta_init = "")

  vmax_init <- ""
  vmax_var <- list(mod = "", eta_init = "")

  km_init <- ""
  km_var <- list(mod = "", eta_init = "")

  tlag_init <- ""
  tlag_var <- list(mod = "", eta_init = "")

  n_init <- ""
  n_var <- list(mod = "", eta_init = "")

  mtt_init <- ""
  mtt_var <- list(mod = "", eta_init = "")

  bio_init <- ""
  bio_var <- list(mod = "", eta_init = "")

  D2_init <- ""
  D2_var <- list(mod = "", eta_init = "") # D2 for zero-order

  F1_init <- ""
  F1_var <- list(mod = "", eta_init = "")

  Fr_init <- ""
  Fr_var <- list(mod = "", eta_init = "")

  betainits <- character(0)
  cov_var  <- character(0)

  if (route == "oral" || route == "mixed_iv_oral") {
    if (abs.type == 1) {
      ka_init <- initialize_param("ka", param_table)
      ka_var <- add_variability("ka", eta.ka, param_table)
    }
    if (abs.type == 2) {
      D2_init <- initialize_param("D2", param_table)
      D2_var <- add_variability("D2", eta.D2, param_table)
    }
    if (abs.type == 3) {
      ka_init <- initialize_param("ka", param_table)
      ka_var <- add_variability("ka", eta.ka, param_table)
      D2_init <- initialize_param("D2", param_table)
      D2_var <- add_variability("D2", eta.D2, param_table)
    }

    if (abs.type == 4) {
      ka_init <- initialize_param("ka", param_table)
      ka_var <- add_variability("ka", eta.ka, param_table)
      D2_init <- initialize_param("D2", param_table)
      D2_var <- add_variability("D2", eta.D2, param_table)
      Fr_init <- initialize_param("Fr", param_table)
      Fr_var <-
        add_variability("Fr", eta.Fr, param_table, param.type = 2)
    }

    if (abs.delay == 1) {
      tlag_init <- initialize_param("tlag", param_table)
      tlag_var <- add_variability("tlag", eta.tlag, param_table)
    }
    if (abs.delay == 2) {
      n_init <- initialize_param("n", param_table)
      n_var <- add_variability("n", eta.n, param_table)
      bio_init <- initialize_param("bio", param_table)
      bio_var <- add_variability("bio", eta.bio, param_table)
      mtt_init <- initialize_param("mtt", param_table)
      mtt_var <- add_variability("mtt", eta.mtt, param_table)
    }
    if (abs.bio == 1 & abs.type != 4) {
      F1_init <- initialize_param("F1", param_table)
      F1_var <-
        add_variability("F1", eta.F1, param_table, param.type = 2)
    }
  }
  cl_init <- initialize_param("cl", param_table)
  cl_var <- add_variability("cl", eta.cl, param_table)

  if (mm == 1) {
    cl_init <- ""
    cl_var <- list(mod = "", eta_init = "")

    vmax_init <- initialize_param("vmax", param_table)
    vmax_var <- add_variability("vmax", eta.vmax, param_table)

    km_init <- initialize_param("km", param_table)
    km_var <- add_variability("km", eta.km, param_table)
  }

  if (no.cmpt==1){
    vc_init <- initialize_param("vc1cmpt", param_table)
    vc_var <- add_variability("vc", eta.vc, param_table)
  }

  if (no.cmpt ==2) {
    vc_init <- initialize_param("vc2cmpt", param_table)
    vc_var <- add_variability("vc", eta.vc, param_table)

    vp_init <- initialize_param("vp2cmpt", param_table)
    vp_var <- add_variability("vp", eta.vp, param_table)

    q_init <- initialize_param("q2cmpt", param_table)
    q_var <- add_variability("q", eta.q, param_table)
  }

  if (no.cmpt == 3) {
    vc_init <- initialize_param("vc3cmpt", param_table)
    vc_var <- add_variability("vc", eta.vc, param_table)

    vp_init <- initialize_param("vp3cmpt", param_table)
    vp_var <- add_variability("vp", eta.vp, param_table)

    vp2_init <- initialize_param("vp2", param_table)
    vp2_var <- add_variability("vp2", eta.vp2, param_table)

    q_init <- initialize_param("q3cmpt", param_table)
    q_var <- add_variability("q", eta.q, param_table)

    q2_init <- initialize_param("q2", param_table)
    q2_var <- add_variability("q2", eta.q2, param_table)
  }

  correlation <- ""
  correlated_params  <- ""
  # Initialize omega blocks and correlated parameters

  omega_block_1 <- ""

  # Define omega_block_1_params
  block1_names <- c("vmax", "km")[c(eta.vmax, eta.km) == 1]
  if (length(block1_names) > 0) {
    omega_block_1_params <- paste0("eta.", block1_names)
  } else {
    omega_block_1_params <- NULL
  }

  omega_block_2 <- ""
  # Define omega_block_2_params
  block2_names <-
    c("cl", "vc", "vp", "q", "vp2", "q2")[c(eta.cl, eta.vc, eta.vp, eta.q, eta.vp2, eta.q2) == 1]
  if (length(block2_names) > 0) {
    omega_block_2_params <- paste0("eta.", block2_names)
  } else {
    omega_block_2_params <- NULL
  }

  # Handle omega_block_1 only if mm == 1
  if (mm == 1) {
    if (length(omega_block_1_params) > 0) {
      if (length(omega_block_1_params) == 1) {
        omega_block_1 <-
          omega_block(param_list = omega_block_1_params,
                      mcorr = 0,
                      eta_table = param_table)
      }
      if (length(omega_block_1_params) == 2) {
        omega_block_1 <-
          omega_block(param_list = omega_block_1_params,
                      mcorr = mcorr,
                      eta_table = param_table)
        if (mcorr == 1) {
          # Combine correlation blocks
          correlation <- c(correlation, omega_block_1)
          correlated_params <- unique(c(omega_block_1_params))
        }
      }
    }
  }

  if (length(omega_block_2_params) > 0) {
    if (length(omega_block_2_params) > 1) {
      omega_block_2 <-
        omega_block(param_list = omega_block_2_params,
                    mcorr = mcorr,
                    eta_table = param_table)
      if (mcorr == 1) {
        # Combine correlation blocks
        correlation <- c(correlation, omega_block_2)
        correlated_params <-
          c(correlated_params , unique(c(omega_block_2_params)))
      }
    }
    if (length(omega_block_2_params) == 1) {
      omega_block_2 <-
        omega_block(param_list = omega_block_2_params,
                    mcorr = 0,
                    eta_table = param_table)
    }
  }

  omega_block_3_params<-NULL
  omega_block_3 <-NULL

  # Identify parameters with eta for parameter in the absorption model
  if (sum(eta.ka,eta.tlag,eta.mtt,eta.n,eta.bio, eta.D2, eta.F1,eta.Fr)>0){
    omega_block_3_params<- paste0(
      "eta.",
      c("ka","tlag","mtt","n","bio", "D2", "F1","Fr")[c(eta.ka,eta.tlag,eta.mtt,eta.n,eta.bio, eta.D2, eta.F1,eta.Fr) == 1]
    )

    omega_block_3 <- omega_block(param_list = omega_block_3_params,mcorr=0, eta_table = param_table)
  }

  # Define a list of all parameter variables
  param_vars <- list(
    ka = ka_var,
    cl = cl_var,
    vc = vc_var,
    vp = vp_var,
    q = q_var,
    vp2 = vp2_var,
    q2 = q2_var,
    vmax = vmax_var,
    km = km_var,
    tlag = tlag_var,
    D2 = D2_var,
    mtt = mtt_var,
    n = n_var,
    bio = bio_var,
    F1 = F1_var,
    Fr = Fr_var
  )

  # Parameter not in the correlation
  exclude_params <- c("ka", "tlag", "mtt", "n", "bio", "D2", "F1", "Fr")

  if (length(correlated_params) > 1) {
    for (param in correlated_params) {
      param_short <- gsub("eta\\.", "", param)
      if (param_short %in% names(param_vars) &&
          !param_short %in% exclude_params) {
        param_vars[[param_short]]$eta_init <- ""
      }
    }
  }

  # Reassign updated variables back to original variables
  ka_var <- param_vars$ka
  cl_var <- param_vars$cl
  vc_var <- param_vars$vc
  vp_var <- param_vars$vp
  q_var <- param_vars$q
  vp2_var <- param_vars$vp2
  q2_var <- param_vars$q2
  vmax_var <- param_vars$vmax
  km_var <- param_vars$km
  tlag_var <- param_vars$tlag
  D2_var <- param_vars$D2
  mtt_var <- param_vars$mtt
  bio_var <- param_vars$bio
  n_var <- param_vars$n
  F1_var <- param_vars$F1
  Fr_var <- param_vars$Fr

 # allometric scaling part
  if (allometric_scaling %in% 1:3) {
    # Determine covariate type based on scaling option
    covar <- switch(
      allometric_scaling,
      "1" = "WT",
      # Allometric scaling 1: Weight
      "2" = "BMI",
      # Allometric scaling 2: BMI
      "3" = "FFM",
      # Allometric scaling 3: Fat-Free Mass
      stop("Invalid allometric_scaling value")
    )

    # Define parameters to process with their activation conditions
    params <- list(
      list(name = "cl",    condition = mm == 0),
      # Clearance (mm=0 condition)
      list(name = "vmax",  condition = mm == 1),
      # Vmax (mm=1 condition)
      list(name = "vc",    condition = TRUE),
      # Central volume (always)
      list(name = "vp",    condition = no.cmpt >= 2),
      # Peripheral volume (2-3 cmpt)
      list(name = "q",     condition = no.cmpt >= 2),
      # Intercomp flow (2-3 cmpt)
      list(name = "vp2",   condition = no.cmpt == 3),
      # 2nd peripheral vol (3 cmpt)
      list(name = "q2",    condition = no.cmpt == 3)  # 2nd intercomp flow (3 cmpt)
    )

    # Process each parameter in sequence
    for (p in params) {
      if (p$condition) {
        # Get corresponding variable object (cl_var, vmax_var, etc.)
        var_name <- paste0(p$name, "_var")
        current_var <- get(var_name)

        # Add covariate relationship using unified interface
        cov_result <- add_covariate(
          param_name = p$name,
          covariate_var = covar,
          param_model = current_var$mod,
          use_fix = TRUE,
          existing_betas = betainits
        )

        # Update model structure and beta initializations
        assign(var_name, `$<-`(current_var, "mod", cov_result$mod))
        betainits <- cov_result$betas
      }
    }
  }

  # Build ode part
  odelines <-
    build_odeline(
      mm = mm,
      no.cmpt = no.cmpt,
      route = route,
      abs.type = abs.type,
      abs.delay = abs.delay,
      abs.bio = abs.bio
    )

  # Residual variability
  res_line <- ""
  sigma_add_init <-""
  sigma_prop_init <-""

  if (rv == 1) {
    sigma_add_init <-
      paste0("sigma_add <-", param_table[param_table$Name == "sigma_add", "init"])
    res_line <- paste0("cp ~ add(sigma_add)")
  } else if (rv == 2) {
    sigma_prop_init <-
      paste0("sigma_prop <-", param_table[param_table$Name == "sigma_prop", "init"])
    res_line <- paste0("cp ~ prop(sigma_prop)")
  } else if (rv == 3) {
    sigma_add_init <-
      paste0("sigma_add <-", param_table[param_table$Name == "sigma_add", "init"])
    sigma_prop_init <-
      paste0("sigma_prop <-", param_table[param_table$Name == "sigma_prop", "init"])
    res_line <- paste0("cp ~ prop(sigma_prop) + add(sigma_add)")
  }
  else if (rv == 4) {
    sigma_add_init <-
      paste0("sigma_add <-", param_table[param_table$Name == "sigma_add", "init"])
    res_line <- paste0("cp ~ lnorm(sigma_add)")
  } else {
    stop("Invalid rv value: ", rv, ". Expected 1-4.")
  }
  # Generate model content
  model_content <- c(
    "f <- function(){",
    "ini({",
    ka_init,
    cl_init,
    vc_init,
    vp_init,
    q_init,
    vp2_init,
    q2_init,
    vmax_init,
    km_init,
    tlag_init,
    D2_init,
    mtt_init,
    n_init,
    bio_init,
    F1_init,
    Fr_init,

    betainits,

    ka_var$eta_init,
    cl_var$eta_init,
    vc_var$eta_init,
    vp_var$eta_init,
    q_var$eta_init,
    vp2_var$eta_init,
    q2_var$eta_init,
    vmax_var$eta_init,
    km_var$eta_init,
    tlag_var$eta_init,
    D2_var$eta_init,
    mtt_var$eta_init,
    n_var$eta_init,
    bio_var$eta_init,
    F1_var$eta_init,
    Fr_var$eta_init,

    correlation,

    sigma_add_init,
    sigma_prop_init,

    "})",
    "model({",

    ka_var$mod,
    cl_var$mod,
    vc_var$mod,
    vp_var$mod,
    q_var$mod,
    vp2_var$mod,
    q2_var$mod,
    vmax_var$mod,
    km_var$mod,
    tlag_var$mod,
    D2_var$mod,
    mtt_var$mod,
    n_var$mod,
    bio_var$mod,
    F1_var$mod,
    Fr_var$mod,

    odelines,
    res_line,
    "})",
    "}"
  )

  model_content <- model_content[model_content != ""]

  if (is.null(out.dir) || !nzchar(out.dir)) {
    out.dir <- getwd()
  }
  # validate user-specified directory
  if (!dir.exists(out.dir)) {
    if (file.exists(out.dir)) {
      stop("out.dir exists but is not a directory: ", out.dir, call. = FALSE)
    } else {
      stop("out.dir does not exist: ", out.dir, call. = FALSE)
    }
  }
  mod_file <- file.path(out.dir, paste0("mod", modi, ".txt"))
  file_conn <- file(mod_file)
  writeLines(model_content, file_conn)
  close(file_conn)

  file_name <- paste0("mod", modi, ".txt")
  if (verbose){
    message(
      "[Success] Model file created:\n",
      normalizePath(mod_file, mustWork = FALSE)
    )
  }
  # Evaluate the function to return it as an object
  if (return.func==T){
  f <- eval(parse(text = paste(model_content, collapse = "\n")))
  return(f)
  }
}

