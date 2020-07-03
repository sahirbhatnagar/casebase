# This is where all functions related to model fitting should appear

#' Fit smooth-in-time parametric hazard functions.
#'
#' Miettinen and Hanley (2009) explained how case-base sampling can be used to
#' estimate smooth-in-time parametric hazard functions. The idea is to sample
#' person-moments, which may or may not correspond to an event, and then fit the
#' hazard using logistic regression.
#'
#' The object \code{data} should either be the output of the function
#' \code{\link{sampleCaseBase}} or the source dataset on which case-base
#' sampling will be performed. In the latter case, it is assumed that
#' \code{data} contains the two columns corresponding to the supplied time and
#' event variables. The variable \code{time} is used for the sampling the base
#' series, and therefore it should represent the time variable on its original
#' (i.e. non transformed) scale. If \code{time} is missing, the function looks
#' for a column named \code{"time"} in the data. Note that the event variable is
#' inferred from \code{formula}, since it is the left hand side.
#'
#' For single-event survival analysis, it is possible to fit the hazard function
#' using \code{glmnet}, \code{gam}, or \code{gbm}. The choice of fitting family
#' is controlled by the parameter \code{family}. The default value is \code{glm},
#' which corresponds to logistic regression. For competing risk analysis, only
#' \code{glm} and \code{glmnet} are allowed.
#'
#' We also provide a matrix interface through \code{fitSmoothHazard.fit}, which
#' mimics \code{glm.fit} and \code{gbm.fit}. This is mostly convenient for
#' \code{family = "glmnet"}, since a formula interface becomes quickly
#' cumbersome as the number of variables increases. In this setting, the matrix
#' \code{y} should have two columns and contain the time and event variables
#' (e.g. like the output of \code{survival::Surv}). We need this linear function
#' of time in order to perform case-base sampling. Therefore, nonlinear
#' functions of time should be specified as a one-sided formula through the
#' argument \code{formula_time} (the left-hand side is always ignored).
#'
#' @param formula an object of class "formula" (or one that can be coerced to
#'   that class): a symbolic description of the model to be fitted. The details
#'   of model specification are given under Details.
#' @param data a data frame, list or environment containing the variables in the
#'   model. If not found in data, the variables are taken from
#'   \code{environment(formula)}, typically the environment from which
#'   \code{fitSmoothHazard} is called.
#' @param time a character string giving the name of the time variable. See
#'   Details.
#' @param family a character string specifying the family of regression models
#'   used to fit the hazard.
#' @param censored.indicator a character string of length 1 indicating which
#'   value in \code{event} is the censored. This function will use
#'   \code{\link[stats]{relevel}} to set \code{censored.indicator} as the
#'   reference level. This argument is ignored if the \code{event} variable is a
#'   numeric.
#' @param ratio integer, giving the ratio of the size of the base series to that
#'   of the case series. Defaults to 100.
#' @param ... Additional parameters passed to fitting functions (e.g.
#'   \code{glm}, \code{glmnet}, \code{gam}).
#' @return An object of \code{glm} and \code{lm} when there is only one event of
#'   interest, or of class \code{\link{CompRisk}}, which inherits from
#'   \code{vglm}, for a competing risk analysis. As such, functions like
#'   \code{summary}, \code{deviance} and \code{coefficients} give familiar
#'   results.
#' @export
#' @examples
#' # Simulate censored survival data for two outcome types from exponential
#' # distributions
#' library(data.table)
#' nobs <- 500
#' tlim <- 20
#'
#' # simulation parameters
#' b1 <- 200
#' b2 <- 50
#'
#' # event type 0-censored, 1-event of interest, 2-competing event
#' # t observed time/endpoint
#' # z is a binary covariate
#' DT <- data.table(z = rbinom(nobs, 1, 0.5))
#' DT[, `:=`(
#'   "t_event" = rweibull(nobs, 1, b1),
#'   "t_comp" = rweibull(nobs, 1, b2)
#' )]
#' DT[, `:=`(
#'   "event" = 1 * (t_event < t_comp) + 2 * (t_event >= t_comp),
#'   "time" = pmin(t_event, t_comp)
#' )]
#' DT[time >= tlim, `:=`("event" = 0, "time" = tlim)]
#'
#' out_linear <- fitSmoothHazard(event ~ time + z, DT, ratio = 10)
#' out_log <- fitSmoothHazard(event ~ log(time) + z, DT, ratio = 10)
#'
#' # Use GAMs
#' library(mgcv)
#' DT[event == 2, event := 1]
#' out_gam <- fitSmoothHazard(event ~ s(time) + z, DT,
#'                            ratio = 10, family = "gam")
#' @importMethodsFrom VGAM summary predict
#' @importFrom VGAM vglm multinomial summaryvglm
#' @importFrom mgcv s te ti t2
fitSmoothHazard <- function(formula, data, time,
                            family = c("glm", "gam", "gbm", "glmnet"),
                            censored.indicator, ratio = 100, ...) {
  family <- match.arg(family)
  if (family == "gbm" && !requireNamespace("gbm", quietly = TRUE)) {
    stop("Pkg gbm needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  if (family == "glmnet" && !requireNamespace("glmnet", quietly = TRUE)) {
    stop("Pkg glmnet needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  formula <- expand_dot_formula(formula, data = data)
  # Infer name of event variable from LHS of formula
  eventVar <- all.vars(formula[[2]])

  if (missing(time)) {
    varNames <- checkArgsTimeEvent(data = data, event = eventVar)
    timeVar <- varNames$time
  } else {
    timeVar <- time
  }

  typeEvents <- sort(unique(data[[eventVar]]))
  # Call sampleCaseBase if class is not cbData
  if (!inherits(data, "cbData")) {
    originalData <- as.data.frame(data)
    if (missing(censored.indicator)) {
      sampleData <- sampleCaseBase(originalData, timeVar, eventVar,
        comprisk = (length(typeEvents) > 2),
        ratio
      )
    } else {
      sampleData <- sampleCaseBase(originalData, timeVar, eventVar,
        comprisk = (length(typeEvents) > 2),
        censored.indicator, ratio
      )
    }
  } else {
    # If class is cbData we no longer have the original data
    originalData <- NULL
    sampleData <- data
  }

  if (family != "glmnet") {
      # Update formula to add offset term
      # glmnet is handled as separate argument
      formula <- update(formula, ~ . + offset(offset))
  }

    # gbm doesn't play nice with interactions and functions of time
    if (family == "gbm") {
        # So warn the user
        if (detect_nonlinear_time(formula, timeVar)) {
          warning(sprintf(paste("You may be using a nonlinear function",
                                "of %s.\ngbm may throw an error."),
                          timeVar),
                  call. = FALSE)
        }
        if (detect_interaction(formula)) {
            warning("gbm may throw an error when using interaction terms",
                    call. = FALSE)
        }
    }

  # Fit a binomial model if there are no competing risks
  if (length(typeEvents) == 2) {
    fittingFunction <- switch(family,
      "glm" = function(formula) glm(formula, data = sampleData,
                                    family = binomial),
      "glmnet" = function(formula) cv.glmnet.formula(formula, sampleData,
                                                     event = eventVar, ...),
      "gam" = function(formula) mgcv::gam(formula, sampleData,
                                          family = "binomial", ...),
      "gbm" = function(formula) gbm::gbm(formula, sampleData,
                                         distribution = "bernoulli", ...)
    )

    out <- fittingFunction(formula)
    out$originalData <- originalData
    out$typeEvents <- typeEvents
    out$timeVar <- timeVar
    out$eventVar <- eventVar
    if (family == "glmnet") out$formula <- formula
    # Reset offset for absolute risk estimation, but keep track of it
    out$offset <- out$data$offset
    out$data$offset <- 0
    # Add new class
    class(out) <- c("singleEventCB", class(out))
  } else {
    # Otherwise fit a multinomial regression
    if (!family %in% c("glm")) {
      stop(sprintf("Competing-risk analysis is not available for family=%s",
                   family),
        .call = FALSE
      )
    }

    fittingFunction <- function(formula) {
      VGAM::vglm(formula,
                 data = sampleData,
                 family = multinomial(refLevel = 1))
      }

    # Turn off warnings from VGAM::vglm.fitter
    withCallingHandlers(model <- fittingFunction(formula),
      warning = handler_fitter
    )

    # The output is an S4 object that extends vglm-class when family='glm'
    # Otherwise it's just an S3 object like above
    out <- switch(family,
      "glm" = new("CompRisk", model,
        originalData = originalData,
        typeEvents = typeEvents,
        timeVar = timeVar,
        eventVar = eventVar
      ),
      "glmnet" = structure(
        c(
          model,
          list(
            "originalData" = originalData,
            "typeEvents" = typeEvents,
            "timeVar" = timeVar,
            "eventVar" = eventVar,
            "formula" = formula
          )
        ),
        class = c("CompRiskGlmnet", class(model))
      )
    )
  }
  return(out)
}

#' @export
#' @rdname fitSmoothHazard
#' @param x Matrix containing covariates.
#' @param y Matrix containing two columns: one corresponding to time, the other
#'   to the event type.
#' @param formula_time A formula describing how the hazard depends on time.
#'   Defaults to linear.
#' @param event a character string giving the name of the event variable.
#' @importFrom stats glm.fit
fitSmoothHazard.fit <- function(x, y, formula_time, time, event,
                                family = c("glm", "gbm", "glmnet"),
                                censored.indicator, ratio = 100, ...) {
  family <- match.arg(family)
  if (family == "gam") stop("The matrix interface is not available for gam")
  if (family == "gbm" && !requireNamespace("gbm", quietly = TRUE)) {
    stop("Pkg gbm needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  if (family == "glmnet" && !requireNamespace("glmnet", quietly = TRUE)) {
    stop("Pkg glmnet needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  # Default to linear term
  if (missing(formula_time)) {
    formula_time <- formula(paste("~", time))
    timeVar <- time
  } else {
    timeVar <- if (length(formula_time) == 3) {
      all.vars(formula_time[[3]])
    } else all.vars(formula_time)
  }
  # There should only be one time variable
  stopifnot(length(timeVar) == 1)

  # Try to infer event from data
  if (missing(event)) {
      varNames <- checkArgsTimeEvent(data = as.data.frame(y), time = timeVar)
      eventVar <- varNames$event
  } else {
      eventVar <- event
  }

    # gbm doesn't play nice with interactions and functions of time
    if (family == "gbm") {
        # So warn the user
        if (detect_nonlinear_time(formula_time, timeVar)) {
            warning(sprintf(paste("You may be using a nonlinear function",
                                  "of %s.\ngbm may throw an error."),
                            timeVar),
                    call. = FALSE)
        }
        if (detect_interaction(formula_time)) {
            warning("gbm may throw an error when using interaction terms",
                    call. = FALSE)
        }
    }

  typeEvents <- sort(unique(y[, eventVar]))
  # Call sampleCaseBase
  originalData <- list(
    "x" = x,
    "y" = y
  )
  class(originalData) <- c(class(originalData), "data.fit")
  if (missing(censored.indicator)) {
    sampleData <- sampleCaseBase(as.data.frame(cbind(y, x)),
      timeVar, eventVar,
      comprisk = (length(typeEvents) > 2),
      ratio
    )
  } else {
    sampleData <- sampleCaseBase(as.data.frame(cbind(y, x)),
      timeVar, eventVar,
      comprisk = (length(typeEvents) > 2),
      censored.indicator, ratio
    )
  }
  # Format everything into matrices and expand variables that need to be
  sample_event <- as.matrix(sampleData[, eventVar])
  sample_time <- if (family %in% c("glmnet", "gbm")) {
    model.matrix(
      update(formula_time, ~ . - 1),
      sampleData
    )
  } else {
    model.matrix(formula_time, sampleData)
  }
  sample_time_x <- cbind(
    sample_time,
    as.matrix(sampleData[, !names(sampleData) %in% c(eventVar, timeVar,
                                                     "offset")])
  )
  sample_offset <- sampleData$offset

  # Fit a binomial model if there are no competing risks
  if (length(typeEvents) == 2) {
    out <- switch(family,
      "glm" = glm.fit(sample_time_x, sample_event,
        family = binomial(),
        offset = sample_offset
      ),
      "glmnet" = cv.glmnet_offset_hack(sample_time_x, sample_event,
        family = "binomial",
        offset = sample_offset, ...
      ),
      "gbm" = gbm::gbm.fit(sample_time_x, sample_event,
        distribution = "bernoulli",
        offset = sample_offset,
        verbose = FALSE, ...
      )
    )

    out$originalData <- originalData
    out$typeEvents <- typeEvents
    out$timeVar <- timeVar
    out$eventVar <- eventVar
    out$matrix.fit <- TRUE
    out$formula_time <- formula_time
    out$offset <- sample_offset
    # Add new class
    class(out) <- c("singleEventCB", class(out))
  } else {
    stop("The matrix interface is not available for competing risks")
  }
  return(out)
}
