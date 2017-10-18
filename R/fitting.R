# This is where all functions related to model fitting should appear

#' Fit smooth-in-time parametric hazard functions.
#'
#' Miettinen and Hanley (2009) explained how case-base sampling can be used to estimate
#' smooth-in-time parametric hazard functions. The idea is to sample person-moments, which may or
#' may not correspond to an event, and then fit the hazard using logistic regression.
#'
#' The object \code{data} should either be the output of the function \code{\link{sampleCaseBase}}
#' or the source dataset on which case-base sampling will be performed. In the latter case, it is
#' assumed that \code{data} contains the two columns corresponding to the supplied time and event
#' variables. If \code{time} is missing, the function looks for a column named \code{"time"} in the
#' data. Note that the event variable is inferred from \code{formula}, since it is the left hand
#' side.
#'
#' For single-event survival analysis, it is possible to fit the hazard function using
#' \code{glmnet}, \code{gam}, or \code{gbm}. The choice of fitting family is controled by the
#' parameter \code{family}. The default value is \code{glm}, which corresponds to logistic
#' regression.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a
#'   symbolic description of the model to be fitted. The details of model specification are given
#'   under Details.
#' @param data a data frame, list or environment containing the variables in the model. If not found
#'   in data, the variables are taken from \code{environment(formula)}, typically the environment
#'   from which \code{fitSmoothHazard} is called.
#' @param time a character string giving the name of the time variable. See Details.
#' @param family a character string specifying the family of regression models used to fit the
#'   hazard.
#' @param censored.indicator a character string of length 1 indicating which value in \code{event}
#'   is the censored. This function will use \code{\link[stats]{relevel}} to set
#'   \code{censored.indicator} as the reference level. This argument is ignored if the \code{event}
#'   variable is a numeric.
#' @param ratio nteger, giving the ratio of the size of the base series to that of the case series.
#'   Defaults to 100.
#' @param ... Additional parameters passed to fitting functions (e.g. \code{glm}, \code{glmnet},
#'   \code{gam}).
#' @return An object of \code{glm} and \code{lm} when there is only one event of interest, or of
#'   class \code{\link{CompRisk}}, which inherits from \code{vglm}, for a competing risk analysis.
#'   As such, functions like \code{summary}, \code{deviance} and \code{coefficients} give familiar
#'   results.
#' @export
#' @examples
#' # Simulate censored survival data for two outcome types from exponential distributions
#' library(data.table)
#' set.seed(12345)
#' nobs <- 5000
#' tlim <- 20
#'
#' # simulation parameters
#' b1 <- 200
#' b2 <- 50
#'
#' # event type 0-censored, 1-event of interest, 2-competing event
#' # t observed time/endpoint
#' # z is a binary covariate
#' DT <- data.table(z=rbinom(nobs, 1, 0.5))
#' DT[,`:=` ("t_event" = rweibull(nobs, 1, b1),
#'           "t_comp" = rweibull(nobs, 1, b2))]
#' DT[,`:=`("event" = 1 * (t_event < t_comp) + 2 * (t_event >= t_comp),
#'          "time" = pmin(t_event, t_comp))]
#' DT[time >= tlim, `:=`("event" = 0, "time" = tlim)]
#'
#' out_linear <- fitSmoothHazard(event ~ time + z, DT, ratio = 10)
#' out_log <- fitSmoothHazard(event ~ log(time) + z, DT, ratio = 10)
#' @importMethodsFrom VGAM summary predict
#' @importFrom VGAM vglm multinomial summaryvglm
#' @importFrom mgcv s te ti t2
fitSmoothHazard <- function(formula, data, time,
                            family = c("glm", "gam", "gbm", "glmnet"),
                            censored.indicator, ratio = 100, ...) {
    family <- match.arg(family)
    if (family == "gbm" && !requireNamespace("gbm", quietly = TRUE)) {
        stop("Pkg gbm needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (family == "glmnet" && !requireNamespace("glmnet", quietly = TRUE)) {
        stop("Pkg glmnet needed for this function to work. Please install it.",
             call. = FALSE)
    }
    # Infer name of event variable from LHS of formula
    eventVar <- as.character(attr(terms(formula), "variables")[[2]])

    if (missing(time)) {
        varNames <- checkArgsTimeEvent(data = data, event = eventVar)
        timeVar <- varNames$time
    } else timeVar <- time


    typeEvents <- sort(unique(data[[eventVar]]))
    # Call sampleCaseBase
    if (!inherits(data, "cbData")) {
        originalData <- as.data.frame(data)
        if (missing(censored.indicator)) {
            sampleData <- sampleCaseBase(originalData, timeVar, eventVar,
                                         comprisk = (length(typeEvents) > 2),
                                         ratio)
        } else {
            sampleData <- sampleCaseBase(originalData, timeVar, eventVar,
                                         comprisk = (length(typeEvents) > 2),
                                         censored.indicator, ratio)
        }
    } else {
        originalData <- NULL
        sampleData <- data
    }

    if (family != "glmnet") {
        # Update formula to add offset term
        formula <- update(formula, ~ . + offset(offset))
    }

    # Fit a binomial model if there are no competing risks
    if (length(typeEvents) == 2) {
        fittingFunction <- switch(family,
                                  "glm" = function(formula) glm(formula, data = sampleData, family = binomial),
                                  "glmnet" = function(formula) cv.glmnet.formula(formula, sampleData, event = eventVar, ...),
                                  "gam" = function(formula) mgcv::gam(formula, sampleData, family = "binomial", ...),
                                  "gbm" = function(formula) gbm::gbm(formula, sampleData, distribution = "bernoulli", ...))

        out <- fittingFunction(formula)
        out$originalData <- originalData
        out$typeEvents <- typeEvents
        out$timeVar <- timeVar
        out$eventVar <- eventVar
        if (family == "glmnet") out$formula <- formula

    } else {
        # Otherwise fit a multinomial regression
        withCallingHandlers(model <- vglm(formula, family = multinomial(refLevel = 1),
                                          data = sampleData),
                            warning = handler_fitter)

        out <- new("CompRisk", model,
                   originalData = originalData,
                   typeEvents = typeEvents,
                   timeVar = timeVar,
                   eventVar = eventVar)
    }
    return(out)
}

#' @export
#' @rdname fitSmoothHazard
#' @param x Matrix containing covariates.
#' @param y Matrix containing two columns: one corresponding to time, the other to the event type.
#' @param event a character string giving the name of the event variable.
#' @importFrom stats glm.fit
fitSmoothHazard.fit <- function(x, y, time, event, family = c("glm", "gbm", "glmnet"),
                                censored.indicator, ratio = 100, ...) {
    family <- match.arg(family)
    if (family == "gam") stop("The matrix interface is not available for gam")
    if (family == "gbm" && !requireNamespace("gbm", quietly = TRUE)) {
        stop("Pkg gbm needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (family == "glmnet" && !requireNamespace("glmnet", quietly = TRUE)) {
        stop("Pkg glmnet needed for this function to work. Please install it.",
             call. = FALSE)
    }
    # Infer name of event variable from LHS of formula
    # eventVar <- as.character(attr(terms(formula), "variables")[[2]])

    eventVar <- event
    if (missing(time)) {
        varNames <- checkArgsTimeEvent(data = as.data.frame(y), event = eventVar)
        timeVar <- varNames$time
    } else timeVar <- time


    typeEvents <- sort(unique(y[,eventVar]))
    # Call sampleCaseBase
    originalData <- cbind(y[, timeVar, drop = FALSE], x)
    if (missing(censored.indicator)) {
        sampleData <- sampleCaseBase(as.data.frame(cbind(y, x)),
                                     timeVar, eventVar,
                                     comprisk = (length(typeEvents) > 2),
                                     ratio)

    } else {
        sampleData <- sampleCaseBase(as.data.frame(cbind(y, x)),
                                     timeVar, eventVar,
                                     comprisk = (length(typeEvents) > 2),
                                     censored.indicator, ratio)
    }
    sample_event <- as.matrix(sampleData[,eventVar])
    sample_time_x <- as.matrix(sampleData[,!names(sampleData) %in% c(eventVar, "offset")])
    sample_offset <- sampleData$offset

    # Fit a binomial model if there are no competing risks
    if (length(typeEvents) == 2) {
        out <- switch(family,
                      "glm" = glm.fit(sample_time_x, sample_event,
                                      family = binomial(),
                                      offset = sample_offset),
                      "glmnet" = glmnet::cv.glmnet(sample_time_x, sample_event,
                                                   family = "binomial",
                                                   offset = sample_offset, ...),
                      "gbm" = gbm::gbm.fit(sample_time_x, sample_event,
                                           distribution = "bernoulli",
                                           offset = sample_offset,
                                           verbose = FALSE, ...))

        out$originalData <- originalData
        out$typeEvents <- typeEvents
        out$timeVar <- timeVar
        out$eventVar <- eventVar
        out$matrix.fit <- TRUE

    } else {
        stop("Not implemented yet")
        # Otherwise fit a multinomial regression
        # withCallingHandlers(model <- vglm(formula, family = multinomial(refLevel = 1),
        #                                   data = sampleData),
        #                     warning = handler_fitter)
        #
        # out <- new("CompRisk", model,
        #            originalData = originalData,
        #            typeEvents = typeEvents,
        #            timeVar = timeVar,
        #            eventVar = eventVar)
    }
    return(out)
}

