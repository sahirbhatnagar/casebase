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
#' event variables. If \code{time} is missing, the function looks for a column
#' named \code{"time"} in the data. Note that the event variable is inferred
#' from \code{formula}, since it is the left hand side.
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
#' @param censored.indicator a character string of length 1 indicating which
#'   value in \code{event} is the censored. This function will use
#'   \code{\link[stats]{relevel}} to set \code{censored.indicator} as the
#'   reference level. This argument is ignored if the \code{event} variable is a
#'   numeric
#' @param ... Additional parameters passed to \code{\link{sampleCaseBase}}. If
#'   \code{data} inherits from the class \code{cbData}, then these parameters
#'   are ignored.
#' @return An object of \code{glm} and \code{lm} when there is only one event of
#'   interest, or of class \code{\link{CompRisk}}, which inherits from
#'   \code{vglm}, for a competing risk analysis. As such, functions like
#'   \code{summary}, \code{deviance} and \code{coefficients} give familiar
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
#' out_linear <- fitSmoothHazard(event ~ time + z, DT)
#' out_log <- fitSmoothHazard(event ~ log(time) + z, DT)
#' @importMethodsFrom VGAM summary predict
#' @importFrom VGAM vglm multinomial summaryvglm
fitSmoothHazard <- function(formula, data, time, censored.indicator, ...) {
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
                                         comprisk = (length(typeEvents) > 2), ...)
        } else {
            sampleData <- sampleCaseBase(originalData, timeVar, eventVar,
                                         comprisk = (length(typeEvents) > 2),
                                         censored.indicator, ...)
        }
        # if (length(list(...)) != 2) {
        #     message("sampleCaseBase is using some default values; see documentation for more details.")
        # }
    } else {
        originalData <- NULL
        sampleData <- data
    }

    # Update formula to add offset term
    formula <- update(formula, ~ . + offset(offset))

    # Fit a binomial model if there are no competing risks
    if (length(typeEvents) == 2) {
        out <- glm(formula, data = sampleData, family = binomial)
        out$originalData <- originalData
        out$typeEvents <- typeEvents
        out$timeVar <- timeVar
        out$eventVar <- eventVar

    } else {
        # # If we have competing risks, we need to reformat the response
        # multiData_mat <- c()
        # for (type in typeEvents[typeEvents != 0]) {
        #     multiData_mat <- cbind(multiData_mat, as.numeric(sampleData[[eventVar]] == type))
        # }
        # # Base series should correspond to last column
        # multiData_mat <- cbind(multiData_mat, 1- rowSums(multiData_mat))
        # multiData_mat <- as.data.frame(multiData_mat)
        # colnames(multiData_mat) <- paste0("EventType", c(typeEvents[typeEvents != 0], 0))
        #
        # formula <- do.call(update,
        #                    list(formula,
        #                         as.formula(paste(paste0("cbind(",
        #                                                 paste(names(multiData_mat),
        #                                                       collapse = ", "), ")"),
        #                                          "~ ."))))
        #
        # combData <- cbind(sampleData, multiData_mat)
        # model <- VGAM::vglm(formula, family = VGAM::multinomial,
        #                     data = combData)
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
