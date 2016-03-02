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
#' @param link A character string, which gives the specification for the model
#'   link function. Default is the \code{logit} link.
#' @param ... Additional parameters passed to \code{\link{sampleCaseBase}}. If
#'   \code{data} inherits from the class \code{cbData}, then these parameters
#'   are ignored.
#' @return An object of class \code{caseBase}, which inherits from the classes
#'   \code{glm} and \code{lm}. As such, functions like \code{summary} and
#'   \code{coefficients} give familiar results.
#' @export
fitSmoothHazard <- function(formula, data, time, link = "logit", ...) {
    # Infer name of event variable from LHS of formula
    eventVar <- as.character(attr(terms(formula), "variables")[[2]])

    varNames <- checkArgsTimeEvent(data = data, time = time, event = eventVar)
    timeVar <- varNames$time

    typeEvents <- sort(unique(subset(data, select=(names(data) == eventVar), drop = TRUE)))
    # Call sampleCaseBase
    if (!inherits(data, "cbData")) {
        originalData <- as.data.frame(data)
        sampleData <- sampleCaseBase(originalData, timeVar, eventVar,
                                     cmprisk = (length(typeEvents) > 2), ...)
        if (length(list(...)) != 2) {
            warning("sampleCaseBase is using some default values; see documentation for more details.")
        }
    } else {
        originalData <- NULL
        sampleData <- data
    }

    # Update formula to add offset term
    formula <- update(formula, ~ . + offset(offset))

    # Fit a binomial model if there are no competing risks
    if (length(typeEvents) == 2) {
        out <- glm(formula, data = sampleData, family = binomial(link=link))
        out$originalData <- originalData
        out$typeEvents <- typeEvents
        out$timeVar <- timeVar
        out$eventVar <- eventVar

    } else {
        # If we have competing risks, we need to reformat the response
        multiData_mat <- c()
        for (type in typeEvents[typeEvents != 0]) {
            multiData_mat <- cbind(multiData_mat, as.numeric(subset(sampleData, select=(names(sampleData) == eventVar), drop = TRUE) == type))
        }
        # Base series should correspond to last column
        multiData_mat <- cbind(multiData_mat, 1- rowSums(multiData_mat))
        multiData_mat <- as.data.frame(multiData_mat)
        colnames(multiData_mat) <- paste0("EventType", c(typeEvents[typeEvents != 0], 0))

        formula <- do.call(update,
                           list(formula,
                                as.formula(paste(paste0("cbind(",
                                                        paste(names(multiData_mat),
                                                              collapse = ", "), ")"),
                                                 "~ ."))))

        combData <- cbind(sampleData, multiData_mat)
        model <- VGAM::vglm(formula, family = VGAM::multinomial,
                            data = combData)
        # Output of vglm is an S4 object
        # model@originalData <- originalData
        # model@typeEvents <- typeEvents
        out <- new("CompRisk", model,
                   originalData = originalData,
                   typeEvents = typeEvents,
                   timeVar = timeVar,
                   eventVar = eventVar)
        # class(out) <- c("compRisk", class(model))
    }

    return(out)
}
