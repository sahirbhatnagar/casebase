# This is where all functions related to model fitting should appear

#' Fit smooth-in-time parametric hazard functions.
#'
#' Miettinen and Hanley (2009) explained how case-base sampling can be used to
#' estimate smooth-in-time parametric hazard functions. The idea is to sample
#' person-moments, which may or may not correspond to an event, and then fit the
#' hazard using logistic regression.
#'
#' The object \code{data} should either be the output of the function
#' \link{\code{sampleCaseBase}} or the source dataset on which case-base
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
    event <- as.character(attr(terms(formula), "variables")[[2]])

    if (missing(time)) {
        if (any(grepl("^time", names(data), ignore.case = TRUE))) {
            time <- grep("^time", names(data), ignore.case = TRUE, value = TRUE)
        } else {
            stop("data does not contain time variable")
        }
    }
    # Call sampleCaseBase
    if (!inherits(data, "cbData")) {
        originalData <- as.data.frame(data)
        data <- sampleCaseBase(originalData, time, event, ...)
        if (length(list(...)) != 2) {
            warning("sampleCaseBase is using some default values; see documentation for more details.")
        }
    } else {
        originalData <- NULL
    }

    # Update formula to add offset term
    formula <- update(formula, ~ . + offset(offset))

    # Fit a binomial model is there are no competing risks
    typeEvents <- unique(subset(data, select=(names(data) == event), drop = TRUE))
    if (length(typeEvents) == 2) {
        model <- glm(formula, data = data, family = binomial(link=link))
        model$originalData <- originalData
    } else {
        # If we have competing risks, we need to reformat the response
        multiData_mat <- c()
        for (type in typeEvents[typeEvents != 0]) {
            multiData_mat <- cbind(multiData_mat, as.numeric(subset(data, select=(names(data) == event), drop = TRUE) == type))
        }
        # Base series should correspond to last column
        multiData_mat <- cbind(multiData_mat, 1- rowSums(multiData_mat))

        formula <- update(formula, multiData ~ .)

        model <- vglm(formula, data = data, family = multinomial)
        model$originalData <- originalData
    }

    class(model) <- c("caseBase", class(model))

    return(model)
}
