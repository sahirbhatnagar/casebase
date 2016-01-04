# This is where all functions related to model fitting should appear

#' Fit smooth-in-time parametric hazard functions.
#'
#' @param formula an object of class "formula" (or one that can be coerced to
#'   that class): a symbolic description of the model to be fitted. The details
#'   of model specification are given under ‘Details’.
#' @param data an optional data frame, list or environment (or object coercible
#'   by as.data.frame to a data frame) containing the variables in the model. If
#'   not found in data, the variables are taken from environment(formula),
#'   typically the environment from which glm is called.
#' @param link A character string, giving the specification for the model link
#'   function.
#' @return An object of class \code{caseBase}, which inherits from the classes
#'   \code{glm} and \code{lm}. As such, functions like \code{summary} and
#'   \code{coefficients} give familiar results.
fitSmoothHazard <- function(formula, data, link) {
    # Update formula to add offset term
    formula <- update(formula, ~ . + offset(offset))
    model <- glm(formula, data = data, family = binomial(link=link))

    class(model) <- c("caseBase", class(model))

    return(model)
}
