# This is where all utility functions should appear
# These functions are not exported
expit <- function(x) 1/(1 + exp(-x))
logit <- function(p) log(p) - log(1 - p)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
roundUp <- function(x) 10^ceiling(log10(x))

# Handling warning messages coming from predictvglm when offset = 0
handler_offset <- function(msg) {
    if (any(grepl("offset", msg))) invokeRestart("muffleWarning")
}
# Handling warning messages coming from predictvglm when using b-splines
handler_bsplines <- function(msg) {
    if (any(grepl("ill-conditioned bases", msg))) invokeRestart("muffleWarning")
}
# Handling warning messages coming from vglm.fitter
handler_fitter <- function(msg) {
    if (any(grepl("vglm.fitter", msg))) invokeRestart("muffleWarning")
}

# Check if provided time and event variables are in the dataset
# and also check for any good substitute
#' @rdname popTime
#' @export
checkArgsTimeEvent <- function(data, time, event) {

    if (missing(time)) {
        if (any(grepl("[\\s\\W_]+time|^time\\b", names(data),
                      ignore.case = TRUE, perl = TRUE))) {
            time <- grep("[\\s\\W_]+time|^time\\b", names(data),
                         ignore.case = TRUE, value = TRUE, perl = TRUE)
            if (length(time) > 1)
                warning(paste0("The following variables for time were found in
                              the data: ",paste0(time, collapse = ", "),". '", time[1],
                               "' will be used as the time variable" )) else
                                   message(paste0("'",time,"'",
                                                  " will be used as the time variable"))
        } else {
            stop("data does not contain time variable")
        }
    }

    if (missing(event)) {
        if (any(grepl("[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b",
                      names(data)[-which(colnames(data) == time[1])],
                      ignore.case = TRUE, perl = TRUE))) {
            event <- grep("[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b",
                          names(data)[-which(colnames(data) == time[1])],
                          ignore.case = TRUE, value = TRUE, perl = TRUE)
            if (length(event) > 1)
                warning(paste0("The following variables for event were found in
                              the data: ",paste0(event, collapse = ", "),". '", event[1],
                               "' will be used as the event variable" )) else
                                   message(paste0("'",event,"'",
                                                  " will be used as the event variable"))
        } else {
            stop("data does not contain event or status variable")
        }
    }

    if (!all(c(time, event) %in% colnames(data))) {
        stop("data does not contain supplied time and/or event variables")
    }

    return(list(time = time[1], event = event[1]))
}


#' Check that Event is in Correct Format
#'
#' Checks for event categories and gives a warning message indicating which
#' level is assumed to be the reference level.
#'
#' @inheritParams popTime
#' @return A list of length two. The first element is the factored event, and
#'   the second element is the numeric representation of the event
#'
#' @export
#' @examples
#'
#' \dontrun{
#' library(survival) # for veteran data
#' checkArgsEventIndicator(data = veteran, event = "celltype", censored.indicator = "smallcell")
#' checkArgsEventIndicator(data = veteran, event = "status")
#' checkArgsEventIndicator(data = veteran, event = "trt") # returns error
#'
#' url <- "https://raw.githubusercontent.com/sahirbhatnagar/casebase/master/inst/extdata/bmtcrr.csv"
#' bmt <- read.csv(url)
#' checkArgsEventIndicator(data = bmt, event = "Sex", censored.indicator = "M")
#' checkArgsEventIndicator(data = bmt, event = "D", censored.indicator = "AML")
#' checkArgsEventIndicator(data = bmt, event = "D", censored.indicator = "AMLL") #returns error
#' checkArgsEventIndicator(data = bmt, event = "Source")
#' checkArgsEventIndicator(data = bmt, event = "Status")
#' checkArgsEventIndicator(data = bmt, event = "Status", censored.indicator = 3)
#' }
#'
checkArgsEventIndicator <- function(data, event, censored.indicator) {

    isFactor <- is.factor(data[,event])
    isNumeric <- is.numeric(data[,event])
    isCharacter <- is.character(data[,event])

    if (!any(isFactor, isNumeric, isCharacter))
        stop(strwrap("event variable must be either a factor,
                     numeric or character variable", width = 60))

    nLevels <- nlevels(factor(data[,event]))
    if (nLevels < 2) stop(strwrap("event variable must have at least two unique values"))

    if (missing(censored.indicator) || is.null(censored.indicator)) {

        if (isFactor) {
            slev <- levels(data[,event])
            warning(paste0("censor.indicator not specified. assuming ",
                           slev[1], " represents a censored observation and ",
                           slev[2], " is the event of interest"))
            event.factored <- data[,event]
        }

        if (isCharacter) {
            event.factored <- factor(data[,event])
            slev <- levels(event.factored)
            warning(paste0("censor.indicator not specified. assuming ",
                           slev[1], " represents a censored observation and ",
                           slev[2], " is the event of interest"))
        }

        if (isNumeric) {

            slev <- sort(unique(data[,event]))
            if (!any(slev %in% 0)) stop(strwrap("event is a numeric variable that
                                                doesn't contain 0. if event is a numeric
                                                it must contain some 0's
                                                to indicate censored observations"))
            event.factored <- if (nLevels == 2) factor(data[,event],
                                                       labels = c("censored","event")) else
                                                           factor(data[,event],
                                                                  labels = c("censored","event",
                                                                             paste0("competing event",
                                                                                    if (nLevels >= 4) 1:(nLevels-2))))
        }

    } else {

        if (!(censored.indicator %in% data[,event]) & any(isCharacter, isFactor))
            stop(strwrap("censored.indicator not found in event variable of data"))

        if (isNumeric) {
            warning(strwrap("censored.indicator specified but ignored because
                                event is a numeric variable"))
            slev <- sort(unique(data[,event]))
            if (!any(slev %in% 0)) stop(strwrap("event is a numeric variable that
                                        doesn't contain 0. if event is a numeric
                                        it must contain some 0's
                                        to indicate censored observations"))
            event.factored <- if (nLevels == 2) factor(data[,event],
                                                       labels = c("censored","event")) else
                                                           factor(data[,event],
                                                                  labels = c("censored","event",
                                                                             paste0("competing event",
                                                                                    if (nLevels >= 4) 1:(nLevels-2))))

        }

        if (isFactor | isCharacter) {

            event.factored <- relevel(factor(data[,event]), censored.indicator)
            slev <- levels(event.factored)
            message(paste0("assuming ",
                           slev[1], " represents a censored observation and ",
                           slev[2], " is the event of interest"))
        }
    }

    return(list(event.factored = event.factored,
                event.numeric = as.numeric((event.factored)) - 1,
                nLevels = nLevels))

}


# Fill in a templated function with default parameter values
# This is pryr::partial almost verbatim
partialize <- function(`_f`, ...) {
    stopifnot(is.function(`_f`))
    fcall <- as.call(c(substitute(`_f`), list(...)))
    fcall[[length(fcall) + 1]] <- quote(...)
    args <- as.pairlist(list(... = quote(expr = )))
    stopifnot(is.language(fcall))
    eval(call("function", args, fcall), parent.frame())
}
