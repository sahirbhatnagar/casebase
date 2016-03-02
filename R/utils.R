# This is where all utility functions should appear
# These functions are not exported
expit <- function(x) 1/(1 + exp(-x))
logit <- function(p) log(p)-log(1-p)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# Handling warning messages coming from predictvglm when offset = 0
handler_offset <- function(msg) {
    if( any(grepl("offset", msg))) invokeRestart("muffleWarning")
}

# Check if provided time and event variables are in the dataset
# and also check for any good substitute
checkArgsTimeEvent <- function(data, time, event) {

    if (missing(time)) {
        if (any(grepl("[\\s\\W_]+time|^time\\b", names(data),
                      ignore.case = TRUE, perl = TRUE))) {
            time <- grep("[\\s\\W_]+time|^time\\b", names(data),
                         ignore.case = TRUE, value = TRUE, perl = TRUE)
            if (length(time)>1)
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
                      names(data)[-which(colnames(data)==time[1])],
                      ignore.case = TRUE, perl = TRUE))) {
            event <- grep("[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b",
                          names(data)[-which(colnames(data)==time[1])],
                          ignore.case = TRUE, value = TRUE, perl = TRUE)
            if (length(event)>1)
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
        stop("data does not contain supplied time and event variables")
    }

    return(list(time = time[1], event = event[1]))
}

# Fill in a templated function with default parameter values
# This is pryr::partial almost verbatim
partialize <- function (`_f`, ...) {
    stopifnot(is.function(`_f`))
    fcall <- as.call(c(substitute(`_f`), list(...)))
    fcall[[length(fcall) + 1]] <- quote(...)
    args <- as.pairlist(list(... = quote(expr = )))
    stopifnot(is.language(fcall))
    eval(call("function", args, fcall), parent.frame())
}
