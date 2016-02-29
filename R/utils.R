# This is where all utility functions should appear
# These functions are not exported
expit <- function(x) 1/(1 + exp(-x))
logit <- function(p) log(p)-log(1-p)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# Handling warning messages coming from predictvglm when offset = 0
handler_offset <- function(msg) {
    if( any(grepl("offset", msg))) invokeRestart("muffleWarning")
}
