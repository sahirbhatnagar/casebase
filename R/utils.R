# This is where all utility functions should appear
# These functions are not exported

`%ni%` <- Negate("%in%")

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
      ignore.case = TRUE, perl = TRUE
    ))) {
      time <- grep("[\\s\\W_]+time|^time\\b", names(data),
        ignore.case = TRUE, value = TRUE, perl = TRUE
      )
      if (length(time) > 1) {
        warning(paste0(
          "The following variables for time were found in the data: ",
          paste0(time, collapse = ", "), ". '", time[1],
          "' will be used as the time variable"
        ))
      } else {
        message(paste0(
          "'", time, "'",
          " will be used as the time variable"
        ))
      }
    } else {
      stop("data does not contain time variable")
    }
  }

  if (missing(event)) {
    if (any(grepl("[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b",
      names(data)[-which(colnames(data) == time[1])],
      ignore.case = TRUE, perl = TRUE
    ))) {
      event <- grep("[\\s\\W_]+event|^event\\b|[\\s\\W_]+status|^status\\b",
        names(data)[-which(colnames(data) == time[1])],
        ignore.case = TRUE, value = TRUE, perl = TRUE
      )
      if (length(event) > 1) {
        warning(paste0(
          "The following variables for event were found in the data: ",
          paste0(event, collapse = ", "), ". '", event[1],
          "' will be used as the event variable"
        ))
      } else {
        message(paste0(
          "'", event, "'",
          " will be used as the event variable"
        ))
      }
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
#' if (requireNamespace("survival", quietly = TRUE)) {
#' library(survival) # for veteran data
#' checkArgsEventIndicator(data = veteran, event = "celltype",
#'                         censored.indicator = "smallcell")
#' checkArgsEventIndicator(data = veteran, event = "status")
#' }
#' data("bmtcrr") # from casebase
#' checkArgsEventIndicator(data = bmtcrr, event = "Sex",
#'                         censored.indicator = "M")
#' checkArgsEventIndicator(data = bmtcrr, event = "D",
#'                         censored.indicator = "AML")
#' checkArgsEventIndicator(data = bmtcrr, event = "Status")
checkArgsEventIndicator <- function(data, event, censored.indicator) {
  isFactor <- is.factor(data[[event]])
  isNumeric <- is.numeric(data[[event]])
  isCharacter <- is.character(data[[event]])

  if (!any(isFactor, isNumeric, isCharacter)) {
    stop(strwrap("event variable must be either a factor,
                     numeric or character variable", width = 60))
  }

  nLevels <- nlevels(factor(data[[event]]))
  if (nLevels < 2) stop(paste("event variable must have",
                              "at least two unique values"))

  if (missing(censored.indicator) || is.null(censored.indicator)) {
    if (isFactor) {
      slev <- levels(data[[event]])
      warning(paste0(
        "censor.indicator not specified. assuming ",
        slev[1], " represents a censored observation and ",
        slev[2], " is the event of interest"
      ))
      event.factored <- data[[event]]
    }

    if (isCharacter) {
      event.factored <- factor(data[[event]])
      slev <- levels(event.factored)
      warning(paste0(
        "censor.indicator not specified. assuming ",
        slev[1], " represents a censored observation and ",
        slev[2], " is the event of interest"
      ))
    }

    if (isNumeric) {
      slev <- sort(unique(data[[event]]))
      if (!any(slev %in% 0)) stop(paste("event is a numeric variable that",
                                        "doesn't contain 0. if event is a",
                                        "numericit must contain some 0's",
                                        "to indicate censored observations"))
      event.factored <- if (nLevels == 2) {
        factor(data[[event]],
          labels = c("censored", "event")
        )
      } else {
        factor(data[[event]],
          labels = c(
            "censored", "event",
            paste0(
              "competing event",
              if (nLevels >= 4) 1:(nLevels - 2)
            )
          )
        )
      }
    }
  } else {
    if (!(censored.indicator %in% data[[event]]) & any(isCharacter, isFactor)) {
      stop(strwrap("censored.indicator not found in event variable of data"))
    }

    if (isNumeric) {
      warning(strwrap("censored.indicator specified but ignored because
                                event is a numeric variable"))
      slev <- sort(unique(data[[event]]))
      if (!any(slev %in% 0)) stop(strwrap("event is a numeric variable that
                                        doesn't contain 0. if event is a numeric
                                        it must contain some 0's
                                        to indicate censored observations"))
      event.factored <- if (nLevels == 2) {
        factor(data[[event]],
          labels = c("censored", "event")
        )
      } else {
        factor(data[[event]],
          labels = c(
            "censored", "event",
            paste0(
              "competing event",
              if (nLevels >= 4) 1:(nLevels - 2)
            )
          )
        )
      }
    }

    if (isFactor | isCharacter) {
      event.factored <- relevel(factor(data[[event]]), censored.indicator)
      slev <- levels(event.factored)
      message(paste0(
        "assuming ",
        slev[1], " represents a censored observation and ",
        slev[2], " is the event of interest"
      ))
    }
  }

  return(list(
    event.factored = event.factored,
    event.numeric = as.numeric((event.factored)) - 1,
    nLevels = nLevels
  ))
}

# Remove offset from formula
# https://stackoverflow.com/a/40313732/2836971


#' @importFrom stats model.matrix
#' @importFrom stats contrasts
#' @details `prepareX` is a slightly modified version of the same function from
#'   the `glmnet` package. It can be used to convert a data.frame to a matrix
#'   with categorical variables converted to dummy variables using one-hot
#'   encoding
#' @rdname fitSmoothHazard
#' @export
prepareX <- function(formula, data) {
  whichfac <- sapply(data, inherits, "factor")
  ctr <- if (any(whichfac)) {
    lapply(subset(data, select = whichfac),
           contrasts, contrast = FALSE)
  } else NULL
  X <- model.matrix(update(formula, ~ . - 1), data = data, contrasts.arg = ctr)
  if (any(whichfac))
    attr(X, "contrasts") <- NULL
  attr(X, "assign") <- NULL
  X
}

cv.glmnet.formula <- function(formula, data, event,
                              competingRisk = FALSE, ...) {
  X <- prepareX(formula, data)
  Y <- data[, event]
  if (competingRisk) {
    fam <- "multinomial"
    offset <- NULL
  } else {
    fam <- "binomial"
    offset <- data[, "offset"]
  }
  cv.glmnet_offset_hack(X, Y, offset = offset, family = fam,
                        type.multinomial = "grouped", ...)
}

cv.glmnet_offset_hack <- function(x, y, offset, ...) {
  # For some values of the offset, cv.glmnet does not converge
  # For constant offset, we can use the hack below
  if (diff(range(offset)) > 1e-06) {
    stop("Glmnet is only available with constant offset",
      call. = FALSE
    )
  }

  offset_value <- unique(offset)[1]
  # 1. Fit without offset
  out <- glmnet::cv.glmnet(x, y, ...)
  # 2. Fix the intercept
  out$glmnet.fit$a0 <- out$glmnet.fit$a0 - offset_value

  return(out)
}

# Montecarlo Integration
# Mimic the interface of integrate
integrate_mc <- function(f, lower, upper, ..., subdivisions = 100L) {
  sampledPoints <- runif(subdivisions,
    min = lower,
    max = upper
  )
  return((upper - lower) * mean(f(sampledPoints, ...)))
}

# Taken from brms package
expand_dot_formula <- function(formula, data = NULL) {
  if (isTRUE("." %in% all.vars(formula))) {
    att <- attributes(formula)
    try_terms <- try(
      stats::terms(formula, data = data),
      silent = TRUE
    )
    if (!is(try_terms, "try-error")) {
      formula <- formula(try_terms)
    }
    attributes(formula) <- att
  }
  formula
}

# Streamlined version of pracma::cumtrapz
trap_int <- function(x, y) {
  x <- as.matrix(c(x))
  m <- length(x)
  y <- as.matrix(y)
  n <- ncol(y)
  dt <- kronecker(matrix(1, 1, n), 0.5 * diff(x))
  ct <- apply(dt * (y[1:(m - 1), ] + y[2:m, ]), 2, cumsum)
  return(rbind(0, ct))
}

# Detect if formula contains a function of time or interaction----
count_matches <- function(pat, vec) sapply(regmatches(vec, gregexpr(pat, vec)),
                                           length)

balance_parentheses <- function(str) {
  num_left <- count_matches("\\(", str)
  num_right <- count_matches("\\)", str)

  str[num_left > num_right] <- sub("\\(", "", str[num_left > num_right])
  str[num_left < num_right] <- sub("\\)", "", str[num_left < num_right])

  return(str)
}

detect_nonlinear_time <- function(formula, timeVar) {
  # Two regular expressions
  # 1. Find function arguments
  pattern_args <- "\\(\\s*([^)]+?)\\s*\\)"
  # 2. Find exactly time as the clean string
  time_regex <- paste0("^", timeVar, "$")
  # Extract variables in RHS of formula
  terms <- attr(terms(formula), "term.labels")
  # Then extract the arguments of any function
  matches <- regmatches(terms, regexpr(pattern_args, terms))
  # Next, detect time within nested calls
  matches <- balance_parentheses(matches)
  while (any(matches != regmatches(matches, regexpr(pattern_args, matches)))) {
    matches <- regmatches(matches, regexpr(pattern_args, matches))
    matches <- balance_parentheses(matches)
  }
  # Check if one of these arguments is timeVar
  contain_time <- lapply(
    strsplit(matches, ","),
    function(str) {
      clean_str <- gsub(
        ".*=", "", # Remove equal signs if they exist
        gsub("(\\(\\s*|\\s*\\))", "", str)
      ) # Remove parentheses
      any(grepl(time_regex, trimws(clean_str)))
    }
  )
  any(unlist(contain_time))
}

detect_interaction <- function(formula) {
  # Extract the order of the terms
  orders <- attr(terms(formula), "order")
  # Check if terms of order > 1
  any(orders > 1)
}

# Get typical covariate profile from dataset
#' @importFrom stats median
get_typical <- function(data) {
  data.frame(lapply(data, function(col) {
    if (is.numeric(col) || inherits(col, "Date")) {
      # For numeric or dates, take median
      median(col, na.rm = TRUE)
    } else {
      # If character string or factor, take most common value
      mode <- names(sort(-table(col)))[1]
      factor(mode, levels = levels(factor(col)))
    }
  }))
}

#' @rdname plot.singleEventCB
incrVar <- function(var, increment = 1) {
  n <- length(var)
  if (n > 1 && length(increment) == 1) {
    increment <- rep(increment, n)
  }
  function(data) {
    for (i in 1:n) {
      if (is.factor(data[[var[i]]])) {
        data[[var[i]]] <- fct_shift_ord(data[[var[i]]],
                                        increment = increment[i])
      } else {
        data[[var[i]]] <- data[[var[i]]] + increment[i]
      }
    }
    data
  }
}


fct_shift_ord <- function(x, increment = 1, cap = TRUE, .fun = `+`) {
  x_nlevel <- nlevels(x)
  x_lables <- levels(x)

  # apply function .fun to the numeric of the ordered vector
  erg <- .fun(as.numeric(x), increment)

  # cap to 1 and x_nlevel if the increment was larger
  # than the original range of the factor levels
  if (cap) {
    erg[erg < 1] <- 1
    erg[erg > x_nlevel] <- x_nlevel
  }
  ordered(erg, levels = 1:x_nlevel, labels = x_lables)
}


hrJacobian <- function(object, newdata, newdata2, term) {

  # Set offset to zero
  newdata$offset <- 0
  newdata2$offset <- 0


  m1 <- stats::model.frame(term,
    data = newdata2,
    na.action = stats::na.pass,
    xlev = object$xlevels
  )
  m0 <- stats::model.frame(term,
    data = newdata,
    na.action = stats::na.pass,
    xlev = object$xlevels
  )

  X1 <- stats::model.matrix(term, m1, contrasts.arg = object$contrasts)
  X0 <- stats::model.matrix(term, m0, contrasts.arg = object$contrasts)

  # this is the jacobian!!
  X1 - X0
}
