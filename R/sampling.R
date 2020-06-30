#' Create case-base dataset for use in fitting parametric hazard functions
#'
#' This function implements the case-base sampling approach described in Hanley
#' and Miettinen (2009). It can be used to fit smooth-in-time parametric
#' functions easily via logistic regression.
#'
#' The base series is sampled using a multinomial scheme: individuals are
#' sampled proportionally to their follow-up time.
#'
#' It is assumed that \code{data} contains the two columns corresponding to the
#' supplied time and event variables. If either the \code{time} or \code{event}
#' argument is missing, the function looks for columns with appropriate-looking
#' names (see \code{\link{checkArgsTimeEvent}}).
#'
#' @section Warning: The offset is calculated using the total follow-up time for
#'   all individuals in the study. Therefore, we need \code{time} to be on the
#'   original scale, not a transformed scale (e.g. logarithmic). Otherwise, the
#'   offset and the estimation will be wrong.
#'
#' @param data a data.frame or data.table containing the source dataset.
#' @param time a character string giving the name of the time variable. See
#'   Details.
#' @param event a character string giving the name of the event variable. See
#'   Details.
#' @param ratio Integer, giving the ratio of the size of the base series to that
#'   of the case series. Defaults to 10.
#' @param comprisk Logical. Indicates whether we have multiple event types and
#'   that we want to consider some of them as competing risks.
#' @param censored.indicator a character string of length 1 indicating which
#'   value in \code{event} is the censored. This function will use
#'   \code{\link[stats]{relevel}} to set \code{censored.indicator} as the
#'   reference level. This argument is ignored if the \code{event} variable is a
#'   numeric
#' @return The function returns a dataset, with the same format as the source
#'   dataset, and where each row corresponds to a person-moment sampled from the
#'   case or the base series.
#' @export
#' @examples
#' # Simulate censored survival data for two outcome types from exponential
#' library(data.table)
#' set.seed(12345)
#' nobs <- 500
#' tlim <- 10
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
#' out <- sampleCaseBase(DT, time = "time", event = "event", comprisk = TRUE)
sampleCaseBase <- function(data, time, event, ratio = 10, comprisk = FALSE,
                           censored.indicator) {
  # Get the variables names for time and event
  varNames <- checkArgsTimeEvent(data = data, time = time, event = event)
  timeVar <- varNames$time
  eventName <- varNames$event
  # Check the event categories
  modifiedEvent <- checkArgsEventIndicator(data, eventName, censored.indicator)
  eventVar <- modifiedEvent$event.numeric
  # Check if we have competing events
  if (!comprisk && modifiedEvent$nLevels > 2) {
    stop(paste(
      "For more than one event type, you should have compRisk=TRUE,",
      "or reformat your data so that there is only one event of interest."
    ),
    call. = FALSE
    )
  }
  # Create survival object from dataset
  if (comprisk) surv_type <- "mstate" else surv_type <- "right"
  survObj <- survival::Surv(data[[timeVar]],
    eventVar,
    type = surv_type
  )

  n <- nrow(survObj) # no. of subjects
  B <- sum(survObj[, "time"]) # total person-time in base
  c <- sum(survObj[, "status"] != 0) # no. of cases (events)
  b <- ratio * c # size of base series
  offset <- log(B / b) # offset so intercept = log(ID | x, t = 0 )

  # We select person-moments from individual proportional
  # to their total follow-up time
  prob_select <- survObj[, "time"] / B
  which_pm <- sample(n, b, replace = TRUE, prob = prob_select)
  bSeries <- as.matrix(survObj[which_pm, ])
  bSeries[, "status"] <- 0
  bSeries[, "time"] <- runif(b) * bSeries[, "time"]

  # Combine base series with covariate data
  selectTimeEvent <- !(colnames(data) %in% c(timeVar, eventName))
  bSeries <- cbind(bSeries,
                   subset(data, select = selectTimeEvent)[which_pm, ,
                                                          drop = FALSE])
  # Rename columns appropriately
  names(bSeries)[names(bSeries) == "status"] <- eventName
  names(bSeries)[names(bSeries) == "time"] <- timeVar

  cSeries <- data[which(subset(data,
                               select = (names(data) == eventName)) != 0), ]

  # Combine case and base series
  cbSeries <- rbind(cSeries, bSeries)
  # Add offset to dataset
  cbSeries <- cbind(cbSeries, rep_len(offset, nrow(cbSeries)))
  names(cbSeries)[ncol(cbSeries)] <- "offset"

  class(cbSeries) <- c("cbData", class(cbSeries))
  return(cbSeries)
}
