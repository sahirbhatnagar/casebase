#' currently not being used
#'
#'
#' #' @export
#' #' @rdname stat_poptime
#'
#' StatPoptime <- ggproto("StatPoptime", Stat,
#'                    required_aes = c("time", "event"), ## biomarker, binary outcome
#'                    default_aes = aes(x = ..xstart.., xend = ..xend2..,
#'                                      y = ..ystart.., yend = ..yend2..),
#'
#'                    # setup_data = function(data, params){
#'                    #
#'                    #     varNames <- checkArgsTimeEvent(data = data,
#'                    #                                    time = data$time,
#'                    #                                    event = data$event)
#'                    #
#'                    #     data[,varNames$event] <- verify_d(data[,varNames$event])
#'                    #     data <- data.table::as.data.table(data)
#'                    #     data.table::setnames(data,varNames$time, "time")
#'                    #     data.table::setnames(data,varNames$event, "event")
#'                    #
#'                    #     data
#'                    # },
#'                    compute_group = function(data, scales, na.rm = TRUE){
#'
#'                        if(na.rm){
#'                            data <- subset(data, !is.na(event) & !is.na(time))
#'                        }
#'
#'                        #varNames <- checkArgsTimeEvent(data = data)#,
#'                                                       #time = data$time,
#'                                                       #event = data$event)
#'
#'                        #data$event <- verify_d(data$event)
#'
#'                        nobs <- nrow(data)
#'                        t <- data$time
#'                        e <- data$event
#'
#'                        idx <- order(t)
#'                        ycoord <- rep(NA, nobs)
#'                        for (i in 1:nobs) {
#'                            ycoord[idx[i]] <- (nobs:1)[i]
#'                        }
#'                        for (i in 1:nobs) {
#'                            if (e[i] == 1) {
#'                                yc <- sample(ycoord[t >= t[i]], 1)
#'                                ycoord[ycoord >= yc & ycoord < ycoord[i]] <- ycoord[ycoord >= yc & ycoord < ycoord[i]] + 1
#'                                ycoord[i] <- yc
#'
#'                            }
#'                        }
#'
#'                        # data[,varNames$event] <- verify_d(data[,varNames$event])
#'                        # data <- data.table::as.data.table(data)
#'                        # data.table::setnames(data,varNames$time, "time")
#'                        # data.table::setnames(data,varNames$event, "event")
#'                        #
#'                        # # people with
#'                        # # short values of t at the top
#'                        # nobs <- nrow(data)
#'                        # data[ data[,order(time)], ycoord := (nobs:1)]
#'                        #
#'                        #
#'                        # # sample y coordinates for each event, so that we can
#'                        # # see the incidence density
#'                        # # on population-time plots. Sampling from people who
#'                        # # have an observed time t greater than that of a
#'                        # # fixed individual who had the event
#'                        # data[, yc := 0L]
#'                        #
#'                        # # need to find out how many controls are left to sample from
#'                        # data[, n_available := 0L]
#'                        #
#'                        # data[event == 1,
#'                        #    n_available := sapply(time, function(i) data[time >=i, .N ] )]
#'                        #
#'                        # data[event == 1 & n_available > 0,
#'                        #    yc := sapply(time,
#'                        #                 function(i)
#'                        #                     sample(data[time >= i,ycoord],1) )]
#'
#'                        data.frame(xstart = rep(0,nobs), xend2 = t,
#'                                   ystart = ycoord, yend2 = ycoord)
#'
#'
#'                    })
#'
#' #' Calculate the empirical Receiver Operating Characteristic curve
#' #'
#' #' Given a binary outcome d and continous measurement m, computes the empirical
#' #' ROC curve for assessing the classification accuracy of m
#' #'
#' #' @inheritParams ggplot2::stat_identity
#' #' @param na.rm Remove missing observations
#' #' @section Aesthetics:
#' #' \code{stat_roc} understands the following aesthetics (required aesthetics
#' #' are in bold):
#' #' \itemize{
#' #'   \item \strong{\code{time}} The time variable
#' #'   \item \strong{\code{event}} The binary event variable
#' #'   \item \code{alpha}
#' #'   \item \code{color}
#' #'   \item \code{linetype}
#' #'   \item \code{size}
#' #' }
#' #' @section Computed variables:
#' #' \describe{
#' #'   \item{false_positive_fraction}{estimate of false positive fraction}
#' #'   \item{true_positive_fraction}{estimate of true positive fraction}
#' #'   \item{cutoffs}{values of m at which estimates are calculated}
#' #' }
#' #' @export
#' #' @rdname stat_roc
#' #' @examples
#' #' D.ex <- rbinom(50, 1, .5)
#' #' rocdata <- data.frame(D = c(D.ex, D.ex),
#' #'                    M = c(rnorm(50, mean = D.ex, sd = .4), rnorm(50, mean = D.ex, sd = 1)),
#' #'                    Z = c(rep("A", 50), rep("B", 50)))
#' #'
#' #' ggplot(rocdata, aes(m = M, d = D)) + stat_roc()
#'
#' stat_poptime <- function(mapping = NULL, data = NULL, geom = "poptime",
#'                      position = "identity", show.legend = NA, inherit.aes = TRUE,
#'                      na.rm = TRUE, ...) {
#'     ggplot2::layer(
#'         stat = StatPoptime,
#'         data = data,
#'         mapping = mapping,
#'         geom = geom,
#'         position = position,
#'         show.legend = show.legend,
#'         inherit.aes = inherit.aes,
#'         params = list(na.rm = na.rm, ...)
#'     )
#'
#' }
#'
#'
#' #' @param n.cuts Number of cutpoints to display along each curve
#' #' @param lineend Line end style (round, butt, square)
#' #' @param linejoin Line join style (round, mitre, bevel)
#' #' @param linemitre Line mitre limit (number greater than 1)
#' #' @param arrow Arrow specification, as created by \code{\link[grid]{arrow}}
#' #' @param alpha.line Alpha level for the lines
#' #' @param alpha.point Alpha level for the cutoff points
#' #' @param size.point Size of cutoff points
#' #' @param labels Logical, display cutoff text labels
#' #' @param labelsize Size of cutoff text labels
#' #' @param labelround Integer, number of significant digits to round cutoff labels
#' #' @param na.rm Remove missing values from curve
#' #' @section Computed variables:
#' #' \describe{
#' #'   \item{false_positive_fraction}{estimate of false positive fraction}
#' #'   \item{true_positive_fraction}{estimate of true positive fraction}
#' #'   \item{cutoffs}{values of m at which estimates are calculated}
#' #' }
#' #' @export
#' #' @rdname geom_roc
#' #' @examples
#' #' D.ex <- rbinom(50, 1, .5)
#' #' rocdata <- data.frame(D = c(D.ex, D.ex),
#' #'                    M = c(rnorm(50, mean = D.ex, sd = .4), rnorm(50, mean = D.ex, sd = 1)),
#' #'                    Z = c(rep("A", 50), rep("B", 50)))
#' #'
#' #' ggplot(rocdata, aes(m = M, d = D)) + geom_roc()
#' #' \donttest{
#' #' ggplot(rocdata, aes(m = M, d = D, color = Z)) + geom_roc()
#' #' ggplot(rocdata, aes(m = M, d = D)) + geom_roc() + facet_wrap(~ Z)
#' #' ggplot(rocdata, aes(m = M, d = D)) + geom_roc(n.cuts = 20)
#' #' ggplot(rocdata, aes(m = M, d = D)) + geom_roc(labels = FALSE)
#' #' }
#'
#' GeomPoptime <- ggproto("GeomPoptime", Geom,
#'                        required_aes = c("x", "y", "xend", "yend"),
#'                        non_missing_aes = c("linetype", "size", "shape"),
#'                        default_aes = aes(colour = "grey90", size = 1, linetype = 1,
#'                                          alpha = NA),
#'
#'                        draw_panel = function(data, panel_scales, coord, arrow = NULL,
#'                                              lineend = "butt", na.rm = FALSE) {
#'
#'                            data <- ggplot2::remove_missing(data, na.rm = na.rm,
#'                                                   c("x", "y", "xend", "yend", "linetype", "size", "shape"),
#'                                                   name = "geom_segment")
#'                            if (empty(data)) return(zeroGrob())
#'
#'                            if (coord$is_linear()) {
#'                                coord <- coord$transform(data, panel_scales)
#'                                return(segmentsGrob(coord$x, coord$y, coord$xend, coord$yend,
#'                                                    default.units = "native",
#'                                                    gp = gpar(
#'                                                        col = alpha(coord$colour, coord$alpha),
#'                                                        fill = alpha(coord$colour, coord$alpha),
#'                                                        lwd = coord$size * .pt,
#'                                                        lty = coord$linetype,
#'                                                        lineend = lineend
#'                                                    ),
#'                                                    arrow = arrow
#'                                ))
#'                            }
#'
#'                            data$group <- 1:nrow(data)
#'                            starts <- subset(data, select = c(-xend, -yend))
#'                            ends <- plyr::rename(subset(data, select = c(-x, -y)), c("xend" = "x", "yend" = "y"),
#'                                                 warn_missing = FALSE)
#'
#'                            pieces <- rbind(starts, ends)
#'                            pieces <- pieces[order(pieces$group),]
#'
#'                            GeomPath$draw_panel(pieces, panel_scales, coord, arrow = arrow,
#'                                                lineend = lineend)
#'                        },
#'
#'                        draw_key = draw_key_path
#' )
#'
#' #' Empirical Receiver Operating Characteristic Curve
#' #'
#' #' Display the empirical ROC curve. Useful for characterizing the classification
#' #' accuracy of continuous measurements for predicting binary states
#' #'
#' #' @section Aesthetics:
#' #' \code{geom_roc} understands the following aesthetics (required aesthetics
#' #' are in bold):
#' #' \itemize{
#' #'   \item \strong{\code{x}} The FPF estimate. This is automatically mapped by \link{stat_roc}
#' #'   \item \strong{\code{y}} The TPF estimate. This is automatically mapped by \link{stat_roc}
#' #'   smallest level in sort order is assumed to be 0, with a warning
#' #'   \item \code{alpha}
#' #'   \item \code{color}
#' #'   \item \code{fill}
#' #'   \item \code{linetype}
#' #'   \item \code{size}
#' #' }
#' #'
#' #' @param stat Use to override the default connection between
#' #'   \code{geom_roc} and \code{stat_roc}.
#' #' @seealso See \code{\link{geom_rocci}} for
#' #'   displaying rectangular confidence regions for the empirical ROC curve, \code{\link{style_roc}} for
#' #'   adding guidelines and labels, and \code{\link{direct_label}} for adding direct labels to the
#' #'   curves. Also \link{export_interactive_roc} for creating interactive ROC curve plots for use in a web browser.
#' #' @inheritParams ggplot2::geom_point
#' #' @export
#' #'
#'
#' geom_poptime <- function(mapping = NULL, data = NULL, stat = "poptime",
#'                          na.rm = TRUE, position = "identity", show.legend = NA,
#'                          inherit.aes = TRUE, ...) {
#'     layer(
#'         geom = GeomPoptime, mapping = mapping, data = data, stat = stat,
#'         position = position, show.legend = show.legend, inherit.aes = inherit.aes,
#'         params = list(na.rm = na.rm, ...)
#'     )
#' }
#'
#'
#'
#' head(veteran)
#' names(veteran)
#' dev.off()
#'
#' ggplot(veteran, aes(time = time, event = status)) +
#'     stat_poptime(na.rm = FALSE)
#'
#' veteran$status
#'
