# # devtools::load_all()
# # library(casebase)
# library(ggplot2)
# data(ERSPC)
# `%ni%` <- Negate('%in%')
# ERSPC$ScrArm <- factor(ERSPC$ScrArm,
#                        levels = c(0, 1),
#                        labels = c("Control group", "Screening group"))
# ERSPC$id <- rev(1:nrow(ERSPC))
#
# dt <- sampleCaseBase(data = ERSPC, event = "DeadOfPrCa", ratio = 1)
# pt_object <- casebase::popTime(ERSPC, event = "DeadOfPrCa")
#
# p1 <- ggplot(pt_object) +
#     geom_segment(color = "grey", aes(x = 0, xend = time, y = ycoord, yend = ycoord)) +
#     theme_minimal()
#
# p1 + geom_point(mapping = aes(x = time, y = yc),
#                 color = "red",
#                 size = 0.8,
#                 data = pt_object[event == 1]) +
#     geom_point(mapping = aes(x = Follow.Up.Time, y = id),
#                color = "blue",
#                size = 0.2,
#                data = dt[dt$DeadOfPrCa==0,])
#
#
#
#
# library(casebase)
# library(tidyverse)
# data("bmtcrr")
# popTimeData <- popTime(data = ERSPC, time = "Follow.Up.Time",
#                        event = "DeadOfPrCa")
# sampleData <- ERSPC %>%
#     arrange(desc(Follow.Up.Time)) %>%
#     mutate(id = row_number()) %>%
#     sampleCaseBase(time = "Follow.Up.Time",
#                    event = "DeadOfPrCa",
#                    comprisk = FALSE,
#                    ratio = 1)
# plot(popTimeData) + geom_point(aes(x = Follow.Up.Time, y = id),
#                                data = filter(sampleData,
#                                              DeadOfPrCa == 0),
#                                colour = "black",
#                                size = 0.5, inherit.aes = FALSE)
#
#
#
# dt[dt$DeadOfPrCa==0,] %>% dim
#
# head(dt)
#
# head(dt)
# base_series <- dt[dt$DeadOfPrCa==0,]
# head(base_series)
# base_series$DeadOfPrCa <- 1
#
# popTime(base_series, event = "DeadOfPrCa")
#
# p1 <- plot(pt_object, legend = TRUE) +
#     scale_colour_manual(values = c("event" = "red"),
#                         labels = c("Death from Prostate Cancer"))
# #
# # dt <- sampleCaseBase(data = ERSPC, event = "DeadOfPrCa", ratio = 10)
# #
# base_series_index <- data.frame(pop = dt,
#                            time = dt[dt$DeadOfPrCa==0,"Follow.Up.Time"])
# #
# # head(base_series_index)
# #
# p1 + geom_point(data = base_series_index, mapping = aes(x = time, y = pop))
# #
# # boxplot(base_series_index$pop)
# #
# # dt %>% head
# #
# # pt_object %>% head
# #
# # which(base_series_index %in% pt_object$yc )
# #
# # which(base_series_index %ni% pt_object$ycoord)
# #
# # table(pt_object$event, pt_object$`event status`)
# # pt_object$`event status` %>% table
# #
# # pt_object[which(tt %in% ycoord)]
# # head(pt_object)
# # which(tt %ni% pt_object$ycoord)
# #
# # pt_object$yc
# #
# #
# # head(dt)
# # ERSPC[42,]
# # head(pt_object)
# # pt_object$n_available
# #
# #     # scale_color_discrete(labels = c("Death from Prostate Cancer"),
# #     #                      color = "blue")
# #
# # head(pt_object)
# #
# #
# # dt <- sampleCaseBase(data = ERSPC, event = "DeadOfPrCa", ratio = 10)
# # dt$DeadOfPrCa <- factor(dt$DeadOfPrCa, levels = 0:1, labels = c("Base series","Case series"))
# # head(dt)
# # # dt <- dt[order(dt$Follow.Up.Time, decreasing = F),]
# # dt$ycoord <- as.numeric(as.character(rownames(dt)))
# # head(dt)
# # plot(dt$ycoord)
# # dt <- dt[order(dt$ycoord, decreasing = T),]
# # dt$ycoord2 <- 1:nrow(dt)
# # # dt$ycoord <- 1:nrow(dt)
# #
# # str(dt)
# # table(dt$DeadOfPrCa)
# # head(pt_object)
# #
# # # p1 <- ggplot(dt, aes(x = 0, xend = Follow.Up.Time, y = ycoord, yend = ycoord))
# #
# # p1 <- ggplot(dt, aes(x = Follow.Up.Time, y = ycoord2, color = DeadOfPrCa))
# #
# # p1 + geom_point() +
# #     # xlab(xlab) +
# #     # ylab(ylab) +
# #     theme_minimal()
# # plot(dt$Follow.Up.Time,1:nrow(dt))
# #
# #
# #
