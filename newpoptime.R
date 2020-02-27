library(casebase)
library(ggplot2)

ERSPC$DeadOfPrCa <- factor(ERSPC$DeadOfPrCa,
                           levels = c(0, 1),
                           labels = c("Censored", "PrCa Death"))
popTimeData <- popTime(data = ERSPC,
                       time = "Follow.Up.Time",
                       event = "DeadOfPrCa")
head(popTimeData)
p1 <- ggplot()

# doesnt work
nogo <- list(data = popTimeData[event == 1], mapping = aes(x = time, y = yc, colour = `event status`), colour = "green")
p1 + do.call("geom_point", nogo)

# works
works <- list(data = popTimeData[event == 1], mapping = aes(x = time, y = yc, colour = `event status`))
p1 + do.call("geom_point", works) + scale_color_manual(values = c("PrCa Death" = "blue"))


dim(popTimeData)
sampleCaseBase(popTimeData)
sampleData <- ERSPC %>%
    arrange(desc(Follow.Up.Time)) %>%
    mutate(id = row_number()) %>%
    sampleCaseBase(time = "Follow.Up.Time",
                   event = "DeadOfPrCa",
                   comprisk = FALSE,
                   ratio = 1)
plot(popTimeData) + geom_point(aes(x = Follow.Up.Time, y = id),
                               data = filter(sampleData,
                                             DeadOfPrCa == 0),
                               colour = "black",
                               size = 0.5, inherit.aes = FALSE)


huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
h <- ggplot(huron, aes(year))

h + geom_ribbon(aes(ymin=0, ymax=level))


popTimeData %>% head

huron %>% head


library(casebase)
library(ggplot2)


devtools::load_all()
ERSPC$ScrArm <- factor(ERSPC$ScrArm,
                       levels = c(0, 1),
                       labels = c("Control group", "Screening group"))
ERSPC$DeadOfPrCa <- factor(ERSPC$DeadOfPrCa,
                       levels = c(0, 1),
                       labels = c("Censored", "PrCa Death"))
popTimeData <- popTime(data = ERSPC,
                       time = "Follow.Up.Time",
                       event = "DeadOfPrCa")

plot(popTimeData)

devtools::load_all()
cols <- c("Case series" = "#D55E00", "Competing event" = "#009E73", "Base series" = "#0072B2")
plot(popTimeData,
     casebase.theme = TRUE,
     add.case.series = TRUE,
     ratio = 1,
     add.base.series = TRUE,
     legend = TRUE,
     base.params = list(show.legend = TRUE))#,
     # legend.params = list(name = element_blank(),
     #                 breaks = c("Case series", "Competing event", "Base series"),
     #                 values = cols))
     # case.params = list(mapping = aes(x = time, y = yc, colour = "Relapse")),
     # base.params = list(mapping = aes(x = time, y = ycoord, colour = "controls")),
     # legend.params = list(breaks = c("Relapse", "controls"), values = c("Relapse" = "green","controls" = "red")))


popTimeData <- popTime(data = bmtcrr, time = "ftime", event = "Status")
cols <- c("Case series" = "black", "Competing event" = "#009E73", "Base series" = "#0072B2")
p1 <- plot(popTimeData,
     casebase.theme = TRUE,
     add.case.series = TRUE,
     ratio = 1,
     add.base.series = TRUE,
     legend = TRUE,
     comprisk = TRUE,
     add.competing.event = FALSE,
     legend.params = list(name = element_blank(),
                     breaks = c("Case series", "Competing event", "Base series"),
                     values = cols),
     theme.params = list(legend.position = "none"))

dev.off()
p1

p1 + labs(caption = c("hellow worls"), title = "Poptime", x = "nine")
     # case.params = list(mapping = aes(x = time, y = yc, colour = "Relapse")),
     # base.params = list(mapping = aes(x = time, y = ycoord, colour = "controls")),
     # legend.params = list(breaks = c("Relapse", "controls"), values = c("Relapse" = "green","controls" = "red")))


h <- ggplot(popTimeData)


devtools::load_all()
# library(casebase)
popTimeData <- popTime(data = bmtcrr, time = "ftime", event = "Status")
plot(popTimeData,
     add.case.series = TRUE,
     add.base.series = TRUE,
     add.competing.event = FALSE,
     comprisk = TRUE,
     # case.params = list(show.legend = TRUE),
     ratio = 1,
     legend = T,
     casebase.theme = TRUE,
     theme.params = list(legend.position = "bottom"))



bt <- data.table::as.data.table(bmtcrr)
popTimeData <- popTime(data = bt, time = "ftime", event = "Status")



head(popTimeData)


rlang::last_error()
?popTime

head(popTimeData)
popTimeData$comprisk.event <- 0
popTimeData$comprisk.event[popTimeData$event == 1] <- 2
popTimeData$comprisk.event[popTimeData$event == 2] <- 1

table(popTimeData$comprisk.event, popTimeData$event)

str(popTimeData)
popTime(data = popTimeData, time = "time", event = "comprisk.event")




# Ribbon ------------------------------------------------------------------

h + geom_ribbon(aes(x = time, ymin=0, ymax=ycoord)) +
    theme_minimal()

pdf("geom_ribbon.pdf")
h + geom_ribbon(aes(x = time, ymin=0, ymax=ycoord)) +
    theme_minimal()
dev.off()



# Segment -----------------------------------------------------------------

h + geom_segment(aes(x = 0, xend = time, y = ycoord, yend = ycoord)) +
    theme_minimal()

pdf("geom_segment.pdf")
h + geom_segment(aes(x = 0, xend = time, y = ycoord, yend = ycoord)) +
    theme_minimal()
dev.off()



# File size in kB ---------------------------------------------------------

file.info("geom_ribbon.pdf")$size/1000
file.info("geom_segment.pdf")$size/1000



x_var <- "displ"
aes(x_var)


aes_(quote(displ), quote(ty), quote(uu))
#> Aesthetic mapping:
#> * `x` -> `displ`
aes_(as.name(x_var))
#> Aesthetic mapping:
#> * `x` -> `displ`
aes_(parse(text = x_var)[[1]])
#> Aesthetic mapping:
#> * `x` -> `displ`

f <- function(x_var) {
    aes_(substitute(x_var))
}
f(displ)
#> Aesthetic mapping:
#> * `x` -> `displ`



library(ggplot2)

mtcars$wt2 <- mtcars$wt*0.9
mtcars$mpg2 <- mtcars$mpg*0.9

mp <- list(size = 3, color = "red", pch = 10, aes(x = wt2, y = mpg2))

add <- FALSE

tt <- list(
    if (add)
        do.call("geom_point",mp)
)

p <- ggplot(mtcars, aes(wt, mpg))
p + geom_point() + tt



library(magrittr)
do.call("geom_point",mp)


cols <- c("LINE1"="#f04546","LINE2"="#3591d1","BAR"="#62c76b")
a <-c("S1","S2","S3","S4","S5","S6","S7","S8","S9") #names
b <-c(0.23,0.26,0.55,0.56,0.36,0.23,0.18,0.06,0.04) #mean t0
c <-c(0.64,0.6,0.81,1.4,0.89,0.55,0.48,0.22,0.09) #mean t1
d <-c(0.20,0.23,0.52,0.53,0.33,0.20,0.15,0.04,0.03) #SD low t0
e <-c(0.26,0.29,0.58,.59,0.39,0.26,0.21,0.08,0.05) #SD high t0
f <-c(0.67,0.63,0.86,1.44,0.93,0.59,0.51,0.25,0.10) #SD high t1
g <-c(0.61,0.57,0.78,1.36,0.85,0.53,0.45,0.19,0.08) #SD low t1
h <-c(0.41,0.34,0.26,0.84,0.53,0.32,0.30,0.16,0.05) #absolute change

data <- data.frame(a,b,c,d,e,f,g,h)
ggplot(data=data,aes(x=a)) +
    geom_bar(stat="identity", aes(y=h, fill = "BAR"),colour="#333333")+ #green
    geom_line(aes(y=b,group=1, colour="LINE1"),size=1.0) +   #red
    geom_point(aes(y=b, colour="LINE1"),size=3) +           #red
    geom_errorbar(aes(ymin=d, ymax=e, colour="LINE1"), width=0.1, size=.8) +
    geom_line(aes(y=c,group=1,colour="LINE2"),size=1.0) +   #blue
    geom_point(aes(y=c,colour="LINE2"),size=3) +           #blue
    geom_errorbar(aes(ymin=f, ymax=g,colour="LINE2"), width=0.1, size=.8) +
    scale_colour_manual(name="Error Bars",values=cols) + scale_fill_manual(name="Bar",values=cols) +
    ylab("Symptom severity") + xlab("PHQ-9 symptoms") +
    ylim(0,1.6) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
    theme(axis.title.y = element_text(size = 15, vjust=0.3))

set.seed(7)
df <- data.frame(
    x = 1:20,
    y.bar = rpois(20, lambda = 5),
    y.line = rpois(20, lambda = 10)
)
ggplot(data = df,
       mapping = aes(x = x)) +

    # specify fill for bar / color for line inside aes(); you can use
    # whatever label you wish to appear in the legend
    geom_col(aes(y = y.bar, fill = "bar.label")) +
    geom_line(aes(y = y.line, color = "line.label")) +

    xlab("Month of year") +
    scale_y_continuous(name = "Daily classifications per Spotter") +

    # the labels must match what you specified above
    scale_fill_manual(name = "", values = c("bar.label" = "grey")) +
    scale_color_manual(name = "", values = c("line.label" = "black")) +

    theme_bw()
