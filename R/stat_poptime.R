p <- ggplot(DT, aes(x=0, xend=time, y=ycoord, yend=ycoord))

p + theme_classic() +
    geom_segment(size=3, colour="grey90") +
    geom_point(aes(x=time, y=yc), data = DT[event==1], size=1, colour="red") +
    theme(axis.text=element_text(size=12, face='bold'))
