pacman::p_load(TH.data)
pacman::p_load(usethis)

data("GBSG2")

GBSG2 <- transform(GBSG2, hormon=as.numeric(GBSG2$horTh)-1)
GBSG2 <- transform(GBSG2, meno=as.numeric(GBSG2$menostat)-1)
brcancer <- GBSG2

usethis::use_data(brcancer, overwrite = TRUE)
