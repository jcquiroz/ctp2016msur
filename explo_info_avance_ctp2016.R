rm(list = ls())

library(ggplot2)
library(dplyr)
library(reshape2)
library(rlist)
library(pipeR)
library(ggthemes)
library(PBSmodelling)


source("read.admb.R")
.VIEWLAS    <- 1
.FIGUREDIR  <- "figures_avance/"

A  <- read.admb('results/outcome_mod_1')      ## modelo Data 1
B  <- read.admb('results/outcome_mod_2')      ## modelo Data 2
C  <- read.admb('results/outcome_mod_3')      ## modelo Data 3
D  <- read.admb('results/outcome_mod_4')      ## modelo Data 4
E  <- read.admb('results/outcome_mod_3a')     ## modelo Data 3 - usando peso empirico 
Ee <- read.admb('results/outcome_mod_2a')     ## modelo Data 2 - usando peso empirico
JJ <- read.rep('data-control/to_R_data4.jcq') ## report Data 4


## ---------- Actualiza data 2014 ------------------------
dat2014 <- read.csv('compilacion.csv', head=TRUE)

# Catage
JJ$catTrawl[38,] <- as.matrix(dat2014[5,-c(1:3)])
JJ$catLong[38,]  <- as.matrix(dat2014[10,-c(1:3)])
JJ$catArti[38,]  <- as.matrix(dat2014[15,-c(1:3)])
JJ$catArtiDat2[38,]  <- as.matrix(dat2014[15,-c(1:3)])
JJ$catSurv[38,]  <- as.matrix(dat2014[18,-c(1:3)])

# Index & landing
JJ$index[38,]     <- as.matrix(c(0.58,0.33,0.21))
JJ$indexDat2[38,] <- as.matrix(c(0.73,0.34,0.28))

p1 = JJ$Ygear[38,2]/(JJ$Ygear[38,2]+JJ$Ygear[38,3])
p2 = (1 - p1)
JJ$Ygear[38,]       <- as.matrix(c(7145,p1*5251,p2*5251))
JJ$Ycorrected[38,]  <- JJ$Ygear[38,] * c(1.619,1.01,1.109)

# catage proportion
JJ$catTrawlp <- prop.table(JJ$catTrawl, 1)

## ---------- Graficos informe 2016 ------------------------
source("layout_figures.R")
.pngPlots(JJ)









