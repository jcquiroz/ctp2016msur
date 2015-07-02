library(ggplot2)
library(dplyr)
library(reshape2)
library(rlist)
library(pipeR)
library(ggthemes)

#setwd("~/Documents/southern_hake/assessment_2015/model-ctp-2016")

rm(list = ls())
source("read.admb.R")
A  <- read.admb('results/outcome_mod_1')
B  <- read.admb('results/outcome_mod_2')
C  <- read.admb('results/outcome_mod_3')
D  <- read.admb('results/outcome_mod_4')
E  <- read.admb('results/outcome_mod_3a')
Ee <- read.admb('results/outcome_mod_2a')

#tmp <- list(A,B,C,D,Ee)
tmp <- list(B,Ee)

#matplot(data.frame(t(list.rbind(list.map(tmp, SB)))), type='l') ## using rlist and pipeR

years <- 1977:2014
years <- melt( years %*% t(rep(1,8)))  # 8  --  20

varest <- melt(list.select(tmp, B.Desovante=SB/1e3, B.Total=BT/1e3, B.Juvenil=B6/1e3, Reclutamiento=R/1e6)) ## using reshape and rlist
varest$years <- years$value
varest$L1 <- paste('Data',varest$L1,sep=" ")


k <- ggplot(varest) 
k <- k + geom_line(aes(x=years, y=value, colour= factor(L1)), size=1.2)
k <- k + facet_wrap( ~ L2, nrow=2, scales='free')
k <- k + theme_calc(12) + scale_colour_hc("darkunica", guide = guide_legend(title = "Modelo"),
#        labels = c("Data 1","Data 2","Data 3","Data 4","Empirico\n(Data2)"))
       labels = c("Data 2","Teorico"))
k <- k + labs(x = "Año", y=" ")
print(k)


H = 1680
W = 1980
Res = 200
gfn  <- paste(.FIGUREDIR, "trendEmpirico2.png", sep="")
png(gfn, width=W*1.3, height=H, res=Res)
print(k)
dev.off()


varestDiff <- read.csv('varest.csv',head=TRUE)

h <- ggplot(varestDiff,aes(x = (years), y = diffPer, fill=L1)) + geom_bar(position="dodge",stat = "identity") + facet_wrap( ~ L2, nrow=2, scales='free')
#h <- h + scale_fill_manual(values = c("springgreen3","orangered3","thistle2","tan"), name=" ")
h <- h + theme_hc(14) + theme(legend.key = element_rect(colour = 'white', fill = 'white'))
h <- h + scale_fill_hc( name=" ")
h <- h + labs(x='Año', y='Diferencia entre modelos [%]')
print(h)


gfn  <- paste(.FIGUREDIR, "DiffEmpirico.png", sep="")
png(gfn, width=W, height=H*1.4, res=200)
print(h)
dev.off()



