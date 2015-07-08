library(ggplot2)
library(reshape2)
library(SpatialTools)



source("read.admb.R")
.VIEWLAS    <- 1


A  <- read.admb('mod_2')      ## modelo Data 1

deple <- A$RPRp  # Yproy/1000 # BDp/1000  # RPRp
x <- seq(2015,2015+49)
y <- seq(0,1.2,0.05)


deple[,c(1,5,9,13,17,21)]


filled.contour(x, y, deple, color.palette = terrain.colors)
cL = contourLines(x,y,deple, levels=0.4)
plot.contourLines(t(cL))


#typical plot
filled.contour(x, y, deple, color.palette = terrain.colors)

#try
cont <- contourLines(y, x, t(deple), levels=c(0.2,0.3,0.4,0.5))
fun <- function(x) x$level
LEVS <- sort(unique(unlist(lapply(cont, fun))))
COLS <- terrain.colors(length(LEVS))
contour(y, x, t(deple), levels=c(0.2,0.3,0.5), lty = "solid", xlim=c(0,1), ylim=c(2015,2040),
        ylab="Periodo recuperacion", xlab="Ponderadores de Fmrs", main="Niveles de recuperacion de BD")
for(i in seq(cont)){
  COLNUM <- match(cont[[i]]$level, LEVS)
  polygon(cont[[i]], col=COLS[COLNUM], border="NA")
}
contour(y, x, t(deple), levels=c(0.4), col=2, add=TRUE, lwd=2) 


grid()


## -----------------lee proyeciones MCMC ------------------
## --------------------------------------------------------
## --------------------------------------------------------

mcmc.py <- read.csv('proyecciones.mcmc.out', header=F, sep=' ')
ft      <- c('vec','year','0.0Fmrs','0.2Fmsr','0.4Fmrs','0.6Fmrs','0.8Fmrs','1.0Fmrs','1.2Fmrs')
colnames(mcmc.py) <- ft
mcmc.py <- mcmc.py[,-9]


median.quartile <- function(x){
  out <- c(mean(x), quantile(x, probs = c(0.5)))
  names(out) <- c("mean","median")
  return(out) 
}

stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, colour="red", geom=geom, size = 3, ...)
}


# Captura --------------------------------
pos.captura <- subset(mcmc.py , vec == 'Captura')
pos.cap.mel <- melt(pos.captura, id.vars=c('vec','year'))

h <- ggplot(pos.cap.mel, aes(as.factor(year),value/1000)) + geom_violin(fill = "orange") 
h <- h + facet_wrap( ~ variable) + 
  stat_summary(fun.y=median.quartile, geom='point', colour='red', size=2)
#h <- h + facet_wrap( ~ variable) + stat_sum_single(mean, geom="line")
h <- h + ggtitle ("Proyeccion Captura") +
  xlab("Ahno") +  ylab ("Toneladas [1000]")
print(h)

# Captura --------------------------------
pos.depleti <- subset(mcmc.py , vec == 'depletion')
pos.dep.mel <- melt(pos.depleti, id.vars=c('vec','year'))

h <- ggplot(pos.dep.mel, aes(as.factor(year),value)) + geom_violin(fill = "orange") 
h <- h + facet_wrap( ~ variable) + geom_hline(yintercept=0.4, linetype=2) + 
  stat_summary(fun.y=median.quartile, geom='point', colour='red', size=2)
#h <- h + facet_wrap( ~ variable) + stat_sum_single(mean, geom="line")
h <- h + ggtitle ("Reduccion respecto BDo") +
  xlab("Año") +  ylab ("Proporcion [BD/Bo]") 
h <- h + scale_x_discrete(breaks=seq(2015,2040,5))
print(h)















