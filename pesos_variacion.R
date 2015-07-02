## exploring length-weight realtionship
library(reshape2)
library(ggplot2)
library(ggthemes)

# parameters
Linf = 121
k    = 0.080
to   = -1.4571

ana = 0.0023892
bna = 3.2409166
asa = 0.0025665
bsa = 3.2250394
asp = 0.0038372
bsp = 3.1400165

rangel = 0:120
age    = 1:24


# curves

lmed = Linf*(1-exp(-k*(age-to)))

lp.an = ana*lmed^bna
lp.as = asa*lmed^bsa
lp.ps = asp*lmed^bsp

tmp <- as.data.frame(cbind(age,lmed,lp.an,lp.as,lp.ps))
tmp <- melt(tmp, id=c("age","lmed"))

h <- ggplot(tmp, aes(lmed,value/1000, colour=variable))
h <- h + geom_line()
h <- h + theme_calc(12) + scale_colour_hc("darkunica", guide = guide_legend(title = "Data"))
h <- h + labs(x = "Talla (cm)", y="Peso (kg)")
h

p_em <- cbind(age=1:24, read.csv('data-control/peso_empirico.ooo', header=TRUE, sep=" "))

k <- ggplot()
k <- k + geom_line(data=tmp, aes(x=age, y=value/1000, colour=variable))
k <- k + theme_calc(16) + scale_colour_hc("darkunica", guide = guide_legend(title = "Data"),
                                          labels = c("Arrastre Norte","Arrastre Sur","Palangre Sur"))
k <- k + labs(x = "Edad", y="Peso (kg)")
k <- k + geom_line(data=p_em, aes(x=age, y=peso/1000), colour=1, size=1, linetype=4)
print(k)

H = 1680
W = 1980
Res = 200
gfn  <- paste(.FIGUREDIR, "pesomedio2.png", sep="")
png(gfn, width=W*1.3, height=H, res=Res)
print(k)
dev.off()

