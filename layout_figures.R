.pngPlots <- function( repObj )
{
  
  fileprefix <- paste(.FIGUREDIR, "Msur_2016.", sep="")
  H = 1680
  W = 1980
  Res = 200
  
  gfn  <- paste(.FIGUREDIR, "AgeCompsTrawl.png", sep="")
  png(gfn, width=W, height=H, res=Res)
  .plotAgecomps.Trawl    ( repObj, meanAge = TRUE, cprox = TRUE  )
  dev.off()
 
  gfn  <- paste(.FIGUREDIR, "AgeCompsLongline.png", sep="")
  png(gfn, width=W, height=H, res=Res)
  .plotAgecomps.Longline    ( repObj, meanAge = TRUE, cprox = TRUE  )
  dev.off()
  
  gfn  <- paste(.FIGUREDIR, "AgeCompsArtisanalDat2.png", sep="")
  png(gfn, width=W, height=H, res=Res)
  .plotAgecomps.Artisanal    ( repObj, meanAge = TRUE, cprox = TRUE  )
  dev.off()

  gfn  <- paste(.FIGUREDIR, "AgeCompsSurvey.png", sep="")
  png(gfn, width=W, height=H, res=Res)
  .plotAgecomps.Survey    ( repObj, meanAge = TRUE, cprox = TRUE  )
  dev.off()
  
  gfn  <- paste(.FIGUREDIR, "index.png", sep="")
  png(gfn, width=W*1.2, height=H, res=Res)
  .plotIndex    ( repObj, annotate=TRUE )
  dev.off()

  gfn  <- paste(.FIGUREDIR, "index2.png", sep="")
  png(gfn, width=W*1.2, height=H, res=Res)
  .plotIndex2    ( repObj, annotate=TRUE )
  dev.off()  
  
  gfn  <- paste(.FIGUREDIR, "index.png", sep="")
  png(gfn, width=W*1.2, height=H, res=Res)
  .plotIndex    ( repObj, annotate=TRUE )
  dev.off()
  
  gfn  <- paste(.FIGUREDIR, "Ygear.png", sep="")
  png(gfn, width=W*1.2, height=H, res=Res)
  .plotY( repObj, annotate=TRUE )
  dev.off()
  
  gfn  <- paste(.FIGUREDIR, "Descarte.png", sep="")
  png(gfn, width=W*1.3, height=H, res=Res+50)
  .plotD( repObj, annotate=TRUE )
  dev.off()
  
  gfn  <- paste(.FIGUREDIR, "wmedadPop.png", sep="")
  png(gfn, width=W*1.2, height=H, res=Res)
  par(mfrow=c(1,2))
  .plotMeanwt( repObj, wm="WmPopu" )
  .plotMeanwt( repObj, wm="WiniPopu" )
  dev.off()
  
  gfn  <- paste(.FIGUREDIR, "wmedadTeo.png", sep="")
  png(gfn, width=W, height=H, res=Res)
  par(mfrow=c(1,1))
  .plotMeanEmp( repObj, wm="WmTrawl" )
  dev.off()
    
  gfn  <- paste(.FIGUREDIR, "wmedad.png", sep="")
  png(gfn, width=W, height=H*1.3, res=Res)
  par(mfrow=c(2,2))
  .plotMeanwt( repObj, wm="WmTrawl" )
  .plotMeanwt( repObj, wm="WmLong" )
  .plotMeanwt( repObj, wm="WmArti" )
  .plotMeanwt( repObj, wm="WmSurv" )  
  dev.off()
  
}

.plotMeanEmp  <- function( repObj, wm="WmTrawl" )
{
  #plot mean weight-at-age by cohort
  with(repObj, {
    xx = yr  	## xaxis labels
    yy = age	## yaxis labels
    nage=length(age)
    p_em <- read.csv('data-control/PesoTeorico.csv', header=TRUE, sep=",")
    library(pracma)
    
    WmPopu <- repmat(p_em[,2],length(xx),1)/1000
    
    if(sum(par("mfcol"))==2)
    {
      xl = "Cohorte"; xlm=""
      yl = "Peso - edad (kg)"; ylm=""
    }
    else
    {
      xlm = "Cohorte" ;xl=""
      ylm = "Peso - edad (kg)"; yl=""
    }
    
    plot(range(xx), range(WmPopu), type="n", axes=FALSE,
         xlab=xl, ylab=yl, main= paste("Peso Medio: ", "Teorico", sep=""))
    axis( side=1 )
    axis( side=2, las=.VIEWLAS )
    box()
    grid()
    
    for(i in 1:dim(WmPopu)[1]-1)
    {
      yy = (diag(as.matrix(WmPopu[0:-i, ]))) 
      xx = 1:length(yy)+yr[i]-min(age)+1
      
      yy[yy==0]=NA;xx[yy==NA]=NA
      lines(xx,yy,col="slategray4")
      
      points(xx[1],yy[1],pch=20,col="steelblue",cex=1.4)
      points(xx[nage],yy[nage],pch=20,col="salmon",cex=1.4)      
      points(xx[10],yy[10],pch=20,col="seagreen4",cex=1.4)  
      points(xx[15],yy[15],pch=20,col="slategray",cex=1.4)
    }
    for(i in 1:dim(WmPopu)[2]-1)
    {
      yy = diag(as.matrix(WmPopu[,-1:-i]))
      n = length(yy)
      xx = yr[1]:(yr[1]+n-1)
      lines(xx, yy,col="slategray4")
      points(xx[n], yy[n], pch=20, col="salmon", cex=1.4) 
      if (length(xx)>=15)
      {
        ubi1 = xx[n-14]
        if (ubi1 >= 1977)
        {
          points(xx[n-14], yy[n-14], pch=20, col="seagreen4", cex=1.4)  
        }
      }
      if (length(xx)>=10)
      {
        ubi1 = xx[n-9]
        if (ubi1 >= 1977)
        {
          points(xx[n-9], yy[n-9], pch=20, col="slategray", cex=1.4)  
        }
      }
    }
    
    mtext(xlm, side=1, outer=FALSE, line=2)
    mtext(ylm, side=2, outer=FALSE, line=2)
  })
}

.plotMeanwt  <- function( repObj, wm="WmTrawl" )
{
  #plot mean weight-at-age by cohort
  with(repObj, {
    posi = which(names(JJ)==wm)
    WmPopu = as.data.frame(repObj[posi])/1000    
    #WmPopu[WmPopu == 0] <- NA
    xx = yr		## xaxis labels
    yy = age	## yaxis labels
    nage=length(age)
    
    if(sum(par("mfcol"))==2)
    {
      xl = "Cohorte"; xlm=""
      yl = "Peso - edad (kg)"; ylm=""
    }
    else
    {
      xlm = "Cohorte" ;xl=""
      ylm = "Peso - edad (kg)"; yl=""
    }
    
    plot(range(xx), range(WmPopu), type="n", axes=FALSE,
         xlab=xl, ylab=yl, main= paste("Peso Medio: ", wm, sep=""))
    axis( side=1 )
    axis( side=2, las=.VIEWLAS )
    box()
    grid()
    
    for(i in 1:dim(WmPopu)[1]-1)
    {
      yy = (diag(as.matrix(WmPopu[0:-i, ]))) 
      xx = 1:length(yy)+yr[i]-min(age)+1
      
      yy[yy==0]=NA;xx[yy==NA]=NA
      lines(xx,yy,col="slategray4")
      
      points(xx[1],yy[1],pch=20,col="steelblue",cex=1.4)
      points(xx[nage],yy[nage],pch=20,col="salmon",cex=1.4)      
      points(xx[10],yy[10],pch=20,col="seagreen4",cex=1.4)  
      points(xx[15],yy[15],pch=20,col="slategray",cex=1.4)
    }
    for(i in 1:dim(WmPopu)[2]-1)
    {
      yy = diag(as.matrix(WmPopu[,-1:-i]))
      n = length(yy)
      xx = yr[1]:(yr[1]+n-1)
      lines(xx, yy,col="slategray4")
      points(xx[n], yy[n], pch=20, col="salmon", cex=1.4) 
      if (length(xx)>=15)
      {
        ubi1 = xx[n-14]
        if (ubi1 >= 1977)
        {
          points(xx[n-14], yy[n-14], pch=20, col="seagreen4", cex=1.4)  
        }
      }
      if (length(xx)>=10)
      {
        ubi1 = xx[n-9]
        if (ubi1 >= 1977)
        {
          points(xx[n-9], yy[n-9], pch=20, col="slategray", cex=1.4)  
        }
      }
    }
          
    mtext(xlm, side=1, outer=FALSE, line=2)
    mtext(ylm, side=2, outer=FALSE, line=2)
  })
}

.plotY   <- function( repObj, annotate=FALSE )
{
  #line plot for relative abundance indices
  des            <- as.data.frame(repObj$Ygear)
  colnames(des)  <- c('Arrastre','Palangre','Artesanal')
  des$Total      <- rowSums(des)
  des$type       <- 'Desembarque'
  des$yr         <- repObj$yr
  
  catch            <- as.data.frame(repObj$Ycorrected)
  colnames(catch)  <- c('Arrastre','Palangre','Artesanal')
  catch$Total      <- rowSums(catch)
  catch$type       <- 'Captura'
  catch$yr         <- repObj$yr
  
  matY             <- rbind(des,catch)
  
  matYY            <- melt(matY, id=c(5,6))
  
#   type        yr    variable     value
#   Desembarque 1977  Arrastre     2250.000
  
 h <- ggplot(matYY) + geom_line(aes(x=yr,y=value/1000,colour=type,group=type), size=1) 
 h <- h + facet_wrap( ~ variable, scales = "free") + theme_bw(15) + labs(x='Year', y=' ')
 h <- h + scale_colour_manual(values = c("palegreen3","orangered3"), name=" ")
 #h <- h + scale_x_discrete(breaks = round((seq(2010,2050,5)),digits = 3))
 h <- h + theme_calc(14) + theme(legend.key = element_rect(colour = 'white', fill = 'white'))
 print(h)
  
}

.plotD   <- function( repObj, annotate=FALSE )
{
  #line plot for relative abundance indices
  des            <- as.data.frame(repObj$Ygear)
  colnames(des)  <- c('Arrastre','Palangre','Artesanal')
  des$Total      <- rowSums(des)
  des$type       <- 'Desembarque'
  des$yr         <- repObj$yr
  
  catch            <- as.data.frame(repObj$Ycorrected)
  colnames(catch)  <- c('Arrastre','Palangre','Artesanal')
  catch$Total      <- rowSums(catch)
  catch$type       <- 'Captura'
  catch$yr         <- repObj$yr
  
  diffC            <- catch[,1:4] - des[,1:4]
  diffC$year       <- repObj$yr
  
  diffCY          <- melt(diffC, id=c(5))
  
  h <- ggplot(subset(diffCY, year>=2000),aes(x = factor(year), y = value, fill=variable)) + geom_bar(position="dodge",stat = "identity")
  h <- h + scale_fill_manual(values = c("springgreen3","orangered3","thistle2","tan"), name=" ")
  h <- h + theme_hc(14) + theme(legend.key = element_rect(colour = 'white', fill = 'white'))
  h <- h + labs(x='A?o', y='Descarte / Subreporte [toneladas]')
  print(h)
  
}

.plotIndex2  <- function( repObj, annotate=FALSE )
{
  #line plot for relative abundance indices
  with(repObj, {
    index <- cbind(indexDat2,survey/max(survey))
    index[index == 0] <- NA
    if(is.matrix(repObj$index)){
      xx=yr
      yy=index
    }else{
      xx=yr
      yy=index
    }
    n=nrow(t(as.matrix(yy)))
    yrange=c(0, max(yy, na.rm=TRUE))
    
    matplot(xx, yy, type="n", axes=FALSE,
            xlab="A?o", ylab="Abundacia Relativa", 
            ylim=yrange , main="Series Temporales")
    
    colores = c("steelblue2","tomato2","peachpuff3","orchid4")
    matlines(xx, yy, col= colores,type="o", pch=1:n, cex=1.5, lwd=2)
    
    axis( side=1 )
    axis( side=2, las=.VIEWLAS )
    box()
    grid()
    
    if( annotate )
    {
      txt= c('Arrastre','Palangre','Artesanal','Acustico')
      legend("top", txt, lty=1:n, pch=1:n, bty="n", col= colores)
    }
  })
}

.plotIndex  <- function( repObj, annotate=FALSE )
{
  #line plot for relative abundance indices
  with(repObj, {
    index <- cbind(index,survey/max(survey))
    index[index == 0] <- NA
    if(is.matrix(repObj$index)){
      xx=yr
      yy=index
    }else{
      xx=yr
      yy=index
    }
    n=nrow(t(as.matrix(yy)))
    yrange=c(0, max(yy, na.rm=TRUE))
    
    matplot(xx, yy, type="n", axes=FALSE,
            xlab="A?o", ylab="Abundacia Relativa", 
            ylim=yrange , main="Series Temporales")
    
    colores = c("steelblue2","tomato2","peachpuff3","orchid4")
    matlines(xx, yy, col= colores,type="o", pch=1:n, cex=1.5, lwd=2)
    
    axis( side=1 )
    axis( side=2, las=.VIEWLAS )
    box()
    grid()
    
    if( annotate )
    {
      txt= c('Arrastre','Palangre','Artesanal','Acustico')
      legend("top", txt, lty=1:n, pch=1:n, bty="n", col= colores)
    }
  })
}

.plotAgecomps.Trawl  <- function(repObj, meanAge = FALSE, cprox = TRUE )
{
  # Bubble plot of age-composition data
  with( repObj, {
    if(!is.null(repObj$catTrawl))
      {
      xx = yr
      zz = t(catTrawl)
      xrange = range(yr)
        
        
        # plot proportions-at-age (cpro=TRUE)
        plotBubbles(zz, xval = xx, yval = age, cpro=cprox, hide0=TRUE,  
                    las=.VIEWLAS, xlab="A?o", ylab="Edad", frange=0.05, size=0.1, 
                    bg=colr("steelblue", 0.5), main="Captura Edad Arrastre", 
                    xlim=xrange)

        grid()
        
        if( meanAge )
        {
          tz = t(zz)
          p = t(tz/rowSums(tz))
          abar = colSums(t(tz/rowSums(tz))*age)
          sbar = sqrt(colSums(p*(1-p)*age))
          sbar = 1.96*colSums(sqrt(p*(1-p))/sqrt(age))
          
          lines( xx, abar, col=colr("springgreen4", 0.85), lwd=3 )
          
          yy = c(exp(log(abar)+log(sbar)), rev(exp(log(abar)-log(sbar))))
          #polygon(c(xx, rev(xx)),yy,border=NA,col=colr("steelblue",0.25))
        }
        
      }
    else
    {
      print("There is no age-composition data")
    }
  })
}

.plotAgecomps.Longline  <- function(repObj, meanAge = FALSE, cprox = TRUE )
{
  # Bubble plot of age-composition data
  with( repObj, {
    if(!is.null(repObj$catLong))
    {
      xx = yr
      zz = t(catLong)
      xrange = range(yr)
      
      
      # plot proportions-at-age (cpro=TRUE)
      plotBubbles(zz, xval = xx, yval = age, cpro=cprox, hide0=TRUE,  
                  las=.VIEWLAS, xlab="A?o", ylab="Edad", frange=0.05, size=0.1, 
                  bg=colr("steelblue", 0.5), main="Captura Edad Palangre", 
                  xlim=xrange)
      
      grid()
      
      if( meanAge )
      {
        tz = t(zz)
        p = t(tz/rowSums(tz))
        abar = colSums(t(tz/rowSums(tz))*age)
        sbar = sqrt(colSums(p*(1-p)*age))
        sbar = 1.96*colSums(sqrt(p*(1-p))/sqrt(age))
        
        lines( xx, abar, col=colr("springgreen4", 0.85), lwd=3 )
        
        yy = c(exp(log(abar)+log(sbar)), rev(exp(log(abar)-log(sbar))))
        #polygon(c(xx, rev(xx)),yy,border=NA,col=colr("steelblue",0.25))
      }
      
    }
    else
    {
      print("There is no age-composition data")
    }
  })
}

.plotAgecomps.Artisanal  <- function(repObj, meanAge = FALSE, cprox = TRUE )
{
  # Bubble plot of age-composition data
  with( repObj, {
    if(!is.null(repObj$catArtiDat2))
    {
      xx = yr
      zz = t(catArtiDat2)
      xrange = range(yr)
      
      
      # plot proportions-at-age (cpro=TRUE)
      plotBubbles(zz, xval = xx, yval = age, cpro=cprox, hide0=TRUE,  
                  las=.VIEWLAS, xlab="A?o", ylab="Edad", frange=0.05, size=0.1, 
                  bg=colr("steelblue", 0.5), main="Captura Edad Artesanal", 
                  xlim=xrange)
      
      grid()
      
      if( meanAge )
      {
        tz = t(zz)
        p = t(tz/rowSums(tz))
        abar = colSums(t(tz/rowSums(tz))*age)
        sbar = sqrt(colSums(p*(1-p)*age))
        sbar = 1.96*colSums(sqrt(p*(1-p))/sqrt(age))
        
        lines( xx, abar, col=colr("springgreen4", 0.85), lwd=3 )
        
        yy = c(exp(log(abar)+log(sbar)), rev(exp(log(abar)-log(sbar))))
        #polygon(c(xx, rev(xx)),yy,border=NA,col=colr("steelblue",0.25))
      }
      
    }
    else
    {
      print("There is no age-composition data")
    }
  })
}

.plotAgecomps.Survey  <- function(repObj, meanAge = FALSE, cprox = TRUE )
{
  # Bubble plot of age-composition data
  with( repObj, {
    if(!is.null(repObj$catSurv))
    {
      xx = yr
      zz = t(catSurv)
      xrange = range(yr)
      
      
      # plot proportions-at-age (cpro=TRUE)
      plotBubbles(zz, xval = xx, yval = age, cpro=cprox, hide0=TRUE,  
                  las=.VIEWLAS, xlab="A?o", ylab="Edad", frange=0.05, size=0.1, 
                  bg=colr("steelblue", 0.5), main="Abundancia Edad Crucero", 
                  xlim=xrange)
      
      grid()
      
      if( meanAge )
      {
        tz = t(zz)
        p = t(tz/rowSums(tz))
        abar = colSums(t(tz/rowSums(tz))*age)
        sbar = sqrt(colSums(p*(1-p)*age))
        sbar = 1.96*colSums(sqrt(p*(1-p))/sqrt(age))
        
        lines( xx, abar, col=colr("springgreen4", 0.85), lwd=3 )
        
        yy = c(exp(log(abar)+log(sbar)), rev(exp(log(abar)-log(sbar))))
        #polygon(c(xx, rev(xx)),yy,border=NA,col=colr("steelblue",0.25))
      }
      
    }
    else
    {
      print("There is no age-composition data")
    }
  })
}
