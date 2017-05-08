cor.data=
  function(predictors , plot = TRUE){
    backg <- randomPoints(predictors , n = 1000 , extf = 1.25)
    colnames(backg) = c( 'long' ,  'lat' )
    backvalues = extract(predictors,backg)
    
    if(plot==T){
      panel.cor <- function(x , y , digits = 2, prefix = "" , cex.cor , ...){
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste0(prefix, txt)
        if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
      }
    
    pairs(backvalues,lower.panel = panel.smooth, upper.panel = panel.cor)
    }
    return(round(cor(backvalues) , 2))
  }

