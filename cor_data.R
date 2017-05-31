cor.data = function(predictors, plot = TRUE) {
    require(raster)
    require(dismo)
    if(class(predictors)=="RasterStack"|class(predictors)=="RasterLayer"){
      backg <- randomPoints(predictors, n = 1000)
      colnames(backg) = c("long", "lat")
      backvalues = extract(predictors, backg)
    }else(backvalues=predictors)
    
    
    if (plot == T) {
        panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
            usr <- par("usr")
            on.exit(par(usr))
            par(usr = c(0, 1, 0, 1))
            r <- abs(cor(x, y))
            txt <- format(c(r, 0.123456789), digits = digits)[1]
            txt <- paste0(prefix, txt)
            if (missing(cex.cor)) 
                cex.cor <- 0.8/strwidth(txt)
            text(0.5, 0.5, txt, cex = cex.cor * r)
        }
        
        panel.hist <- function(x, ...){
            usr <- par("usr"); on.exit(par(usr))
            par(usr = c(usr[1:2], 0, 1.5) )
            h <- hist(x, plot = FALSE)
            breaks <- h$breaks; nB <- length(breaks)
            y <- h$counts; y <- y/max(y)
            rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
        }
        
        pairs(backvalues, lower.panel = panel.smooth, diag.panel= panel.hist, upper.panel = panel.cor)
    }
    return(round(cor(backvalues), 2))
}