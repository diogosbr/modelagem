#' @title Exibe um gráfico com as correlações entre as variáveis
#' @name cor.data
#'
#' @description Uma funcao para exibir um gráfico com as correlações entre as variáveis ambientais informadas.
#'
#' @param abio os rasters. Objeto do tipo _stack_
#' @param plot lógico. Plota um dos rasters cortados.
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed: "pearson" (default), "kendall", or "spearman".
#' @param rep númerico. Número de pontos gerados para extrair os valores dos raters e utilizar na correlação. O padrão é 1000. 
#'
#' @details O índice de correlação utilizado é spearman
#'
#' @return Retorna uma tabela com os valores de correlação entre as variáveis. 
#'
#' @author Diogo S. B. Rocha
#'
#' @seealso 
#'
#' @examples
#' fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE )
#' predictors <- stack(fnames)
#' cor.data(abio = predictors)
#' 
#' @import raster
#' @import dismo
#'
#' @export
cor.data = function(abio, plot = TRUE, method = "pearson", rep = 1000) {

    if(class(abio)=="RasterStack"|class(abio)=="RasterLayer"){
      backg <- randomPoints(abio, n = rep)
      colnames(backg) = c("long", "lat")
      backvalues = extract(abio, backg)
    }else(backvalues=abio)
    
    
    if (plot == T) {
        panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
            usr <- par("usr")
            on.exit(par(usr))
            par(usr = c(0, 1, 0, 1))
            r <- abs(cor(x, y, method = method))
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
            rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
        }
        
        pairs(backvalues, lower.panel = panel.smooth, diag.panel= panel.hist, upper.panel = panel.cor)
    }
    return(round(cor(backvalues, method = method), 2))
}