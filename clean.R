clean = function(coord, predictors) {
    if (dim(coord)[2] == 2) {
        if (exists("predictors")) {
            # selecionar os pontos �nicos e sem NA
            pts = coord
            mask = predictors[[1]]
            # Selecionar pontos espacialmente �nicos #
            cell <- raster::cellFromXY(mask, pts)  # get the cell number for each point
            dup <- duplicated(cell)
            pts1 <- pts[!dup, ]  # select the records that are not duplicated
            pts1 <- pts1[!is.na(extract(mask, pts1)), ]  #selecionando apenas pontos que tem valor de raster
            
            cat(dim(pts)[1] - dim(pts1)[1], "pontos retirados\n")
            cat(dim(pts1)[1], "pontos espacialmente �nicos\n")
            # pts1
            names(pts1) = c("long", "lat")#
            return(pts1)
        } else (cat("Indique o objeto com as vari�veis preditoras"))
    } else (stop("Tabela de coordenadas tem mais de duas colunas.\nEsta tabela deve ter apenas long e lat, nesta ordem."))
}
