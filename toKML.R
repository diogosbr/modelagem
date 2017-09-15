#' @title Transforma o raster em arquivo KML
#' @name toKML
#'
#' @description Uma função para transformar o modelo de nicho de formato raster em arquivo KML.
#'
#' @param modelo raster com o modelo de distribuição
#' @param name nome do arquivo KML gerado.
#' @param zeros lógico. Se FALSE, os zeros do modelo são convertivos em NA. 
#' @param open lógico. Se TURE, abre o arquivo KML no Google Earth.
#' 
#' @details aaa
#'
#' @return data.frame contendo longitude e latitude.
#'
#' @author Diogo S. B. Rocha
#'
#' @seealso \code{\link[plotKML]{plotKML}}
#'
#' @examples
#' fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE )
#' predictors <- stack(fnames)
#' toKML(modelo = predictors, name = "exemplo")
#' 
#' @import plotKML
#'
#' @export
toKML = function(modelo, name = "meuKML", zeros = FALSE, open = FALSE){
  if(missing(modelo)){stop("Informe o nome do objeto que contem o raster do modelo")}
  if(zeros==FALSE){
    values(modelo)[values(modelo)==0]=NA
    plotKML(modelo,folder.name=name,open.kml = open)
  }else(plotKML(modelo,folder.name=name,open.kml = open))
}
