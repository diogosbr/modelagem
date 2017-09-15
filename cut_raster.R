#' @title Corta raster com base em um shape
#' @name cut.raster
#'
#' @description Uma funcao para cortar um raster com base em um shapefile informado.
#'
#' @param abio os rasters a serem cortados. Aceita um objeto do tipo _stack_
#' @param shape.dir aqui tem duas opções. A primeira é informar um ojeto no formato shapefile. A segunda opção é informar um diretório que contem o shapefile.
#' @param extension extensão de saída dos rasters cortados. O padrão é .asc, veja \code{\link[raster]{write.raster}} para mais possibilidades de formatos de saída.
#' @param br lógico. Se TRUE, utiliza o shape do brasil formecido por \code{\link[maptools]{wrld_simpl}}. Se FALSE (padrão), utliza o shape informado.
#' @param plot lógico. Plota um dos rasters cortados.
#' @param trim lógico. Se TRUE, exclui os pixels com NA (mais demorado). Se FALSE (padrão), mantém os NAs.
#'
#' @details Quando argumento shapedir for uma pasta onde contém o shape de referência, esta pasta deve conter somente o shape e seus arquivos auxiliares.
#'
#' @return Arquivos raster em uma pasta chamada "Cortados".
#'
#' @author Diogo S. B. Rocha
#'
#' @seealso \code{\link[raster]{crop}}, \code{\link[raster]{mask}}
#'
#' @examples
#' fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE )
#' predictors <- stack(fnames)
#' cut.raster(abio = predictors, br=T, extension = ".tif")
#' 
#' @import raster
#'
#' @export
cut.raster = 
  function(abio, shape.dir, extension = ".asc", br = FALSE, plot = TRUE, trim = FALSE){
    
    require(raster)
    # original=getwd()
    
    if (dir.exists("Cortados") == F) {
        dir.create("Cortados")
    }
    
    # Definir shape para cortar
    if(missing(shape.dir)){
      # Possibilidade de cortar para o Brasil
      #importando shape do brasil
      if(br==T){
        data(wrld_simpl, package = "maptools")
        br=subset(wrld_simpl, wrld_simpl$NAME=="Brazil")
      }else(stop("Não selecionou a pasta contendo o shape de corte"))
      }
    
    if(class(shape.dir)== "character") {
      shape = rgdal::readOGR(list.files(shape.dir, pattern = ".shp", full.names = T)[1])
        stop("Não selecionou a pasta contendo o shape de corte")
    } else (shape = shape.dir)
    

    
   
    crs(shape) = crs(abio)
    crs(br) = crs(abio)
    
    # loop para cortar todos os rasters
    
    # sem trim
    if (trim == F) {
        ini = Sys.time()
        for (i in 1:length(names(abio))) {
            ini1 = Sys.time()
            
            mask(crop(abio[[i]], extent(shape)), shape, filename = paste("./Cortados", "/", names(abio)[i], 
                sep = "", extension), overwrite = TRUE)
            
            print(Sys.time())
            cat("\n", paste("Tá indo", i))
            
            fim1 = Sys.time()
            cat(paste("\n", round(as.numeric(fim1 - ini1), 2), units(fim1 - ini1)))
            if (i == length(names(abio))) {
                cat("\n", "Acabou!")
                fim = Sys.time()
                fim - ini
            }
        }
    }
    
    
    # com trim
    if (trim == T) {
        ini = Sys.time()
        for (i in 1:length(names(abio))) {
            ini1 = Sys.time()
            
            writeRaster(trim(mask(crop(abio[[i]], extent(shape)), shape)), 
                filename = paste(".\\Cortados", "\\", names(abio)[i], ".tif", sep = ""), 
                format = "GTiff", overwrite = TRUE, NAflag = -9999)
            cat("\r")
            print(Sys.time())
            cat("\n", paste("Tá indo", i))
            fim1 = Sys.time()
            cat(paste("\n", round(as.numeric(fim1 - ini1), 2), units(fim1 - ini1)))
            if (i == length(names(abio))) {
                cat("\n", "Acabou!", "\n")
                fim = Sys.time()
                fim - ini
            }
        }
    }
    
    # plotando a primeira variavel cortada
    if (plot == T) {
        plot(raster(list.files("./cortados", pattern = "", full.names = TRUE)[1], native = TRUE))
        plot(shape, add = T)
    }
}
