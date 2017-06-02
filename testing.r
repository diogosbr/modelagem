source("https://raw.githubusercontent.com/diogosbr/modelagem/master/cut_raster.R")  #carrega fun��o para cortar raster

wc = raster::getData("worldclim", var = "bio", res = 10)
unlink("./wc10/bio_10m_bil.zip")
download.file("http://biogeo.ucdavis.edu/data/gadm2.8/shp/BRA_adm_shp.zip", destfile = "bra.zip")
unzip("bra.zip", exdir = "./brasil")
unlink("bra.zip")
cut.raster(raster.dir = "wc10", shape.dir = "brasil", extension = ".bil")


source("https://raw.githubusercontent.com/diogosbr/modelagem/master/cor_data.R")  #carrega fun��o para ver correla��o entre rasters

wc = raster::getData("worldclim", var = "bio", res = 10)
cor.data(predictors = wc)


source("https://raw.githubusercontent.com/diogosbr/modelagem/master/clean.R")  #carrega fun��o para selecionar apenas pontos espacialmente �nicos e sem NA

library(dismo)
data(acaule)
pontos = as.data.frame(na.omit(cbind(acaule$lon, acaule$lat)))
wc = raster::getData("worldclim", var = "bio", res = 10)
clean(coord = pontos, predictors = wc)


source("https://raw.githubusercontent.com/diogosbr/modelagem/master/modelos.R")  #carrega fun��o para gerar os modelos

library(dismo)
data(acaule)

# Roda somente o Bioclim
pontos = as.data.frame(na.omit(cbind(acaule$lon, acaule$lat)))
wc = raster::getData("worldclim", var = "bio", res = 10)
predictors = wc
modelos(coord = pontos, diretorio = "solanum")

source("https://raw.githubusercontent.com/diogosbr/modelagem/master/toKML.R")  #carrega fun��o para gerar os modelos

library(dismo)

fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
                     pattern='grd', full.names=TRUE )
predictors <- stack(fnames)
toKML(predictors[[4]],open=T)

