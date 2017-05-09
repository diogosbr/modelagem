source("cut_raster.r")  #carrega função para cortar raster

wc = raster::getData("worldclim", var = "bio", res = 10)
unlink("./wc10/bio_10m_bil.zip")
download.file("http://biogeo.ucdavis.edu/data/gadm2.8/shp/BRA_adm_shp.zip", destfile = "bra.zip")
unzip("bra.zip", exdir = "./brasil")
unlink("bra.zip")
cut.raster(raster.dir = "wc10", shape.dir = "brasil", extension = ".bil")


source("cor_data.r")  #carrega função para ver correlação entre rasters

wc = raster::getData("worldclim", var = "bio", res = 10)
cor.data(predictors = wc)


source("clean.r")  #carrega função para selecionar apenas pontos espacialmente únicos e sem NA

library(dismo)
data(acaule)
pontos = as.data.frame(na.omit(cbind(acaule$lon, acaule$lat)))
wc = raster::getData("worldclim", var = "bio", res = 10)
clean(coord = pontos, predictors = wc)


source("modelos.r")  #carrega função para gerar os modelos

# Roda somente o Bioclim
pontos = as.data.frame(na.omit(cbind(acaule$lon, acaule$lat)))
wc = raster::getData("worldclim", var = "bio", res = 10)
predictors = wc
modelos(coord = pontos, diretorio = "solanum")


