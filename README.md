# Modelagem de nicho
Funções para auxiliar a modelagem de distribuição de nicho

##Funções disponíveis
* cut.raster
* cor.data
* clean
* modelos 
* toKML

No arquivo *testing.R* tem exemplos os exemplos para testar as funções. 

## **Funções** 
### *cut.raster()*
Função para recortar as variáveis ambientais (*raster*) a partir de um *shapefile*

cut.raster(raster.dir , shape.dir , extension = ".asc" , plot = TRUE , trim = FALSE)

**Argumentos:**

 * raster.dir: diretório que contém os rasters a serem cortados. Se não for informado, então vai procurar os rasters na pasta local.
 * shape.dir: diretório que está o shape que será usado como máscara para cortar os rasters. Obrigatório.
 * extension: conjuntos de caracteres com a extensão dos *ratsers*. ".asc" é o padrão.
 * plot: plota os rasters. TRUE é o padrão.
 * trim: se for TRUE os NAs gerados após o corte dos rasters são removidos. FALSE é o padrão.

**Exemplo:**

    wc=raster::getData('worldclim', var = 'bio' , res = 10 )
    unlink("./wc10/bio_10m_bil.zip")
    download.file("http://biogeo.ucdavis.edu/data/gadm2.8/shp/BRA_adm_shp.zip", destfile = "bra.zip") 
    unzip("bra.zip",exdir = "./brasil")
    unlink("bra.zip")
    cut.raster(raster.dir = "wc10" , shape.dir = "brasil" , extension = ".bil")

---

### *cor.data()*
Função para verificar a correlação entre as variáveis ambientais

cor.data(predictors , plot = TRUE)

**Argumentos:**

* predictors: objeto com as variáveis ambientais (*raster*)
* plot: se for TRUE (padrão), plota gráfico e valores de correlação

**Exemplo:** 

    wc=raster::getData('worldclim', var = 'bio' , res = 10 )
	cor.data(predictors=wc)

---
### *clean*
Função para selecionar os pontos espacialmente únicos e retirar os pontos que estão fora dos limites do *raster*.

clean(coord , predictors)

**Argumentos:**

* coord: matriz ou *dataframe* contendo duas colunas com longitude e latitude em graus decimais nesta ordem.
* predictors: objeto com as variáveis ambientais (*raster*)

**Exemplo:** 

    library(dismo)
	data(acaule)
	pontos = as.data.frame(na.omit(cbind(acaule$lon,acaule$lat)))
	wc = raster::getData('worldclim', var = 'bio' , res = 10 )
	clean(coord = pontos , predictors = wc)

---
### *modelos()*
Função que roda os algoritmos

modelos (coord , k=3 , diretorio = "teste" , plot = T , bc = T , mx = T , GLM = T , RF = T , SVM = T , dm = F , mah = F)

**Argumentos:**

 * coord: matriz ou *dataframe* contendo duas colunas com longitude e latitude em graus decimais nesta ordem. 
 * k: número de partições.
 * diretorio: nome do diretório que serão salvos  os resultados da modelagem.
 * plot: se TRUE (padrão), então plota o modelo *ensemble* final
 * bc, mx, GLM, RF, SVM, dm, mah: algorítmos disponíveis para modelagem. 

**Importante:<p>**
Necessário criar um objeto chamado *predictors* que contém todas as variáveis preditoras, antes de rodar a função.

Exemplo de uso:

    #Roda somente o Bioclim 
	pontos = as.data.frame(na.omit(cbind(acaule$lon,acaule$lat)))
	wc = raster::getData('worldclim', var = 'bio' , res = 10 )
    predictors=wc
	modelos(coord = pontos , diretorio = "solanum")
    
---
### *toKML()*
Função que roda os algoritmos

toKML (modelo, name = "meuKML", zeros = FALSE, open = FALSE) 

**Argumentos:**

 * modelo: *raster* do modelo. 
 * name: nome do arquivo KML a ser gerado.
 * zeros: se for *TRUE*, os valores 0 são mantidos no KML. Caso contrário, são substituídos por NA.
 * open: se *TRUE* o arquivo KML é aberto no *Google Earth*. 


Exemplo de uso:

    library(dismo)
    fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE )
    predictors <- stack(fnames)
    toKML(predictors[[4]],open=T)