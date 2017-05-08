# Modelagem de nicho
Funções para auxiliar a modelagem de distribuição de nicho

##Índice
* cut.raster
* cor.data
* clean
* modelos 

## **Funções** 
### *cut.raster()*
Função para recortar as variáveis ambientais (*raster*) a partir de um *shapefile*

cut.raster(raster.dir , shape.dir , extension = ".asc" , plot = TRUE , trim = FALSE)

**Argumentos:**

 * raster.dir: diretório que contém os rasters a serem cortados. Se não for informado, então vai procurar os rasters na pasta local.
 * shape.dir: diretório que stá o shape que será usado como máscara para cortar os rasters. Obrigatório.
 * extension: conjuntos de caracteres com a extensão dos *ratsers*. ".asc" é o padrão.
 * plot: plota os rasters. TRUE é o padrão.
 * trim: se for TRUE os NAs gerados após o corte dos rasters são removidos. FALSE é o padrão.

**Exemplo:**

    wc=raster::getData('worldclim', var = 'bio' , res = 10 )
    for(x in c("shx","shp",'cpg','dbf','csv','prj'))
{download.file(paste0("https://github.com/diogosbr/modelagem/blob/master/Dataset/shape/BRA_adm0.",x),destfile = paste0("br.",x),quiet = T)}
    cut.raster(raster.dir = "wc10" , shape.dir = "./" , extension = ".tif")

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
    
