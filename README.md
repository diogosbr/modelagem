# Modelagem de nicho
Funções para auxiliar a modelagem de distribuição de nicho

##Funções disponíveis
* clean()
* cor.data()
* cut.raster()
* dist_euc()
* modelos()
* toKML()

No arquivo *testing.R* tem exemplos os exemplos para testar as funções. 

## **Funções** 
### *clean()*
Função para selecionar os pontos espacialmente únicos e retirar os pontos que estão fora dos limites do *raster*.

clean(coord , predictors)

**Argumentos:**

* coord: matriz ou *dataframe* contendo duas colunas com longitude e latitude em graus decimais nesta ordem.
* predictors: objeto com as variáveis ambientais (*raster*)

**Exemplo:** 

    library(dismo)
	data(acaule)
	
	filenames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
              pattern='grd', full.names=TRUE )	
	wc <- stack(filenames[1:8])
	pontos = as.data.frame(na.omit(cbind(acaule$lon,acaule$lat)))
	clean(coord = pontos , predictors = wc)

---
### *cor.data()*
Função para verificar a correlação entre as variáveis ambientais

cor.data(predictors , plot = TRUE, rep = 1000)

**Argumentos:**

* predictors: objeto com as variáveis ambientais (*raster*)
* plot: se for TRUE (padrão), plota gráfico e valores de correlação

**Exemplo:** 

    filenames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
              pattern='grd', full.names=TRUE )
	wc <- stack(filenames[1:8])
	cor.data(predictors=wc)


---
### *cut.raster()*
Função para recortar as variáveis ambientais (*raster*) a partir de um *shapefile*

cut.raster(raster.dir , shape.dir , extension = ".asc" , plot = TRUE , trim = FALSE)

**Argumentos:**

 * raster.dir: diretório que contém os rasters a serem cortados. Se não for informado, então vai procurar os rasters na pasta local.
 * shape.dir: diretório que está o shape que será usado como máscara para cortar os rasters. Obrigatório.
 * extension: conjuntos de caracteres com a extensão dos *rasters*. ".asc" é o padrão.
 * plot: plota os rasters. TRUE é o padrão.
 * trim: se for TRUE os NAs gerados após o corte dos rasters são removidos. FALSE é o padrão.

**Exemplo:**

    filenames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
              pattern='grd', full.names=TRUE )
	
	wc <- stack(filenames[1:8])
    download.file("http://biogeo.ucdavis.edu/data/gadm2.8/shp/BRA_adm_shp.zip", destfile = "bra.zip") 
    unzip("bra.zip",exdir = "./brasil")
    unlink("bra.zip")
    cut.raster(raster.dir = "wc10" , shape.dir = "brasil" , extension = ".bil")

---
### *dist_euc()*
Função para gerar modelo de distância ambiental, utilizando a distância euclidiana. Pode ser utilizado apenas um ponto. 

Função que calcula modelos com base na distância euclidiana.

dist_euc(occ, env, method = "mean", decostand.method = "standardize", suitability = FALSE)

**Argumentos:**

 * occ: 
 * env:
 * method:
 * decostand.method:
 * suitability:

Exemplo de uso:

    #Roda somente o Bioclim 
	library(dismo)
	
	# file with presence points
	occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
	occ <- read.table(occurence, header=TRUE, sep=',')[,-1]

	fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
              pattern='grd', full.names=TRUE )	
	result <- dist_euc(occ, predictors, method = "mean", suitability = FALSE,
                   decostand.method = "standardize")
	predictors <- stack(fnames)
	plot(results)


---
### *modelos()*
Função que roda os algoritmos
modelos(coord, abio, k = 3, diretorio = "teste", plot = T, bc = T, mx = F, GLM = F, RF = F, SVM = F, dm = F, mah = F, proj, buffer, geo.filt = T, br, mod = 'before', tss)


**Argumentos:**

 * coord: matriz ou *dataframe* contendo duas colunas com longitude e latitude em graus decimais nesta ordem. 
 * abio: objeto com as variáveis preditoras empilhadas.
 * k: número de partições.
 * diretorio: nome do diretório que serão salvos  os resultados da modelagem.
 * plot: se TRUE (padrão), então plota o modelo *ensemble* final
 * bc, mx, GLM, RF, SVM, dm, mah: algoritmos disponíveis para modelagem.
 * proj: se não for informado, então os modelos são projetados utilizando os mesmos *rasters* informados em *abio*.
 * buffer: pode não ser informado e neste caso, os pontos de pseudo ausência são criados dentro dos limites de *abio*. Pode ser informado "mean" e são criados os pontos dentro de um buffer com a distância média entre os pontos. Caso seja "max", então o buffer criado é a distância máxima entre os pontos de ocorrência. Note que neste último caso, pode acontecer de a área de criação das pseudo ausências serem iguais ou superior às informadas em *abio*.
 * geo.filt: se TRUE (*padrão*), aplica um filtro geográfico retirando pontos muito próximos. Atualmente, são mantidos pontos de ocorrência com mais de 20 Km (ou 10 min).
 * br: argumento obrigatório (por enquanto), com o shape do Brasil.
 * mod:pode ser "before" e neste caso cada partição é cortada pelo ser próprio TSSth e depois são feitos os ensembles. Se for "after", então é feito o ensemble de cada algoritmo e então cortado pelo TSSth médio das partições do algoritmo. 
 * tss: númerico. Serão utilizados apenas os modelos com TSS maior do que o informado para gerar o ensemble geral.  Se não for informado, são utilizados todos os modelos de todas as partições para gerar o ensemble geral.


Exemplo de uso:

    #Roda somente o Bioclim 
	
	library(dismo)
	data(acaule)
	pontos = as.data.frame(na.omit(cbind(acaule$lon,acaule$lat)))
	filenames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
              pattern='grd', full.names=TRUE )
	
	wc <- stack(filenames[1:8])
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
    filenames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), pattern='grd', full.names=TRUE )
    predictors <- stack(filenames)
    toKML(predictors[[4]],open=T)