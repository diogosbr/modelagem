# Modelagem de nicho
Funções para auxiliar a modelagem de distribuição de nicho

## **Funções** <p>
### *cut.raster*
Função para recortar as variáveis ambientais

### *cor.data*
Função para verificar a correlação entre as variáveis ambientais

### *clean*
Função para selecionar os pontos espacialmente únicos e retirar os NAs

### *modelos*
Função que roda os algoritmos

---
# Como usar
###cut.raster (raster.dir,shape.dir,extension=".asc",plot=F,trim=F)

Argumentos:

 * raster.dir: diretório que contém os rasters a serem cortados. Se não for informado, então vai procurar os rasters na pasta local.
 * shape.dir: diretório que stá o shape que será usado como máscara para cortar os rasters. Obrigatório.
 * extension: conjuntos de caracteres com a extensão dos ratsers. ".asc" é o padrão.
 * plot: plota os rasters. TRUE é o padrão.
 * trim: se for TRUE os NAs gerados após o corte dos rasters são removidos. FALSE é o padrão.

Exemplo de uso:


###modelos (coord,k=3,diretorio="teste",plot=T,
                 bc=T,mx=T,GLM=T,RF=T,SVM=T,dm=F,mah=F)

Argumentos:

 * coord: matriz ou dataframe contendo duas colunas com longitude e latitude em graus decimais nesta ordem. 
 * k: número de partições.
 * diretorio: nome do diretório que serão salvos  os resultados da modelagem.
 * plot: se TRUE (padrão), então plota o modelo *ensemble* final
 * bc, mx, GLM, RF, SVM, dm, mah: algorítmos disponíveis para modelagem. 

Importante:
necessário criar um objeto chamado *predictors* que contém todas as variáveis preditoras

Exemplo de uso:


And here's some code! :+1:

```Rscript
modelos()
```

