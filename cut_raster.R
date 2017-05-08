
cut.raster=function(raster.dir , shape.dir , extension = ".asc" , plot=TRUE , trim=FALSE ,teste = FALSE){
  
  require(raster)
  #original=getwd()
  
  if(dir.exists("Mask_temp")==F){dir.create("Mask_temp")}
  if(dir.exists("Mask_temp2")==F){dir.create("Mask_temp2")}
  if(dir.exists("Cortados")==F){dir.create("Cortados")}
  
  #Definir shape para cortar
  if(missing(shape.dir)){stop("Não selecionou a pasta contendo o shape de corte")
  } else(shape=rgdal::readOGR(list.files(shape.dir,pattern = ".shp",full.names = T)[1]))
  
  #shape= maptools::readShapeSpatial (shape.dir)
  #
  #subsetar
  #MA_CAA=brasil2[brasil2$CD_LEGENDA=="MATA ATL?NTICA"|brasil2$CD_LEGENDA=="CAATINGA",]
  
  #Lendo os rasters
  if(missing(raster.dir)){
    wrdclim <- list.files(pattern=extension, full.names=TRUE )
    predictors <- raster::stack(wrdclim); #predictors
  }else(predictors=raster::stack(list.files(raster.dir,pattern=extension,full.names = TRUE)))
  

  #loop para cortar todos os rasters
  
  #sem trim
  if(trim==F){
    ini=Sys.time()
    for(i in 1:length(names(predictors)))
    {
      ini1=Sys.time()
      crop(predictors[[i]],extent(shape),filename=paste("./Mask_temp", '/',names(predictors)[i],sep="",extension),overwrite=TRUE) 
      croped=list.files("./Mask_temp",pattern=extension, full.names=TRUE )
      mask(raster(croped[i]),shape,filename= paste("./Cortados","/", names(predictors)[i],sep="",extension),overwrite=TRUE)
      
      #if(i!="1"){cat("\r")}
      print(Sys.time())
      cat("\n",paste("Tá indo",i))
      
      fim1= Sys.time()
      cat(paste("\n",round(as.numeric(fim1-ini1),2),units(fim1-ini1)))
      if(i==length(names(predictors))){cat("\n","Acabou!")
        fim= Sys.time()
        fim-ini}
    }
  }
  
  
  #com trim
  if(trim==T){
    ini=Sys.time()
    for(i in 1:length(names(predictors)))
    {
      ini1=Sys.time()
      crop(predictors[[i]],extent(shape),filename=paste(".\\Mask_temp", '\\',names(predictors)[i],sep="",extension),overwrite=TRUE) 
      croped=list.files(".\\Mask_temp",pattern=extension, full.names=TRUE )
      mask(raster(croped[i]),shape,filename= paste(".\\Mask_temp2","\\", names(predictors)[i],sep="",extension),overwrite=TRUE)
      
      writeRaster(
        trim(
          raster(
            list.files(".\\Mask_temp2", pattern=extension, full.names=TRUE )[i]
          )
        )
        ,filename=paste(".\\Cortados","\\", names(predictors)[i],".tif",sep='')
        , format="GTiff", overwrite=TRUE, NAflag=-9999)
      cat("\r")
      print(Sys.time())
      cat("\n",paste("TÃ¡ indo",i))
      fim1= Sys.time()
      cat(paste("\n",round(as.numeric(fim1-ini1),2),units(fim1-ini1)))
      if(i==length(names(predictors))){cat("\n","Acabou!","\n")
        fim= Sys.time()
        fim-ini}
    }
  }
  
  #plotando a primeira variavel cortada
  if(plot==T){
    plot(raster(list.files("./cortados", pattern="", full.names=TRUE )[1],native=TRUE))
    plot(shape,add=T)
  }
  unlink("Mask_temp",recursive = T,force = T)
  unlink("Mask_temp2",recursive = T,force = T)
}
