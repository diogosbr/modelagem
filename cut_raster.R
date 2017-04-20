
cut=function(raster.dir,shape.dir,ext=".asc",plot=F,trim=F){
  
  require(raster)
  #original=getwd()
  
  if(dir.exists("Mask_temp")==F){dir.create("Mask_temp")}
  if(dir.exists("Mask_temp2")==F){dir.create("Mask_temp2")}
  if(dir.exists("Cortados")==F){dir.create("Cortados")}
  
  #Definir shape para cortar
  if(missing(shape.dir)){stop("N�o selecionou o shape de corte")
    } else(shape= maptools::readShapeSpatial (shape.dir))
  
  #subsetar
  #MA_CAA=brasil2[brasil2$CD_LEGENDA=="MATA ATL?NTICA"|brasil2$CD_LEGENDA=="CAATINGA",]
  
  #Lendo os rasters
  if(missing(raster.dir)){
    wrdclim <- list.files(pattern=ext, full.names=TRUE )
    predictors <- raster::stack(wrdclim); #predictors
  }else(predictors=raster::stack(list.files(raster.dir,pattern=ext,full.names = TRUE)))
  
  #plotando a primeira variavel 
  if(plot=T){
    plot(predictors[[1]])
    plot(shape,add=T)
  }

  #loop para cortar todos os rasters
  
  #sem trim
  if(trim=F){
    ini=Sys.time()
    for(i in 1:length(names(predictors)))
    {
      ini1=Sys.time()
      crop(predictors[[i]],extent(shape),filename=paste("./Mask_temp", '/',names(predictors)[i],sep="",".tif"),overwrite=TRUE) 
      croped=list.files("./Mask_temp",pattern='.tif', full.names=TRUE )
      mask(raster(croped[i]),shape,filename= paste("./Cortados","/", names(predictors)[i],sep="",".tif"),overwrite=TRUE)
      
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
      crop(predictors[[i]],extent(shape),filename=paste(".\\Mask_temp", '\\',names(predictors)[i],sep="",".tif"),overwrite=TRUE) 
      croped=list.files(".\\Mask_temp",pattern='.tif', full.names=TRUE )
      mask(raster(croped[i]),shape,filename= paste(".\\Mask_temp2","\\", names(predictors)[i],sep="",".tif"),overwrite=TRUE)
      
      writeRaster(
        trim(
          raster(
            list.files(".\\Mask_temp2", pattern='.tif', full.names=TRUE )[i]
          )
        )
        ,filename=paste(".\\Cortados","\\", names(predictors)[i],".tif",sep='')
        , format="GTiff", overwrite=TRUE, NAflag=-9999)
      cat("\r")
      print(Sys.time())
      cat("\n",paste("Tá indo",i))
      fim1= Sys.time()
      cat(paste("\n",round(as.numeric(fim1-ini1),2),units(fim1-ini1)))
      if(i==length(names(predictors))){cat("\n","Acabou!","\n")
        fim= Sys.time()
        fim-ini}
    }
  }
  
}
