modelos=function(coord,predictors,k=3,diretorio="teste",pasta,plot=T,
                 bc=T,mx=T,dm=T,GLM=T,RF=T,SVM=T,mah=F){
  original=getwd()
  #escolha da pasta 
  dir.create(paste0("./",diretorio))
  setwd(paste0("./",diretorio))
  
  #instalando pacotes
  packages=c('maps','mapdata',  'maptools', 'dismo', 'rgdal', 'raster', 'jsonlite',"rJava")
  for (p in setdiff(packages, installed.packages()[,"Package"])){
    install.packages(p)
  }
  
  #MaxEnt .jar####
  #baixa e descompacta o maxent java 
  jar <- paste0(system.file(package="dismo"), "/java/maxent.jar")
  if(file.exists(jar)!=T){
    url="http://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download"
    download.file(url, dest="maxent.zip", mode="wb") 
    unzip ("maxent.zip",files="maxent.jar", exdir = system.file("java", package="dismo"))
    unlink("maxent.zip")
  }else(cat("O maxent.jar está na pasta java do dismo\n"))
  
  #Abrindo bibliotecas necessárias####
  library(maps)
  library(mapdata)
  library(maptools)
  library(dismo)
  library(rgdal)
  library(raster)
  
  #-----------#
  # Raster ####
  #-----------#
  
  #Escolhento pasta das variáveis ambientais
  if(exists('pasta')){
    pasta=pasta
  }else(c(pasta="D:/modelagem/Asc",
        warning("A pasta das variáveis ambientais não foi inoformada e foi alterada automaticamente para 'D:/modelagem/Asc'")))
  
  if(exists("predictors")!=T){
    lista <- list.files(pasta,pattern='.asc', full.names=TRUE )
    predictors <- stack(lista)
  }
  
  
  ##------------------##
  #Pontos de ocorrência
  ##------------------##
  
  #Extrair os valores ambientais das localidades onde há registros de ocorrência
  if(exists('pts')){
    pts=pts
    if(dim(pts)[2]==2){pts=pts}else(stop("Verique o número de colunas de planilha com as coordenadas"))
  }else(stop("Não existe objeto com os pontos de ocorrência.","Verifique o nome do objeto"))
  
  source("https://raw.githubusercontent.com/diogosbr/modelagem/master/clean.R")
  pts1=clean(pts,predictors)
  
  #presvals=extract(predictors,pts1)
  
  #--------------#
  # Modelando ####
  #--------------#
  
  aval=as.data.frame(matrix(NA,k*6,7))
  
  backg <- randomPoints(predictors, n=1000, extf = 1.25)
  colnames(backg) = c( 'long' ,  'lat' )
  backvalues=extract(predictors,backg)
  
  group.p <- kfold(pts1, k)
  group.a <- kfold(backg, k)
  
  dir.create("./modelos")
  dir.create("./modelos/bin")
  dir.create("./ensembles")
  dir.create("./final")
  dir.create("./png")
  
  if(bc==T){
    #Bioclim #####
    
    for(i in 1:k){
      cat(c("\r","Começou a partição", i,"Bioclim"))
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != 1, ]
      backg_test <- backg[group.a == 1, ]
      bc <- bioclim(predictors, pres_train)
      e=evaluate( pres_test, backg_test,bc, predictors)
      tr=threshold(e,"spec_sens")
      aval[i,]=threshold(e)
      aval[i,7]="Bioclim"
      bc.mod=predict(predictors,bc)
      writeRaster(bc.mod,paste0("./modelos/","bc_",i,"_con.tif"),format="GTiff",overwrite=T)
      writeRaster(bc.mod>tr,paste0("./modelos/bin/","bc_",i,"_bin.tif"),format="GTiff",overwrite=T)
      png(paste0("./png/",'bc_',i,'con.png'))
      plot(bc.mod,main=paste0("Bioclim con part_",i))
      dev.off()
      png(paste0("./png/",'bc_',i,'bin.png'))
      plot(bc.mod>tr,main=paste0("Bioclim bin part_",i))
      dev.off()
      if(i==k){
        bc.stack=stack(list.files("./modelos",pattern = c("bc",".tif"),full.names = T))
        bc.ens=mean(bc.stack)
        writeRaster(bc.ens,paste0("./ensembles/","bc_",'ensemble ',".tif"),format="GTiff",overwrite=T)
        png(paste0("./png/",'bc_','ensemble ','con.png'))
        plot(bc.ens,main="Bioclim ensemble")
        dev.off()
        cat(c("\r","Terminou! Verifique seus modelos"))
      }
    }
  }
  
  if(mx==T){
    #Maxent #####
    
    for(i in 1:k){
      cat(c("\r","Começou a partição", i,"Maxent"))
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != 1, ]
      backg_test <- backg[group.a == 1, ]
      mx <- maxent(predictors, pres_train)
      e=evaluate( pres_test, backg_test,mx, predictors)
      tr=threshold(e,"spec_sens")
      aval[i+3,]=threshold(e)
      aval[i+3,7]="Maxent"
      mx.mod=predict(predictors,mx)
      writeRaster(mx.mod,paste0("./modelos/","Maxent_",i,"_con.tif"),format="GTiff",overwrite=T)
      writeRaster(mx.mod>tr,paste0("./modelos/bin/","Maxent_",i,"_bin.tif"),format="GTiff",overwrite=T)
      png(paste0("./png/",'Maxent_',i,'con.png'))
      plot(mx.mod,main=paste0("Maxent con part_",i))
      dev.off()
      png(paste0("./png/",'Maxent_',i,'bin.png'))
      plot(mx.mod>tr,main=paste0("Maxent bin part_",i))
      dev.off()
      if(i==k){
        mx.stack=stack(list.files("./modelos",pattern = c("Maxent",".tif"),full.names = T))
        mx.ens=mean(mx.stack)
        writeRaster(mx.ens,paste0("./ensembles/","Maxent_",'ensemble',".tif"),format="GTiff",overwrite=T)
        png(paste0("./png/",'Maxent_','ensemble','con.png'))
        plot(mx.ens,main="Maxent ensemble")
        dev.off()
        
        cat(c("\r","Terminou! Verifique seus modelos"))
      }
    }
  }
  
  if(dm==T){
    #Domain #####
    
    for(i in 1:k){
      cat(c("\r","Começou a partição", i,"Domain"))
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != 1, ]
      backg_test <- backg[group.a == 1, ]
      dm <- domain(predictors, pres_train)
      e=evaluate( pres_test, backg_test,dm, predictors)
      tr=threshold(e,"spec_sens")
      aval[i+6,]=threshold(e)
      aval[i+6,7]="Domain"
      dm.mod=predict(predictors,dm)
      writeRaster(dm.mod,paste0("./modelos/","Domain_",i,"_con.tif"),format="GTiff",overwrite=T)
      writeRaster(dm.mod>tr,paste0("./modelos/bin/","Domain_",i,"_bin.tif"),format="GTiff",overwrite=T)
      png(paste0("./png/",'Domain_',i,'con.png'))
      plot(dm.mod,main=paste0("Domain con part_",i))
      dev.off()
      png(paste0("./png/",'Domain_',i,'bin.png'))
      plot(dm.mod>tr,main=paste0("Domain bin part_",i))
      dev.off()
      if(i==k){
        dm.stack=stack(list.files("./modelos",pattern = c("Domain",".tif"),full.names = T))
        dm.ens=mean(dm.stack)
        writeRaster(dm.ens,paste0("./ensembles/","Domain_",'ensemble',".tif"),format="GTiff",overwrite=T)
        png(paste0("./png/",'domain_','ensemble','con.png'))
        plot(dm.ens,main="Domain ensemble")
        dev.off()
        
        cat(c("\r","Terminou! Verifique seus modelos"))
      }
    }
  }
  
  if(mah==T){
    #Mahalanobis #####
    
    for(i in 1:k){
      cat(c("\r","Começou a partição", i,"Mahalanobis"))
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != 1, ]
      backg_test <- backg[group.a == 1, ]
      mah <- mahal(predictors, pres_train)
      e=evaluate( pres_test, backg_test,mah, predictors)
      tr=threshold(e,"spec_sens")
      mah.mod=predict(predictors,mah)
      writeRaster(mah.mod,paste0("./modelos/","Mahalanobis_",i,"_con.tif"),format="GTiff",overwrite=T)
      writeRaster(mah.mod>tr,paste0("./modelos/bin/","Mahalanobis_",i,"_bin.tif"),format="GTiff",overwrite=T)
      png(paste0("./png/",'Mahalanobis_',i,'con.png'))
      plot(mah.mod,main=paste0("Mahalanobis con part_",i))
      dev.off()
      png(paste0("./png/",'Mahalanobis_',i,'bin.png'))
      plot(mah.mod>tr,main=paste0("Mahalanobis bin part_",i))
      dev.off()
      if(i==k){
        mah.stack=stack(list.files("./modelos",pattern = c("mahalanobis",".tif"),full.names = T))
        mah.ens=mean(mah.stack)
        writeRaster(mah.ens,paste0("./ensembles/","mahalanobis_",'ensemble',".tif"),format="GTiff",overwrite=T)
        png(paste0("./png/",'mahalanobis_','ensemble','con.png'))
        plot(mah.ens,main="Mahalanobis ensemble")
        dev.off()
        
        cat(c("\r","Terminou! Verifique seus modelos"))
      }
    }
  }
  
  if(GLM==T){
    #GLM ####
    for(i in 1:k){
      cat(c("\r","Começou a partição", i,"GLM"))
      
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != 1, ]
      backg_test <- backg[group.a == 1, ]
      
      train <- rbind(pres_train, backg_train)
      pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
      envtrain <- extract(predictors, train)
      envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
      
      GLM <- glm(pa ~.,family = gaussian(link = "identity"), data=envtrain)
      e=evaluate( pres_test, backg_test,GLM, predictors)
      tr=threshold(e,"spec_sens")
      aval[i+9,]=threshold(e)
      aval[i+9,7]="GLM"
      GLM.mod=predict(predictors,GLM)
      writeRaster(GLM.mod,paste0("./modelos/","GLM_",i,"_con.tif"),format="GTiff",overwrite=T)
      writeRaster(GLM.mod>tr,paste0("./modelos/bin/","GLM_",i,"_bin.tif"),format="GTiff",overwrite=T)
      png(paste0("./png/",'GLM_',i,'con.png'))
      plot(GLM.mod,main=paste0("GLM con part_",i))
      dev.off()
      png(paste0("./png/",'GLM_',i,'bin.png'))
      plot(GLM.mod>tr,main=paste0("GLM bin part_",i))
      dev.off()
      if(i==k){
        GLM.stack=stack(list.files("./modelos",pattern = c("GLM",".tif"),full.names = T))
        GLM.ens=mean(GLM.stack)
        writeRaster(GLM.ens,paste0("./ensembles/","GLM_",'ensemble',"bin.tif"),format="GTiff",overwrite=T)
        png(paste0("./png/",'GLM_','ensemble','.png'))
        plot(GLM.ens,main="GLM ensemble")
        dev.off()
        
        cat(c("\r","Terminou! Verifique seus modelos"))
      }
    }
  }
  
  if(RF==T){
    #Random Forest ####
    library(randomForest)
    for(i in 1:k){
      cat(c("\r","Começou a partição", i,"RandomForest"))
      
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != 1, ]
      backg_test <- backg[group.a == 1, ]
      
      train <- rbind(pres_train, backg_train)
      pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
      envtrain <- extract(predictors, train)
      envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
      
      RF <- randomForest(pa ~., data=envtrain)
      e=evaluate( pres_test, backg_test,RF, predictors)
      tr=threshold(e,"spec_sens")
      aval[i+12,]=threshold(e)
      aval[i+12,7]="Random Forest"
      RF.mod=predict(predictors,RF)
      writeRaster(RF.mod,paste0("./modelos/","RF_",i,"_con.tif"),format="GTiff",overwrite=T)
      writeRaster(RF.mod>tr,paste0("./modelos/bin/","RF_",i,"_bin.tif"),format="GTiff",overwrite=T)
      png(paste0("./png/",'RF_',i,'con.png'))
      plot(RF.mod,main=paste0("RF con part_",i))
      dev.off()
      png(paste0("./png/",'RF_',i,'bin.png'))
      plot(RF.mod>tr,main=paste0("RF bin part_",i))
      dev.off()
      if(i==k){
        RF.stack=stack(list.files("./modelos",pattern = c("RF",".tif"),full.names = T))
        RF.ens=mean(RF.stack)
        writeRaster(RF.ens,paste0("./ensembles/","RF_",'ensemble',".tif"),format="GTiff",overwrite=T)
        png(paste0("./png/",'RF_','ensemble','.png'))
        plot(RF.ens,main="Random Forest ensemble")
        dev.off()
        
        cat(c("\r","Terminou! Verifique seus modelos"))
      }
    }
  }
  
  if(SVM==T){
    #SVM ####
    library(kernlab)
    for(i in 1:k){
      cat(c("\r","Começou a partição", i,"SVM"))
      
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != 1, ]
      backg_test <- backg[group.a == 1, ]
      
      train <- rbind(pres_train, backg_train)
      pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
      envtrain <- extract(predictors, train)
      envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
      
      SVM <- ksvm(pa ~., data=envtrain)
      e=evaluate( pres_test, backg_test,SVM, predictors)
      tr=threshold(e,"spec_sens")
      aval[i+15,]=threshold(e)
      aval[i+15,7]="SVM"
      SVM.mod=predict(predictors,SVM)
      writeRaster(SVM.mod,paste0("./modelos/","SVM_",i,"_con.tif"),format="GTiff",overwrite=T)
      writeRaster(SVM.mod>tr,paste0("./modelos/bin/","SVM_",i,"_bin.tif"),format="GTiff",overwrite=T)
      png(paste0("./png/",'SVM_',i,'con.png'))
      plot(SVM.mod,main=paste0("SVM con part_",i))
      dev.off()
      png(paste0("./png/",'SVM_',i,'bin.png'))
      plot(SVM.mod>tr,main=paste0("SVM bin part_",i))
      dev.off()
      if(i==k){
        SVM.stack=stack(list.files("./modelos",pattern = c("SVM",".tif"),full.names = T))
        SVM.ens=mean(SVM.stack)
        writeRaster(SVM.ens,paste0("./ensembles/","SVM_",'ensemble',".tif"),format="GTiff",overwrite=T)
        png(paste0("./png/",'SVM_','ensemble','.png'))
        plot(SVM.ens,main="Random Forest ensemble")
        dev.off()
        
        cat(c("\r","Terminou! Verifique seus modelos"))
      }
    }
  }
  
  
  #Ensemble final ####
  
  final=stack(list.files("./ensembles",pattern = ".tif",full.names = T))
  #final.ens=mean(final)
  
  bc.cut=final[[grep("bc",names(final))]]
  dm.cut=final[[grep("Domain",names(final))]]
  GLM.cut=final[[grep("GLM",names(final))]]
  mx.cut=final[[grep("Maxent",names(final))]]
  RF.cut=final[[grep("RF",names(final))]]
  SVM.cut=final[[grep("SVM",names(final))]]
  
  #recorte com TSSth
  values(bc.cut)[values(bc.cut)<mean(aval[grep("Bioclim",aval[,7]),2])]=0
  #padronizando de 0 a 1
  values(bc.cut)=values(bc.cut)/bc.cut@data@max
  
  #recorte com TSSth
  values(dm.cut)[values(dm.cut)<mean(aval[grep("Domain",aval[,7]),2])]=0
  #padronizando de 0 a 1
  values(dm.cut)=values(dm.cut)/dm.cut@data@max
  
  #recorte com TSSth
  values(GLM.cut)[values(GLM.cut)<mean(aval[grep("GLM",aval[,7]),2])]=0
  #padronizando de 0 a 1
  values(GLM.cut)=values(GLM.cut)/GLM.cut@data@max
  
  #recorte com TSSth
  values(mx.cut)[values(mx.cut)<mean(aval[grep("Maxent",aval[,7]),2])]=0
  #padronizando de 0 a 1
  values(mx.cut)=values(mx.cut)/mx.cut@data@max
  
  #recorte com TSSth
  values(RF.cut)[values(RF.cut)<mean(aval[grep("Random Forest",aval[,7]),2])]=0
  #padronizando de 0 a 1
  values(RF.cut)=values(RF.cut)/RF.cut@data@max
  
  #recorte com TSSth
  values(SVM.cut)[values(SVM.cut)<mean(aval[grep("SVM",aval[,7]),2])]=0
  #padronizando de 0 a 1
  values(SVM.cut)=values(SVM.cut)/SVM.cut@data@max
  
  mm=mean(bc.cut,
          dm.cut,
          GLM.cut,
          mx.cut,
          RF.cut,
          SVM.cut)
  
  values(mm)=values(mm)/mm@data@max
  
  writeRaster(mm,paste0("./final/","Geral",'ensemble',".tif"),format="GTiff",overwrite=T)
  write.table(aval,"Aval.csv",sep=";",dec=".")
  
  png(paste0("./png/",'Geral_','ensemble','.png'))
  plot(mm,main="Geral ensemble")
  dev.off()
  
  png(paste0("./png/",'Geral_','ensemble',"pontos",'.png'))
  plot(mm,main="Geral ensemble")
  points(pts1)
  dev.off()
  
  ### Modelagem até aqui ###
  
  if(plot==T){
    plot(mm)
    points(pts1,pch=16)
  }
  setwd=original
  cat('Veja o modelo final')
}



