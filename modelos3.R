modelos = function(coord, abio, k = 3, diretorio = "teste", plot = T, bc = T, mx = F, GLM = F, RF = F, 
                   SVM = F, dm = F, mah = F, proj, buffer, geo.filt = T, br, mod = 'before', tss) {
  
  if(missing(abio)){stop("Informe as variÃ¡veis abióticas")}else(predictors=abio)
  original = getwd()
  # escolha da pasta
  dir.create(paste0("./", diretorio))
  setwd(paste0("./", diretorio))
  
  # instalando pacotes
  packages = c("dismo", "rgdal", "raster", "randomForest", "kernlab")
  for (p in setdiff(packages, installed.packages()[, "Package"])) {
    install.packages(p, dependencies = T)
  }
  
  # MaxEnt .jar#### baixa e descompacta o maxent java
  jar <- paste0(system.file(package = "dismo"), "/java/maxent.jar")
  if (file.exists(jar) != T) {
    url = "http://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download"
    download.file(url, dest = "maxent.zip", mode = "wb")
    unzip("maxent.zip", files = "maxent.jar", exdir = system.file("java", package = "dismo"))
    unlink("maxent.zip")
    warning("Maxent foi colocado no diretório")
  } 
  
  # Abrindo bibliotecas necessárias####
  library(dismo)
  library(rgdal)
  
  ##--------------------------##
  # Pontos de ocorrência####
  ##------------------------##
  
  # Extrair os valores ambientais das localidades onde há registros de ocorrência
  if (exists("coord")) {
    pts = coord
    if (dim(pts)[2] == 2) {
      pts = pts
    } else (stop("Verique o número de colunas de planilha com as coordenadas"))
  } else (stop("Não existe objeto com os pontos de ocorrência.", "Verifique o nome do objeto"))
  
  source("https://raw.githubusercontent.com/diogosbr/modelagem/master/clean.R")
  pts1 = clean(pts, predictors = predictors)
  names(pts1) = c("long", "lat")
  aa=data.frame(Originais=dim(pts)[1], Unicos=dim(pts1)[1], Retirados=dim(pts)[1]-dim(pts1)[1])
  
  #Filtros geográficos####
  if(geo.filt==T){
    res=0.1666667#10min - 20km
    r=raster(extent(range(pts1[,1]), range(pts1[,2])) + res)
    res(r)=res
    pts1=gridSample(pts1,r, n=1)
    cat(paste0(dim(pts1)[1], ' após o filtro geográfico de 20Km'))
    aa$geo.filt=dim(pts1)[1]
  }
  write.table(aa, "Nocc.csv", row.names = F, sep = ";")
  
  #--------------#
  # Modelando #####
  #--------------#
  
  #criando tabela de saída para armazenar valores de desempenho dos modelos
  aval = as.data.frame(matrix(NA, k * 7, 11))
  
  #Buffer####
  if( exists(buffer) ){
    pts2=pts1
    names(pts2)=c("lon",'lat')
    coordinates(pts2) <- ~lon + lat
    if(buffer=="mean"){
      dist.buf <- mean(spDists(x = pts1, longlat = FALSE, segments = TRUE))
    } 
    if(missing(buffer)){dist.buf <- max(spDists(x = pts1, longlat = FALSE, segments = TRUE))}
    
    buffer <- raster::buffer(pts2, width = dist.buf, dissolve = TRUE)
    buffer <- SpatialPolygonsDataFrame(buffer, data = as.data.frame(buffer@plotOrder), 
                                       match.ID = FALSE)
    crs(buffer) <- crs(predictors)
    crs(br) = crs(predictors)
    buffer=rgeos::gIntersection(br, buffer, byid = T)
    backg = spsample(buffer, 1000, type="random")
    backg = as.data.frame(backg)
    rm(buffer)
    rm(pts2)
  } else(backg <- randomPoints(predictors, n = 1000, extf = 1))
  
  colnames(backg) = c("long", "lat")
  backvalues = extract(predictors, backg)
  
  group.p <- kfold(pts1, k)
  group.a <- kfold(backg, k)
  
  dir.create("./modelos")
  dir.create("./modelos/bin")
  dir.create("./ensembles")
  dir.create("./final")
  dir.create("./png")
  #dir.create("./temporario")
  
  if (bc == T) {
    # Bioclim #####
    
    for (i in 1:k) {
      cat(c("\n", "Começou a partição", i, "Bioclim"))
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != i, ]
      backg_test <- backg[group.a == i, ]
      bc <- bioclim(predictors, pres_train)
      e = evaluate(pres_test, backg_test, bc, predictors)
      tr = e@t[which.max(e@TPR + e@TNR)]
      TSS.calc=max(e@TPR + e@TNR) - 1
      aval[i, ] = threshold(e)
      aval[i, 7] = "Bioclim"
      aval[i, 8] = e@auc
      aval[i, 9] = max(e@TPR + e@TNR) - 1
      aval[i, 10] = tr
      aval[i,11] = i
      if(missing(proj)){bc.mod = predict(predictors, bc)
      }else(bc.mod = predict(proj, bc))
      if(mod=="before"){values(bc.mod)[values(bc.mod) < tr] = 0}
      if(missing(tss)){
        writeRaster(bc.mod, paste0("./modelos/", "bc_", i, "_con.tif"), format = "GTiff", 
                    overwrite = T)
        writeRaster(bc.mod > tr, paste0("./modelos/bin/", "bc_", i, "_bin.tif"), format = "GTiff", 
                    overwrite = T)
        png(paste0("./png/", "bc_", i, "con.png"))
        plot(bc.mod, main = paste0("Bioclim con part_", i))
        dev.off()
        png(paste0("./png/", "bc_", i, "bin.png"))
        plot(bc.mod > tr, main = paste0("Bioclim bin part_", i))
        dev.off()
      }
      if(missing(tss)==FALSE){
        if(TSS.calc>tss){
          writeRaster(bc.mod, paste0("./modelos/", "bc_", i, "_con.tif"), format = "GTiff", 
                      overwrite = T)
          writeRaster(bc.mod > tr, paste0("./modelos/bin/", "bc_", i, "_bin.tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "bc_", i, "con.png"))
          plot(bc.mod, main = paste0("Bioclim con part_", i))
          dev.off()
          png(paste0("./png/", "bc_", i, "bin.png"))
          plot(bc.mod > tr, main = paste0("Bioclim bin part_", i))
          dev.off()
        }
      }
      #Ensemble do algoritmo
      if (i == k) {
        if(length(list.files("./modelos", pattern = c("bc", ".tif")))!=0){
          bc.stack = stack(list.files("./modelos", pattern = c("bc", ".tif"), full.names = T))
          bc.ens = mean(bc.stack)
          
          # padronizando de 0 a 1
          values(bc.ens) = values(bc.ens)/bc.ens@data@max
          
          # recorte com TSSth
          if(mod=="after"){values(bc.ens)[values(bc.ens) < mean(aval[grep("Bioclim", aval[, 7]), 2])] = 0}
          writeRaster(bc.ens, paste0("./ensembles/", "bc_", "ensemble ", ".tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "bc_", "ensemble ", "con.png"))
          plot(bc.ens, main = "Bioclim ensemble")
          dev.off()
          rm(bc.ens)
        }
        cat(c("\n", "Terminou Bioclim"))
      }
    }
  }
  
  if (mx == T) {
    # Maxent #####
    
    for (i in 1:k) {
      cat(c("\n", "Começou a partição", i, "Maxent"))
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != i, ]
      backg_test <- backg[group.a == i, ]
      mx <- maxent(predictors, pres_train)
      e = evaluate(pres_test, backg_test, mx, predictors)
      tr = e@t[which.max(e@TPR + e@TNR)]
      TSS.calc=max(e@TPR + e@TNR) - 1
      aval[i + 3, ] = threshold(e)
      aval[i + 3, 7] = "Maxent"
      aval[i + 3, 8] = e@auc
      aval[i + 3, 9] = max(e@TPR + e@TNR) - 1
      aval[i + 3, 10] = tr
      aval[i + 3,11] = i
      if(missing(proj)){mx.mod = predict(predictors, mx)
      }else(mx.mod = predict(proj, mx))
      if(mod=="before"){values(mx.mod)[values(mx.mod) < tr] = 0 }
      if(missing(tss)){
        writeRaster(mx.mod, paste0("./modelos/", "Maxent_", i, "_con.tif"), format = "GTiff", 
                    overwrite = T)
        writeRaster(mx.mod > tr, paste0("./modelos/bin/", "Maxent_", i, "_bin.tif"), format = "GTiff", 
                    overwrite = T)
        png(paste0("./png/", "Maxent_", i, "con.png"))
        plot(mx.mod, main = paste0("Maxent con part_", i))
        dev.off()
        png(paste0("./png/", "Maxent_", i, "bin.png"))
        plot(mx.mod > tr, main = paste0("Maxent bin part_", i))
        dev.off()
      }
      if(missing(tss)==FALSE){
        if(TSS.calc>tss){
          writeRaster(mx.mod, paste0("./modelos/", "Maxent_", i, "_con.tif"), format = "GTiff", 
                      overwrite = T)
          writeRaster(mx.mod > tr, paste0("./modelos/bin/", "Maxent_", i, "_bin.tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "Maxent_", i, "con.png"))
          plot(mx.mod, main = paste0("Maxent con part_", i))
          dev.off()
          png(paste0("./png/", "Maxent_", i, "bin.png"))
          plot(mx.mod > tr, main = paste0("Maxent bin part_", i))
          dev.off()
        }
      }
      
      if (i == k) {
        if(length(list.files("./modelos", pattern = c("Maxent", ".tif"), full.names = T))!=0){
          mx.stack = stack(list.files("./modelos", pattern = c("Maxent", ".tif"), full.names = T))
          mx.ens = mean(mx.stack)
          
          # padronizando de 0 a 1
          values(mx.ens) = values(mx.ens)/mx.ens@data@max
          
          # recorte com TSSth
          if(mod=="after"){values(mx.ens)[values(mx.ens) < mean(aval[grep("Maxent", aval[, 7]), 2])] = 0}
          writeRaster(mx.ens, paste0("./ensembles/", "Maxent_", "ensemble", ".tif"), 
                      format = "GTiff", overwrite = T)
          png(paste0("./png/", "Maxent_", "ensemble", "con.png"))
          plot(mx.ens, main = "Maxent ensemble")
          dev.off()
        }
        
        cat(c("\n", "Terminou Maxent"))
      }
    }
  }
  
  if (dm == T) {
    # Domain #####
    for (i in 1:k) {
      cat(c("\n", "Começou a partição", i, "Domain"))
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != i, ]
      backg_test <- backg[group.a == i, ]
      dm <- domain(predictors, pres_train)
      e = evaluate(pres_test, backg_test, dm, predictors)
      tr = e@t[which.max(e@TPR + e@TNR)]
      TSS.calc=max(e@TPR + e@TNR) - 1
      aval[i+6, ] = threshold(e)
      aval[i+6, 7] = "Domain"
      aval[i+6, 8] = e@auc
      aval[i+6, 9] = max(e@TPR + e@TNR) - 1
      aval[i+6, 10] = tr
      aval[i+6,11] = i
      if(missing(proj)){dm.mod = predict(predictors, dm)
      }else(dm.mod = predict(proj, dm))
      if(mod=="before"){values(dm.mod)[values(dm.mod) < tr] = 0}
      if(missing(tss)){
        writeRaster(dm.mod, paste0("./modelos/", "dm_", i, "_con.tif"), format = "GTiff", 
                    overwrite = T)
        writeRaster(dm.mod > tr, paste0("./modelos/bin/", "dm_", i, "_bin.tif"), format = "GTiff", 
                    overwrite = T)
        png(paste0("./png/", "dm_", i, "con.png"))
        plot(dm.mod, main = paste0("Domain con part_", i))
        dev.off()
        png(paste0("./png/", "dm_", i, "bin.png"))
        plot(dm.mod > tr, main = paste0("Domain bin part_", i))
        dev.off()
      }
      if(missing(tss)==FALSE){
        if(TSS.calc>tss){
          writeRaster(dm.mod, paste0("./modelos/", "dm_", i, "_con.tif"), format = "GTiff", 
                      overwrite = T)
          writeRaster(dm.mod > tr, paste0("./modelos/bin/", "dm_", i, "_bin.tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "dm_", i, "con.png"))
          plot(dm.mod, main = paste0("Domain con part_", i))
          dev.off()
          png(paste0("./png/", "dm_", i, "bin.png"))
          plot(dm.mod > tr, main = paste0("Domain bin part_", i))
          dev.off()
        }
      }
      #Ensemble do algoritmo
      if (i == k) {
        if(length(list.files("./modelos", pattern = c("dm", ".tif")))!=0){
          dm.stack = stack(list.files("./modelos", pattern = c("dm", ".tif"), full.names = T))
          dm.ens = mean(dm.stack)
          
          # padronizando de 0 a 1
          values(dm.ens) = values(dm.ens)/dm.ens@data@max
          
          # recorte com TSSth
          if(mod=="after"){values(dm.ens)[values(dm.ens) < mean(aval[grep("Domain", aval[, 7]), 2])] = 0}
          writeRaster(dm.ens, paste0("./ensembles/", "dm_", "ensemble ", ".tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "dm_", "ensemble ", "con.png"))
          plot(dm.ens, main = "Domain ensemble")
          dev.off()
          rm(dm.ens)
        }
        cat(c("\n", "Terminou Domain"))
      }
    }
  }
  
  if (mah == T) {
    # Mahalanobis #####
    
    for (i in 1:k) {
      cat(c("\n", "Começou a partição", i, "Mahalanobis"))
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != i, ]
      backg_test <- backg[group.a == i, ]
      mah <- mahal(predictors, pres_train)
      e = evaluate(pres_test, backg_test, mah, predictors)
      tr = e@t[which.max(e@TPR + e@TNR)]
      TSS.calc=max(e@TPR + e@TNR) - 1
      aval[i+18, ] = threshold(e)
      aval[i+18, 7] = "Mahalanobis"
      aval[i+18, 8] = e@auc
      aval[i+18, 9] = max(e@TPR + e@TNR) - 1
      aval[i+18, 10] = tr
      aval[i+18,11] = i
      if(missing(proj)){mah.mod = predict(predictors, mah)
      }else(mah.mod = predict(proj, mah))
      if(mod=="before"){values(mah.mod)[values(mah.mod) < tr] = 0}
      if(missing(tss)){
        writeRaster(mah.mod, paste0("./modelos/", "mah_", i, "_con.tif"), format = "GTiff", 
                    overwrite = T)
        writeRaster(mah.mod > tr, paste0("./modelos/bin/", "mah_", i, "_bin.tif"), format = "GTiff", 
                    overwrite = T)
        png(paste0("./png/", "mah_", i, "con.png"))
        plot(mah.mod, main = paste0("Mahalanobis con part_", i))
        dev.off()
        png(paste0("./png/", "mah_", i, "bin.png"))
        plot(mah.mod > tr, main = paste0("Mahalanobis bin part_", i))
        dev.off()
      }
      if(missing(tss)==FALSE){
        if(TSS.calc>tss){
          writeRaster(mah.mod, paste0("./modelos/", "mah_", i, "_con.tif"), format = "GTiff", 
                      overwrite = T)
          writeRaster(mah.mod > tr, paste0("./modelos/bin/", "mah_", i, "_bin.tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "mah_", i, "con.png"))
          plot(mah.mod, main = paste0("Mahalanobis con part_", i))
          dev.off()
          png(paste0("./png/", "mah_", i, "bin.png"))
          plot(mah.mod > tr, main = paste0("Mahalanobis bin part_", i))
          dev.off()
        }
      }
      #Ensemble do algoritmo
      if (i == k) {
        if(length(list.files("./modelos", pattern = c("mah", ".tif")))!=0){
          mah.stack = stack(list.files("./modelos", pattern = c("mah", ".tif"), full.names = T))
          mah.ens = mean(mah.stack)
          
          # padronizando de 0 a 1
          values(mah.ens) = values(mah.ens)/mah.ens@data@max
          
          # recorte com TSSth
          if(mod=="after"){values(mah.ens)[values(mah.ens) < mean(aval[grep("Mahalanobis", aval[, 7]), 2])] = 0}
          writeRaster(mah.ens, paste0("./ensembles/", "mah_", "ensemble ", ".tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "mah_", "ensemble ", "con.png"))
          plot(mah.ens, main = "Mahalanobis ensemble")
          dev.off()
          rm(mah.ens)
        }
        cat(c("\n", "Terminou Mahalanobis"))
      }
    }
  }
  
  if (GLM == T) {
    # GLM ####
    for (i in 1:k) {
      cat(c("\n", "Começou a partição", i, "GLM"))
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != i, ]
      backg_test <- backg[group.a == i, ]
      
      train <- rbind(pres_train, backg_train)
      pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
      envtrain <- extract(predictors, train)
      envtrain <- data.frame(na.omit(cbind(pa = pb_train, envtrain)))
      
      envtest_p=envtrain[envtrain[,1]==1,]
      envtest_a=envtrain[envtrain[,1]==0,]
      
      #adaptado model-R
      null.model <- glm(pa ~ 1, data = envtrain, family = "binomial")
      full.model <- glm(pa ~ ., data = envtrain, family = "binomial")
      GLM <- step(object = null.model, scope = formula(full.model), direction = "both", 
                  trace = F)
      e <- evaluate(envtest_p, envtest_p, model = GLM, type = "response") 
      #e <- evaluate(pres_test, backg_test, predictors, type = "response") 
      
      #Exemplo do dismo
      #GLM <- glm(pa ~ ., family = gaussian(link = "identity"), data = envtrain)
      #e = evaluate(pres_test, backg_test, GLM, predictors)
      
      tr = e@t[which.max(e@TPR + e@TNR)]
      TSS.calc=max(e@TPR + e@TNR) - 1
      aval[i+9, ] = threshold(e)
      aval[i+9, 7] = "GLM"
      aval[i+9, 8] = e@auc
      aval[i+9, 9] = max(e@TPR + e@TNR) - 1
      aval[i+9, 10] = tr
      aval[i+9,11] = i
      if(missing(proj)){GLM.mod = predict(predictors, GLM)
      }else(GLM.mod = predict(proj, GLM))
      if(mod=="before"){values(GLM.mod)[values(GLM.mod) < tr] = 0}
      if(missing(tss)){
        writeRaster(GLM.mod, paste0("./modelos/", "GLM_", i, "_con.tif"), format = "GTiff", 
                    overwrite = T)
        writeRaster(GLM.mod > tr, paste0("./modelos/bin/", "GLM_", i, "_bin.tif"), format = "GTiff", 
                    overwrite = T)
        png(paste0("./png/", "GLM_", i, "con.png"))
        plot(GLM.mod, main = paste0("GLM con part_", i))
        dev.off()
        png(paste0("./png/", "GLM_", i, "bin.png"))
        plot(GLM.mod > tr, main = paste0("GLM bin part_", i))
        dev.off()
      }
      if(missing(tss)==FALSE){
        if(TSS.calc>tss){
          writeRaster(GLM.mod, paste0("./modelos/", "GLM_", i, "_con.tif"), format = "GTiff", 
                      overwrite = T)
          writeRaster(GLM.mod > tr, paste0("./modelos/bin/", "GLM_", i, "_bin.tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "GLM_", i, "con.png"))
          plot(GLM.mod, main = paste0("GLM con part_", i))
          dev.off()
          png(paste0("./png/", "GLM_", i, "bin.png"))
          plot(GLM.mod > tr, main = paste0("GLM bin part_", i))
          dev.off()
        }
      }
      #Ensemble do algoritmo
      if (i == k) {
        if(length(list.files("./modelos", pattern = c("GLM", ".tif")))!=0){
          GLM.stack = stack(list.files("./modelos", pattern = c("GLM", ".tif"), full.names = T))
          GLM.ens = mean(GLM.stack)
          
          # padronizando de 0 a 1
          values(GLM.ens) = values(GLM.ens)/GLM.ens@data@max
          
          # recorte com TSSth
          if(mod=="after"){values(GLM.ens)[values(GLM.ens) < mean(aval[grep("GLM", aval[, 7]), 2])] = 0}
          writeRaster(GLM.ens, paste0("./ensembles/", "GLM_", "ensemble ", ".tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "GLM_", "ensemble ", "con.png"))
          plot(GLM.ens, main = "GLM ensemble")
          dev.off()
          rm(GLM.ens)
        }
        cat(c("\n", "Terminou GLM"))
      }
    }
  }
  
  if (RF == T) {
    # Random Forest ####
    library(randomForest)
    for (i in 1:k) {
      cat(c("\n", "Começou a partição", i, "RandomForest"))
      
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != i, ]
      backg_test <- backg[group.a == i, ]
      
      train <- rbind(pres_train, backg_train)
      pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
      envtrain <- extract(predictors, train)
      envtrain <- data.frame(na.omit(cbind(pa = pb_train, envtrain)))
      
      RF <- randomForest(pa ~ ., data = envtrain)
      e = evaluate(pres_test, backg_test, RF, predictors)
      tr = e@t[which.max(e@TPR + e@TNR)]
      TSS.calc=max(e@TPR + e@TNR) - 1
      aval[i + 12, ] = threshold(e)
      aval[i + 12, 7] = "Random Forest"
      aval[i + 12, 8] = e@auc
      aval[i + 12, 9] = max(e@TPR + e@TNR) - 1
      aval[i + 12, 10] = tr
      aval[i + 12, 11] = i
      if(missing(proj)){RF.mod = predict(predictors, RF)
      }else(RF.mod = predict(proj, RF))
      if(missing(tss)){
        writeRaster(RF.mod, paste0("./modelos/", "RF_", i, "_con.tif"), format = "GTiff", 
                    overwrite = T)
        writeRaster(RF.mod > tr, paste0("./modelos/bin/", "RF_", i, "_bin.tif"), format = "GTiff", 
                    overwrite = T)
        png(paste0("./png/", "RF_", i, "con.png"))
        plot(RF.mod, main = paste0("RF con part_", i))
        dev.off()
        png(paste0("./png/", "RF_", i, "bin.png"))
        plot(RF.mod > tr, main = paste0("RF bin part_", i))
        dev.off()
      }
      if(missing(tss)==FALSE){
        if(TSS.calc>tss){
          writeRaster(RF.mod, paste0("./modelos/", "RF_", i, "_con.tif"), format = "GTiff", 
                      overwrite = T)
          writeRaster(RF.mod > tr, paste0("./modelos/bin/", "RF_", i, "_bin.tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "RF_", i, "con.png"))
          plot(RF.mod, main = paste0("RF con part_", i))
          dev.off()
          png(paste0("./png/", "RF_", i, "bin.png"))
          plot(RF.mod > tr, main = paste0("RF bin part_", i))
          dev.off()
        }
      }
      
      if (i == k) {
        if(length(list.files("./modelos", pattern = c("RF", ".tif"), full.names = T))!=0){
          RF.stack = stack(list.files("./modelos", pattern = c("RF", ".tif"), full.names = T))
          RF.ens = mean(RF.stack)
          
          # padronizando de 0 a 1
          values(RF.ens) = values(RF.ens)/RF.ens@data@max
          
          # recorte com TSSth
          if(mod=="after"){values(RF.ens)[values(RF.ens) < mean(aval[grep("Random Forest", aval[, 7]), 2])] = 0}
          
          writeRaster(RF.ens, paste0("./ensembles/", "RF_", "ensemble", ".tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "RF_", "ensemble", ".png"))
          plot(RF.ens, main = "Random Forest ensemble")
          dev.off()
          rm(RF.ens)
        }
        cat(c("\n", "Terminou Random Forest"))
      }
    }
  }
  
  if (SVM == T) {
    # SVM ####
    library(kernlab)
    for (i in 1:k) {
      cat(c("\n", "Começou a partição", i, "SVM"))
      
      pres_train <- pts1[group.p != i, ]
      pres_test <- pts1[group.p == i, ]
      backg_train <- backg[group.a != i, ]
      backg_test <- backg[group.a == i, ]
      
      train <- rbind(pres_train, backg_train)
      pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
      envtrain <- extract(predictors, train)
      envtrain <- data.frame(na.omit(cbind(pa = pb_train, envtrain)))
      
      #SVM <- ksvm(pa ~ ., data = envtrain)
      SVM <- ksvm(pa ~ ., data = envtrain, cross = k)
      
      e = evaluate(pres_test, backg_test, SVM, predictors)
      tr = e@t[which.max(e@TPR + e@TNR)]
      TSS.calc=max(e@TPR + e@TNR) - 1
      aval[i + 15, ] = threshold(e)
      aval[i + 15, 7] = "SVM"
      aval[i + 15, 8] = e@auc
      aval[i + 15, 9] = max(e@TPR + e@TNR) - 1
      aval[i + 15, 10] = tr
      aval[i + 15, 11] = i
      if(missing(proj)){SVM.mod = predict(predictors, SVM)
      }else(SVM.mod = predict(proj, SVM))
      
      if(mod=="before"){values(SVM.mod)[values(SVM.mod) < tr] = 0}
      if(missing(tss)){
        writeRaster(SVM.mod, paste0("./modelos/", "SVM_", i, "_con.tif"), format = "GTiff", 
                    overwrite = T)
        writeRaster(SVM.mod > tr, paste0("./modelos/bin/", "SVM_", i, "_bin.tif"), format = "GTiff", 
                    overwrite = T)
        png(paste0("./png/", "SVM_", i, "con.png"))
        plot(SVM.mod, main = paste0("SVM con part_", i))
        dev.off()
        png(paste0("./png/", "SVM_", i, "bin.png"))
        plot(SVM.mod > tr, main = paste0("SVM bin part_", i))
        dev.off()
      }
      
      if(missing(tss)==FALSE){
        if(TSS.calc>tss){
          writeRaster(SVM.mod, paste0("./modelos/", "SVM_", i, "_con.tif"), format = "GTiff", 
                      overwrite = T)
          writeRaster(SVM.mod > tr, paste0("./modelos/bin/", "SVM_", i, "_bin.tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "SVM_", i, "con.png"))
          plot(SVM.mod, main = paste0("SVM con part_", i))
          dev.off()
          png(paste0("./png/", "SVM_", i, "bin.png"))
          plot(SVM.mod > tr, main = paste0("SVM bin part_", i))
          dev.off()
        }
      }
      
      #Ensemble do algoritmo
      if (i == k) {
        if(length(list.files("./modelos", pattern = c("SVM", ".tif")))!=0){
          SVM.stack = stack(list.files("./modelos", pattern = c("SVM", ".tif"), full.names = T))
          SVM.ens = mean(SVM.stack)
          
          # padronizando de 0 a 1
          values(SVM.ens) = values(SVM.ens)/SVM.ens@data@max
          
          # recorte com TSSth
          if(mod=="after"){values(bc.ens)[values(SVM.ens) < mean(aval[grep("SVM", aval[, 7]), 2])] = 0}
          writeRaster(SVM.ens, paste0("./ensembles/", "SVM_", "ensemble ", ".tif"), format = "GTiff", 
                      overwrite = T)
          png(paste0("./png/", "SVM_", "ensemble ", "con.png"))
          plot(SVM.ens, main = "SVM ensemble")
          dev.off()
          rm(SVM.ens)
        }
        cat(c("\n", "Terminou SVM"))
      }
    }
  }
  
  
  # Ensemble final ####
  
  names(aval)[1:6] = names(threshold(e))
  names(aval)[7:11] = c("Algoritmo", "AUC", "TSS", "TSSth", "Partição")
  write.table(na.omit(aval), "Avaliação.csv", sep = ";", dec = ".",row.names = F)
  
  if(length(list.files("./ensembles", pattern = ".tif", full.names = T))!=0){
    final = stack(list.files("./ensembles", pattern = ".tif", full.names = T))
    mm = mean(final)
    values(mm) = values(mm)/mm@data@max
    
    writeRaster(mm, paste0("./final/", "Geral_", "ensemble", ".tif"), format = "GTiff", overwrite = T)
    
    png(paste0("./png/", "_Geral_", "ensemble", ".png"))
    plot(mm, main = paste0("Ensemble ", "geral"))
    dev.off()
    
    png(paste0("./png/", "_Geral_", "ensemble", "pontos", ".png"))
    plot(mm, main = paste0("Ensemble ", "geral"))
    points(pts1)
    dev.off()
    
    ### Modelagem até aqui ###
    
    if (plot == T) {
      plot(mm)
      points(pts1)
    }
    setwd(original)
    cat("\nTerminou! Verifique seus modelos.")
  }else(cat(paste0("\nNenhum modelo apresentou TSS maior que ",tss)))
  setwd(original)
  #stop(paste0("Nenhum modelo apresentou TSS maior que ",tss),call. = F)
}
