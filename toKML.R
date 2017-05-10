toKML = function(modelo, name = "meuKML", zeros = FALSE, open = FALSE){
  # instalando pacote
  packages = c("plotKML")
  for (p in setdiff(packages, installed.packages()[, "Package"])) {
    install.packages(p,dependencies = T)
  }
  if(missing(modelo)){stop("Informe o nome do objeto que contem o raster do modelo")}
  if(zeros==FALSE){
    values(modelo)[values(modelo)==0]=NA
    plotKML(modelo,folder.name=name,open.kml = open)
  }else(plotKML(modelo,folder.name=name,open.kml = open))
}
