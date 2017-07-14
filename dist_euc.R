# Function to calculate enviromental distance
# by Bruno Vilela (bvilela@wustl.edu)
# occ a two column matrix or a data.frame with the first column being the 
#  longitude and the second the latitude
# env a raster with the enviromental variables, croped to the projection area
#  method if mean, the "mean" distance to all points is calculated, if "min"
#  the minimun distance to any point is calculated.
# decostand.method What method should be applied to standardize the enviromental data.
#  see the vegan funtion decostand, argument method.
# suitability Whether the function should return a suitability (TRUE) 
#  or the actual distance calculated (FALSE).
#
# The function returns a raster with the distance (for suitability = FALSE) 
# or the suitability values (for suitability = TRUE).

dist_euc <- function (occ, env, method = "mean", 
                      decostand.method = "standardize",
                      suitability = FALSE) {
  require(raster)
  require(vegan)
  if (class(env) != "raster" & class(env) != "RasterStack"
      & class(env) != "RasterBrick") {
    stop("env has to be a raster")
  }
  if (class(occ) != "matrix" & class(occ) != "data.frame") {
    stop("occ has to be a matrix or data.frame")
  }
  if (ncol(occ) != 2) {
    stop("occ has to be a matrix or data.frame of 2 columns (x and y)")
  }
  values <- values(env)
  values <- apply(values, 2, decostand, method = decostand.method,
                  na.rm = TRUE)
  values(env) <- values
  values_occ <- extract(env, occ)
  pos <- is.na(values[, 1])
  values2 <- values[!pos, ]
  n <- nrow(values2)
  eu <- numeric(n)
  values_occ <- na.omit(values_occ)
  
  for (i in 1:n) {
    for (j in 1:ncol(values_occ)) {
      temp <- eu[i] + ((values2[i, j] - values_occ[, j]) ^ 2)
      if (method == "mean") {
        eu[i] <- mean(temp, na.rm = TRUE)
      }
      if (method == "min") {
        eu[i] <-  min(temp, na.rm = TRUE)  
      }
    }
    eu[i] <- sqrt(eu[i])
  }
  env <- raster(env, 1)
  if (!suitability) {
    values(env)[!pos] <- eu
    return(env)
  } 
  if (suitability) {
    values(env)[!pos] <- decostand(-eu, method = "range", na.rm = TRUE)
    return(env)
  }
}

#Elaborated by Bruno Vilela
#5/23/2017
#http://rpubs.com/Bruno_Vilela/279257
