#' Function to perform ALRA on a Seurat2 object
#' 
#' This function performs ALRA on the "data" slot in a given Seurat object, 
#' does the neccessary pre- and post-processing to fit with the object, use 
#' the ALRA imputed data to REPLACE the "data" slot in the Seurat object, and 
#' returns that updated Seurat object.
#' 
#' @param obj the Seurat object to run ALRA
#' @param ... parameters passing to alra()
#' 
#' @return an updated Seurat object with ALRA imputed data in the "data" slot
#' 
alraSeurat2 <- function(obj, ...){
  data_alra <- t(alra(t(as.matrix(obj@data)), ...)[[3]])
  colnames(data_alra) <- obj@cell.names
  data_alra <- Matrix(data_alra, sparse = T)
  obj@data <- data_alra
  return(obj)
}

