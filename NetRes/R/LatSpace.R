#' @include NetRes.R
NULL

# @description The top-level function to calculate latent variables from residuals
NetRes$set("private", "calculateLatVars", function(residuals, method='pca', scale=FALSE, algorithm.args=NULL, BPPARAM = BiocParallel::DoparParam()) {
  if (is.null(algorithm.args)) {
    maxRank = ncol(residuals) - 1
  } else {
    maxRank = algorithm.args$max.rank
  }
  maxRank = min(round(ncol(residuals)*0.1), maxRank)
  print(paste('Maximum latent space rank set to', maxRank))
  startTime = proc.time()
  ##how many PCs are we dealing with?  Use Horn's parallel analysis to find out
  resPpca = parallelPCA(
    as.matrix(residuals), 
    max.rank = maxRank,
    niters = 500,
    threshold = 0.05/ncol(residuals), #bonferroni-corrected p-vaoue cut-off
    transposed = TRUE, #variables in columns
    scale = scale,
    BPPARAM = BPPARAM
  )
  print("Done with Horn's PFA")

  print(paste('Took', proc.time()[1] - startTime[1], 'seconds'))
  if (resPpca$n == 0) {
    stop('Latent space not identified - aborting')
  }  
  
  if (method == 'pca') {
    latVars = private$calculateLatVarsPCA(residuals, resPpca = resPpca, scale=scale, algorithm.args=algorithm.args)  
  } else if (method == 'sparse.pca') {
    latVars = private$calculateLatVarsSparsePCA(residuals, resPpca = resPpca, scale=scale, algorithm.args=algorithm.args)  
  } else if (method == 'robust.sparse.pca') {
    latVars = private$calculateLatVarsRobustSparsePCA(residuals, resPpca = resPpca, scale=scale, algorithm.args=algorithm.args)  
  } else if (method == 'NNMF') {
    latVars = private$calculateLatVarsNNMF(residuals, resPpca = resPpca, scale=scale, algorithm.args=algorithm.args)      
  } else if (method == 'VAE') {
    ##Note: the argument "scale" is neither supported nor necessary
    latVars = private$calculateLatVarsVAE(residuals, resPpca = resPpca, algorithm.args=algorithm.args)      
  } else {
    stop(paste('Method', method, 'for latent space analysis has not been implemented yet'))
  }

  ##The methods above return a common output list with three fields; format it as a common R6 object
  LatentResult = R6Class('LatentResult',
                         public=list(
                           n = NULL,
                           v = NULL,
                           generator = NULL,
                           initialize = function(lvList) {
                             self$n = lvList$nLatVars
                             self$v = lvList$latVars
                             colnames(self$v) = paste('U_', colnames(self$v), sep='')
                             self$generator = lvList$lvPredictor
                           },
                           predict = function(newdata = NULL) {
                             if (is.null(newdata)) {
                               return(predict(self$generator))
                             }
                             return(predict(self$generator, newdata=newdata))
                           }
                         )
  )
  latRes = LatentResult$new(latVars)

  return(latRes)
})


# @description Calculate latent variables using NNMF
# residuals samples x vars matrix of residuals
NetRes$set("private", "calculateLatVarsNNMF", function(residuals, resPpca, scale=F, algorithm.args=NULL) {
  ##no native scaling in NNMF, so scale the residuals matrix if needed
  if (scale) {
    residuals = scale(residuals)
  }
  resPc <- nnmf(as.matrix(residuals), k=resPpca$n)
  
  scores = resPc$W #left matrix is rows by eigenspace
  ##call these PCs for compatibilty
  colnames(scores) = paste('PC', 1:ncol(scores), sep='')
  resPc$W = scores
  return(list(nLatVars = resPpca$n, 
              latVars = resPc$W,
              lvPredictor = resPc))
})
