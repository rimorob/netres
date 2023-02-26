#' @include NetRes.R

# @description The top-level function to calculate latent variables from residuals
NetRes$set("private", "calculateLatVars", function(residuals, method='pca', scale=FALSE, algorithm.args=NULL) {
  if (is.null(algorithm.args)) {
    maxRank = ncol(residuals) - 1
  } else {
    maxRank = algorithm.args$max.rank
  }
  maxRank = min(round(ncol(residuals)*0.1), maxRank)
  print(paste('Maximum latent space rank set to', maxRank))
  ##how many PCs are we dealing with?  Use Horn's parallel analysis to find out
  resPpca = parallelPCA(
    as.matrix(residuals), #rows have to be variables, so transpose
    max.rank = maxRank,
    niters = 500,
    threshold = 0.05/ncol(residuals), #bonferroni-corrected p-vaoue cut-off
    transposed = TRUE, #variables in columns
    scale = scale,
    ##BSPARAM = IrlbaParam(),
    BPPARAM = BiocParallel::MulticoreParam()
  )
  
  if (resPpca$n == 0) {
    stop('Latent space not identified - aborting')
  }  
  
  if (method == 'pca') {
    latVars = private$calculateLatVarsPCA(residuals, resPpca = resPpca, scale=scale, algorithm.args=algorithm.args)  
  } else if (method == 'sparsepca') {
    latVars = private$calculateLatVarsSparsePCA(residuals, resPpca = resPpca, scale=scale, algorithm.args=algorithm.args)  
  } else if (method == 'autoencoder') {
    stop('Autoencoder is not ported to the new code version yet')  
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

# @description Calculate latent variables using linear PCA
# residuals samples x vars matrix of residuals
NetRes$set("private", "calculateLatVarsPCA", function(residuals, resPpca, scale=F, algorithm.args=NULL) {
  resPc = prcomp(residuals, rank=resPpca$n, scale=scale)
  
  return(list(nLatVars = resPpca$n, 
              latVars = resPc$x,
              lvPredictor = resPc))
})

# @description Calculate latent variables using linear PCA
# residuals samples x vars matrix of residuals
NetRes$set("private", "calculateLatVarsSparsePCA", function(residuals, resPpca, scale=F, algorithm.args=NULL) {
  ##resPc = prcomp(residuals, rank=resPpca$n, scale=scale)
  resPc = rspca(residuals, k=resPpca$n, scale=scale)

  scores = resPc$scores
  colnames(scores) = paste('PC', 1:ncol(scores), sep='')
  resPc$scores = scores
  return(list(nLatVars = resPpca$n, 
              latVars = resPc$scores,
              lvPredictor = resPc))
})
