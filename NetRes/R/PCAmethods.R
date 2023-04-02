#' @include LatSpace.R
#' 
# @description Calculate latent variables using linear PCA
# residuals samples x vars matrix of residuals
NetRes$set("private", "calculateLatVarsPCA", function(residuals, resPpca, scale=F, algorithm.args=NULL) {
  resPc = prcomp(residuals, rank=resPpca$n, scale=scale)
  
  return(list(nLatVars = resPpca$n, 
              latVars = resPc$x,
              lvPredictor = resPc))
})

# @description Calculate latent variables using sparse PCA
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

# @description Calculate latent variables using sparse, robust PCA
# residuals samples x vars matrix of residuals
NetRes$set("private", "calculateLatVarsRobustSparsePCA", function(residuals, resPpca, scale=F, algorithm.args=NULL) {
  ##resPc = prcomp(residuals, rank=resPpca$n, scale=scale)
  resPc = robspca(residuals, k=resPpca$n, scale=scale)
  
  scores = resPc$scores
  colnames(scores) = paste('PC', 1:ncol(scores), sep='')
  resPc$scores = scores
  return(list(nLatVars = resPpca$n, 
              latVars = resPc$scores,
              lvPredictor = resPc))
})