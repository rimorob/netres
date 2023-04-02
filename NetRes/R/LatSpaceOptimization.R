#' @include LatSpace.R

##helper function to calcualte new lat vars given a linear combination of coefficients
mappedLatVars = function(coefs, latVars) {
  n = (1 + sqrt(1 + 4 * length(coefs)))/2
  coefs = t(matrix(c(rep(1, n), coefs), n, n)) #restore the parameters' shape; first row is unity (i.e., reference)
  ##compute the new linear combination
  newLatVars = scale(latVars %*% coefs)
  
  dim(newLatVars) = dim(latVars)
  colnames(newLatVars) = colnames(latVars)
  
  return(newLatVars)
}

##optimize latent variables to minimize the Bayes score
##do this by calculating a linear combination to optimize the ensemble BIC
##note that the sign of returned BIC is flipped for minimization by default
##also note that nBoot overrides the function default unless otherwise specified -
##this is done because optimization over a grid doesn't require oversampling and therefore
##we want to cash in on a speedup by using a smaller ensemble size
NetRes$set("private", "latVarLinComb", function(coefs, latVars, algorithm, algorithm.args, lvPrefix, dframe, cluster,
                         maximize=TRUE, nBoot=detectCores() - 2, returnEnsemble = FALSE) {
  newLatVars = mappedLatVars(coefs, latVars)
  
  ##drop the latent variables from the data frame
  newDframe = private$exciseLatVarsFromDataFrame(dframe, lvPrefix = '^U\\_')
  newDframe = cbind(newDframe, newLatVars)
  
  ##update the blacklist for new variable names
  new.algorithm.args = algorithm.args                        
  new.algorithm.args$blacklist = rbind(algorithm.args$blacklist,
                                       private$makeBlacklist(newDframe, lvPrefix)) 
  ##add the new latent space back
  if (returnEnsemble) { # then run a bootstrap with the new lat vars
    ##newEns = bn.boot(newDframe, statistic = dud, R=nBoot, algorithm = algorithm, algorithm.args = new.algorithm.args, cluster=cluster)
    newEns = bn.boot(newDframe, statistic = dud, R=nBoot, algorithm = algorithm, algorithm.args = new.algorithm.args, cluster=cluster)
    
    ##newBIC = mean(sapply(newEns, function(net, data, algorithm.args) {
    ##  score(net, data, type = new.algorithm.args$score, prior = new.algorithm.args$prior)
    ##}, newDframe, new.algorithm.args))
    newBIC = mean(parSapply(cluster, newEns, function(net, data, algorithm.args) {
      ##note: latent space is regularized during Horn's PFA separately and is excluded from score computation
      scores = score(net, data, type = algorithm.args$score, prior = algorithm.args$prior, by.node=TRUE)
      idx = grep(lvPrefix, names(scores))
      if (length(idx) > 0) {
        return(sum(scores[-idx]))
      } else {
        return(sum(scores))
      }
    }, newDframe, new.algorithm.args))
    
  } else { #calculate a single network to be used as part of GA optimization
    ##learn the network structure once rather than via a bootstrap
    net = do.call(algorithm, c(list(x = newDframe), new.algorithm.args))
    ##newBIC = score(net, newDframe, type = new.algorithm.args$score, prior = new.algorithm.args$prior)
    ##note: latent space is regularized during Horn's PFA separately and is excluded from score computation
    newBIC = {
      scores = score(net, newDframe, type = algorithm.args$score, prior = algorithm.args$prior, by.node=TRUE)
      idx = grep(lvPrefix, names(scores))
      if (length(idx) > 0) {
        sum(scores[-idx])
      } else {
        sum(scores)
      }}
    
    newEns = list(net)
  }
  
  newBIC = ifelse(maximize, newBIC, -newBIC)
  
  if (returnEnsemble) {
    return(list(BIC = newBIC, ens = newEns, mappedLVs = newLatVars))
  } else {
    return(newBIC)
  }
}) 