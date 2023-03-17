#' @title NetRes NetRes class definition
#' @import R6 dplyr PCAtools
#' @importFrom bnlearn as.bn as.igraph parents drop.arc as.graphNEL bn.boot nodes remove.node bn.fit
#' @importFrom parallel makeCluster stopCluster detectCores parLapply parSapply 
#' @importFrom doParallel registerDoParallel
#' @importFrom BiocSingular IrlbaParam
#' @importFrom corrplot corrplot
#' @importFrom paran paran
#' @importFrom GA ga
#' @importFrom GenSA GenSA
#' @import sparsepca
#' @importFrom NNLM nnmf
#' 
##An R6 class for a generated scale-free network
#' @export
NetRes <- R6Class("NetRes", 
                  public = list(
                    #' @field ensemble The final ensemble as inferred after some iterations
                    ensemble = NULL,
                    #' @field latent.space The list of latent spaces after a given iteration
                    latent.space = NULL,
                    #' @field latent.space.transform The function to map PCs to the final basis vector
                    latent.space.transform = NULL,
                    #' @field latent.space.transform.coefs The transform that maps PCs to the final basis vector form s.t. the PCs are conditionally de-correlated given their descendants (i.e., the network is optimized)
                    latent.space.transform.coefs = NULL,
                    #' @field train.data The initial training data
                    train.data = NULL,
                    #' @field latent.data The true data for the latent space (if provided)
                    latent.data = NULL,
                    #' @field BIC The list of mean ensemble BICs per run - used as the stopping condition
                    BIC = list(),
                    #' @field true.graph The true graph (if any)
                    true.graph = NULL,
                    #' @description Constructor that fits a NetRes object from a data fram
                    #' @param dframe Training data frame, optionally with true latent variables (these will be used for performance reporting)
                    #' @param nIter Number of iterations (NULL to autodetect - not yet supported)
                    #' @param nBoot Number of bootstraps (default = 50)
                    #' @param algorithm Algorithm to use (default = 'tabu')
                    #' @param algorithm.args Arguments to the algorithm (must be a list)
                    #' @param lvPrefix Prefix of variables to be treated as latent (if any).  If not present, defautls to "U\\_".  If no latent vars found, 
                    #' no performance assessment of latent discovery will (or indeed can) be made, but whatever latent space is discovered will be duly reported
                    #' @param latentSpaceParentOnly If true (default), latent space variables can only be parents, never children                    #' 
                    #' @return A NetRes object
                    initialize = function(dframe, true.graph = NULL, nIter, nBoot=50, algorithm='tabu', algorithm.args, 
                                          lvPrefix = "^U\\_", mode=NULL,
                                          weightedResiduals = FALSE, scale=FALSE, debug=FALSE,
                                          latentSpaceParentOnly = TRUE,
                                          latentSpaceMethod = 'pca',
                                          optimizeLatentSpace=FALSE,
                                          nCores=NULL,
                                          BPPARAM=BiocParallel::MulticoreParam()) {
                      if(debug)
                        browser()
                      self$latent.data = dframe %>% select_if(grepl(lvPrefix, names(.)))
                      self$train.data = dframe %>% select_if(!grepl(lvPrefix, names(.)))
                      self$true.graph = true.graph
                      require(parallel)
                      if (mode == 'oracular') {
                        warning('!!!RUNNING IN ORACULAR MODE - USING TRUE LATENT VARIABLES FOR INFERENCE OF STRUCTURE!!!')
                        train = cbind(self$train.data, self$latent.data)
                      } else { #normal run
                        train = self$train.data
                      }
                      stopifnot(is.list(algorithm.args))
                      if(is.null(nCores)){
                          nCores = detectCores() - 1
                      }else{
                          nCores=min(detectCores()-1,nCores)
                      }
                      for (ni in 1:nIter) {
                        if(debug){
                          message("In interation ", ni)
                          browser()
                        }                      
                        ##cluster creation moved inside the loop in order to help clean up what seems like a memory leak in cluster creation
                        cluster = makeCluster(nCores)
                        registerDoParallel(cluster)
                        setDefaultCluster(cl=cluster)
                        clusterEvalQ(cluster, library(bnlearn))

                        curRes = private$runOneIteration(train, nBoot, algorithm, algorithm.args, cluster, lvPrefix = lvPrefix, 
                                                         weightedResiduals, scale=scale, optimizeLatentSpace=optimizeLatentSpace,
                                                         learnLatentSpace = ifelse(ni == 1, FALSE, TRUE),
                                                         latentSpaceParentOnly=latentSpaceParentOnly,
                                                         latentSpaceMethod = latentSpaceMethod,debug=debug,BPPARAM=BPPARAM)
                        if(debug){
                          message("After runOneIteration ")
                          browser()
                        }                      
                        self$ensemble[[ni]] = curRes$ensemble
                        self$latent.space[[ni]] = curRes$latent.space                        
                        ##self$latent.space.transform[[ni]] = curRes$latent.space.transform
                        ##self$latent.space.transform.coefs[[ni]] = curRes$latent.space.transform.coefs
                        print(paste('ni:', ni))
                        self$BIC[[ni]] = curRes$BIC #for PREVIOUS latent space - this one hasn't been evaluated yet - that's in the next iteration!
                        print(paste('BIC history:', as.numeric(self$BIC)))
                        
                        if (ni > 1 && self$BIC[[ni]] <= self$BIC[[ni-1]]) {
                          warning('Stopping early due to convergence')
                          break  
                        }

                        if (ni == 1) { #then skip to the next iteration
                          next
                        }                        
                        
                        if (is.null(mode) || mode != 'oracular') {
                          train = cbind(self$train.data, curRes$latent.space$v)
                        } else {
                          warning('!!!RUNNING IN ORACULAR MODE - USING TRUE LATENT VARIABLES FOR INFERENCE OF STRUCTURE!!!')
                          train = cbind(self$train.data, self$latent.data)
                        }
                        corrplot(cor(cbind(self$latent.data, curRes$latent.space$v)), method='ellipse', order='AOE', diag=F)
                        print('pausing to admire the corrplot')
                        Sys.sleep(5)
                        stopCluster(cluster)
                      }
                      if (ni == nIter) {
                        print('Maximum number of iterations reached - stopping')
                      }

                    },
                    #' @description assess Assess the inferred ensemble against the true graph
                    #' @param true.graph The true graph to use; defaults to the one provided at initialization (if any)
                    #' @param lvPrefix The latent variable-identifying regular expression, as elsewhere; defaults to "^U\\_"
                    assess = function(true.graph = self$true.graph, lvPrefix = "^U\\_",nCores=NULL) {
                      if (is.null(true.graph)) {
                        stop('Cannot assess performance without the true graph')
                      }
                      if(is.null(nCores)){
                          nCores = detectCores() - 2
                      }else{
                          nCores=min(detectCores()-2,nCores)
                      }
                      cluster = makeCluster(nCores)
                      registerDoParallel(cluster)                      
                      ##excise the latent space from the true graph
                      true.graph = private$exciseLatVarsFromEnsemble(list(true.graph), cluster, lvPrefix)[[1]]
                      true.graph.ig=as.igraph(true.graph)
                      aucs = c()
                      prAucs = c()
                      f1maxes = c()
                      ##test performance at every iteration
                      for (ni in 1:length(self$ensemble)) { 
                          print(paste('step', ni))
                          browser()
                          curEnsemble = private$exciseLatVarsFromEnsemble(self$ensemble[[ni]], cluster, lvPrefix)
                        curStrength = bnlearn::custom.strength(curEnsemble, bnlearn::nodes(true.graph))
                          curStrength=curStrength %>%
                              as.data.frame() %>%
                              mutate(freq=strength*direction)
                          perf=network_performance(true.graph.ig,curStrength)
                          
                        ##pred = as.prediction(curStrength, true.graph)
                        ##perf = ROCR::performance(pred, "tpr", "fpr")
                          ##auc = round(ROCR::performance(pred, "auc")@y.values[[1]], 3)
                          auc=minet::auc.roc(perf)
                          ##plot(perf, main = paste("AUC:", auc), colorize=TRUE)
                          aucpr=minet::auc.pr(perf)
                        ##aucpr = round(ROCR::performance(pred, "aucpr")@y.values[[1]], 3)                        
                        perf = ROCR::performance(pred, "prec", "sens")
                        plot(perf, main = paste("PR-AUC:", aucpr), colorize=TRUE)    
                        perf = ROCR::performance(pred, 'f')
                        plot(perf, main = paste('F1-AUC', mean(perf@y.values[[1]][-1])))
                        f1max = round(max(perf@y.values[[1]][-1]), 3)
                        aucs[[ni]] = auc
                        prAucs[[ni]] = aucpr
                        f1maxes[[ni]] = f1max
                      }
                      plot(1:length(self$ensemble), aucs, main='AUCs over iterations', xlab='Iteration', ylab='AUC')
                      plot(1:length(self$ensemble), prAucs, main='PR-AUCs over iterations', xlab='Iteration', ylab='PR-AUC')                      
                      plot(1:length(self$ensemble), f1maxes, main='F1max values over iterations', xlab='Iteration', ylab='F1max')
                      stopCluster(cluster)
                    },
                    #' @description plot The function to plot (some) networks in the ensemble
                    #' @param networks Indexes into the ensemble
                    #' @param lvPrefix The latent variable prefix (default = 'U_'; use perl regexp)
                    #' @returns Nothing
                    plot = function(networks = NULL, lvPrefix = "^U\\_") {
                      if (!is.null(networks)) {
                        stopifnot(!prod(networks %in% 1:length(self$ensemble)))
                      }
                      if (is.null(networks)) {#pick a single random network to plot
                        networks = sample(1:length(self$ensemble), 1)
                      }
                      latVars = grep(lvPrefix, colnames(self$train.data), perl=TRUE, value=TRUE)
                      for (ni in networks) {
                        graphviz.plot(self$ensemble[[ni]], main = paste("DAG with ", length(self$outDegree), " vertices", sep=''), layout='dot',
                                      highlight = list(nodes = latVars, fill="orange"))
                      }
                      browser()
                    }
                  ),
                  private = list(
                    # @description A utility function to create a blacklist from a data frame so that no variable drives a latent variable
                    makeBlacklist = function(dframe, lvPrefix = "^U\\_") {
                      latVars = grep(lvPrefix, colnames(dframe), value=T)
                      nonLatVars = setdiff(colnames(dframe), latVars)
                      blacklist = data.frame(from=NULL, to=NULL)
                      
                      for (nlv in nonLatVars) {
                        for (lv in latVars) {
                            blacklist = rbind(blacklist,
                                              data.frame(from=nlv,to=lv))
                                        # need to convert to data.frame:
                                        # - otherwise you lose column names
                                        # - then you can't merge to provided 
                        }
                      }
                      if (nrow(blacklist) == 0) {
                        blacklist=NULL
                      }
                      return(blacklist)
                    },
                    # @description The top-level function to build an ensemble and calculate residuals
                    # @param dframe The training data frame
                    # @param nBoot Number of bootstraps
                    # @param algorithm The network inference algorithm, as in bnlearn
                    # @param algorithm.args The list of arguments to the algorithm, as in bnlearn
                    # @param cluster The cluster object, as returned by makeCluster from package "parallel"
                    # @param lvPrefix The latent variable prefix (default = "U_"; use perl regexp)
                    # @param latentSpaceParentOnly If true (default), latent space variables can only be parents, never children
                    # @param learnLatentSpace If true (default), learn latent space, else return an ensemble without inferring latent space
                    # @param optimizeLatentSpace If false (default), return first-pass estimate of the latent space without conditional BIC-based optimization (which is very slow)
                    # @return TBD
                    runOneIteration  = function(dframe, nBoot, algorithm, algorithm.args, cluster, lvPrefix = "^U\\_", 
                                                weightedResiduals=FALSE, scale=FALSE, learnLatentSpace = TRUE,
                                                optimizeLatentSpace=FALSE, 
                                                latentSpaceParentOnly=TRUE, latentSpaceMethod = 'pca',
                                                BPPARAM=BiocParallel::MulticoreParam(),
                                                debug=FALSE) {
                      dud = function(x) x
                      
                      new.algorithm.args = algorithm.args
                      new.algorithm.args$blacklist = rbind(algorithm.args$blacklist,
                                                       private$makeBlacklist(dframe, lvPrefix))
                      ens = bn.boot(dframe, statistic = dud, R=nBoot, algorithm = algorithm, algorithm.args = new.algorithm.args, cluster=cluster)
                      ##ensBIC =  mean(parSapply(cluster, ens, function(net, data, algorithm.args) {
                      ##  score(net, data, type = algorithm.args$score, prior = algorithm.args$prior)
                      ##}, dframe, new.algorithm.args))
                      ensBIC = mean(parSapply(cluster, ens, function(net, data, algorithm.args) {
                        ##note: latent space is regularized during Horn's PFA separately and is excluded from score computation
                        scores = score(net, data, type = algorithm.args$score, prior = algorithm.args$prior, by.node=TRUE)
                        idx = grep(lvPrefix, names(scores))
                        if (length(idx) > 0) {
                          return(sum(scores[-idx]))
                        } else {
                          return(sum(scores))
                        }
                      }, dframe, new.algorithm.args))
                      
                      netWeights = private$calcBayesFactors(ens, dframe, cluster, algorithm.args)
                      ens2 = private$exciseLatVarsFromEnsemble(ens, cluster, lvPrefix)                      
                      if (!learnLatentSpace) { #then return the ensemble-specific BIC (w/o latent space learning)
                        return(list(ensemble=ens2,
                                    latent.space = NULL,
                                    BIC = ensBIC))
                      }

                      res = private$calculateResiduals(ens2, netWeights, weightedResiduals, cluster)
                      if(debug){
                          message("before latent variables estimation")
                          browser()                          
                      }
                      ##paran(res)
                      latvars = private$calculateLatVars(res, method=latentSpaceMethod, scale=scale, algorithm.args=new.algorithm.args,BPPARAM=BPPARAM) 

                      ##debugging step
                      trueCoefs = NULL
                      for (ci in 1:ncol(self$latent.data)) {
                        curDf = cbind(data.frame(out = self$latent.data[, ci]), latvars$v)
                        m = lm(out ~ ., curDf)
                        coefs = as.numeric(coef(m)) 
                        coefs = coefs[2:length(coefs)]
                        if (is.null(trueCoefs)) {
                          trueCoefs = coefs
                          dim(trueCoefs) = c(1, length(coefs))
                        } else {
                          trueCoefs = rbind(trueCoefs, coefs)
                        }
                      }

                      rownames(trueCoefs) = colnames(self$latent.data)
                      colnames(trueCoefs) = colnames(latvars$v)
                      print('True coefficients:')
                      print(trueCoefs)
                      
                      ##Note that the number of parameters being optimized is one less than rows 
                      ##Otherwise, we'd be optimizing weights where the maximum lies on a manifold
                      ##In other words, the first coefficient of each weight vector will be arbitrarily set to one
                      latCoefs = array(data = 1, dim=c(ncol(latvars$v) - 1, ncol(latvars$v)))

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
                      latVarLinComb = function(coefs, latVars, algorithm, algorithm.args, lvPrefix, dframe, cluster,
                                               maximize=TRUE, nBoot=detectCores() - 2, returnEnsemble = FALSE) {
                        print(paste('nBoot inside function:', nBoot))
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
                     } 

                      if (ncol(latCoefs) > 1 && optimizeLatentSpace) { 
                        print('optimizing the basis vector of the latent space')                        
                        startTime = Sys.time()
                        optimRes = ga(type = 'real-valued', fitness=latVarLinComb, latVars = latvars$v, 
                                      algorithm = algorithm,
                                      algorithm.args = algorithm.args,
                                      lvPrefix = lvPrefix,
                                      dframe = dframe,
                                      cluster = cluster,
                                      nBoot = 1, #"fast" regime
                                      lower = rep(-100, length(as.numeric(latCoefs))),
                                      upper = rep(100, length(as.numeric(latCoefs))),
                                      run = 10, ##increase this eventually, or leave at default (maxiter)
                                      monitor=plot,
                                      parallel = cluster,
                                      optim = FALSE, #use optim for local optimization
                                      popSize = 500 #and elitism defaults to 5 winners
                        )
                        print(Sys.time() - startTime)
                        mappedLVs = mappedLatVars(optimRes@solution, latvars$v)
                        ensBICoptimized = optimRes@fitnessValue                          

                        ##FIX THE BELOW - CURRENTLY OVERRIDING LATVARS$V AND UNABLE TO PREDICT CORRECTLY OOS, NEED TO PASS ON THE MAPPING BASIS
                        ##FOR NOW, JUST RETURN THE CORRECT VECTOR OVER TRAINING DATA                          

                        ##latvars$v = mappedLVs
                        remappedRes = latVarLinComb(optimRes@solution, latvars$v, algorithm = algorithm,
                                                    algorithm.args = algorithm.args,
                                                    lvPrefix = lvPrefix,
                                                    dframe = dframe,
                                                    cluster = cluster,
                                                    nBoot = nBoot, #for recalculating the ensemble
                                                    returnEnsemble = TRUE)

                        print(Sys.time())
                        print('debug:')
                        print(coefs)
                        print(remappedRes$newBIC)
                        print('correlation structure BEFORE optimization:')
                        print(cor(cbind(self$latent.data, latvars$v)))
                        print('correlation structure AFTER optimization:')                        
                        print(cor(cbind(self$latent.data, remappedRes$mappedLVs)))
                        print('---')
                        
                        latvars$v = remappedRes$mappedLVs #hack for now - eventually don't overwrite the original latent space basis
                        ens = remappedRes$ens
                        ensBIC = remappedRes$BIC
                      } else { #recalculate the network with the identified latent space
                        newDframe = private$exciseLatVarsFromDataFrame(dframe, lvPrefix = '^U\\_')
                        newDframe = cbind(newDframe, latvars$v)                        
                        new.algorithm.args = algorithm.args                        
                        new.algorithm.args$blacklist = rbind(algorithm.args$blacklist,
                                                             private$makeBlacklist(newDframe, lvPrefix))                         
                        ens = bn.boot(newDframe, statistic = dud, R=nBoot, algorithm = algorithm, algorithm.args = new.algorithm.args, cluster=cluster)
                        #ensBIC =  mean(parSapply(cluster, ens, function(net, data, algorithm.args) {
                        #  score(net, data, type = algorithm.args$score, prior = algorithm.args$prior)
                        #}, newDframe, new.algorithm.args))
                        ensBIC = mean(parSapply(cluster, ens, function(net, data, algorithm.args) {
                            ##note: latent space is regularized during Horn's PFA separately and is excluded from score computation
                            scores = score(net, data, type = algorithm.args$score, prior = algorithm.args$prior, by.node=TRUE)
                            idx = grep(lvPrefix, names(scores))
                            if (length(idx) > 0) {
                                return(sum(scores[-idx]))
                            } else {
                                return(sum(scores))
                            }
                        }, newDframe, new.algorithm.args))
                      }
                      
                      return(list(ensemble=ens,
                                  latent.space = latvars, 
                                  ##latent.space.transform = mappedLatVars,
                                  ##latent.space.transform.coefs = optimRes@solution,
                                  BIC = ensBIC))                      
                    },
                    
                    ##
                    # @description calcBayesFactors Calculate bayes factors for every (unexcised) network in the ensemble
                    # Use same score and prior as used to infer the ensemble to calculate the bayes factors as well
                    # @param ens ensemble
                    # @param dframe data frame
                    # @param cluster cluster to distribute to
                    # @param algorithm.args algorithm arguments to get score from
                    calcBayesFactors = function(ens, dframe, cluster, algorithm.args) {
                      ##compute all scores relative to the first network
                      BFs = parSapply(cluster, ens, function(net, data, algorithm.args) {
                        bnlearn::BF(net, ens[[1]], dframe, score = algorithm.args$score, prior = algorithm.args$prior, log=F)
                      }, dframe, algorithm.args)
                      BFs[BFs > 1e20] = 1e20 #cap max values

                      ##renormalize, with some floor, roughly O(1/(2*length(nets))) 
                      ##That is, the hoi polloi get half the weight, the best nets the other half
                      scores = (BFs/sum(BFs) + 1/length(BFs))
                      scores = scores/sum(scores)
                      return(scores)  
                    },

                    # @description The private helper function to calculate residuals for a given ensemble
                    calculateResiduals = function(ens, netWeights, weightedResiduals, cluster) {
                      residuals = parLapply(cluster, ens, function(net) {
                        ##fit a single network using training data
                        fitted = bn.fit(net, self$train.data[, bnlearn::nodes(net)])
                        netRes = residuals(fitted)

                        return(as.data.frame(netRes))
                      })
                      
                      if (weightedResiduals) {
                        for (ri in 1:length(residuals)) {
                          residuals[[ri]] = residuals[[ri]] * netWeights[ri]
                        }
                        
                        ##average the residuals across all networks (not clear if this is the best way, but likely it is)                      
                        residuals = Reduce(`+`, residuals) / sum(netWeights)                      
                      } else {
                        residuals = Reduce(`+`, residuals) / length(netWeights)                                              
                      }
                      return(residuals)
                    },
                    # @description The private helper function to remove latent vars from the ensemble
                    # @param ens The bnlearn network ensemble
                    # @param cluster The cluster object, as returned by makeCluster from package "parallel"
                    # @param lvPrefix The latent variable prefix (default = 'U_'; use perl regexp)
                    # @return The ensemble with the latent vars excised 
                    exciseLatVarsFromEnsemble = function(ens, cluster, lvPrefix = '^U\\_') {
                      ens = parLapply(cluster, ens, function(net) { 
                          vars = bnlearn::nodes(net)
                          latVars = grep(lvPrefix, vars, perl=TRUE, value=TRUE)
                          for (lv in latVars) {
                            net = bnlearn::remove.node(net, lv)
                          }
                          return(net)
                        })
                      return(ens)
                    },
                    # @description The private helper function to remove latent vars from a data frame
                    # @param dframe data frame
                    # @param lvPrefix The latent variable prefix (default = 'U_'; use perl regexp)
                    # @return The data frame with the latent vars excised 
                    exciseLatVarsFromDataFrame = function(dframe, lvPrefix = '^U\\_') {                    
                      latIdx = grep(lvPrefix, colnames(dframe))
                      if (length(latIdx) > 0) {
                        newDframe = dframe[, -latIdx]
                      } else {
                        newDframe = dframe
                      }
                      return(newDframe)
                    }
                  )
)

plotIgraph = function(g,
                      edgelabels=F,
                      edgeweights = FALSE,
                      edgelabelsFilter = 0,
                      edgelabelsFilter_useabs = TRUE,
                      lwdmin = 0.5,
                      lwdmax = 3,
                      damping = 0.2,
                      overlap = F,
                      splines = TRUE,
                      nodesep = 1,
                      pad = .5,
                      sep = 1,
                      ranksep = 0.05,
                      start = 123,
                      layoutType = 'dot',
                      saveToFile = F,
                      filename = 'net.pdf',
                      width = 1000 / 100,
                      height = 1000 / 100,
                      other_pdf_options = list(),
                      nodeThickness = 1,
                      nodeThickness_important = 2,
                      fill = NULL,
                      edge_color = NULL,
                      edge_labels = NULL,
                      label_pad = 2
                      ){
    allnodes = names(igraph::V(g))
    if(!is.null(fill)){
	fill = expandListRegex(fill, allnodes)
    }
    if(!is.null(edge_color)){
	edge_color = expandDfRegex(edge_color, allnodes)
    }
    ##library(Rgraphviz)
    node_width = NULL
    node_height = NULL
    node_fixedSize = FALSE
    gr = igraph:::as_graphnel(g)
    eAttrs <- list()
					#w <- w[setdiff(seq(along=w), removedEdges(gr))]
    if(length(igraph::get.edge.attribute(g)) == 0){
        edgelabels = FALSE
        edgeweights = FALSE
    }
    if(edgeweights | edgelabels){
        w = signif(igraph::get.edge.attribute(g)[[1]], 2)
        names(w) = sub("\\|", "~", attributes(igraph::E(g))$vnames)
        ##names(w) <- edgeNames(gr, recipEdges="distinct")
        ##names(eAttrs$lwd) = edgeNames(gr, recipEdges="distinct")
    }
    if(edgelabels){
        if(!is.null(edge_labels)){
            edgs = paste(edge_labels[[1]], edge_labels[[2]], sep = '~')
            alledges = graph::edgeNames(gr)
            alledgeslab= rep("", length(alledges))
            names(alledgeslab) = alledges
            alledgeslab[edgs] = as.character(edge_labels[[3]])
            namesw = names(alledgeslab)
            alledgeslab = paste0(paste0(rep(" ",label_pad),collapse=""),alledgeslab)
            names(alledgeslab) = namesw
            eAttrs$label = alledgeslab
            
        }else{
            if(edgelabelsFilter_useabs)
                iiw = which(abs(w) > edgelabelsFilter)
            else
                iiw = which(w > edgelabelsFilter)
            eAttrs$label = rep("", length(w))
            names(eAttrs$label) = names(w)
            wval = w[iiw]
            wval = paste0(paste0(rep(" ",label_pad),collapse=""),wval)
            eAttrs$label[names(w[iiw])] = wval
        }
    }
    if(edgeweights){
        wn = as.numeric(w)
        eAttrs$lwd <- (lwdmax - lwdmin) * (wn - min(wn)) / (max(wn) - min(wn)) + lwdmin
        names(eAttrs$lwd) = sub("\\|", "~", attributes(igraph::E(g))$vnames)
    }
    eAttrs$direction = rep("forward", length(graph::edgeNames(gr, recipEdges="distinct")))
    names(eAttrs$direction) = graph::edgeNames(gr, recipEdges="distinct")
    ## edge colors
    if(!is.null(edge_color)){
        alledges = graph::edgeNames(gr)
        alledgescol = rep("black", length(alledges))
        names(alledgescol) = alledges
        edgs = paste(edge_color[[1]], edge_color[[2]], sep = '~')
        alledgescol[edgs] = as.character(edge_color[[3]])
        eAttrs$color = alledgescol
    }
    attrs = list(node = list(shape = "ellipse",
                             fixedsize = node_fixedSize,
                             width = node_width,
                             height = node_height,
                             lwd = nodeThickness,
                             color = 'black'
                             ),
                 edge=list(
                           direction = 'forward',
                           concentrate = F
                           ),
                 graph = list(damping = damping,
                              nodesep = nodesep,
                              pad = pad,
                              ranksep = ranksep,
                              splines = splines,
                              start = start,
                              dim = 2,
                              sep = sep,
                              concentrate = FALSE,
                              overlap = overlap))
    nAttrs = list()
    gr = Rgraphviz::layoutGraph(gr, layoutType = layoutType,
                                attrs = attrs,
                                nodeAttrs=nAttrs,
                                edgeAttrs=eAttrs,
                                recipEdges="distinct"
                                )
    ## add filler
    nAttrs$fill = fill
    ## add node edge width
    if(!is.null(node_width)){
        nAttrs$width = node_width
        ## graph::nodeRenderInfo(gr) = list(
        ##     width = node_width)
    }
    if(!is.null(node_height)){
        nAttrs$height = node_height
        ## graph::nodeRenderInfo(gr) = list(
        ##     height = node_height)
    }
    nAttrs$fixedSize = node_fixedSize
    graph::nodeRenderInfo(gr) = nAttrs
    ##graph::nodeRenderInfo(gr) =list(fixedSize = node_fixedSize)
    graph::edgeRenderInfo(gr) = eAttrs
    if(saveToFile){
        other_pdf_options = c(other_pdf_options,
                              file = filename,
                              width = width,
                              height = height)
        do.call(pdf, args = other_pdf_options)
					#pdf(filename, width = width, height = height)
        Rgraphviz::renderGraph(gr)
        dev.off()
    }else
        Rgraphviz::renderGraph(gr)
}


expandListRegex = function(mylist, allnames){
    newlist = list()
    for(ll in 1:length(mylist)){
	rg = names(mylist)[ll]
	rg = paste0("^", rg, "$")
	val = mylist[[ll]]
	mv = grep(rg, allnames, value = T)
	tmp = rep(list(val), length(mv))
	names(tmp) = mv
	newlist = c(newlist, tmp)
    }
    return(newlist)
}

expandDfRegex = function(mydf, allnames){
    newlist = map_df(1:nrow(mydf),
		     function(ii){
			 val = mydf[ii, ]
			 inputs = grep(paste0("^", val$inp, "$"),
				       allnames,
				       value = T)
			 outputs = grep(paste0("^", val$out, "$"),
					allnames,
					value = T)
			 expand.grid(inputs, outputs, stringsAsFactors=F) %>%
			     mutate(color = val[[3]])
		     })
    return(newlist)
}

##' @description
##' Given a true graph in igraph format and a estimated edge frequency create a data.frame with performance metrics
##'
##' @details
##' This function convert true and estimated network to adjacency matrixed and use the minet package and validate function to generate data.frame with four columns named  thrsh, tp, fp, fn  estimated at different threhsolds.
##' You can then use function in minet to estimate roc auc pr auc etc. See ?minet::vis.res
##' @title network_performance: estimate
##' @param true_igraph 
##' @param edges 
##' @return data.frame generated by validate function in minet package. 
##' @author Fred Gruber
network_performance = function(true_igraph, edges){
    require(minet)
    checkmate::assertClass(true_igraph, 'igraph')
    checkmate::assertDataFrame(edges)
    est_ig=igraph::graph_from_data_frame(edges)
    est_adj= est_ig %>% igraph::get.adjacency(attr="freq") %>% as.matrix()

    true_adj=igraph::get.adjacency(true_igraph) %>% as.matrix
    ## make sure nodes are ordered the same way
    ## need to complete tis
    est_adj=est_adj[colnames(true_adj),colnames(true_adj)]
    ##get edgelist
    true_edges=reshape2::melt(true_adj)  %>%
        filter(Var1!=Var2) ## remove self edges
    prc_random=sum(true_edges$value)/length(true_edges$value)
    compa=minet::validate(est_adj,true_adj)
    attributes(compa)$auc_pr_random=prc_random
    return(compa)
}


  plot_auc = function(perf, title=""){
    rocgg = minet::rates(perf) %>%
      ggplot(aes(x=fpr, y=tpr)) +
      geom_line() +
      geom_abline(intercept=0,
		  slope=1,
		  colour='red',
		  linetype="dashed",
		  alpha=0.5)+ xlab("FP Rate") + ylab("TP Rate") +
      ggtitle(sprintf("ROC AUC=%0.2g",
		      minet::auc.roc(perf)
		      ))
    prgg = minet::pr(perf) %>%
      ggplot(aes(x=r, y=p)) +
      geom_point()+
      geom_hline(yintercept=attributes(perf)$auc_pr_random,
		  colour='red',
		  linetype="dashed",
		 alpha=0.5)+
      xlab("Recall (TP Rate)") + ylab("Precision") +
      ggtitle(sprintf("PR AUC=%0.2g (%0.2g)",
		      minet::auc.pr(perf),
		      attributes(perf)$auc_pr_random
		      ))
    require(patchwork)
    allplot = (rocgg / prgg) +
      plot_annotation(title=paste0("Network Inference Performance. ", title))
    allplot & theme_light()
  }
