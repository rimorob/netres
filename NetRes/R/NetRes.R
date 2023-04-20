#' @title NetRes NetRes class definition
#' @import R6 dplyr PCAtools
#' @importFrom bnlearn as.bn as.igraph parents drop.arc as.graphNEL bn.boot nodes remove.node bn.fit
#' @importFrom parallel makeCluster stopCluster detectCores parLapply parSapply 
#' @importFrom doParallel registerDoParallel
#' @importFrom BiocSingular IrlbaParam
#' @importFrom corrplot corrplot corrplot.mixed
#' @importFrom ggplot2 qplot
#' @importFrom paran paran
#' @importFrom GA ga
#' @importFrom GenSA GenSA
#' @import sparsepca
#' @importFrom NNLM nnmf
#' @importFrom torch nn_module torch_tensor nn_linear nn_relu nn_sequential torch_exp torch_randn optim_adam nn_bce_loss
#' @importFrom rlang ns_env
#' @importFrom minet rates auc.roc auc.pr pr
#' @import ggplot2

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
                    #' @param mode If set to "oracular", uses the (presumably provided) latent space to infer a "perfect" solution - used in benchmarking 
                    #' performance and should be considered a deprecated argument
                    #' @param latentSpaceParentOnly If true (default), latent space variables can only be parents, never children                    #' 
                    #' @return A NetRes object
                    initialize = function(dframe, true.graph = NULL, nIter, nBoot=50, algorithm='tabu', algorithm.args, 
                                          lvPrefix = "^U\\_", mode=NULL,
                                          weightedResiduals = FALSE, scale=FALSE, debug=FALSE,
                                          latentSpaceParentOnly = TRUE,
                                          latentSpaceMethod = 'sparse.pca',
                                          optimizeLatentSpace=FALSE,
                                          nCores=NULL,
                                          BPPARAM=BiocParallel::DoparParam()) {
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
                        if(debug) {
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
                                                         latentSpaceMethod = latentSpaceMethod,debug=debug, BPPARAM=BPPARAM)
                        ##clean up
                        gc()
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
                          print('Stopping early due to convergence')
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
                        
                        self$assess(lvPrefix=lvPrefix,cluster=cluster,ci=TRUE)

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

                    assess = function(true.graph = self$true.graph, lvPrefix = self$lvPrefix, nCores=NULL, return_roc=FALSE, save_to_pdf=NULL, ci=FALSE,Nboot=200, iteration = NULL,oracle=NULL,cluster=NULL) {
                        require(patchwork)
                        require(corrplot)
                        cutoff=0.5

                      ##plot corrplot of inferred vs true latent space, assuming true latent space exists
                       if (!is.null(self$latent.data)) {
                        if (is.null(iteration)) ##plot the last one by default
                          iteration = length(self$latent.space)
                        
                        if ("v" %in% names(self$latent.space[[iteration]])) {
                            ##dev.set(1)
                            message("plotting correlation")
                            corrplot.mixed(cor(cbind(self$latent.data, self$latent.space[[iteration]]$v)), upper='ellipse', order='AOE', insig='blank')
                        }
                      }
                      if (is.null(true.graph)) { #then can't plot other metrics
                        return  
                      }

                      if(is.null(nCores)){
                          nCores = detectCores() - 2
                      }else{
                          nCores=min(detectCores()-2,nCores)
                      }
                      if (is.null(cluster)) {
                        cluster = makeCluster(nCores)
                        registerDoParallel(cluster)                      
                        cleanUpCluster = TRUE
                      } else {
                        cleanUpCluster = FALSE
                      }
                      ##excise the latent space from the true graph
                      true.graph = private$exciseLatVarsFromEnsemble(list(true.graph), cluster, lvPrefix)[[1]]
                      true.graph.ig=as.igraph(true.graph)

                      ##aucs = c()
                      ##prAucs = c()
                      ##f1maxes = c()
                      ##sids=c()
                      ##test performance at every iteration
                      permetrics=NULL
                      allplots=list()
                        Ne=length(self$ensemble)
                        if(iteration=="all")
                            inters=1:Ne
                        else
                            inters=iteration
                      for (ni in inters) { 
                          print(paste('step', ni))
                          curEnsemble = private$exciseLatVarsFromEnsemble(self$ensemble[[ni]], cluster, lvPrefix)
                          curStrength = bnlearn::custom.strength(curEnsemble, bnlearn::nodes(true.graph))
                          curStrengthdf=curStrength %>%
                              as.data.frame() %>%
                              mutate(freq=strength*direction)
                          if(!is.null(oracle)){
                              oracleEnsemble = private$exciseLatVarsFromEnsemble(oracle$ensemble[[1]], cluster, lvPrefix)
                              oracleStrength = bnlearn::custom.strength(oracleEnsemble, bnlearn::nodes(true.graph))
                              oracleStrengthdf=oracleStrength %>%
                                  as.data.frame() %>%
                                  mutate(freq=strength*direction)
                              perf=private$network_performance(true.graph.ig,curStrengthdf,ci=ci,cutoff=cutoff,oracle=oracleStrengthdf,Nboot=Nboot)
                          }else{
                              perf=private$network_performance(true.graph.ig,curStrengthdf,ci=ci,cutoff=cutoff,Nboot=Nboot)
                          }
                          ##perf=network_performance(true.graph.ig,curStrengthdf,ci=ci,cutoff=0.5 )
                          allplots[[ni]]=perf+plot_annotation(title=sprintf("Iteration %d",ni))
                          ## get performance metrics
                          allmetr=grep("perf_",names(attributes(perf)),value=T)
                           for(pp in allmetr){
                              val=attributes(perf)[[pp]]
                              permetrics=bind_rows(
                                  permetrics,
                                  tibble(Metric=val$names,Value=val$value,th=cutoff,iteration=ni)
                              )
                          }
                          ## sids[ni]=attributes(perf)$perf_sid$sid
                          ## auc=attributes(perf)$perf_auc_roc
                          ## aucpr=attributes(perf)$perf_aurc_pr
                          ## aucs[ni] = auc
                          ## prAucs[ni] = aucpr
                          ## allmes=attributes(perf)$other %>% as.data.frame()
                          ## f1maxes[ni] = pracma::interp1(
                          ##                           filter(allmes,type=='fscore')$x,
                          ##                           filter(allmes,type=='fscore')$y,0.5,method='linear')
                      }
                      plotstats=permetrics %>%
                              ggplot(aes(x=iteration,y=Value,colour=Metric))+
                              geom_line()+geom_point()+theme_light()+
                              ggtitle(sprintf("Cutoff: %0.2g",permetrics$th[1]))+
                              facet_wrap(~Metric,ncol=2)
                      ## ggp1=qplot(1:length(self$ensemble), aucs, main='AUCs over iterations', xlab='Iteration', ylab='AUC')+geom_line()+theme_light()
                      ## ggp2=qplot(1:length(self$ensemble), prAucs, main='PR-AUCs over iterations', xlab='Iteration', ylab='PR-AUC') +geom_line()+theme_light()                 
                      ##   ggp3=qplot(1:length(self$ensemble), f1maxes, main='F1max values over iterations', xlab='Iteration', ylab='F1max')+geom_line()+theme_light()
                      ##   ggp4=qplot(1:length(self$ensemble), sids, main='Structural Intervention Distance over iterations', xlab='Iteration', ylab='SID')+geom_line()+theme_light()
                      ##   ggpcomb=ggp1/ggp2/ggp3/ggp4
                        if(!is.null(save_to_pdf)){
                            message("saving to file ",save_to_pdf)
                            pdf(save_to_pdf)
                            ##print(ggpcomb)
                            print(plotstats)
                            for(ni in 1:length(allplots)){
                                if(!is.null(self$latent.space[[ni]]$v)){
                                    corrplot.mixed(cor(cbind(self$latent.data, self$latent.space[[ni]]$v)), upper='ellipse', order='AOE', insig='blank') %>% print()
                                }
                                print(allplots[[ni]])
                            }
                            dev.off()
                        }else if(return_roc){
                            allplots
                        }else{
                            message("Plotting metrics")
                            print(
                                allplots[[length(allplots)]] | ~corrplot.mixed(cor(cbind(self$latent.data, self$latent.space[[length(allplots)]]$v)), upper='ellipse', order='AOE', insig='blank')
                            )
                        }
                        if (cleanUpCluster)
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
                    # @param optimizeLatentSpace If false (default), return first-ass estimate of the latent space without conditional BIC-based optimization (which is very slow)
                    # @return TBD
                    runOneIteration  = function(dframe, nBoot, algorithm, algorithm.args, cluster, lvPrefix = "^U\\_", 
                                                weightedResiduals=FALSE, scale=FALSE, learnLatentSpace = TRUE,
                                                optimizeLatentSpace=FALSE, 
                                                latentSpaceParentOnly=TRUE, latentSpaceMethod = 'pca',
                                                BPPARAM=BiocParallel::DoparParam(),
                                                debug=FALSE) {
                        if(debug){
                            message("start of runOneIteration")
                            browser()
                        }
                      dud = function(x) x
                      
                      new.algorithm.args = algorithm.args
                      new.algorithm.args$blacklist = rbind(algorithm.args$blacklist,
                                                       private$makeBlacklist(dframe, lvPrefix))

                      ##This part seems to be slow, at least on small networks.  Parallelizes to only a couple of cores at a time
                      ##It would be nice to figure out why
                      ens = bn.boot(dframe, statistic = dud, R=nBoot, algorithm = algorithm, algorithm.args = new.algorithm.args, cluster=cluster, debug=TRUE)
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

                      ##if latent space was present, then let's determine the residuals
                      res = private$calculateResiduals(ens2, netWeights, weightedResiduals, cluster)
                      
                      ##paran(res)
                      latvars = private$calculateLatVars(res, method=latentSpaceMethod, scale=scale, algorithm.args=new.algorithm.args, BPPARAM=BPPARAM) 

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

                      if (ncol(latCoefs) > 1 && optimizeLatentSpace) { 
                          print('optimizing the basis vector of the latent space')

                          ##learn the markov boundary of the latent vars to speed up the process
                          mb = lapply(as.data.frame(latvars$v), function(lv, dframe) {
                              lv = data.frame(lv)
                              colnames(lv) = 'LV'
                              dframe.tmp = cbind(dframe, lv)
                              learn.nbr(dframe.tmp, colnames(lv), 'mmpc')
                          }, dframe)
                          mb = unique(unlist(mb))

                        startTime = Sys.time()
                        optimRes = ga(type = 'real-valued', fitness=private$latVarLinComb, latVars = latvars$v, 
                                      algorithm = algorithm,
                                      algorithm.args = algorithm.args,
                                      lvPrefix = lvPrefix,
                                      ##use markov boundary for optimization
                                      dframe = dframe[, mb],
                                      cluster = cluster,
                                      nBoot = 1, #"fast" regime
                                      lower = rep(-100, length(as.numeric(latCoefs))),
                                      upper = rep(100, length(as.numeric(latCoefs))),
                                      run = algorithm.args$max.iter/5, 
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

                          ##optimization may have been done on the Markov boundary, but now fit the whole ensemble
                        remappedRes = private$latVarLinComb(optimRes@solution, latvars$v, algorithm = algorithm,
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

