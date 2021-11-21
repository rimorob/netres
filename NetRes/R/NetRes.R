#' @title NetRes NetRes class definition
#' @import R6 dplyr PCAtools
#' @importFrom bnlearn as.bn as.igraph parents drop.arc as.graphNEL bn.boot nodes remove.node bn.fit
#' @importFrom parallel makeCluster stopCluster detectCores parLapply parSapply
#' @importFrom BiocSingular IrlbaParam
#' @importFrom corrplot corrplot
#' 
##An R6 class for a generated scale-free network
#' @export
NetRes <- R6Class("NetRes", 
                  public = list(
                    #' @field ensemble The final ensemble as inferred after some iterations
                    ensemble = NULL,
                    #' @field latent.space The list of latent spaces after a given iteration
                    latent.space = NULL,
                    #' @field train.data The initial training data
                    train.data = NULL,
                    #' @field latent.data The true data for the latent space (if provided)
                    latent.data = NULL,
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
                    #' @return A NetRes object
                    initialize = function(dframe, true.graph = NULL, nIter, nBoot=50, algorithm='tabu', algorithm.args, 
                                          lvPrefix = "^U\\_", mode=NULL,
                                          debug=FALSE) {
                        if(debug)
                            browser()
                      self$latent.data = dframe %>% select_if(grepl(lvPrefix, names(.)))
                      self$train.data = dframe %>% select_if(!grepl(lvPrefix, names(.)))
                      self$true.graph = true.graph
                      require(parallel)
                      train = self$train.data
                      stopifnot(is.list(algorithm.args))
                      nCores = detectCores() - 1
                      cluster = makeCluster(nCores)
                      clusterEvalQ(cluster,library(bnlearn))
                        for (ni in 1:nIter) {
                            if(debug){
                                message("In interation ",ni)
                                browser()
                            }
                            curRes = private$runOneIteration(train, nBoot, algorithm, algorithm.args, cluster, lvPrefix = lvPrefix)
                        self$ensemble[[ni]] = curRes$ensemble
                        self$latent.space[[ni]] = curRes$latent.space
                        if (is.null(mode) || mode != 'oracular') {
                          train = cbind(self$train.data, curRes$latent.space$v)
                        } else {
                          warning('!!!RUNNING IN ORACULAR MODE - USING TRUE LATENT VARIABLES FOR INFERENCE OF STRUCTURE!!!')
                          train = cbind(self$train.data, self$latent.data)
                        }
                        corrplot(cor(cbind(self$latent.data, curRes$latent.space$v)), method='ellipse', order='AOE', diag=F)
                        print('pausing to admire the corrplot')
                        Sys.sleep(5)
                      }
                      stopCluster(cluster)
                    },
                    #' @description assess Assess the inferred ensemble against the true graph
                    #' @param true.graph The true graph to use; defaults to the one provided at initialization (if any)
                    #' @param lvPrefix The latent variable-identifying regular expression, as elsewhere; defaults to "^U\\_"
                    assess = function(true.graph = self$true.graph, lvPrefix = "^U\\_") {
                      if (is.null(true.graph)) {
                        stop('Cannot assess performance without the true graph')
                      }
                      
                      nCores = detectCores() - 1
                      cluster = makeCluster(nCores)
                      ##excise the latent space from the true graph
                      true.graph = private$exciseLatVarsFromEnsemble(list(true.graph), cluster, lvPrefix)[[1]]
                      
                      aucs = c()
                      prAucs = c()
                      ##test performance at every iteration
                      for (ni in 1:length(self$ensemble)) {
                        curEnsemble = private$exciseLatVarsFromEnsemble(self$ensemble[[ni]], cluster, lvPrefix)
                        curStrength = custom.strength(curEnsemble, bnlearn::nodes(true.graph))  
                        pred = as.prediction(curStrength, true.graph)
                        perf = ROCR::performance(pred, "tpr", "fpr")
                        auc = round(ROCR::performance(pred, "auc")@y.values[[1]], 2)
                        plot(perf, main = paste("AUC:", auc), colorize=TRUE)
                        aucpr = round(ROCR::performance(pred, "aucpr")@y.values[[1]], 2)                        
                        perf = ROCR::performance(pred, "prec", "sens")
                        plot(perf, main = paste("PR-AUC:", aucpr), colorize=TRUE)     
                        aucs[ni] = auc
                        prAucs[ni] = aucpr
                      }
                      plot(1:length(self$ensemble), aucs, main='AUCs over iterations', xlab='Iteration', ylab='AUC')
                      plot(1:length(self$ensemble), prAucs, main='PR-AUCs over iterations', xlab='Iteration', ylab='PR-AUC')                      
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
                    # @return TBD
                    runOneIteration  = function(dframe, nBoot, algorithm, algorithm.args, cluster, lvPrefix = "^U\\_") {
                      dud = function(x) x
                      
                      algorithm.args$blacklist = rbind(algorithm.args$blacklist,
                                                       private$makeBlacklist(dframe, lvPrefix))
                      ens = bn.boot(dframe, statistic = dud, R=nBoot, algorithm = algorithm, algorithm.args = algorithm.args, cluster=cluster)
                      netWeights = private$calcBayesFactors(ens, dframe, cluster, algorithm.args)
                      ens2 = private$exciseLatVarsFromEnsemble(ens, cluster, lvPrefix)
                      res = private$calculateResiduals(ens2, netWeights, cluster)
                      latvars = private$calculateLatVars(res, method='pca') 
                      return(list(ensemble=ens, latent.space = latvars))
                    },
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
                    # @description The top-level function to calculate latent variables from residuals
                    calculateLatVars = function(residuals, method='pca') {
                      if (method == 'pca') {
                        latVars = private$calculateLatVarsPCA(residuals)  
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
                    },
                    # @description Calculate latent variables using linear PCA
                    # residuals samples x vars matrix of residuals
                    calculateLatVarsPCA = function(residuals) {
                      ##how many PCs are we dealing with?  Use Horn's parallel analysis to find out
                      resPpca = parallelPCA(
                        as.matrix(residuals),
                        max.rank = round(ncol(residuals)*0.1),
                        niters = 50,
                        threshold = 0.1,
                        transposed = TRUE, #variables in columns
                        ##BSPARAM = IrlbaParam(),
                        BPPARAM = BiocParallel::MulticoreParam()
                      )
                      
                      if (resPpca$n == 0) {
                        stop('Latent space not identified - aborting')
                      }
                      
                      resPc = prcomp(residuals, rank=resPpca$n, scale=TRUE)

                      return(list(nLatVars = resPpca$n, 
                                  latVars = resPc$x,
                                  lvPredictor = resPc))
                    },
                    # @description The private helper function to calculate residuals for a given ensemble
                    calculateResiduals = function(ens, netWeights, cluster) {
                      residuals = parLapply(cluster, ens, function(net) {
                        ##fit a single network using training data
                        fitted = bn.fit(net, self$train.data[, bnlearn::nodes(net)])
                        netRes = residuals(fitted)

                        return(as.data.frame(netRes))
                      })
                      for (ri in 1:length(residuals)) {
                        residuals[[ri]] = residuals[[ri]] * netWeights[ri]
                      }

                      ##average the residuals across all networks (not clear if this is the best way, but likely it is)                      
                      residuals = Reduce(`+`, residuals) / sum(netWeights)                      
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
                            net = remove.node(net, lv)
                          }
                          return(net)
                        })
                      return(ens)
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
