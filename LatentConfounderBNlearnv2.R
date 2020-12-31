##create a private global environment for keeping track of shared variables that are
##not being passed around
##more broadly, this whole show needs to be upgraded to S4 or R6 to allow for class variables
.pkgGlobalEnv <- new.env(parent=emptyenv())

##only variables that are disallowed to be driven by each known variable are fixed
inferFixedVars <- function(data=NULL, blacklist=NULL) {
    if (is.null(data)) {
        data = .pkgGlobalEnv$data
        blacklist = .pkgGlobalEnv$blacklist
    }
    allVars = colnames(data)
    fixed = c()
    for (v in unique(blacklist$to)) {
        otherVars = setdiff(allVars, v)
        drivers = subset(blacklist, to == v)$from
        if (length(setdiff(otherVars, drivers)) == 0) {
            fixed = c(fixed, v)
        }
    }
    return(fixed)
}

isDiscrete = function (x, th = sqrt(length(x)))
{
    if (length(unique(x)) < th) {
	return(TRUE)
    }
    return(FALSE)
}

library(bnlearn)
library(data.table)
library(igraph)
library(Rgraphviz)
library(tidyverse)
getEnsemble = function(train, Nboot = 50, algorithm = "hc", cluster = NULL, output = "Z", ...){
    message("Bootstrapping bnlearn ", Nboot, " times.")
    pdaboot = bn.boot(train, statistic = function(x) x,
		      algorithm.args = list(...),
		      algorithm = algorithm, R = Nboot, cluster = cluster)
    message("Calculating edges")
    alledge = purrr::map_df(1:length(pdaboot),
			    function(ff){
				arc.strength(pdaboot[[ff]],train) %>% as.data.frame %>%
				    mutate(Boot = ff)
			    })
    alledge = alledge %>%
	dplyr::group_by(from, to) %>%
	dplyr::summarise(freq = n() / Nboot) %>%
	arrange( -freq)
    message("Fit all networks")
    ## first marginalize latent variables if present
    latvars = grep("LV_.*", colnames(train), value = TRUE)
    for(ll in latvars){
	train[[ll]] = mean(train[[ll]])
    }
    allmod = purrr::map(pdaboot,
			function(ff){
			    bn.fit(ff, train)
			})
    message("calculate coefficients")
    allcoef = purrr::map_df(1:length(allmod),
			    function(ii){
				allmod[[ii]][[output]]$coefficients %>% t() %>%
				    as.data.frame %>% mutate(Boot = ii)
			    })
    allcoef[is.na(allcoef)] = 0
    structure(list(
	Nboot = Nboot,
	algorithm = "hc",
	other_params = list(...),
	edges = alledge,
	coef = allcoef,
	fitmodels = allmod,
	allnet = pdaboot,
	data = as.tibble(train)),
	class = 'bnlearn_ens')
}


## print.bnlearn_ens = function(obj){
##     message("Bnlearn Ensemble")
##     message("\nNboot: ", obj$Nboot)
##     message("\nAlgorithn: ", obj$algorithm)
##     message("\nTraining Data: ")
##     print(preview(obj$data))
##     if(length(obj$other_params) > 0){
##         message("\nOther Paramters: ")
##         for(oo in 1:length(obj$other_params)){
##             message("\n", names(obj$other_params)[oo], ":")
##             head(obj$other_params[[oo]]) %>% print
##         }
##     }
##     message("Other fields:")
##     print(setdiff(names(obj), c("Nboot", "algorithm", "other_params")))
## }

## plot.bnlearn_ens = function(obj, output, ensid = 0, cutoff = 0.5, maxpath = 3,direction = 'upstream', ...){
##     if(ensid == 0)
##         ig = igraph::graph_from_data_frame(
##                          dplyr::filter(obj$edges, freq >= cutoff))
##     else
##         ig = igraph::igraph.from.graphNEL(bnlearn::as.graphNEL(obj$allnet[[ensid]]))
##     if(direction == 'upstream')
##         mode = 'in'
##     else if(direction == "downsteam")
##         mode = 'out'
##     else if(direction == "both")
##         model = "all"
##     else
##         stop("direction not recognized. Only: 'upstream','downsteam', and 'both'.")
##     parents = names(igraph::neighborhood(ig, order = maxpath, output, mode = mode)[[1]])
##     sig = igraph::induced.subgraph(ig, vids = parents)

##     ppPlotIgraph(ig, ...)
## }

library(bnlearn)
  library(foreach)
  library(igraph)
getEnsemble2 = function(train, Nboot = 50, algorithm = "hc",parallel = FALSE, output = NULL, ...){
    if(!is.null(output)){
        if(!output %in% colnames(train)){
            stop("output ", output, " is not in train")
        }
    }
    if(parallel){
            print(paste('Distributing ensemble learning'))
	  `%op%` <-  `%dopar%`
    }else{
        print('Doing ensemble learning locally')
	  `%op%` <-  `%do%`
      }
    message("Bootstrapping bnlearn ", Nboot, " times.")
    
    mycomb = function(xx, yy) {
        list(allnet = c(xx[[1]], yy[[1]]), fitmodels = c(xx[[2]], yy[[2]]),
             edges = rbind(xx[[3]], yy[[3]]),
             boot_order = c(xx[[4]], yy[[4]]))
    }
    pdaboot = foreach(ii = 1:Nboot,
                      .combine = mycomb, .packages = c("bnlearn", "foreach", "igraph","dplyr")) %op% {
                          bootorder = sample(1:nrow(train), nrow(train), replace = TRUE)
                          trainb = train[bootorder, ]
                          if(algorithm == "hc") {
                              net = bnlearn::hc(trainb, ...)
                          } else if (algorithm == 'tabu') {
                              net = bnlearn::tabu(trainb, ...)
                          } else {
                              stop("Unknown Algorithm. Only hc and tabu.")
                          }
                          mod = bn.fit(net, train[bootorder, ])
                          ## fit edges
                          edgedf = arc.strength(net,trainb) %>% as.data.frame %>%
                              mutate(Boot = ii)
                          list(allnet=list(net),
                               fitmodels=list(mod),
                               edges=edgedf,
                               boot_order=list(bootorder))
                      }

      pdaboot$edges = pdaboot$edges %>%
	  dplyr::group_by(from, to) %>%
	  dplyr::summarise(freq = n() / Nboot) %>%
	  arrange( -freq) %>%
	  ungroup()
	## pdaboot = bn.boot(train, statistic = function(x) x,
	##   		algorithm.args = list(...),
	##   		algorithm = algorithm, R = Nboot, cluster = cluster)
	## message("Fit all networks")
	## allmod = purrr::map(pdaboot,
	##   	   function(ff){
	##   	       bn.fit(ff, train)
	##   	   })
    message("calculate coefficients")
    if(!is.null(output)){
        allcoef = purrr::map_df(1:length(pdaboot$fitmodels),
                                function(ii){
                                  pdaboot$fitmodels[[ii]][[output]]$coefficients %>% t() %>%
				 as.data.frame %>% mutate(Boot = ii)
			 })
      allcoef[is.na(allcoef)] = 0
    }else{
        allcoef = NULL
      }
      structure(list(
	  Nboot = Nboot,
	  boot_orders = pdaboot$boot_order,
	  algorithm = algorithm,
	  other_params = list(...),
	  edges = pdaboot$edges,
	  coef = allcoef,
	  fitmodels = pdaboot$fitmodels,
	  allnet = pdaboot$allnet,
	  data = as.tibble(train)),
	  class = 'bnlearn_ens')
  }

  print.bnlearn_ens = function(obj){
	message("Bnlearn Ensemble")
	message("\nNboot: ", obj$Nboot)
	message("\nAlgorithn: ", obj$algorithm)
	message("\nTraining Data: ")
	head(as.tibble(obj$data))
	if(length(obj$other_params) > 0){
	    message("\nOther Paramters: ")
	    for(oo in 1:length(obj$other_params)){
		message("\n", names(obj$other_params)[oo], ":")
		head(obj$other_params[[oo]]) %>% print
	    }
	}
	message("Other fields:")
	print(setdiff(names(obj), c("Nboot", "algorithm", "other_params")))
    }

plot.bnlearn_ens = function(obj,
                            output,
                            ensid = 0,
                            freqth = 0.5,
                            cutoff = 0.5,
                            maxpath = 3,
                            coef_filter = 0.05, 
                            direction = 'upstream',
                            edge_labels = "none",
                            labels_regex = output, ...){
	if(ensid == 0)
	    ig = igraph::graph_from_data_frame(
			     dplyr::filter(obj$edges, freq >= freqth))
	else
	    ig = igraph::igraph.from.graphNEL(bnlearn::as.graphNEL(obj$allnet[[ensid]]))
        allnodes = names(igraph::V(ig))
	if(direction == 'upstream')
	    mode = 'in'
	else if(direction == "downstream")
	    mode = 'out'
	else if(direction == "both")
	    mode = "all"
	else
	    stop("direction not recognized. Only: 'upstream','downstream', and 'both'.")
	##parents = names(igraph::neighborhood(ig, order = maxpath, output, mode = mode)[[1]])
        if (!missing(output)) {
            parents = getDrivers(obj,
                                 output = output,
                                 maxpath = maxpath,
                                 cutoff = cutoff,
                                 direction = direction)$Drivers
            subnode = c(output, parents)
            subnode = subnode[subnode %in% allnodes]
        } else {
            subnode = allnodes
        }
        if(length(subnode > 1))
            sig = igraph::induced.subgraph(ig, vids = subnode)
        else
            sig = ig
        if(edge_labels == 'none'){
            edgelabels = FALSE
            edge_labels = NULL
        }else if (edge_labels == 'frequency'){
            edgelabels = TRUE
            edge_labels = NULL
        }else if(edge_labels == "coefficients" | edge_labels == 'coef'){
            edgelabels = TRUE
            edge_labels = getCoef(obj,labels_regex, as.regex = TRUE) %>%
                filter(wcoef>coef_filter) %>%
                dplyr::select(input,output,wcoef) %>%
                mutate_if(is.numeric, signif, 2)
        }else{
            stop("wrong edge_labels. Only 'none', 'frequency','coefficients'")
        }
	plotIgraph(sig,edgelabels = edgelabels,edge_labels = edge_labels, ...)
    }

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

getGraphResiduals = function(obj, variables, data = NULL, Nsamples = 1){
    Nboot = length(obj$fitmodels)
    ## first marginalize latent variables if present
    if(is.null(data))
        trainmarg = obj$data
    else
        trainmarg = data
    latvars = grep("LV_.*", colnames(obj$data), value = TRUE)
    if(length(latvars) > 0){
        message("Marginalizing latent variables")
        for(ll in latvars){
            if(ll %in% colnames(trainmarg))
                trainmarg[[ll]] = mean(trainmarg[[ll]])
            else
                trainmarg[[ll]] = mean(obj$data[[ll]])
        }
    }
    Ndatasamp = nrow(trainmarg)
    res = array(rep(0, Nboot * Ndatasamp * Nsamples),
                dim = c(Ndatasamp, Nboot, Nsamples),
                dimnames = list(paste0("samp", 1:Ndatasamp),
                                paste0('Boot', 1:Nboot),
                                paste0("Sample", 1:Nsamples)))
    resall = list()
    ## bootstraps
    for(bb in 1:Nboot){
	bootdata = obj$fitmodels[[bb]]
	## remove input only
	net = obj$allnet[[bb]]
	inonly = root.nodes(net)
	vars = names(bootdata)
	vars = vars[vars %in% variables]
	vars = setdiff(vars,
		       setdiff(
			   inonly,
			   latvars
		       ))
	for(vv in vars){
	    tmp = bootdata[[vv]]
            mod = as.lm(tmp, data = trainmarg)
            pred = predict(mod, trainmarg, se.fit = TRUE)
            sds = pred$se.fit*sqrt(nrow(trainmarg))
            ##pred = predict(bootdata,vv, trainmarg)
	    cres = resall[[vv]]
	    if(is.null(cres))
		cres= res
            if(Nsamples > 1){
                for(ss in Nsamples){
                    samp = rnorm(nrow(trainmarg), pred$fit, sds)
                    cres[, paste0("Boot", bb), paste0("Sample", ss)] = samp
                }
            }else{
                cres[, paste0("Boot", bb), "Sample1"] = pred$fit
            }
	    resall[[vv]] = cres
	}
    }
    list(variables = names(resall), predProb = resall)
}

##trainingData: raw training data
##modelPrediction: for now, only the maximum likelihood prediction;
##  eventually should implement sampled posterior implementation
##at present, discretizes all non-discrete data
##dim(X) = samples x variables
##isOrdinal: TRUE if ordinal data, FALSE if categorical; this is not an optional argument
##  Ultimately, we need a better way of figuring out what's ordinal (per variable)
residualDeviance <- function(trainingData, modelPrediction, isOrdinal, missing_code = -1) {
    ##iterate over all variables
    devX = trainingData[, modelPrediction$variables]
    for (vi in 1:length(modelPrediction$variables)) {
        ep = modelPrediction$variables[[vi]]
        Xbar = modelPrediction$predProb[[vi]] #model

        if (is.null(dim(Xbar))) { #just one state
            dim(Xbar) = c(length(Xbar), 1)
        }
        print(paste('vi:', vi))
        X = trainingData[[ep]] #data

        ##For discrete ordinal data:
        ##compute distance from this paper (Shepherd et al):
        ##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5364820/
        ##these are probability-scale residuals:
        ##predicted mass above the true value minus predicted mass
        ##below the true value
        ##For discrete data:
        ##Old method:
        ##compute KL divergence between 1 (training data) and p(training data|model)
        ##New method:
        ##compute PSR based on class/out of class designation
        dX = X
        for (ri in 1:nrow(Xbar)) {
            if (isOrdinal) {
                ## for binary outcomes
                if(length(unique(X)) == 2){
                    ## According to section 4.1 of Shepherd et al for binary cases the PSR reduces to
                    ## y-pr(Y=1).
                    ## It is not clear how to best get the probability of one:
                    ## - should I use the 2nd column as I'm doing here
                    ## - or should  I use column called "1"
                    ##
                    dX[ri] = X[ri] - Xbar[ri, 2, 1] ##note the extra column from GNS REFS usage
                } else if (isDiscrete(X)) {
                    ##probably broken because of posterior dimensionality; needs to be fixed
                     ##will crash as-is
                    ## assuming the level on the first column is the lowest
                    ## FG: fixed now.
                    ## FG: assuming

                    ##BH (1/31/20): changing in the following ways:
                    ##1. ordinal variable is assumed to be represented as integer levels
                    ##(this should change into as.ordered() representation later)
                    ##2. results are presumed to be numeric

                    if ((is.na(missing_code) & is.na(X[ri])) ||
                        (!is.na(missing_code) & X[ri] == missing_code)) {
                        dX[ri] = NA
                    }else{
                        ## remove -1 column if present (BH: what is the -1 level???)
                        selcol = colnames(Xbar)
                        selcol = selcol[selcol != "-1"]
                        Xbar = Xbar[, selcol, 1, drop=F]

                        massBelow = length(which(Xbar[ri,,] < X[ri]))
                        massAbove = length(which(Xbar[ri,,] > X[ri]))
                        dX[ri] = massBelow - massAbove #as per definition of PSR in Shepherd et al
                    }

                } else { #continuous
                    ##currently compute probabilities based on counts
                    ##but consider using densities to be more precise
                    ##on the other hand, that may not be robust in the
                    ##presence of outlie""rs
                    ##flatten xbar, which is data samples x network samples x networks, but may want to eventually compute this per network to check which is better
                    massAbove = length(which(as.numeric(Xbar[ri,,]) > X[ri]))/length(as.numeric(Xbar[ri,,]))
                    massBelow = length(which(as.numeric(Xbar[ri,,]) < X[ri]))/length(as.numeric(Xbar[ri,,]))
                    dX[ri] = massBelow - massAbove #as per definition of PSR in Shepherd et al
                }



            } else {
                if (FALSE) { #old method - KLD
                    ##for each row, we assume that probability of observed data was 1
                    ##(this is different than saying that every pX(X_ri == k) = p(k) because it allows for
                    ##(unobserved) context
                    ##thus deviance happens from uniform distribution
                    ##That is, P(data) * log (P(data)/P(data|model)) =
                    ##1 * log (1 / P(data = k|model))
                    dX[ri] = 1*log(1/Xbar[ri, as.character(X[ri])])
                } else { #new method - related to PSR (for consistency)
                    ##different from PSR in that it is not signed
                    ##signed method allows "tuning" of sensitivity vs specificity to be treated as
                    ##bias, but for nominal variables we only know "right" vs "wrong"

                    rightCt = length(which(Xbar[ri,,] == X[ri]))
                    allCt = prod(dim(Xbar[ri,,]))

                    dX[ri] = (allCt - rightCt)/allCt
                }
            }
        }

        ##drop infinity states and warn
        if (length(which(is.infinite(dX))) > 0) {
            dX[is.infinite(dX)] = 0
        }
        devX[, vi] = dX
    }
    ##order variables by name
    allvars = sort(colnames(devX))
    devX = devX[, allvars]
    return(devX)
}

##Identify latent variables using an autoencoder
##Accepts:
## data - an m samples by n features matrix/data frame from which to fit the latenspace
## architecture - the list of encoder layer geometries, left to right
##   This list will be mirrored in the decoder.
##   Last layer is PC layer and should be set from PCA or similar
##Returns a predict-enabled object that can predict the values for new data
##  We do not try to sample w/variance - not using the variational approach
##  Variational autoencoders seem to work worse
##For now, train and validate on the same data, but consider doing cross-validation later
fitAeLatent <- function(data, architecture, lossFunction = 'mean_squared_error',
                        optimizer = 'RMSprop', metrics = 'mean_squared_error', drRate = 0.2,
                        summarize=TRUE, fname='aeRunSummary.pdf', valData=NULL,
                        epochs=100, batch_size = NULL, activation = 'sigmoid', activation_coding = "sigmoid", activation_output = "linear", use_batch_norm = TRUE, learning_rate=0.001, validation_split = 0.1, callback = NULL) {
    ## reticulate::conda_list()
    ## reticulate::use_condaenv("tensorflow_p36")
    stopifnot(keras::is_keras_available())
    library(ggplot2)
    library(keras)
    library(tidyverse)
    library(caret)
    library(reticulate)
    x_train = as.matrix(data)
    if (is.null(valData)) {
        x_test = as.matrix(data)
    } else {
        x_test = as.matrix(valData)
    }
    if (summarize) {
        print(paste('Analyzing a matrix of', nrow(x_train), 'samples x', ncol(x_train), 'features'))
    }
    input_dim = ncol(x_train)
    input_layer <- layer_input(shape = c(input_dim))
    ## dense1 = layer_dense(units = architecture[1])
    ## dense2 = layer_dense(units = architecture[2])
    ## denselat = layer_dense(units = architecture[length(architecture)])
    ## dense1_a = layer_tied_dense(master_layer = dense1)
    ## dense2_a = layer_tied_dense(master_layer = dense2)
    ## ouput = input_layer %>%
    ##     dense1() %>%
    ##     layer_activation(activation) %>%
    ##     dense2() %>%
    ##     layer_activation(activation) %>%
    ##     denselat() %>%
    ##     layer_activation(activation) %>%
    ##     dense2_a() %>%
    ##     layer_activation(activation) %>%
    ##     dense1_a()
    ##encoder <- input_layer
    ##construct the encoder
    encoder = Sequential()
    mylayers = list()
    if(length(architecture) > 1){
        for (li in 1:(length(architecture)-1)) {
            curL = architecture[li]
            curLs = as.character(curL)
            if(li == 1){
                linput = as.integer(input_dim)
                mylayers[[as.character(li)]] = Dense(as.integer(curL),
                                       activation = activation,
                                       input_shape = list(linput),
                                       use_bias = TRUE
                                       )
            }else{
                linput = as.integer(architecture[li - 1])
                mylayers[[as.character(li)]] = Dense(as.integer(curL),
                                       activation = activation,
                                       use_bias = TRUE
                                       )
            }
            encoder$add(mylayers[[as.character(li)]])
            ## layer_dense(units = curL,
            ##                                 kernel_constraint=keras::constraint_unitnorm())
            ## encoder <-
            ##     encoder %>%
            ##     mylayers[[curLs]]() %>%
            ##                     layer_activation(activation)
            if(use_batch_norm)
                encoder$add(BatchNormalization())
            ##     encoder = encoder %>% layer_batch_normalization()
            encoder$add(Dropout(rate = drRate))
        }
    }
    ##add middle layer or coding layer
    ##This should probably be PCA-constrained or at worst constrained from Donoho:
    ##https://arxiv.org/abs/1305.5870
    ii = length(architecture)
    middla = architecture[length(architecture)] %>% as.integer
    middlal = as.character(middla)
    linput = as.integer(architecture[length(architecture) - 1])
    mylayers[[as.character(ii)]] = Dense(middla,
                                         activation = activation_coding,
                                         kernel_regularizer = WeightsOrthogonalityConstraint(middla, weightage=1., axis=0L),
                                         kernel_constraint=UnitNorm(axis=0L),
                                         activity_regularizer=UncorrelatedFeaturesConstraint(middla, weightage = 1.)
                                         )
    encoder$add(mylayers[[as.character(ii)]])
    if(use_batch_norm)
        encoder$add(BatchNormalization())
    Nlayenc = length(encoder$layers)
    ## start decoder
    message("build decoder")
    mylayersdec = list()
    ## mylayersdec[[middlal]] = DenseTied(linput,
    ##                         tied_to = mylayers[[middla]],
    ##                         activation = activation,
    ##                         use_bias = TRUE
    ##                         )
    ## encoder$add(mylayers[[middladt]])
    ## encoder = encoder %>%
    ##     layer_dropout(rate = drRate)
    ##browser()
    ##construct the decoder
    decoder <- encoder
    if(length(architecture) > 1){
        myarch = rev(architecture)
        alliids = rev(1:length(architecture))
        for (li in 1:(length(myarch) - 1)) {
            ##curL = myarch[li + 1]
            ii = as.character(alliids[li])
            ##curLn = as.character(curL)
            clin = as.integer(myarch[li + 1])
            myl = DenseTied(clin,
                            tied_to = mylayers[[ii]],
                            activation = activation,
                            use_bias = TRUE
                            )
            ##myl = layer_dense(units = curL,
            ##                              kernel_constraint=keras::constraint_unitnorm())
            ##myl = layer_tied_dense(master_layer = mylayers[[curLn]])
            mylayersdec[[li]] = myl
            decoder$add(myl)
            if(use_batch_norm)
                decoder$add(BatchNormalization())
            ##     encoder = encoder %>% layer_batch_normalization()
            decoder$add(Dropout(rate = drRate))
        }
        ## output layer: for the last layer use linear activation
        ii = as.character(alliids[length(alliids)])
        ##curL = myarch[length(myarch)]
        ##curLn = as.character(curL)
        clin = as.integer(input_dim)
        myl = DenseTied(clin,
                    tied_to = mylayers[[ii]],
                    activation = activation_output,
                    use_bias = TRUE
                    )
        mylayersdec[[ii]] = myl
        decoder$add(myl)
    }
    ##decoder$add(Dense(as.integer(input_dim), activation = "linear"))
    message("compile")
    decoder$compile(
                loss=lossFunction,
                optimizer=optimizer,
                metrics = list(metrics)
                )
    if (summarize) {
        print(summary(decoder))
    }
    message("Fit")
    if(is.null(callback)){
            history <- decoder$fit(
                    x_train,
                    x_train,
                    epochs = as.integer(epochs),
                    batch_size=as.integer(ifelse(is.null(batch_size),
                                                 nrow(x_train),
                                                 batch_size)),
                    shuffle=TRUE,
                    validation_split = validation_split
                )
    }else{
            history <- decoder$fit(
                    x_train,
                    x_train,
                    epochs = as.integer(epochs),
                    batch_size=as.integer(ifelse(is.null(batch_size),
                                                 nrow(x_train),
                                                 batch_size)),
                    shuffle=TRUE,
                    validation_split = validation_split,
                    callbacks = callback
                )
    }

    if (summarize) {
        message('Saving training history in ', fname)
        ##pdf(file = fname, width=11 )
        ggsave(fname, plot_history(history, metrics), width=11)
        ##dev.off()
    }
    ##self.encoder = Model(inputs=self.autoencoder.input,
    ##  outputs=self.autoencoder.get_layer('encoder').output)
    ##return only the encoder!
    message("building encoder with trained layers")
    encoder = Sequential()
    for(li in 1:Nlayenc){
        ll = decoder$layers[[li]]
        encoder$add(ll)
    }
    message("Finished autoencoder fitting. returning.")
    encoder$history = history
    return(encoder)
}

##Find latent PCs in the deviance residuals matrix (formally, this is a type of ICA on residuals)
##resDev: samples x vars matrix
##scale should probably be FALSE by default because, by extension of Renyi theorem, KL-divergence
##  should be insensitive to reparameterization (discretization).  This needs to be finalized
##method: kernel or linear

findLatentVars <- function(resDev, scale. = FALSE, nIter = 100, method = 'kernel', kernel='rbfdot', alpha = 0.05, pval_method = 'beta',
                           missing_encoding = NA, ## how to encode missing at the end.
                           discretize_confounders = FALSE,
                           method_disc = 'cluster',
                           number_of_groups = 3,
                           maxLatentVars = Inf,
                           multiple_comparison_correction = TRUE,
                           architecture = NULL, activation = 'sigmoid',
                           activation_coding = "sigmoid",
                           activation_output ="linear",
                           drRate=0.2,
                           use_batch_norm = TRUE, batch_size = 32,
                           optimizer = 'RMSprop',
                           metrics = 'mean_squared_error', learning_rate = 0.1,validation_split = 0.1,
                           fname = "aeRunSummary.pdf",
                           callback = NULL
                           ) {
    ## check for missingness
    missing = which(!complete.cases(resDev))
    resDev = resDev[complete.cases(resDev), ]
    stopifnot(method %in% c('kernel', 'linear', 'robustLinear', "ica", "autoencoder"))
    print(paste('Starting decomposition at time', date()))
    prVar = rep(0, ncol(resDev))
    if (method == 'linear') {
        pr = prcomp(resDev, scale. = scale.)
        ##$x contains principal components in columns and the rotation matrix
        ##specifies the recipe, data %*% rotation, to obtain them,
        ##subject to scale(data) if scale. is TRUE
        prVar[1:length(prVar)] = stats:::summary.prcomp(pr)$importance['Proportion of Variance',]
    } else if (method == 'robustLinear') {
        library(rpca)
        pr = rpca(as.matrix(resDev))
        ##prVar[1:length(prVar)] = prop.table(svd(pr$L)$d^2)
        sdv = svd(scale(resDev))$d / sqrt(nrow(resDev) - 1)
        prVar[1:length(prVar)] = sdv^2/sum(sdv)
    } else if(method == "ica"){
        message("Using ICA for reconstruction of latent variables")
        library(ica)
        ##do a quick run of linear to see how many components to expect
        tmp = findLatentVars(resDev, scale., nIter, method = 'linear', multiple_comparison_correction = multiple_comparison_correction)
        nLatent = min(maxLatentVars, ncol(tmp$confounders))
        pr = icafast(resDev, nLatent)
        prVar[1:length(prVar)] = pr$vafs
    } else if(method == "pica"){
        message("Using ICA for reconstruction of latent variables")
        library(ica)
        ##do a quick run of linear to see how many components to expect
        tmp = findLatentVars(resDev, scale., nIter, method = 'linear', multiple_comparison_correction = multiple_comparison_correction)
        nLatent = min(maxLatentVars, ncol(tmp$confounders))
        pr = icafast(resDev, nLatent)
        prVar[1:length(prVar)] = pr$vafs
    } else if (method == 'autoencoder') {
        ##do a quick run of linear to see how many components to expect
        tmp = findLatentVars(resDev, scale., nIter, method = 'linear', multiple_comparison_correction = multiple_comparison_correction)
        nLatent = ncol(tmp$confounders)
        if (nLatent == 0) {
            stop('There seems to be no evidence for latent variables - stopping')
        }
        ##hacked outer dimension - min of 200 or half # of vars
        ##hacked # of layers
        if(is.null(architecture))
            ##be more parsimonious than PCA, use nLatent/2 for the middle layer
            ##PCA tends to be very broad
            architecture = rev(round(exp(seq(log(max(nLatent/2, 2)), log(min(100, ncol(resDev)/2)), length.out=4))))
        encoder <- fitAeLatent(data=resDev, architecture=architecture, valData=resDev, epochs=nIter, activation = activation, drRate = drRate, use_batch_norm = use_batch_norm, batch_size = batch_size, metrics = metrics, optimizer = optimizer, learning_rate = learning_rate,activation_coding = activation_coding, activation_output = activation_output,validation_split = validation_split,
                               fname = fname, callback = callback)
        encoded = encoder$predict(as.matrix(resDev)) %>% as.data.frame
        pr = encoder
        pr$x = encoded #append the predicted latent space to pr
    } else {
        library(kernlab)
        ##do a quick run of linear to see how many components to expect
        tmp = findLatentVars(resDev, scale., nIter, method = 'linear', multiple_comparison_correction = multiple_comparison_correction)
        nLatent = ncol(tmp$confounders)
        ##pr = kha(as.matrix(resDev), kernel = kernel, features = nLatent, kpar=list())
        pr = kpca(as.matrix(resDev), kernel = kernel, features = nLatent)
        ##pr = kpca(as.matrix(resDev), kernel = 'polydot', kpar=list(degree = 3, scale = scale., offset = 1))
        ##pr = kpca(as.matrix(resDev), kernel = 'splinedot', kpar=list())
        ##prVar = pr@eig/sum(pr@eig) ##override the array size!!!
        prVar = pr@eig ##override the array size!!!  do not normalize!!!
    }
    print(paste('Done at time', date()))
    ##now shuffle and get the background distribution for this
    print('Shuffling')


    if (method == 'autoencoder') {
        print('No shuffling to be done for an autoencoder for now - using dimensionality of shuffled PCA for the time being')
    } else {
        ##for (ii in 1:nIter) {
        prVarShuf = foreach (ii = 1:nIter, .combine='rbind') %do% {
            print(paste('Shuffling iteration', ii, 'of', nIter))
            rdShuf = sapply(as.data.frame(resDev), sample)
            if (method == 'linear') {
                prShuf = prcomp(rdShuf, scale. = scale.)
                prVarShuf = stats:::summary.prcomp(prShuf)$importance['Proportion of Variance',]
                ##prVarShuf = rbind(prVarShuf, stats:::summary.prcomp(prShuf)$importance['Proportion of Variance',])
            } else if (method == 'robustLinear') {
                prShuf = rpca(as.matrix(rdShuf))
                ##prVarShuf = rbind(prVarShuf, prop.table(prShuf$L.svd$d^2))
                ##prVarShuf = rbind(prVarShuf, prop.table(svd(prShuf$L)$d^2))
                sdv = svd(scale(rdShuf))$d / sqrt(nrow(rdShuf) - 1)
                prVarShuf = sdv^2/sum(sdv)
                ##prVarShuf[ii, ] = prop.table(svd(prShuf$L)$d^2)
            }else if (method == 'ica') {
                stop('ICA does not support prediction out of sample easily and has been disabled')
                prShuf = icafast(scale(rdShuf), nLatent)
                prVarShuf = prShuf$vafs
            } else {
                prShuf = kpca(as.matrix(rdShuf), kernel = kernel, features = nLatent)
                ##prShuf = kha(as.matrix(rdShuf), kernel = kernel, features = nLatent)
                eig = prShuf@eig
                if (length(eig) < length(pr@eig)) {
                    eig[length(eig):length(pr@eig)] = pr@eig[length(eig):length(pr@eig)]
                } else if (length(eig) > length(pr@eig)) {
                    eig = eig[1:length(pr@eig)]
                }
                print(eig)
                ##prVarShuf = rbind(prVarShuf, eig/sum(eig))
                ##prVarShuf[ii, ] = eig/sum(eig)
                prVarShuf = eig ##do not normalize
                print(prVarShuf)
            }
            matrix(as.numeric(prVarShuf), ncol = length(prVarShuf))
        }## finish foreach
    }
    if (method == 'kernel') {
        prVarShuf = prVarShuf[, 1:nLatent, drop = F]
    }
    if (method == 'autoencoder') {
        p.val = NA
        sigPc = 1:ncol(pr$x)
    } else {
        ##for each principal component, estimate its p-value using beta distribution
        library(Matrix) #for calculating rank
        p.val = c()
        for (pci in 1:ncol(prVarShuf)) {
            if (pci == rankMatrix(resDev)) { #special case - do not fit!
                p.val[pci] = 1
                next
            }
            ve = prVarShuf[, pci]
            if(pval_method == 'beta'){
                library(fitdistrplus)
                ##param = mledist(ve,"beta")
                param = fitdist(ve, "beta")
                p.val[pci] = pbeta(prVar[pci],
                                   param$estimate['shape1'],
                                   param$estimate['shape2'],
                                   lower.tail=F)
            }else if(pval_method == 'empirical'){
                p.val[pci] = mean(ve >= prVar[pci])
            }else{
                stop("pval_method not recognized. Only valids are: 'beta' and 'empirical'")
            }
        }
        if(any(is.na(p.val)))
            stop("some p-values where NA. Try another pval_method")
        if(multiple_comparison_correction)
            sigPc = which(p.val < alpha/ncol(prVarShuf))
        else
            sigPc = which(p.val < alpha)
    }
    ##alternate = numPcs(resDev) #not used - slightly worse; turn on later?

    if (length(sigPc) == 0) {
        ##stop('No latent confounders detected. Sorry about throwing an exception here - try using tryCatch around this line!!!')
        message("#################################################")
        message("\n\n\n No latent confounders detected. !!\n\n\n")
        message("#################################################")
        return(-1)
    }
    if (method == 'linear') {
        confounders = pr$x[, sigPc, drop=F]
    } else if (method == 'robustLinear') {
        confounders = pr$L.svd$u[, sigPc, drop = F]
        class(pr) = append(class(pr), 'robustLinear')
    } else if(method == 'ica'){
        confounders = pr$S[, sigPc]
    } else if (method == 'autoencoder') {
        confounders = pr$x
        class(pr) = append(class(pr), 'autoencoder')
    }else {
        confounders = pr@pcv[, sigPc, drop=F]
    }
    colnames(confounders) = paste('LV', 1:ncol(confounders), sep='_')
    message("Found ", ncol(confounders), " confounders")
    res = list(confounders = confounders,
        details = list(sigPcIdx = sigPc,
            p.val = ifelse(is.na(p.val), NA, p.val*ncol(prVarShuf)),
            scale. = scale.,
            method = method,
            kernel = kernel,
            maxLatentVars = maxLatentVars,
            alpha = alpha,
            nIter = nIter,
            pval_method = pval_method,
            pcObj = pr,
            originalData = as.matrix(resDev)))
    class(res) <- 'LatentConfounder'

    ##  discretize confounder?
    if(discretize_confounders){ #does this need to be checked for correctness for autoencoders? - NOTE
        message("Discretizing confounders")
        res = as.discrete.LatentConfounder(res,
                                           method = method_disc,
                                           num_clusters = number_of_groups)
    }
    if(length(missing) > 0){
        message("Encoding missing values as ", missing_encoding)
        for(mm in missing){
            res$confounders = insertRow(res$confounders,
                                     rep(missing_encoding, ncol(res$confounders)),
                                     mm)
        }
    }
    return(res)
}

library(foreach)
## Latentconfounder with BNlearn
## The following code does the following:
## - take as input:
##   - run location
##   - otf ensemble location
##   - data
##   - variables to use
##   - nSamples
##   - expand
##   - marginalizeVars
##   - isOrdinal
##   - scale.
##   - method
##   - number of repeats
## - calculates residuals
## - calculate deviance
## - find latent variables
## - add latent variable to dataframe
## - restart otf at last checkpoint temperature
## - repeat
##' .Latent confounder with OTF
##' .. content for \details{} ..
##' @title Latent confounder estimation in OTF
##' @param path path for the OTF run
##' @param ensloc location of ensemble to use
##' @param dataloc location of data to use
##' @param output output variable whose parents will be used for latent variable discovery; if null, all variables will be used (this can be slow for large networks)
##' @param freqCutoff
##' @param maxpath
##' @param expand restrict search to parents of vars, which are presumed to be endpoints, else use just the vars
##' @param nSamples
##' @param isOrdinal
##' @param scale.
##' @param method
##' @param nEpochs number of OTF runs do to
##' @param workdir directory where to store the results
##' @param truelatent data.frame with true latent variables. Enters diagnostic mode where you provide the true latent variable. it will plot a real-time curve of the R^2 estimation of the true latent variable
##' @param truecoef data.frame with true coefficient values: data.frame(input=c("G1","G2"),coef=c(1,2)
##' @param textmode print performance of models in the console.No graphics supports needed
##' @param res results data.frame. in case you are making multiple repetiions of the run.    usually you would take the previous run of the code, say latvars, and then assign res to latvars$details$Diagnostics.
##' @param typemar zero, first, mean. How to marginalized variables
##' @param seed random seed set in R and OTF run
##' @param storeDirs FLAG whether to store each of the directories created during the run. Useful for debugging purposes
##' @param vars variable to include in analysis
##' @return LatentConfounder object

latentDiscovery = function(
			   ens,
			   data,
			   output = NULL,
			   freqCutoff = 0.1,
			   maxpath = 1,
			   alpha = 0.05,
			   expand = F,
			   nSamples = 5,
			   isOrdinal = T,
			   scale. = T,
			   method = "linear",
			   multiple_comparison_correction = TRUE,
			   kernel = 'anovadot',
			   nItera = 5,
			   workpath = "latent_Discovery",
			   filename = 'evolution_latent.pdf',
			   truelatent = NULL,
			   truecoef = NULL,
			   textmode = FALSE,
			   res = NULL,
			   seed = 1234,
			   include_downstream = TRUE,
			   include_upstream = TRUE,
			   include_output = TRUE,
			   make_latent_inonly = TRUE,
			   recompute_vars = TRUE,
			   node = NULL, freq_net = 0.5, maxpath_net = 2,
			   useResiduals = TRUE,
			   maxLatentVars = Inf,
			   latent_iterations = 100,
                           Nsamples = 1, 
			   parallel = FALSE,
			   nSimTasks = 1,
			   varmap = NULL,
			   debug = F,
			   discretize_confounders = F,
			   missing_code = NA,
                           testdata = NULL,
                           height = 10,
                           width = 12,
                           showplots = TRUE,
                           ...
			   ){
    if (!exists('blacklist', envir=.pkgGlobalEnv)) {
        assign("blacklist", ens$other_params$blacklist, envir=.pkgGlobalEnv)
        assign('data', data, envir=.pkgGlobalEnv)
    }
    if(debug){
	message("At beginning...")
	browser()
    }
    if(!is.null(output) && !output %in% colnames(data))
        stop("output ", output, " is not in data")
    message("Create directory")
    dir.create(workpath)
    if(is.null(node) & !is.null(output))
	node = output[1]
    if(is.null(truecoef) & !is.null(truelatent)){
	stop("truelatent is given but truecoef is not.")
    }
    if(!is.null(truelatent)){
	if(!"data.frame" %in% is(truelatent))
	    stop("truelatent should be a data.frame")
    }
    if(!is.null(seed))
	set.seed(seed)
    ##if output is not null, use its drivers to infer the latent variable
    ##if it is null, use all variables
    if (!is.null(output)) {
        vars = getDrivers(ens, output, maxpath = maxpath, cutoff = freqCutoff)$Drivers
    } else {
        vars = colnames(data)
    }
    if(include_output){
	vars = unique(c(vars, output)) #unique in case output is already in vars
    }
    ##remove fixed vars from consideration for residual space since they can never be explained away by the graph
    vars = setdiff(vars, getFixedVars(ens)) #only use non-fixed vars

    nnet = ens$Nboot
    newdata = data
    newens = ens
    newvars = vars
    if(!is.null(truelatent) & !is.null(truecoef)){
	allrs = getScores(ens=ens,truecoef=truecoef, lat_estimates=truelatent,output=output)
	if(!is.null(res)){
	    allrs$Repeat = max(res$Repeat) + 1
	    diffcols=setdiff(colnames(res),colnames(allrs))
	    for(cc in diffcols)
	       allrs[[cc]]=NA
	    allrs=allrs[colnames(res)]
	    allrs = rbind(res, allrs)
	}
    }
    for(ii in 1:nItera){
	if(debug){
	    message("Inside loop...")
	    message("Iteration ", ii, " of ", nItera)
	    browser()
	}
	message("\nIteration ", ii, " of ", nItera)
	###########################################
	## start calculation of latent variables ##
	###########################################
	##  residuals
	message("Calculate Residuals")
	if(useResiduals){
	    resList = getGraphResiduals (newens,
                                         newvars,
                                         newdata,
                                         Nsamples = Nsamples)
	    ## deviance
	    message("Calculates deviances")
	    resDev = residualDeviance(newdata,
				      resList,
				      isOrdinal=isOrdinal,
				      missing_code = missing_code)
	}else{
	    resDev = newdata[, newvars, drop = F]
	}
	message("Find Latent Vars")
	if(debug)
	    browser()
	newlatVars = findLatentVars(resDev,
				 scale. = scale.,
				 method = method,
				 kernel = kernel,
				 maxLatentVars = maxLatentVars,
				 discretize_confounders = discretize_confounders,
				 missing_encoding = missing_code,
				 alpha = alpha,
				 nIter = latent_iterations,
				 multiple_comparison_correction = multiple_comparison_correction,
                                 ...)
	## check for missing latent
	if(length(newlatVars) == 1 && newlatVars == -1){
	    break()
	}else{
	    latVars = newlatVars
	}
	df_toadd = latVars$confounders %>% as.data.frame
	Nlat = ncol(df_toadd)
	newdata = cbind(data,
			df_toadd)
	blacklist = ens$other_params$blacklist
	if(make_latent_inonly){
	    ## add incoming edges to blacklist
	    addblacklist = map_df(colnames(df_toadd),
				  function(vv){
				      data.frame(
					  from = colnames(data),
					  to = vv,
					  stringsAsFactors=F
				      )
				  })
	    blacklist = rbind(blacklist,
			      addblacklist)
	}

        ##also ensure the latent variables never drive the true fixed variables
        ##whether they themselves are defined as fixed or not
        trueFixed = inferFixedVars()
        if (!is.null(trueFixed)) {
            addblacklist = map_df(colnames(df_toadd),
                                  function(vv) {
                                      data.frame(
                                          from = vv,
                                          to = trueFixed,
                                          stringsAsFactors=F
                                      )
                                  })
            blacklist = rbind(blacklist,
                              addblacklist)
        }
	#################
	## run bnlearn ##
	#################
	if(debug){
	    message("\n\nAbout to run bnlearn ...")
	    browser()
	}
        oldens = newens
	message("Run bnlearn")
	    newens = getEnsemble2(
		newdata,
		blacklist = blacklist,
		restart = ens$other_params$restart,
		Nboot = ens$Nboot,
		prior = ens$other_params$prior,
		algorithm = ens$algorithm,
		score = ens$other_params$score,
		parallel = parallel,
                output = output
	    )
	## }else{
	##     newens = getEnsemble(
	##         newdata,
	##         blacklist = blacklist,
	##         restart = ens$other_params$restart,
	##         Nboot = ens$Nboot,
	##         prior = ens$other_params$prior,
	##         algorithm = ens$algorithm,
	##         score = ens$other_params$score
	##     )
	## }
	if(recompute_vars){
	    newvars = getDrivers(newens, output,
				 maxpath = maxpath,
				 cutoff = freqCutoff)$Drivers
	}
	latvars = grep("LV.*", colnames(newdata), value = T)
	if(include_downstream){
	    for(vvv in latvars){
		newvars = c(newvars,
			    getDrivers(newens,
				       vvv,
				       cutoff = freqCutoff,
				       maxpath = maxpath,
				       direction = 'downstream')$Drivers
			    )
	    }
	    newvars = unique(newvars)
	    newvars = setdiff(newvars, output)
	}
	if(include_upstream){
	    for(vvv in latvars){
		newvars = c(newvars,
			    getDrivers(newens,
				       vvv,
				       cutoff = freqCutoff,
				       maxpath = maxpath,
				       direction = 'upstream')$Drivers
			    )
		newvars = unique(newvars)
		newvars = setdiff(newvars, output)
	    }
	}
	if(include_output){
	    newvars = unique(c(newvars, output))
	}
	newvars = setdiff(newvars, latvars)
	####################################
	## if in diagnostic mode run plot ##
	####################################
	if(!is.null(truelatent)){
	    message("Entering diagnostic mode...")
	    if(debug){
		browser()
	    }
            if(!is.null(testdata)){
                ## predict latent variables
                message("predicting latent variables in test set")
                latVars_test = predict(latVars,
                                       testdata,
                                       isOrdinal = isOrdinal,
                                       missing_code = missing_code,
                                       useResiduals = useResiduals,
                                       ens = oldens, allPCs=TRUE)
                allpcs = latVars_test$confounders
                sigpc = latVars_test$details$sigPcIdx
                latVars_test$confounders = latVars_test$confounders[, sigpc, drop = F]
            }else{
                latVars_test = latVars
                allpcs = NULL
                }
	    allrs = getScores(allrs,
			      ens=newens,
			      truecoef=truecoef,
			      truelatent = truelatent,
			      lat_estimates = latVars_test,
			      output=output,
			      ignore = "U.*|LV.*",
                              allpcs = allpcs
                              )
	    if(debug){
		message("About to save plot diagnostics")
		browser()
	    }
	    if(textmode){
		diagnosticPlots_text(allrs,
				     ii,
				     truelatent,
				     truecoef)
	    }else{
		message("Saving plot to disk: ", filename)
		cairo_pdf(file.path(workpath, filename), width = 15, height = 10)
		diagnosticPlots(allrs,
				ii,
				truelatent,
				truecoef,
				node = node,
				freq = freq_net,
				maxpath = maxpath_net,
				varmap = varmap,
                                pvalue=FALSE)
		dev.off()
                if(showplots){
		if(ii == 1)
		    dev.new(width = width, height = height)
		if(debug){
		    message("About to plot...")
		    browser()
		}
		flush.console()
		diagnosticPlots(allrs,
				ii,
				truelatent,
				truecoef,
				node = output,
				freq = freq_net,
				maxpath = maxpath_net,
				varmap = varmap,
                                pvalue=FALSE)
		## ggplot2:::print.ggplot(
		##     ggpall
		## )
                }
	    }
	}
    }
    if(debug){
	message("After loops...")
	browser()
    }
    message("Finished bnlearn Runs")
    if(!is.null(truelatent) & !is.null(truecoef)){
        latVars$details$Diagnostics = allrs
        }
    finaldataloc = file.path(
	    workpath,
	    "train_final.csv"
    )
    fwrite(newdata,
	   file = finaldataloc
	   )
    latVars$details$dataloc = finaldataloc
    ## add parameters
    latVars$details$isOrdinal = isOrdinal
    latVars$details$missing_code = missing_code
    latVars$details$useResiduals = useResiduals
    latVars$details$final_ensemble = newens
    latVars$details$latvar_ensemble = oldens
    saveRDS(latVars,
	    file = file.path(
		workpath,
		"latVars.RDS")
	    )
    ## latvafname = file.path(workpath,
    ##                        paste0(tools::file_path_sans_ext(filename), ".RDS"))
    ## saveRDS(latVars, file = latvafname)
    return(latVars)
}


## plotting functions
pcaPlot = function(ii,obj, true, type = 'Confounders', pvalue = FALSE){
    pcs0 = obj[[type]][[ii]]
    if (pvalue) {
        pvals = obj$Pvalues[[ii]]
        names(pvals) = colnames(pcs0)
        pvals = head(pvals, n)
    }
    n = 10
    pcs0 = pcs0[, head(colnames(pcs0), n)]
    pcs = pcs0 %>%
	as.data.frame %>% cbind(true)

    pcs2 = gather(pcs, var, val, -matches(paste0(colnames(true), collapse = "|")))
    pcs3 = gather(pcs2, var2, val2, -var, -val)
    pcs3s = pcs3 %>%
	group_by(var, var2) %>%
	summarise(Cor = cor(val, val2))
    if(length(unique(pcs3$val2)) < 5){
	if(pvalue){
	    pcs3 %>%
		ggplot() +
		geom_boxplot(aes(x = factor(val2), y = val,
				 fill = pvals[var])) +
		geom_smooth(aes(x = val2, y = val), method = 'lm',
			    colour = 'red') +
		geom_text(data = pcs3s, x = -Inf, y = Inf,
			  hjust = -0.2, vjust = 1.2,
			  aes(label = paste0("Cor=",signif(Cor, 2)))) +
		facet_grid(var2 ~ var, scales = 'free') +
		scale_fill_gradient(low="lightcoral", high="skyblue") +
		ylab("PCA") + ylab("True LV")

	}else{
	    pcs3 %>%
		ggplot() +
		geom_boxplot(aes(x = factor(val2), y = val,
				 )) +
		geom_smooth(aes(x = val2, y = val), method = 'lm',
			    colour = 'red') +
		geom_text(data = pcs3s, x = -Inf, y = Inf,
			  hjust = -0.2, vjust = 1.2,
			  aes(label = paste0("R2=",signif(Cor^2, 2)))) +                 facet_grid(var2 ~ var, scales = 'free')
	}
    }else{
	if(pvalue){
	    pcs3 %>%
		ggplot(aes(x = val, y = val2)) +
		geom_point(aes(colour = pvals[var])) +
		geom_text(data = pcs3s, x = -Inf, y = Inf,
			  hjust = -0.2, vjust = 1.2,
			  aes(label = paste0("Cor=",signif(Cor, 2)))) +
		facet_grid(var2 ~ var, scales = 'free') +
		scale_color_gradient(low="red", high="blue")

	}else{
	    pcs3 %>%
		ggplot(aes(x = val, y = val2)) +
		geom_point() +
		geom_text(data = pcs3s, x = -Inf, y = Inf,
			  hjust = -0.2, vjust = 1.2,
			  aes(label = paste0("R2=",signif(Cor^2, 2)))) +
		facet_grid(var2 ~ var, scales = 'free')
	}
    }
}

corplot = function(pcs, true){
    pcs = pcs %>%
	as.data.frame %>% mutate(True = true)
    tmp =lm(formula = True ~ ., data = pcs)
    pcs %>%
	gather(var, val, -True) %>%
	ggplot(aes(x = `True`, y = val)) +
	geom_point() + facet_wrap( ~ var, scales='free') +
	ggtitle(paste0("Adj R2=", signif(summary(tmp)$adj.r.squared, 2), ". Cor=",
		       paste0(signif(cor(dplyr::select(pcs, -True), pcs$True), 2), collapse = ", ")))
}

library(grid)
define_region <- function(row, col){
    viewport(layout.pos.row = row, layout.pos.col = col)
}
library(grid)
library(gridBase)
library(ggplot2)
library(tidyr)
diagnosticPlots = function(obj, ii, true, truecoef, node = NULL, freq = 0.5, maxpath = 2, varmap = NULL, pvalue=FALSE){
    ## pval sv
    ## pvals = obj$Pvalues[[ii]]
    ## ggppva = qplot(paste0("PC",1:length(pvals)),pvals)
    ## PCA
    ##ggpc = pcaPlot(ii, obj, true, "Confounders")
    if(ii > 0){
	ggp1 = pcaPlot(ii, filter(obj, Iteration > 0),
		       true, "PCAs", pvalue = pvalue)
    }else
	ggp1 = NULL
    ## R2
    if(ii > 0){
	allrsp = obj %>%
	    filter(Iteration > 0) %>%
	    dplyr::select(Repeat, Iteration, starts_with("R2a")) %>%
	    gather(var, val, -Repeat, -Iteration)
	ggp2 = allrsp %>%
	    ggplot(aes(x = Iteration, y = val)) +
	    geom_point() + geom_line() + geom_smooth() +
	    geom_point(data = filter(allrsp, Iteration == ii),
		       aes(x = Iteration, y = val),
		       colour = 'red', size = 3) +
	    facet_wrap( ~ var, nrow = 1)
    }else
	ggp2 = NULL
    ## coefficients
    ## cc = obj$Coef[[ii]]
    ## cc = gather(cc, input, coef)
    ## ggpc = ggplot() +
    ##     geom_boxplot(data = cc, aes(x = input, y = coef)) +
  ##     ggbeeswarm::geom_quasirandom(data = cc, aes(x = input, y = coef)) +
    ##     geom_point(data = truecoef, aes(x = input, y = coef), colour = 'red', size = 3) + ylim(c(-3, 3))
    ## error in coef
    if(ii > 0){
	errall = obj %>%
	    dplyr::select(Repeat, Iteration, Coef) %>%
	    unnest(Coef)
	err = errall %>%
	     filter(True_coef >= 0.02)
	##            input != colnames(true)) %>%
	##     dplyr::select(-coef, -Estimate) %>%
	##     mutate(input=otfname2name(varmap,input))
	ggpe = err %>%
	    ggplot(aes(x = Iteration, y = abs(error), colour = input)) +
	    geom_point() +
	    geom_line() +
	    geom_vline(xintercept = ii, colour = 'red') +
	    ##ylim(0, max(err$error)) +
	    ggtitle("Errors (excluding FP)")##+ geom_smooth()
	ggpet = obj %>%
	    ggplot(aes(x = Iteration, y = Error_rmse)) +
	    geom_point() +
	    geom_line() +
	    geom_point(data = filter(obj, Iteration == ii),
		       aes(x = Iteration, y = Error_rmse),
		       colour = 'red', size = 4) +
	    ggtitle("Total Error (including FP)")
    }
    ## combine everything
    if(is.null(node) & ii > 0){
	cowplot::plot_grid(ggp2, ggpet, ggpe, ncol = 1,
			   rel_heights = c(0.2, 0.2, 0.6))
    }else{
	ens = obj$Ensemble[[ii + 1]]
	plot.new()
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = 4,
						   ncol = 2)))
	pushViewport(define_region(row = 1, col = 2))
	par(fig = gridFIG(), new = TRUE)
	if(!is.null(varmap)){
	    node_mn = name2otfname(varmap, node)
	    node_dn = otfname2name(varmap, node)
	}else{
	    node_mn = node
	    node_dn = node
	}
	drivers = getDrivers(ens, node_mn, maxpath = 1, cutoff = freq)
	if(!is.null(varmap))
	    drivers_dn = otfname2name(varmap,drivers)
	else
	    drivers_dn = drivers
	edge_color = tibble(inp = drivers,
			    out = node_mn,
			    col = "red")
	plot(ens,
	     node_mn,
	     freqth = freq,
	     cutoff=0.4,
	     maxpath = maxpath,
	     nodesep = 0.001,
	     edge_color = edge_color,
	     edge_labels="coefficients"
	     )
	popViewport()
	if(ii > 0){
	    ggplot2:::print.ggplot(ggp2, vp = define_region(row = 1, col = 1), newpage = F)
	    ggplot2:::print.ggplot(ggpe, vp = define_region(row = 2, col = 1), newpage = F)
	    ggplot2:::print.ggplot(ggpet, vp = define_region(row = 2, col = 2), newpage = F)
	    ggplot2:::print.ggplot(ggp1, vp = define_region(row = 3:4, col = 1:2), newpage = F)
	}
    }
}

diagnosticPlots_text = function(obj, ii, true, truecoef){
    library(tidyr)
    df = obj %>% dplyr::select(
			    Repeat,
			    Iteration,
			    starts_with("R2a")
			)
    print(knitr::kable(filter(df, Iteration == ii), format = 'pandoc'))
    message("Latent Variable Estimation R2 using selected confounders")
    for(uu in colnames(true)){
	message("\n", uu)
	txtplot::txtplot(df$Iteration, df[[paste0("R2a_", uu)]])
    }
    message("Latent Variable Estimation R2 using all pcas")
    for(uu in colnames(true)) {
	message("\n", uu)
	txtplot::txtplot(df$Iteration, df[[paste0("R2a_allpca_", uu)]])
    }
    ## coefficients
    dfc = obj %>% dplyr::select(
			     Repeat,
			     Iteration,
			     Coef
			 ) %>%
	unnest(Coef) %>%
	unique %>%
	group_by(Repeat, Iteration)
    message("Error in True Drivers")
    dfc2 = dfc %>%
	filter(abs(True_coef) > 0) %>%
	dplyr::select(Repeat, Iteration, input, error) %>%
	spread(input, error)
    print(knitr::kable(filter(dfc2, Iteration %in% c(0, ii)),
		       format = 'pandoc'))
    message("RMSE in all coefficients")
    txtplot::txtplot(obj$Iteration, obj$Error_rmse)
}



#########################################
## calculate coeficients from ensemble ##
#########################################
## using fsTermsFrequencies
ens2coef = function(ens, variable, summary = F, sumfun = mean, freqCutoff = 0.05){
    out = variable
    ## fixNames = function(vv){
    ##     sub("\\((.*) = .*)", "\\1", vv)
    ## }
    ## input = fixNames(input)
    coefmean = fsTermFrequencies(ens,
				 incParamStats=T) %>%
	filter(output==out,
	       freq > freqCutoff) %>%
	mutate(coef = coef.mean * freq) %>%
	dplyr::select(input, coef)
    if(summary)
	coefmean = df2vector(coefmean)
    return(coefmean)
}
df2vector = function(df, names = 1, values = 2){
    if(names == 'rownames')
	mynames = rownames(df)
    else
	mynames = df[[names]]
    myvalues = df[[values]]
    mylist = myvalues
    names(mylist) = mynames
    return(mylist)
}

ens2coefold = function(ens, variable, summary = FALSE, sumfun = mean){
    library(stringr)
    ## get fragments for the variable
    allfrags = fsGetFrags(ens, variable)
    ## get list of variables present in the
    allvars = allfrags$input %>% unlist %>% unique
    allvars = c("(Intercept)", allvars)
    allforms = allfrags$formula
    allforms2 = str_split(allforms, "~")
    Nnet = nrow(allfrags)
    Nvar = length(allvars)
    matcoef = matrix(rep(0, Nnet * Nvar), ncol = Nvar)
    colnames(matcoef) = allvars
    for(jj in 1:Nnet){
	ff = allforms2[[jj]]
	tmp = str_split(ff[2], "\\+")[[1]]
	tmp = str_trim(tmp)
	tmp = str_split(tmp, " ")
	tmp = sapply(tmp, function(fff) {
	    if(length(fff) == 1){
		val = as.numeric(fff)
		names(val) = '(Intercept)'
	    }else{
		val = as.numeric(fff[1])
		names(val) = fff[2]
	    }
	    return(val)
	})
	newNames = setdiff(names(tmp), colnames(matcoef))
	if (length(newNames) > 0) {
	    matcoef = cbind(matcoef, matrix(rep(0, Nnet * length(newNames)), ncol=length(newNames)))
	    colnames(matcoef)[(ncol(matcoef)-length(newNames)):ncol(matcoef)] = newNames
	}
	matcoef[jj, names(tmp)] = tmp
    }
    if(summary)
	return(
	    apply(matcoef, 2, sumfun)
	)
    return(as.data.frame(matcoef))
}


##  allrsp = allrs %>% dplyr::select(Repeat, Epoch, R2a, R2a_allpca) %>% gather(var, val, -Repeat, -Epoch)
## ggp = ggplot(data = allrsp,
##                       aes(x = Epoch, y = val)) +
##      geom_point() +
##      geom_smooth() +
##      geom_line(aes(group = Repeat,
##                    colour = factor(Repeat))) +
##     facet_wrap( ~ var, nrow = 1)
##  tmp =lm(formula = True ~ ., data = pcs)
##  ggppva = qplot(paste0("PC",1:length(latVars$details$p.val)),latVars$details$p.val)
##  ggppca = pcs %>%
##      gather(var, val, -True) %>%
##      ggplot(aes(x = `True`, y = val)) +
##      geom_point() + facet_wrap( ~ var, scale = "free") +
##      ggtitle(paste0("Adj R2=", signif(summary(tmp)$adj.r.squared, 2), ". Cor=",
##                     paste0(signif(cor(dplyr::select(pcs, -True), pcs$True), 2), collapse = ", "))
##              )


## other stuff

##' Calculate PCA-like latent variables estimates using mutual information estimates instead of covariance matrix
##'
##' Calculates the mutual information matrix between all variables, calculate the eigendecomposition and projects the data into this new coordinate given by the eigenvectors.
##' @title miPCAm
##' @param train data.frame
##' @param mi_type type of mutual information estimate
##' @param scale_data whether to scale the daftda
##' @param center_data whether to center the data
##' @param disctype type of discretization to use
##' @param nbins number of bins to use in B-Spline method
##' @param splineorder order of spline in B-Spline method
##' @return data.framae with latent varibles. Also an attribute with eigenvalues
miPCA = function(train,
		 mi_type = c(
		     "B-Spline MI",
		     "MIC",
		     "MAS",
		     "MEV",
		     "MCN",
		     "MICR2",
		     "GMIC",
		     "TIC",
		     "Empirical",
		     "Miller-Madow",
		     "Shrinkage-Dirichlet",
		     "Schurmann-Grassberger",
		     "Pearson",
		     'Spearman',
		     "Covariance"
		 ),
		 scale_data = TRUE,
		 center_data = TRUE,
		 disctype = c(
		     "equalfreq",
		     "equalwidth",
		     "globalequalwidth"
		 ),
		 nbins = 7,
		 splineorder = 3
		 ){
    disctype = disctype[1]
    if(! mi_type %in% c(
			  "B-Spline MI",
			  "MIC",
			  "MAS",
			  "MEV",
			  "MCN",
			  "MICR2",
			  "GMIC",
			  "TIC",
			  "Empirical",
			  "Miller-Madow",
			  "Shrinkage-Dirichlet",
			  "Schurmann-Grassberger",
			  "Pearson",
			  'Spearman',
			  "Covariance"
		      )){
	stop("mi_type not allowed")
    }else
	mi_type = mi_type[1]
    if(scale_data){
	train = scale(train, center = center_data)
    }
    trainm = as.matrix(train)
    if(mi_type == "Empirical"){
	mi = minet::build.mim(trainm, estimator = "mi.empirical", disc = disctype)
    }else if(mi_type == "Covariance"){
	mi = cov(trainm)
    }else if(mi_type == "Miller-Madow"){
	mi = minet::build.mim(trainm, estimator = "mi.mm", disc = disctype)
    }else if(mi_type == "Shrinkage-Dirichlet"){
	mi = minet::build.mim(trainm, estimator = "mi.shrink", disc = disctype)
    }else if(mi_type == "Schurmann-Grassberger"){
	mi = minet::build.mim(trainm, estimator = "mi.sg", disc = disctype)
    }else if(mi_type == "Pearson"){
	mi = minet::build.mim(trainm, estimator = "pearson")
    }else if(mi_type == "Spearman"){
	mi = minet::build.mim(trainm, estimator = "spearman")
    }else if(mi_type == "B-Spline MI"){
	mi = BSplineMI::calcSplineMI(t(trainm), nbins, splineorder)
    }else{
	mi = minerva::mine(trainm)[[mi_type]]
    }
    ## get spectral decomposition
    mieig = eigen(mi)
    lat = trainm %*% mieig$vectors %>% as.data.frame
    colnames(lat) = paste0 ("LV_", 1:ncol(lat))
    attributes(lat)$eigenvalues = mieig$values
    return(lat)
}


refitLatent = function(latobj,
		       testdf,
		       output,
		       include_output = FALSE,
		       scale.,
		       nIter,
		       method,
		       kernel,
		       alpha,
		       pval_method,
		       maxLatentVars){
    latdet = latvars_simple$details
    if(missing(scale.))
	scale. = latdet$scale.
    if(missing(nIter))
	nIter = latdet$nIter
    if(is.null(nIter))
	nIter = 100
    if(missing(method))
	method = latdet$method
    if(is.null(method))
	method = 'linear'
    if(missing(kernel))
	kernel = latdet$kernel
    if(is.null(kernel))
	kernel = 'anovadot'
    if(missing(alpha))
	alpha = latdet$alpha
    if(is.null(alpha))
	alpha = 0.05
    if(missing(pval_method))
	pval_method = latdet$pval_method
    if(is.null(pval_method))
	pval_method = 'beta'
    if(missing(maxLatentVars))
	maxLatentVars = latdet$maxLatentVars
    if(is.null(maxLatentVars))
	maxLatentVars = Inf
    resid = latdet$originalData
    if(!include_output){
	resid[[output]] = NULL
	message("find latent variables model without outcome")
	latVars = findLatentVars(
	    resid,
	    scale. = scale.,
	    nIter = nIter,
	    method = method,
	    kernel = kernel,
	    alpha = alpha,
	    pval_method = pval_method,
	    maxLatentVars = maxLatentVars
	)
    }else{
	latVars = latobj
    }
    allvars = latVars$details$originalData %>% colnames
    testdfsel = testdf %>% dplyr::select_(.dots = allvars)
    testdfsel = testdfsel[, allvars, drop = F]
    latVars_test = predict(latVars, testdfsel, )
    return(latVars_test)
}

getErrorTable = function(truecoef, estcoef, allvariables, truelatent_names){
    errcoef = truecoef %>%
	unique %>%
	mutate(Estimate = ifelse(is.na(estcoef[input]),
				 0,
				 estcoef[input]),
	       Error = abs(coef - Estimate))
    othervars = setdiff(names(estcoef),
			c("(Intercept)",
			  errcoef$input))
    othervars = othervars[which(!grepl("LV_.*", othervars))]
    errcoef = rbind(
	errcoef,
	tibble(input = othervars,
	       coef = 0,
	       Estimate = estcoef[othervars]
	       ) %>% mutate(Error = abs(coef - Estimate))
    )
    if(!missing(truelatent_names)){
	errcoef = errcoef %>% filter(!input %in%
				     truelatent_names)
    }
    return(errcoef)
}

insertRow <- function(existingDF, newrow, r) {
    existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    existingDF
}



getScores = function(res = NULL, ens, truecoef, truelatent = NULL, lat_estimates, output, ignore = "U.*|LV.*", allpcs = NULL){
    currentRepeat=max(res$Repeat)
    coefmean0= getCoef(ens, output)
    latVars = lat_estimates
    missingvars = setdiff(coefmean0$input, c(truecoef$input, output))
    truecoef = rbind(
	truecoef %>%
	mutate(freq =1,
	       wcoef = coef
	       ),
	tibble(
	    input = missingvars,
	    output = output,
	    freq = 0,
	    coef = 0,
	    wcoef = 0
	)
    )
    truecoef = filter(truecoef, !grepl(ignore, input))
    coefmean0 = filter(coefmean0, !grepl(ignore, input))
    alldf = inner_join(
	dplyr::select(truecoef, input, output, True_coef = wcoef),
	dplyr::select(coefmean0, input, output, Est_coef = wcoef),

        by = c("input", "output")
    ) %>% mutate(error = True_coef - Est_coef)
    rmseval = alldf %>% summarise(rmse=sqrt(mean(error^2))) %>% with(rmse)
    crs = tibble(
	Repeat = currentRepeat,
	Iteration = 0,
	PCAs = list(NULL),
	Confounders = list(NULL),
	Residuals = list(NULL),
	Coef = list(alldf),
	Error_rmse = rmseval,
	Ensemble = list(ens),
	Pvalues = list(latVars$details$p.val)
    )
    if(!is.null(truelatent)){
	## get R2 for true latent variable
	method =latVars$details$method
	tst = list()
	latVars = lat_estimates
	dflat0 = latVars$confounders %>%
	    as.data.frame
	for(uu in colnames(truelatent)){
	    dflat = dflat0 %>%
		cbind(U.lat = truelatent[[uu]])
	    tst[[uu]] = lm(formula = U.lat ~ ., data = dflat)
	    form = paste('~',
			 paste(colnames(latVars$confounders),
			       sep = '', collapse = '+'))
	    tst[[uu]] = stepAIC(tst[[uu]],
				k = log(nrow(truelatent)),
				scope = list(lower = "~1",
					     upper = form),
				trace = 0,
				direction = 'both')
	}
        if(is.null(allpcs)){
            if(method == 'linear'){
            pcs = latVars$details$pcObj$x %>%
                    as.data.frame
	} else if(method == 'robustLinear'){
	    pcs = svd(latVars$details$pcObj$L)$u %>%
					     as.data.frame
	    colnames(pcs) = paste0("PC",1:ncol(pcs))
	} else if (method == 'ica'){
	    pcs = latVars$details$pcObj$S %>% as.data.frame
        } else if (method == 'autoencoder') {
            pcs = latVars$confounders %>% as.data.frame
	} else {
	    pcs = latVars$details$pcObj@pcv %>%
		as.data.frame
	}
        }else{
            pcs = allpcs %>%
                as.data.frame
        }
	crs$PCAs = list(pcs)
	tst2 = list()
	for(uu in colnames(truelatent)){
	    dflat = pcs %>%
		cbind(U.lat = truelatent[[uu]])
	    tst2[[uu]] = lm(data = dflat, formula = U.lat ~ 1)
	    form = paste('~',
			 paste(colnames(pcs),
			       sep = '', collapse = '+')
			 )
	    tst2[[uu]] = stepAIC(tst2[[uu]],
				 k = log(nrow(truelatent)),
				 scope = list(lower = "~1",
					      upper = form
					      ),
				 trace = 0,
				 direction = 'both')
	}            ## if on the first interation also add initial case
	for(uu in colnames(truelatent)){
	    dftmp = tibble(
		R2 = summary(tst[[uu]])$r.squared,
		R2a = summary(tst[[uu]])$adj.r.squared,
		R2_allpca = summary(tst2[[uu]])$r.squared,
		R2a_allpca = summary(tst2[[uu]])$adj.r.squared
	    )
	    colnames(dftmp) = paste0(
		colnames(dftmp),
		"_",
		uu
	    )
	    crs = bind_cols(
		crs,
		dftmp
	    )
	}
    }
    if(!is.null(res)){
	crs$Iteration = max(res$Iteration) + 1
	allrs = bind_rows(res,
		      crs
		      )
    }else
	allrs = crs
    return(allrs)
}

getAvgPredictions = function(ens, variables){
    results = map_df(variables,
		     function(vv){
			 fit = ens$fitmodels
			 map_df(1:length(fit),
				function(nn){
				    resi = fit[[nn]][[vv]]$residuals
				    tibble(
					patid = 1:length(resi),
					Variable = vv,
					Network = nn,
					Residual = resi)
				})
		     })
    results %>% group_by(patid, Variable) %>%
	summarise(Residual = mean(Residual)) %>%
	ungroup() %>%
	spread(Variable, Residual)
}

getFixedVars = function(ens, vars = ".*"){
    if (length(vars) == 1) {
	allvars = grep(vars, colnames(ens$data), value = T)
    } else {
	allvars = vars
    }

    all_inputonly = purrr::map(1:ens$Nboot,
			function(bb){
			    logv = purrr::map_lgl(allvars,
                                           function(vv){
                                               if(length(ens$fitmodels[[bb]][[vv]]$parents) > 0)
                                                   FALSE
                                               else
                                                   TRUE
                                           })
			    allvars[logv]
			})
    counts = table(unlist(all_inputonly))
    ##this only works if there are enough counts, which is a combination of variables and bootstraps
    if (sum(counts) > 60) { #rule of thumb for a gaussian is 30, so 2x30....
        library(mclust)
        res = Mclust(counts)
        ##use the cluster with the highest counts, discard the others
        keepIdx = which(res$classification == max(res$classification))
        return(names(keepIdx))
    }
    
    return(all_inputonly)
        
}

getDrivers = function(ens,output, maxpath =4, cutoff = 0.5, direction = 'upstream'){
    nnet = length(ens$allnet)
    if(direction == 'upstream')
	mode = 'in'
    else if(direction == "downstream")
	mode = 'out'
    else if(direction == "both")
	mode = "all"
    else
	stop("direction not recognized. Only: 'upstream','downstream', and 'both'.")
    allres = purrr::map_df(1:nnet, function(ii){
	ig = igraph::igraph.from.graphNEL(bnlearn::as.graphNEL(ens$allnet[[ii]]))
	parents = igraph::neighborhood(ig,order = maxpath, output, mode = mode)
	parents = names(parents[[1]])
	parents = setdiff(parents, output)
	tibble(network = ii,
	       Drivers = parents)
    })
    allres %>%
	dplyr::group_by(Drivers) %>%
	dplyr::summarise(freq = n() / nnet) %>%
	dplyr::arrange(-freq) %>%
	dplyr::filter(freq >= cutoff)
}

getCoef = function(ens, outcome, as.regex = FALSE){
    if(as.regex){
        alloutcomes = grep(outcome, colnames(ens$data), value = T)
        allcoef = map_df(
            alloutcomes,
            function(oo){
                getCoef(ens, oo, as.regex = FALSE)
            }
        )
    }else{
        nnet = ens$Nboot
        allcoef = map_df(1:nnet,
                         function(ii){
                             tmp = ens$fitmodels[[ii]][[outcome]]$coefficients
                             tmp %>% t() %>%  as.data.frame %>% mutate(Boot = ii)
                         })
        allcoef = gather(allcoef, input, coef, -Boot)
        allcoef = allcoef %>% group_by(input) %>% summarise(
                                                      freq = sum(!is.na(coef) ) / nnet,
                                                      coef = mean(coef[!is.na(coef)])
                                                  ) %>%
            ungroup() %>%
            mutate(output = outcome,
                   wcoef = coef * freq) %>%
            arrange(-freq, abs(coef)) %>%
            dplyr::select(input, output, freq, coef, wcoef)
        allvars = colnames(ens$data)
        missingvars = setdiff(allvars, c(allcoef$input, outcome))
        allcoef = rbind(allcoef,
                        tibble(
                            input = missingvars,
                            output = "output",
                            freq = 0,
                            coef = 0,
                            wcoef = 0
                        ))
    }
    return(allcoef)
}

## tied weights autoencoder:
## https://github.com/dfalbel/deep-autoencoder-netflix/blob/master/tied-dense-layer.R#L4
## TiedDenseLayer <- R6::R6Class(
##   "TiedDenseLayer",
##   inherit = KerasLayer,
##   public = list(
##     master_layer = NULL,
##     kernel = NULL,
##     bias = NULL,
##     output_dim = NULL,
##     initialize = function(output_dim, master_layer) {
##       self$master_layer <- master_layer
##     },
##     build = function(input_shape) {
##         self$kernel <- k_transpose(self$master_layer$kernel)
##         self$output_dim <- self$kernel$shape$as_list()[[2]]
##         self$bias <- self$add_weight(
##                               name = 'bias',
##                               shape = list(self$output_dim),
##                               initializer = initializer_constant(0),
##                               trainable = TRUE
##                           )
##     },
##     call = function(x, mask = NULL) {
##         k_dot(x, self$kernel) + self$bias
##     },
##     compute_output_shape = function(input_shape) {
##         list(input_shape[[1]], self$output_dim)
##     }
##   )
##   )

## layer_tied_dense <- function(object, master_layer, name = NULL, trainable = TRUE) {
##     create_layer(TiedDenseLayer, object, list(
##                                              master_layer = master_layer,
##                                              name = name,
##                                              trainable = trainable
##                                          ))
## }


TiedDenseLayer <- R6::R6Class(
  "TiedDenseLayer",
  inherit = KerasLayer,
  public = list(
    master_layer = NULL,
    kernel = NULL,
    bias = NULL,
    output_dim = NULL,
    initialize = function(master_layer) {
        self$master_layer <- master_layer
    },
    build = function(input_shape) {
        browser()
      self$kernel <- k_transpose(self$master_layer$kernel)
      self$output_dim <- self$kernel$shape$as_list()[[2]]
      ## self$bias <- self$add_weight(
      ##   name = 'bias',
      ##   shape = list(self$output_dim),
      ##   initializer = initializer_constant(0),
      ##   trainable = TRUE
      ##   )
     },
    call = function(x, mask = NULL) {
      k_dot(x, self$kernel) ##+ self$bias
    },
    compute_output_shape = function(input_shape) {
      list(input_shape[[1]], self$output_dim)
    }
  )
)

layer_tied_dense <- function(object, master_layer, name = NULL, trainable = TRUE) {i
  create_layer(TiedDenseLayer, object, list(
    master_layer = master_layer,
    name = name,
    trainable = trainable
  ))
}


## TiedDenseLayer <- R6::R6Class(
##   "TiedDenseLayer",
##   inherit = KerasLayer,
##   public = list(

##     master_layer = NULL,
##     W = NULL,
##     b = NULL,
##     output_dim = NULL,

##     initialize = function(output_dim, master_layer) {
##       self$master_layer <- master_layer
##     },

##     build = function(input_shape) {
##         message("build function")
##         cw = self$master_layer$weights
##         if(length(cw) == 0){
##             self$W = list()
## p            self$output_dim = 0
##         }else if(length(cw > 1)){
##             self$W <- k_transpose(self$master_layer$weights[[1]])
##             self$output_dim <- self$W$shape$as_list()[[2]]
##         }else{
##             self$W = k_transpose(self$master_layer$weights)
##             self$output_dim <- self$W$shape$as_list()[[2]]
##         }

##       self$b <- self$add_weight(
##         name = 'bias',
##         shape = list(self$output_dim),
##         initializer = initializer_constant(0),
##         trainable = TRUE
##       )

##     },

##     call = function(x, mask = NULL) {
##         message("call function")
##       k_dot(x, self$W) + self$b
##     },

##     compute_output_shape = function(input_shape) {
##       list(input_shape[[1]], self$output_dim)
##     }

##   )
## )

## layer_tied_dense <- function(object, master_layer, name = NULL, trainable = TRUE) {
##   create_layer(TiedDenseLayer, object, list(
##     master_layer = master_layer,
##     name = name,
##     trainable = trainable
##   ))
## }


plot_history = function(history, type = "loss"){
    res = history$history
    if(type == "loss"){
        lossdf = tibble(val_loss = as.numeric(res$val_loss),
                loss = as.numeric(res$loss),
                Epoch = 1:length(loss)
                )
    gather(lossdf, key, value, -Epoch) %>% ggplot(aes(x = Epoch, y = value, colour = key)) + geom_point() + geom_line()
    }else{
        val = paste0("val_", type)
        msedf = tibble(
            val_error = as.numeric(res[[val]]),
            error = as.numeric(res[[type]]),
            Epoch = 1:length(error))
        colnames(msedf)[1:2] = c(val, type)
        gather(msedf, key, value, -Epoch) %>% ggplot(aes(x = Epoch, y = value, colour = key)) + geom_point() + geom_line()
    }

}


predict.robustLinear <- function(pr, newdata) {
    as.matrix(newdata %*% t(pr$L.svd$vt))
}

##Note that, for out of sample data, the details section never gets modified,
##keeping the origins of the derivation of the confounder within reach
predict.LatentConfounder <- function(latConf,
                                     newdata = NULL,
                                     ens = NULL,
                                     allPCs = FALSE,
                                     isOrdinal = NULL,
                                     missing_code = NULL,
                                     useResiduals = NULL
                                     ) {
    if (is.null(newdata)) { #use original residuals from training
        newdata = latConf$details$originalData
    }else{
        message("Calculating residuals on new dataset")
        if(is.null(useResiduals))
            useResiduals = latConf$details$useResiduals
        if(useResiduals){
            if(is.null(isOrdinal))
                isOrdinal = latConf$details$isOrdinal
            if(is.null(missing_code))
                missing_code = latConf$details$missing_code
            ## calculate the residuals
            allvars = latConf$details$originalData %>% colnames
            if(is.null(ens))
                ens = latConf$details$latvar_ensemble
            message("get graph residuals")
            resList = getGraphResiduals(
                ens,
                allvars,
                data = newdata
             )
            message("Calculating PSR")
            resDev = residualDeviance(newdata,
                                      resList,
                                      isOrdinal=isOrdinal,
                                      missing_code = missing_code
                                      )
            newdata = resDev
        }
    }
    ##res = predict(latConf$details$pcObj, newdata)
    message("predicting latent variables")
    if(latConf$details$method == 'kernel'){
        varsel = colnames(kernlab::xmatrix(latConf$details$pcObj))
        res = kernlab::predict(latConf$details$pcObj,
                               x = as.matrix(newdata[, varsel, drop = F])
                               )
    }else
        res = predict(latConf$details$pcObj, as.matrix(newdata))
    if (!allPCs) {#then significant PCs only
        res = res[, latConf$details$sigPcIdx, drop = F]
    }
    ##colnames(res) = gsub('PC', 'LV_', colnames(res))
    colnames(res) = paste('LV', 1:ncol(res), sep = '_')
    latConf$confounders = res
    return(latConf)
}
