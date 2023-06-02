#' @title SFNetwork SFNetwork class definition
#' @import R6
#' @importFrom pcalg randDAG rmvDAG
#' @importFrom igraph topo_sort ego_size V
#' @importFrom bnlearn as.bn as.igraph parents drop.arc as.graphNEL graphviz.plot
#' @importFrom brms rskew_normal

##An R6 class for generating a scale-free network and some associated data
#' @export
SFNetwork <- R6Class("SFNetwork", 
                     public = list(
                       #' @field dag The generated directed acyclic graph
                       dag = NULL,
                       #' @field vRank vertex topological rank
                       vRank = NULL,
                       #' @field outDegree out-degree of vertices based on direct descendants
                       outDegree = NULL,
                       #' @description
                       #' Generate a DAG using Barabasi's game and generate matching data
                       #' @param numVertices number of vertices in the DAG (default=20)
                       #' @param topology power/star
                       #' @param numChildren number of children of each vertex, on average (i.e., node denisty, default=2)
                       #' @param nParentOnly how many nodes to make parent-only hubs by removing their high order but low out-degree parents (default=4)
                       #' @param addOneEdge whether to add a single edge not involving a latent variable to star or twin-star topologies (helps with AUC calculation; default=TRUE)
                       #' @param par1 Par1, as in documentation for pcalg::randDAG.  For Power Law nets, has the meaning of power of preferential attachment.  Defaults to 1
                       #' @param supersize.top.hubs If true and nParentOnly > 0, then super-size the top n hubs by supersize.factor (up to the total number of nodes in the network)
                       #' @param supersize.factor Multiply the number of automatically generated nodes by this factor (up to the network size minus the number of top hubs)
                       #' @return A `SFNetwork` object
                       initialize = function(numVertices = 20, topology='power', numChildren = 2, nParentOnly = 4, addOneEdge = TRUE, par1 = 1, 
                                             supersize.top.hubs = FALSE, supersize.factor = 2) {
                         if (topology == 'power') {
                           res = private$makePowerDag(numVertices = numVertices, numChildren = numChildren, nParentOnly = nParentOnly, 
                                                      par1 = par1, supersize.top.hubs = supersize.top.hubs, supersize.factor = supersize.factor)
                         } else if (topology == 'twin.star') {
                           res = private$makeTwinStarDag(numVertices = numVertices, addOneEdge = addOneEdge)                                                      
                         } else {
                           res = private$makeStarDag(numVertices = numVertices, addOneEdge = addOneEdge)                           
                         }
                         return(res)
                       },
                       #' @description 
                       #' prints the details of the DAG
                       print = function() {
                         print(self$dag)
                       },
                       #' @description
                       #' plots the DAG
                       plot = function() {
                         opar = par()
                         par(mfrow = c(2, 1), oma = c(3, 3, 3, 3))
                         ## plot the DAG
                         graphviz.plot(self$dag, main = paste("DAG with ", length(self$outDegree), " vertices", sep=''), layout='dot')
                         plot(length(self$outDegree) - rank(self$outDegree), self$outDegree, xlab = 'rank', ylab = 'direct out-degree')
                         par(opar)
                       },
                       #' @param doPlot whether or not to generate a network plot (default = FALSE)  
                       #' @param numSamples number of samples in the generated data frame (default = 1000)
                       #' @param errDist error distribution of data generated for the network, as in rmvDAG
                       #' @param latIdx indices of hubs (in decreasing order) to make latent; by default 2nd and 3rd largest.  These will be prefixed with "U_" (unobserved).  If none, set to NULL
                       #' @param rescaleNoiseBy If a positive number, scales sd from default parameter for all variables by a constant; 
                       #' default behavior is scale everything through a left-skewed gaussian to generate several proxy variables for latent space; 
                       #' if NULL, no rescaling
                       generateData = function(numSamples = 1000,
                                               doPlot = FALSE,
                                               errDist = 'normal',
                                               latIdx = c(2, 3),
                                               rescaleNoiseBy=function(mu=0.5, sigma=0.5, alpha=-5) {
                                                 return(rskew_normal(1, mu, sigma, alpha))
                                               }) {
                         ##recompute the rDAG
                         rDAG = as.graphNEL(self$dag)
                         dataMat <- rmvDAG(numSamples, rDAG, errDist = errDist)
                         
                         if (!is.null(latIdx)) {
                           latent = self$vRank[latIdx]
                           colIdx = match(latent, colnames(dataMat))
                           colnames(dataMat)[colIdx] = paste('U_', colnames(dataMat)[colIdx], sep = '')
                         }
                         graph::nodes(rDAG) = colnames(dataMat)                                                  

                         if (!is.null(rescaleNoiseBy)) { #then add some noise per variable
                           scales = c()
                           for (ci in 1:ncol(dataMat)) {
                             scaleFactor = ifelse(ci %in% latIdx, 
                                                  0.1, #add little noise to latent space
                                                  max(0.1, rescaleNoiseBy()) #keep scale factor positive
                             )
                             scales[ci] = scaleFactor
                             dataMat[, ci] = dataMat[, ci] + rnorm(nrow(dataMat), 0, scaleFactor * sd(dataMat[, ci]))
                           }
                         }
                         
                         return(list(data = as.data.frame(dataMat), graph = as.bn(rDAG)))
                       }
                     ), 
                     private = list(
                       ##used for debugging                       
                       makeStarDag = function(numVertices = 20, addOneEdge = TRUE) {
                         ##start with an empty DAG
                         dag = empty.graph(paste('v', as.character(1:numVertices), sep=''))
                         
                         ##make the last vertex the center
                         for (vi in 1:(numVertices-1)) {
                           from = paste('v', numVertices, sep='')
                           to = paste('v', vi, sep='')                           
                           dag = set.arc(dag, from, to, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
                         }
                         if (addOneEdge) { #add a "hack" - an edge from vertex 1 to vertex 2 - so that the assess() function in NetRes doesn't crash
                           dag = set.arc(dag, from='v1', 
                                         to = 'v2' , check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
                         }
                         
                         ##and topologically sort it
                         nOrder <- topo_sort(as.igraph(dag))
                         
                         aM = as(as.graphNEL(dag), 'matrix')
                         aM = aM[nOrder, nOrder]
                         rDAG = as(aM, 'graphNEL')
                         dag = as.bn(rDAG)
                         
                         ##compute #s of direct children
                         iGraph = as.igraph(dag)
                         
                         outDegree = ego_size(
                           iGraph,
                           order = 1,
                           nodes = V(iGraph),
                           mode = c("out"),
                           mindist = 1
                         )
                         sorted = order(outDegree, decreasing = T)
                         sorted = V(iGraph)[sorted]
                         
                         self$dag = dag
                         self$vRank = names(sorted)
                         self$outDegree = outDegree
                       },
                       makeTwinStarDag = function(numVertices = 20, addOneEdge = TRUE) {
                         ##start with an empty DAG
                         dag = empty.graph(paste('v', as.character(1:numVertices), sep=''))
                         
                         ##make the last vertex the center
                         for (vi in 1:(numVertices-2)) {
                           ##pick parents from either/both stars with equal probability
                           toss = runif(1)
                           to = paste('v', vi, sep='')      
                           if (toss < 1/3) {
                             from = paste('v', numVertices - 1, sep='')
                           } else { #both for 1/3 < toss < 2/3 and toss > 2/3
                             from = paste('v', numVertices, sep='')                             
                           } 
                           dag = set.arc(dag, from, to, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)                           
                           if (toss > 2/3) { #set an additional edge
                             from = paste('v', numVertices - 1, sep='')
                             dag = set.arc(dag, from, to, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)                                                        
                           }
                         }
                         if (addOneEdge) { #add a "hack" - an edge from vertex 1 to vertex 2 - so that the assess() function in NetRes doesn't crash
                           dag = set.arc(dag, from='v1', 
                                         to = 'v2' , check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
                         }
                         
                         ##and topologically sort it
                         nOrder <- topo_sort(as.igraph(dag))
                         
                         aM = as(as.graphNEL(dag), 'matrix')
                         aM = aM[nOrder, nOrder]
                         rDAG = as(aM, 'graphNEL')
                         dag = as.bn(rDAG)
                         
                         ##compute #s of direct children
                         iGraph = as.igraph(dag)
                         
                         outDegree = ego_size(
                           iGraph,
                           order = 1,
                           nodes = V(iGraph),
                           mode = c("out"),
                           mindist = 1
                         )
                         sorted = order(outDegree, decreasing = T)
                         sorted = V(iGraph)[sorted]
                         
                         self$dag = dag
                         self$vRank = names(sorted)
                         self$outDegree = outDegree
                       },                       
                       makePowerDag = function(numVertices = 20,  numChildren = 2, nParentOnly = 4, par1 = 1.2, 
                                               supersize.top.hubs = FALSE, supersize.factor = 2) {
                         ## generate random DAG
                         ##note that d is numChildren * 2 in this implementation of randDAG (symmetric parent and children couns)
                         rDAG <- randDAG(n = numVertices, d = numChildren*2, par1 = par1, method = 'barabasi')
                         
                         ##iGraph = as.igraph(as.bn(rDAG))
                         
                         ##and topologically sort it
                         nOrder <- as.character(topo_sort(as.igraph(as.bn(rDAG))))
                         
                         aM = as(rDAG, 'matrix')
                         aM = aM[nOrder, nOrder]
                         rownames(aM) = 1:nrow(aM)
                         colnames(aM) = 1:ncol(aM)
                         
                         ##The rows are parents and the columns are children
                         ##For top regulators, add edges if necessary
                         if (supersize.top.hubs && nParentOnly > 0) {
                           nOut = apply(aM, 1, function(x) { length(which(x != 0))})
                           top.hubs = order(nOut, decreasing = TRUE)[1:nParentOnly]
                           ##up to supersize.factor * nOut but no more than ncol(aM) - 1 total
                           n.edges.to.add = pmin(rep(ncol(aM) - 1, nParentOnly), supersize.factor*nOut[top.hubs]) - nOut[top.hubs]
                           ##for each hub
                           aM2 = aM
                           for (thi in 1:length(top.hubs)) {
                             th = top.hubs[thi]
                             ##find empty spaces
                             zIdx = setdiff(which(aM[th,] == 0), th)
                             useIdx = sample(zIdx, n.edges.to.add[thi])
                             ##sample values from a normal distribution of observed edges
                             aM2[th, useIdx] = rnorm(n.edges.to.add[thi], mean(aM[th,]), sd(aM[th,]))
                           }
                         }
                         ##continue here
                         ##AT THIS POINT, am2 is not a dag.  May need to remove parent edges to the hubs first
                         
                         rDAG = as(aM, 'graphNEL')
                         graph::nodes(rDAG) = paste('v', graph::nodes(rDAG), sep = '') #name (a little) better for convenience
                         
                         ##compute #s of direct children
                         iGraph = as.igraph(as.bn(rDAG))
                         
                         outDegree = ego_size(
                           iGraph,
                           order = 1,
                           nodes = V(iGraph),
                           mode = c("out"),
                           mindist = 1
                         )
                         sorted = order(outDegree, decreasing = T)
                         sorted = V(iGraph)[sorted]
                         
                         ##now kill the parents of the top direct drivers of the network - i.e., get rid of nodes with very high order but low out-degree
                         bnDAG = as.bn(rDAG)
                         for (curLat in names(sorted[1:nParentOnly])) {
                           curParents = parents(bnDAG, curLat)
                           for (cp in curParents) {
                             ##remove the edge from current parent to the current latent var
                             bnDAG = drop.arc(bnDAG, cp, curLat)
                           }
                         }
                         
                         self$dag = bnDAG
                         self$vRank = names(sorted)
                         self$outDegree = outDegree
                       }                       
                     )
)

addConfounderIgraph <- function(ig,
                                prob = 0.8,
                                ucoef = 1,
                                uname = "U",
                                exclude = "U.*") {
    allvars <- names(V(ig))
    toexc <- grep(exclude, allvars, value = T)
    allvars <- setdiff(allvars, toexc)
    ii <- rbinom(length(allvars), 1, prob)
    newedg <- map(allvars[ii == 1], function(vv) c(uname, vv)) %>% unlist()
    ucoef <- map(
        allvars[ii == 1],
        function(vv) {
            list(coefficients = set_names(c(ucoef), uname))
        }
    ) %>% set_names(allvars[ii == 1])
    ig_u <- ig %>%
        add_vertices(1, attr = list(name = uname)) %>%
        add_edges(newedg)
    attributes(ig_u)$ucoef <- ucoef
    return(ig_u)
}


##' Generates a simcausal object from a igraph object
##'
##' Take an igraph object and generated a simcausal by generating coefficients and distribution for each of the nodes of the network.
##' @title igraph2simcausal
##' @param ig igraph object
##' @param default_dist default distribution of the nodes in the inside of network. For now rnorm, rbern or rweibull
##' @param default_param list with parameters of the distribution. Depends on the default_dist. Ex, for rnorm the list should provide sd.
##' @param default_source_dist default distribution for source nodes.rnorm or rbern
##' @param default_transformation default transformation of the parent nodes. identiy, sigmoid or square.
##' @param default_coef default value for coefficient of parents of node. it can be a number or characters 'runif' or 'rnorm'
##' @param default_coef_square coefficients for the square term. Only applies if default_transformation is 'squared'. 'runif', 'rnorm', or number
##' @param default_coef_square_param if default_coef_square='runif' then list with min and max. if 'rnorm' then list with mean and sd.
##' @param default_source_param list with default parameter. Parameter depend on distribution. Ex, for rnorm you need to provide meand sd.
##' @param default_coef_param  if default_coef_square='runif' then list with min and max. if 'rnorm' then list with mean and sd.
##' @param specific_parameters list with parameters and distribution for specific nodes and parents.
##' @return simcausal object
##' @author Fred Gruber
igraph2simcausal <- function(ig,
			     default_dist = "rnorm",
			     default_param = list(sd = 1),
			     default_source_dist = "rnorm",
			     default_transformation = "identity",
			     default_coef = "runif",
			     default_coef_square = "runif",
			     default_coef_square_param = list(min = -2, max = 2),
			     default_source_param = list(mean = 0, sd = 1),
			     default_coef_param = list(min = -2, max = 2),
			     specific_parameters = list()) {
  require(simcausal)
  ## check arguments
  checkmate::assertClass(ig, "igraph")
  checkmate::assertChoice(default_dist, c("rnorm", "rbern", "rweibull"))
  checkmate::assertChoice(default_source_dist, c("rnorm", "rbern"))
  checkmate::assertChoice(default_transformation, c("identity", "sigmoid", "square"))
  checkmate::assertMultiClass(default_coef, classes = c("numeric", "character"))
  checkmate::assertList(default_coef_param)
  checkmate::assertList(specific_parameters)
  if(default_coef == "runif"&!all(c("min","max") %in% names(default_coef_param)))
    stop("For default_coef=='runif', default_coef_param should contain min and max")
  if(default_coef == "rnorm"&!all(c("mean","sd") %in% names(default_coef_param)))
    stop("For default_coef=='rnorm', default_coef_param should contain mean and sd")
  ## fix names
  V(ig)$name=make.names(names(V(ig)))
  ## start converting
  allnodes=as.bn(ig) %>% node.ordering()
  ##allnodes <- V(ig) %>% names()
  D <- DAG.empty()
  for (vv in allnodes) {
    parents <- igraph::neighbors(ig, vv, mode = "in") %>% names()
    ## node distributions
    vvdist <- specific_parameters[[vv]][["dist"]]
    ## are there parents?
    if (length(parents) == 0) {
      ## source
      source_params <-  specific_parameters[[vv]][["default_source_param"]]
      if (is.null(vvdist)) {
	vvdist <- default_source_dist
      }
      if(is.null(source_params))
	source_params <- default_source_param
      D <- D + node(vv,
		    distr = vvdist,
		    params = source_params
		    )
    }else{
      ## inner variables
      if (is.null(vvdist)) {
	vvdist <- default_dist
      }
      ## transformation for each parent
      alltrans = rep(default_transformation,length(parents)) %>% set_names(parents)
      providedtrans = specific_parameters[[vv]][["transformations"]]
      trans = specific_parameters[[vv]][["transformations"]][parents]
      ## coefficients for every parent
      coefs = specific_parameters[[vv]][["coefficients"]][parents]
      myform = NULL
      for(pp in parents){
	 ## if missing parameters use default
	if(is.null(trans[pp])||is.na(trans[pp])){
	  trans[pp] = default_transformation
	}
	if(is.null(coefs[pp])||is.na(coefs[pp])){
	  if(default_coef == "runif"){
	    coefs[pp] = runif(1,default_coef_param[["min"]],default_coef_param[["max"]])
	  }else if(is.numeric(default_coef))
	    coefs[pp] = default_coef
	}
	## generate formula
	if(trans[pp] == 'identity'){
	  myform = c(myform,
		     paste0(coefs[pp],"*",pp)
		     )
	}else if(trans[pp] == "sigmoid"){
	  myform = c(myform,
		     sprintf("%g/(1+exp(-%s))",coefs[pp],pp)
		     )
	}else if(trans[pp] == 'square'){
	  coefsq = coefs[paste0(pp,"^2")]
	  if(is.null(coefsq)|is.na(coefsq)){
	    if(default_coef_square == "runif"){
	      coefsq = runif(1,default_coef_square_param[["min"]],
			     default_coef_square_param[["max"]])
	    }else if(default_coef_square == 'rnorm'){
	      coefsq = rnorm(1,default_coef_square[["mean"]],default_coef_square[["sd"]])
	    }else if(is.numeric(default_coef_square)){
	      coefsq = default_coef_square
	    }
	  }
	  myform = c(myform,
		     sprintf("%g*%s+%g*%s^2",coefs[pp],pp,coefsq,pp)
		     )
	}
      }
      myform = paste0(myform,collapse = "+")
      if(vvdist == 'rnorm'){
	mysd = specific_parameters[[vv]][["sd"]][parents]
      if(is.null(mysd)||is.na(mysd))
	mysd = default_param$sd
      myparams = list(mean = myform,sd = mysd)
      }else if(vvdist == 'rbern'){
	myparams = list(prob = sprintf("1/(1+exp(-(%s)))",myform))
      }
      D <- D + node(vv,
		    distr = vvdist,
		    params = myparams
    )      
    }
  }
  return(D)
}
