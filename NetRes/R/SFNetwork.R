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
                                               rescaleNoiseBy=function(mu=0.5, sigma=0.2, alpha=-5) {
                                                 return(rskew_normal(1, mu, sigma, alpha))
                                               }) {
                         ##recompute the rDAG
                         rDAG = self$dag #as.graphNEL(self$dag)
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
                             scaleFactor = max(0.1, rescaleNoiseBy()) #keep scale factor positive
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
                         
                         self$dag = rDAG
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
                         
                         self$dag = rDAG
                         self$vRank = names(sorted)
                         self$outDegree = outDegree
                       },                       
                       makePowerDag = function(numVertices = 20,  numChildren = 2, nParentOnly = 4, par1 = 1.2, 
                                               supersize.top.hubs = FALSE, supersize.factor = 2) {
                         ## generate random DAG
                         ##note that d is numChildren * 2 in this implementation of randDAG (symmetric parent and children couns)
                         rDAG <- randDAG(n = numVertices, d = numChildren*2, par1 = par1, method = 'barabasi')
                         
                         ##The rows are parents and the columns are children
                         ##For top regulators, add edges if necessary
                         if (supersize.top.hubs && nParentOnly > 0) {
                           ##find the top hubs
                           bnDAG = as.bn(rDAG)
                           ##nOrder <- as.character(topo_sort(as.igraph(bnDAG)))
                           nOrder <- RBGL::tsort(rDAG)    
                           allNodes = nodes(bnDAG)                        
                           
                           ##try to add nodes to the top hubs
                           for (thi in 1:nParentOnly) {
                             ##remove edges to parents 
                             curParents = parents(bnDAG, nOrder[thi])
                             for (parent in curParents) {
                               bnDAG = drop.arc(bnDAG, parent, nOrder[thi])
                             } 
                             ##how many children does this hub have?
                             curChildren = children(bnDAG, nOrder[thi])
                             ##figure out which nodes can be added  
                             curNotChildren = setdiff(allNodes, c(curChildren, nOrder[1:nParentOnly]))
                             ##and how many parents they have - use this to compute sampling weights:
                             ##nodes with fewer parents are more likely to be sampled in order to promote orthogonality
                             ctParents = sapply(curNotChildren, function(node) {
                               length(parents(bnDAG, node))
                             })
                             samplingWeights = 1/(1+ctParents)
                             ##add all available not-children or up to supersize*length(children), whichever is smaller
                             nEdgesToAdd = min(length(curNotChildren),
                                               (supersize.factor- 1)*length(curChildren))
                             ##try to add nodes
                             ##at present, if addition fails due to DAG constraint, the node is skipped silently without substitution
                             childrenToAdd = sample(curNotChildren, nEdgesToAdd, prob=samplingWeights, replace=F)
                             for (ci in 1:length(childrenToAdd)) {
                               bnDAG = set.arc(bnDAG, nOrder[thi], childrenToAdd[ci])
                               ##need to set arc strength, right?
                             }
                           }
                         } 

                         rDAG = as.graphNEL(bnDAG)
                         graph::nodes(rDAG) = paste('v', graph::nodes(rDAG), sep = '') #name (a little) better for convenience                         
                         nOrder = RBGL::tsort(rDAG)                         
                         nChildren = sapply(nOrder, function(v) {
                           length(graph::adj(rDAG, v)[[1]])
                         })
                         
                         ##re-sort the graph
                         ##convert to matrix for ease of certain operations
                         aM = as(rDAG, 'matrix')
                         aM = aM[nOrder, nOrder]
                         ##set edge weights to something from a normal distribution
                         aM[which(aM != 0)] = pmin(pmax(rnorm(length(which(aM != 0)), 0.5, 0.1), 0.1), 0.9) 

                         rDAG = as(aM, 'graphNEL')
                         bnDAG = as.bn(rDAG)
                         
                         self$dag = rDAG
                         self$vRank = nOrder
                         self$outDegree = nChildren
                       }                       
                     )
)
