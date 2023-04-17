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
                       #' @param numNeighbors number of neighbors of each vertex, on average (i.e., node denisty, default=4)
                       #' @param nParentOnly how many nodes to make parent-only hubs by removing their high order but low out-degree parents (default=4)
                       #' @param addOneEdge whether to add a single edge not involving a latent variable to star or twin-star topologies (helps with AUC calculation; default=TRUE)
                       #' @return A `SFNetwork` object
                       initialize = function(numVertices = 20, topology='power', numNeighbors = 4, nParentOnly = 4, addOneEdge=TRUE) {
                         if (topology == 'power') {
                           res = private$makePowerDag(numVertices = numVertices, numNeighbors = numNeighbors, nParentOnly = nParentOnly)
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
                                               rescaleNoiseBy=function(mu=1, sigma=0.15, alpha=-5) {
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
                       makePowerDag = function(numVertices = 20,  numNeighbors = 4, nParentOnly = 4) {
                         ## generate random DAG
                         ##rDAG <- randDAG(n = numVertices, method = 'power', d = numNeighbors)
                         rDAG <- randDAG(n = numVertices, method = 'barabasi', d = numNeighbors)
                         
                         ##iGraph = as.igraph(as.bn(rDAG))
                         
                         ##and topologically sort it
                         nOrder <- as.character(topo_sort(as.igraph(as.bn(rDAG))))
                         
                         aM = as(rDAG, 'matrix')
                         aM = aM[nOrder, nOrder]
                         rownames(aM) = 1:nrow(aM)
                         colnames(aM) = 1:ncol(aM)
                         
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
