#' @include NetRes.R

NetRes$set("private", "plotIgraph", function(g,
                                             edgelabels = F,
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
                                             layoutType = "dot",
                                             saveToFile = F,
                                             filename = "net.pdf",
                                             width = 1000 / 100,
                                             height = 1000 / 100,
                                             other_pdf_options = list(),
                                             nodeThickness = 1,
                                             nodeThickness_important = 2,
                                             fill = NULL,
                                             edge_color = NULL,
                                             edge_labels = NULL,
                                             label_pad = 2) {
  allnodes <- names(igraph::V(g))
  if (!is.null(fill)) {
    fill <- expandListRegex(fill, allnodes)
  }
  if (!is.null(edge_color)) {
    edge_color <- expandDfRegex(edge_color, allnodes)
  }
  ## library(Rgraphviz)
  node_width <- NULL
  node_height <- NULL
  node_fixedSize <- FALSE
  gr <- igraph:::as_graphnel(g)
  eAttrs <- list()
  # w <- w[setdiff(seq(along=w), removedEdges(gr))]
  if (length(igraph::get.edge.attribute(g)) == 0) {
    edgelabels <- FALSE
    edgeweights <- FALSE
  }
  if (edgeweights | edgelabels) {
    w <- signif(igraph::get.edge.attribute(g)[[1]], 2)
    names(w) <- sub("\\|", "~", attributes(igraph::E(g))$vnames)
    ## names(w) <- edgeNames(gr, recipEdges="distinct")
    ## names(eAttrs$lwd) = edgeNames(gr, recipEdges="distinct")
  }
  if (edgelabels) {
    if (!is.null(edge_labels)) {
      edgs <- paste(edge_labels[[1]], edge_labels[[2]], sep = "~")
      alledges <- graph::edgeNames(gr)
      alledgeslab <- rep("", length(alledges))
      names(alledgeslab) <- alledges
      alledgeslab[edgs] <- as.character(edge_labels[[3]])
      namesw <- names(alledgeslab)
      alledgeslab <- paste0(paste0(rep(" ", label_pad), collapse = ""), alledgeslab)
      names(alledgeslab) <- namesw
      eAttrs$label <- alledgeslab
    } else {
      if (edgelabelsFilter_useabs) {
        iiw <- which(abs(w) > edgelabelsFilter)
      } else {
        iiw <- which(w > edgelabelsFilter)
      }
      eAttrs$label <- rep("", length(w))
      names(eAttrs$label) <- names(w)
      wval <- w[iiw]
      wval <- paste0(paste0(rep(" ", label_pad), collapse = ""), wval)
      eAttrs$label[names(w[iiw])] <- wval
    }
  }
  if (edgeweights) {
    wn <- as.numeric(w)
    eAttrs$lwd <- (lwdmax - lwdmin) * (wn - min(wn)) / (max(wn) - min(wn)) + lwdmin
    names(eAttrs$lwd) <- sub("\\|", "~", attributes(igraph::E(g))$vnames)
  }
  eAttrs$direction <- rep("forward", length(graph::edgeNames(gr, recipEdges = "distinct")))
  names(eAttrs$direction) <- graph::edgeNames(gr, recipEdges = "distinct")
  ## edge colors
  if (!is.null(edge_color)) {
    alledges <- graph::edgeNames(gr)
    alledgescol <- rep("black", length(alledges))
    names(alledgescol) <- alledges
    edgs <- paste(edge_color[[1]], edge_color[[2]], sep = "~")
    alledgescol[edgs] <- as.character(edge_color[[3]])
    eAttrs$color <- alledgescol
  }
  attrs <- list(
    node = list(
      shape = "ellipse",
      fixedsize = node_fixedSize,
      width = node_width,
      height = node_height,
      lwd = nodeThickness,
      color = "black"
    ),
    edge = list(
      direction = "forward",
      concentrate = F
    ),
    graph = list(
      damping = damping,
      nodesep = nodesep,
      pad = pad,
      ranksep = ranksep,
      splines = splines,
      start = start,
      dim = 2,
      sep = sep,
      concentrate = FALSE,
      overlap = overlap
    )
  )
  nAttrs <- list()
  gr <- Rgraphviz::layoutGraph(gr,
    layoutType = layoutType,
    attrs = attrs,
    nodeAttrs = nAttrs,
    edgeAttrs = eAttrs,
    recipEdges = "distinct"
  )
  ## add filler
  nAttrs$fill <- fill
  ## add node edge width
  if (!is.null(node_width)) {
    nAttrs$width <- node_width
    ## graph::nodeRenderInfo(gr) = list(
    ##     width = node_width)
  }
  if (!is.null(node_height)) {
    nAttrs$height <- node_height
    ## graph::nodeRenderInfo(gr) = list(
    ##     height = node_height)
  }
  nAttrs$fixedSize <- node_fixedSize
  graph::nodeRenderInfo(gr) <- nAttrs
  ## graph::nodeRenderInfo(gr) =list(fixedSize = node_fixedSize)
  graph::edgeRenderInfo(gr) <- eAttrs
  if (saveToFile) {
    other_pdf_options <- c(other_pdf_options,
      file = filename,
      width = width,
      height = height
    )
    do.call(pdf, args = other_pdf_options)
    # pdf(filename, width = width, height = height)
    Rgraphviz::renderGraph(gr)
    dev.off()
  } else {
    Rgraphviz::renderGraph(gr)
  }
})


expandListRegex <- function(mylist, allnames) {
  newlist <- list()
  for (ll in 1:length(mylist)) {
    rg <- names(mylist)[ll]
    rg <- paste0("^", rg, "$")
    val <- mylist[[ll]]
    mv <- grep(rg, allnames, value = T)
    tmp <- rep(list(val), length(mv))
    names(tmp) <- mv
    newlist <- c(newlist, tmp)
  }
  return(newlist)
}

expandDfRegex <- function(mydf, allnames) {
  newlist <- map_df(
    1:nrow(mydf),
    function(ii) {
      val <- mydf[ii, ]
      inputs <- grep(paste0("^", val$inp, "$"),
        allnames,
        value = T
      )
      outputs <- grep(paste0("^", val$out, "$"),
        allnames,
        value = T
      )
      expand.grid(inputs, outputs, stringsAsFactors = F) %>%
        mutate(color = val[[3]])
    }
  )
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
NetRes$set("private", "network_performance", function(true_igraph, edges, Nboot = 200, ci = FALSE, seed = 123,cutoff=0.5,oracle=NULL,debug=FALSE) {
  require(ROCR)
  require(precrec)
  require(SID)
  checkmate::assertClass(true_igraph, "igraph")
  checkmate::assertDataFrame(edges)
  
  ## convert edges to igraph object
  est_ig <- igraph::graph_from_data_frame(edges)
  ## get adjacency matrix
  true_adj <- igraph::get.adjacency(true_igraph) %>% as.matrix()
  est_adj <- est_ig %>%
    igraph::get.adjacency(attr = "freq") %>%
      as.matrix()
    ## make sure nodes are ordered the same way
  est_adj <- est_adj[colnames(true_adj), colnames(true_adj)]
  ## get edgelist with ALL edges
  true_edges <- reshape2::melt(true_adj) %>%
      filter(Var1 != Var2) ## remove self edges
  if(sum(true_edges$value)>0){
      ## not all 0. Calculate SID
      ## get structural interention distance
      sid <- SID::structIntervDist(true_adj, est_adj)
  }
  est_edges <- reshape2::melt(est_adj) %>%
    filter(Var1 != Var2) ## remove self edges
  comb_edges <- full_join(
    true_edges %>%
      rename(True = value),
    est_edges %>%
      rename(Prob = value),
    by = c("Var1", "Var2")
  )
  if(!is.null(oracle)){
      oracle_ig=igraph::graph_from_data_frame(oracle)
      oracle_adj <- oracle_ig %>%
          igraph::get.adjacency(attr = "freq") %>%
          as.matrix()
      oracle_adj <- oracle_adj[colnames(true_adj), colnames(true_adj)]
      oracle_edges=reshape2::melt(oracle_adj) %>%
          filter(Var1 != Var2) ## remove self edges
      comb_edges=full_join(
          comb_edges,
          oracle_edges %>% rename(Oracle=value)
      )
      oracleest=comb_edges$Oracle
      if(debug){
          message("in perormance assesment before roc.test")
          browser()
      }
  }

  trueval <- comb_edges$True
  scores <- comb_edges$Prob
  
  if(length(unique(trueval))==2){
      if(!is.null(oracle))
          roctest=pROC::roc.test(comb_edges$True,comb_edges$Prob,comb_edges$Oracle)
      ## there are positive and negative cases
      pred <- prediction(comb_edges$Prob, comb_edges$True)
      perf <- performance(pred, "lift", "rpp")
      RPP <- perf@x.values[[1]]
      lift <- perf@y.values[[1]]
      ggp_lift <- qplot(RPP, lift, geom = "line") + theme_light() + xlab(perf@x.name) + ylab(perf@y.name) + ggtitle("Lift Graph")
      if (ci) {
          ## generate confidence intervals
          if(is.null(oracle)){
              ## no oracle
              set.seed(seed)
              ids=replicate(Nboot,sample(1:length(scores),replace=TRUE))
              resampled_scores= list(
                  scores=map(1:Nboot,
                             function(vv){
                                 scores[ids[,vv]]
                             }),
                  labels=map(1:Nboot,
                             function(vv){
                                 trueval[ids[,vv]]
                             }),
                  modnames="Inferred",
                  dsids=1:Nboot
              ) %>% set_names(c("scores","labels","modnames","dsids"))
              }else{
                  ## with oracle
                  set.seed(seed)
                  ids=replicate(Nboot,sample(1:length(scores),replace=TRUE))
                  resampled_scores= list(
                      scores=map(1:Nboot,
                                 function(vv){
                                     list(scores[ids[,vv]],
                                          oracleest[ids[,vv]]
                                          )
                                 }),
                      labels=map(1:Nboot,
                                 function(vv){
                                     list(trueval[ids[,vv]],
                                          trueval[ids[,vv]]
                                          )
                                 }),
                      modnames=c("Inferred",'Oracle'),
                      dsids=1:Nboot
                  ) %>% set_names(c("scores","labels","modnames","dsids"))
                  ##ressampled_oracle=replicate(Nboot,sample(oracleest,replace=TRUE))
              }
              cidat <- mmdata(resampled_scores[["scores"]],
                              resampled_scores[["labels"]],
                              modnames = resampled_scores[["modnames"]],
                              dsids = resampled_scores[["dsids"]]
                              )
      } else {
          ## no CI
          if(!is.null(oracle)){
              scores_withora= list(
                  scores=list(
                      list(scores,
                           oracleest
                           )
                  ),
                  labels=trueval,
                  modnames=c("Inferred",'Oracle'),
                  dsids=1
              )
              cidat <- mmdata(scores_withora[["scores"]],
                              scores_withora[["labels"]],
                              modnames=scores_withora[["modnames"]],
                              dsids=scores_withora[["dsids"]]
                              )
          }else{
              cidat <- mmdata(scores, trueval)
          }
      }
      ## get AUC
      sscurves_auc <- precrec::evalmod(cidat)
      if(!is.null(oracle)){
          aucs <- auc(sscurves_auc) %>%
              filter(modnames=="Inferred") %>% 
              group_by(curvetypes) %>%
              summarise(AUC = mean(aucs))
      }else{
          aucs <- auc(sscurves_auc) %>%
              group_by(curvetypes) %>%
              summarise(AUC = mean(aucs))
      }
      auc_roc <- filter(aucs, curvetypes == "ROC")$AUC
      auc_prc <- filter(aucs, curvetypes == "PRC")$AUC
      ## get other metrics
      sscurves_bas <- precrec::evalmod(cidat, mode = "basic")
      ## generate plots
      prc_random <- sum(true_edges$value) / length(true_edges$value)
      ggp_fm <- autoplot(sscurves_bas, "fscore",show_cb=ci)
      ggp_prc <- autoplot(sscurves_auc, "PRC",show_cb=ci) +
          ggtitle(sprintf(
              "PR Curve. AUC = %0.2g.\nRandom: %0.2g",
              signif(auc_prc, 2), prc_random
          ))
      ggp_roc <- autoplot(sscurves_auc, "ROC",show_cb=ci)
      if(!is.null(oracle)){
          ggp_roc=ggp_roc+ggtitle(sprintf(
              "ROC. AUC = %0.2g.\nRandom: %0.2g. \nroc.test pvalue =%s",
              signif(auc_roc, 2), 0.5, roctest$p.value
          ))
      }else{
          ggp_roc=ggp_roc+ggtitle(sprintf(
              "ROC. AUC = %0.2g.\nRandom: %0.2g",
              signif(auc_roc, 2), 0.5
          ))
          }
      final_ggp <- (ggp_prc / ggp_roc) | (ggp_fm / ggp_lift) + plot_layout(width = 1)
      attributes(final_ggp)$auc_pr_random <- prc_random
      attributes(final_ggp)$edges <- comb_edges
      attributes(final_ggp)$other <- sscurves_bas
      attributes(final_ggp)$perf_auc_roc <- list(names="rocAUC",value=auc_roc)
      attributes(final_ggp)$perf_auc_pr <- list(names="prAUC",value=auc_prc)
      attributes(final_ggp)$perf_sid <- list(names="SID",value=sid$sid) 
  }else{
      message("True vector only have 1 level")
      pred <- prediction(comb_edges$Prob, factor(comb_edges$True,levels=c(0,1)))
      perf <- performance(pred, "tnr")
      tnr=perf@y.values[[1]]
      tnr_th=pracma::interp1(
                        perf@x.values[[1]],tnr,
                        cutoff
                    )%>% suppressWarnings()
      ggp_tnr <- qplot(perf@x.values[[1]],tnr, geom = "line") + theme_light() + xlab(perf@x.name) + ylab(perf@y.name) + ggtitle("True Negative Rate")
      perf <- performance(pred, "npv")
      npv <- perf@y.values[[1]]
      npv_th=pracma::interp1(
                         perf@x.values[[1]],npv,
                         cutoff
                     )%>% suppressWarnings()
      ggp_npv <- qplot(perf@x.values[[1]],npv, geom = "line") + theme_light() + xlab(perf@x.name) + ylab(perf@y.name) + ggtitle("NPV")     
      perf=performance(pred,"fpr")
      fpr=perf@y.values[[1]]
      fpr_th=pracma::interp1(
                         perf@x.values[[1]],fpr,
                         cutoff
                     ) %>% suppressWarnings()
      ggp_fpr <- qplot(perf@x.values[[1]],fpr, geom = "line") + theme_light() + xlab(perf@x.name) + ylab(perf@y.name) + ggtitle("FPR")
      final_ggp <- (ggp_fpr+ggp_npv) / ggp_tnr + plot_layout(width = 1)
      attributes(final_ggp)$perf_tnr <- list(names='TNR',value=tnr_th)
      attributes(final_ggp)$edges <- comb_edges
      attributes(final_ggp)$perf_fpr <- list(names="FPR",value=fpr_th)
      attributes(final_ggp)$perf_npv <- list(names="NPV",value=npv_th)
  }
  return(final_ggp)
})


plot_auc <- function(perf, title = "") {
  rocgg <- minet::rates(perf) %>%
    ggplot(aes(x = fpr, y = tpr)) +
    geom_line() +
    geom_abline(
      intercept = 0,
      slope = 1,
      colour = "red",
      linetype = "dashed",
      alpha = 0.5
    ) +
    xlab("FP Rate") +
    ylab("TP Rate") +
    ggtitle(sprintf(
      "ROC AUC=%0.2g",
      minet::auc.roc(perf)
    ))
  prgg <- minet::pr(perf) %>%
    ggplot(aes(x = r, y = p)) +
    geom_point() +
    geom_hline(
      yintercept = attributes(perf)$auc_pr_random,
      colour = "red",
      linetype = "dashed",
      alpha = 0.5
    ) +
    xlab("Recall (TP Rate)") +
    ylab("Precision") +
    ggtitle(sprintf(
      "PR AUC=%0.2g (%0.2g)",
      minet::auc.pr(perf),
      attributes(perf)$auc_pr_random
    ))
  require(patchwork)
  allplot <- (rocgg / prgg) +
    plot_annotation(title = paste0("Network Inference Performance. ", title))
  allplot & theme_light()
}


calculate.f1 <- function(predicted, actual) {
  suppressPackageStartupMessages(library("caret"))
  iinna <- which(!is.na(actual))
  predicted <- predicted[iinna]
  actual <- actual[iinna]
  cm <- confusionMatrix(data = predicted, reference = actual, positive = "1")
  cm$byClass[["F1"]]
}


prcAUC <- function(ytrue, yprob) {
  if (class(ytrue) == "factor") {
    ytrue <- as.numeric(as.vector(ytrue))
  }
  if (class(yprob) == "character" && yprob == "random") {
    return(sum(ytrue) / length(ytrue))
  }
  if (length(unique(ytrue)) == 1) {
    warning("label had only 1 level")
    return(NA)
  }
  iinna <- which(!is.na(ytrue))
  yprob <- yprob[iinna]
  ytrue <- ytrue[iinna]
  ## library(precrec)
  sscurves <- precrec::evalmod(scores = yprob, labels = ytrue)
  filter(precrec::auc(sscurves), curvetypes == "PRC")$aucs
}

rocAUC <- function(ytrue, yprob) {
  if (class(ytrue) == "factor") {
    ytrue <- as.numeric(as.vector(ytrue))
  }
  if (class(yprob) == "character" && yprob == "random") {
    return(0.5)
  }
  if (length(unique(ytrue)) == 1) {
    warning("label had only 1 level")
    return(NA)
  }
  ## library(precrec)
  sscurves <- precrec::evalmod(scores = yprob, labels = ytrue)
  filter(precrec::auc(sscurves), curvetypes == "ROC")$aucs
}

