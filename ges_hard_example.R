source("LatentConfounderBNlearnv2.R")
source('data_generating.R')

if (0) { #what generates the data structure
    p <- 50 # Number of obvserved variables 
    h <- 2 # Number of hidden variables
    n <- 500 # Number of samples
    set.seed(0)
    toy.data <- simulate.latent.DAG.data(nl = h, nv = p, ss = n, sp = 0.05)
    X <- toy.data$data # The observed data
    X <- scale(X)
    Xdf = data.frame(X)
    ##note - Xdf has latent vars in the first columns
}

Xdf = genData(
    N = 1000,
    V = 25, 
    num_out_drivers = 2,
    Npx = 3,
    latent_type = "Gaussian",
    Npu = 0, 
    uparentcoef=0.3, 
    Np = 10,
    parentcoef = 0.3,
    latout_coef = 0.5,
    latdriver_coef1 = 1,
    latdriver_coef2 = 1,
    latother_max = 0.2,
    obs_hub_coef = 0.5,
    addnoisesd = 0.5, 
    coef_driver1 = 1,
    coef_driver2 = 1,
    outbias = 0
)

##reshape our data frame to match the format of LRPS data generating function
inEdges = as.data.frame(attr(Xdf, 'allcoef'))
adjMat = matrix(0L, ncol = ncol(Xdf), nrow = ncol(Xdf),
                dimnames = list(colnames(Xdf), colnames(Xdf)))
for (outIdx in 1:nrow(inEdges)) {
    adjMat[inEdges[outIdx, 2], inEdges[outIdx, 1]] = inEdges[outIdx, 3]
}

##toy.data$true.obs.dag.amat

## a function to compute the precision-recall curve for a
## (custom) BNLearn ensemble
computePrForBnEns <- function(bnEns, trueNet) {
    ##iterate over all ensemble networks to get individual points
    ##for each network, calculate the PR curve based on edge strength
    library(ROCR)

    pred <- prediction(predictions, y);

    RP.perf <- performance(pred, "prec", "rec");
    
}

library(doParallel)
cl <- makeCluster(4) ## for multi-threading
registerDoParallel(cl)

clusterExport(cl = cl, c("estimateGraph_lrps0", "optimalFactorization",
                          "getOptNetwork0", "generate.data.for.GES"))



res_ges = getEnsemble2(Xdf,
                       Nboot = 2, ## used 2 for testing
                       algorithm = "GES",
                       parallel = T
                       )


results_ges = getResultsGES(res_ges$allnet, toy.data = list(true.obs.dag.amat = adjMat))


res_lrps2 = getEnsemble2(Xdf, 
                         Nboot = 2,
                         algorithm = "LRPS",
                         parallel = TRUE
                         )

results_lrps2 = getResultsGES(res_lrps2$allnet, toy.data = list(true.obs.dag.amat = adjMat))

rbind(
      results_ges %>%
      mutate(Method = "LatentConfounder + GES"),
      results_lrps2 %>%
      mutate(Method = "LRPS (bootstrapped) + GES")
      ) %>%
  ggplot(aes(x = rec.sk, y = prec.sk, colour = Method)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 10), se = FALSE) +
  ggpubr::theme_pubclean()
