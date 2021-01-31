source("LatentConfounderBNlearnv2.R")


p <- 50 # Number of obvserved variables 
h <- 2 # Number of hidden variables
n <- 500 # Number of samples
set.seed(0)
toy.data <- simulate.latent.DAG.data(nl = h, nv = p, ss = n, sp = 0.05)
X <- toy.data$data # The observed data
X <- scale(X)
Xdf = data.frame(X)


library(doParallel)
cl <- makeCluster(5) ## for multi-threading
registerDoParallel(cl)

clusterExport(cl = cl, c("estimateGraph_lrps0", "optimalFactorization",
                          "getOptNetwork0", "generate.data.for.GES"))



res_ges = getEnsemble2(Xdf,
                       Nboot = 2,$ = ## used 2 for testing
                       algorithm = "GES",
                       parallel = T
                       )


results_ges = getResultsGES(res_ges$allnet)


res_lrps2 = getEnsemble2(Xdf, 
                         Nboot = 2,
                         algorithm = "LRPS",
                         parallel = F,
                         debug = T
                         )

results_lrps2 = getResultsGES(res_lrps2$allnet)

rbind(
      results_ges %>%
      mutate(Method = "GES"),
      results_lrps2 %>%
      mutate(Method = "LRPS boot")
      ) %>%
  ggplot(aes(x = rec.sk, y = prec.sk, colour = Method)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 10), se = FALSE) +
  ggpubr::theme_pubclean()
