source("LatentConfounderBNlearn.R")
load("final_model_nolvp_novp.RData", verbose = T)

train = datalist$data_noisy
trainlv = train %>% dplyr::select(-U1.out, -U2.out)


blacklistlv = rbind(data.frame(from = "Z", to = colnames(trainlv)))

library(doParallel)
cl <- makeCluster(3) ## for multi-threading
registerDoParallel(cl)



res_missing_small = getEnsemble2(trainlv, blacklist = blacklistlv,
			  restart = 100, Nboot = 10,
			  prior = "vsp",
			  score = "bge",
			  algorithm = 'hc',
			  parallel = TRUE
			  )


test_wrong =  latentDiscovery(
    res_missing_small,
    nItera=5,
    data = trainlv,
    "Z",
    workpath="pca_wrong",
    freqCutoff = 0.01,
    maxpath = 1,
    alpha = 0.01,
    scale. = TRUE,
    method = "linear",
    latent_iterations = 100,
    truecoef = datalist$coef %>% filter(output=="Z"),
    truelatent=train %>% dplyr::select("U1.out","U2.out"),
    include_downstream = TRUE,
    include_output = TRUE,
    multiple_comparison_correction = T,
    debug = F,
    parallel = TRUE
)


source("LatentConfounderBNlearnv2.R")

test_right =  latentDiscovery(
    res_missing_small,
    nItera=10,
    data = trainlv,
    "Z",
    workpath="pca_right",
    freqCutoff = 0.01,
    maxpath = 1,
    alpha = 0.05,
    scale. = TRUE,
    method = "linear",
    latent_iterations = 100,
    truecoef = datalist$coef %>% filter(output=="Z"),
    truelatent=train %>% dplyr::select("U1.out","U2.out"),
    include_downstream = TRUE,
    include_output = TRUE,
    multiple_comparison_correction = T,
    debug = F,
    parallel = TRUE,
    wrongway = FALSE ## this undo the fix in getGraphResiduals
)

latvar_learned = predict(test_right, newdata = trainlv)

latvar_learned$confounders %>% cor(train[, c("U1.out", "U2.out")])

test_right$confounders %>% cor(train[, c("U1.out", "U2.out")])

test_right$confounders %>% cor(latvar_learned$confounders)
