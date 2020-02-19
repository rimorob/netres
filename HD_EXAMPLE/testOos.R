library(tidyverse)
library(data.table)
library(caret)
library(doParallel)
library(reticulate)

reticulate::use_condaenv("TFv1")
source('../LatentConfounderBNlearnv2.R')

source_python("../aupy.py")

##code to clear the GPU session from
##https://github.com/rstudio/keras/issues/739
py <- function(x){
    reticulate::py_run_string(x)
}
##py("from keras import backend as K")
py("from keras.backend.tensorflow_backend import set_session;")
py("from keras.backend.tensorflow_backend import clear_session;")
py("from keras.backend.tensorflow_backend import get_session;")

clear_gpu <- function(){
    ##py("cfg = tf.ConfigProto();")
    #3py("cfg.gpu_options.allow_growth = True;")
    #3k_clear_session()
    ##py("K.set_session(K.tf.compat.v1.Session(config=cfg))")
    py("K.clear_session();")
    print('cleared the session')
    ##py("from numba import cuda;")
    ##py("cuda.select_device(0);")
    ##py("cuda.close();")
}


##start global parameter block
NBOOT = 3
NITER = 10 #10 is the desired size
RAND_SEED = 1234
NVAR_PER_TYPE = 50 ##make a tiny network
NLATENT=200 #above this #, the validation error seems to go up
##end global parameter block


set.seed(RAND_SEED) #for reproducibility

if (!file.exists('trainSmall.csv')) {
    ##Load the CHDI training data
    train = fread('Striatum_naive_GeneExprColumnNames.csv', data.table=F)
    
    ##pre-process
    ##train$CAG = log(train$Input_CAG) #log transform or not?
    ##rank-transform CAG so that it's better-behaved?
    train$CAG = as.numeric(log(train$Input_CAG))
    train$Gender = as.numeric(train$Input_Gender)
    train$Age = as.numeric(train$Input_Age)
    ##drop some columns
    train$Input_Age = NULL
    train$Input_CAG = NULL
    train$Input_Gender = NULL
    train$Mouse.ID = NULL
    train$Input_AgeXCAG = NULL
    train$Input_AgeXCAG2 = NULL
    train$PCWeight = NULL
    miRnaIdx = grep('^ENS', colnames(train))
    train = train[, -miRnaIdx]
    
    psidx = grep('^Psych', colnames(train))
    train = train[, -psidx]
    
    ##find a few top gene expressions and 30 top proteins correlated with CAG
    pidx = grep('^Protein', colnames(train))
    gidx = grep('^GeneExpr', colnames(train))
    
    pCorCAG = sapply(pidx, function(pi) { cor.test(train[, pi], train$CAG, method='spearman')$p.value })
    gCorCAG = sapply(gidx, function(gi) { cor.test(train[, gi], train$CAG, method='spearman')$p.value })
    
    N = NVAR_PER_TYPE
    topProts = pidx[order(pCorCAG)[1:N]]
    topGenes = gidx[order(gCorCAG)[1:N]]
    ##keep one of each (gexpr, prot) so that we can look at gexpr->prot edge recovery
    keepVars = colnames(train)[c(topProts, topGenes)]
    prefs = as.character(sapply(keepVars, function(x) {gsub("^.+?_([^_]+).*?$", "\\1", x, perl=TRUE)}))
    prefs = c(unique(prefs), 'Htt') #add huntingtin
    prefs = paste('_?(', paste(prefs, collapse='|'), ')_?', sep='')
    geneProtCols = grep(prefs, colnames(train))

    trainSmall = data.frame(train[, c('Age', 'CAG', 'Gender', 'Weight')],
                            train[, geneProtCols])
    write.csv(trainSmall, 'trainSmall.csv', row.names=F)    
} else {
    trainSmall = read.csv('trainSmall.csv')
    trainSmall$Age = as.numeric(trainSmall$Age)
    trainSmall$Gender = as.numeric(trainSmall$Gender)    
}
    
fixedVars = c('Gender', 'CAG', 'Age')
expAndProt = setdiff(colnames(trainSmall), fixedVars)
blacklist = rbind(data.frame(from = rep(expAndProt, each=length(fixedVars)),
                             to = fixedVars,
                             stringsAsFactors=F))

##Build a reference model WITH CAG
cl <- makeCluster(4) ## for multi-threading
registerDoParallel(cl)

##learn the true model (with CAG)
resFull = getEnsemble2(trainSmall, blacklist = blacklist,
                      restart = 100, Nboot = NBOOT,
                      prior = "vsp",
                      score = "bge",
                      algorithm = 'tabu',
                      parallel = TRUE
                      )

##helper function to allow fitting of ensembles of length 1
getOneNetwork <- function(ens, nIdx) {
    oneNetEns = ens
    oneNetEns$Nboot = 1
    oneNetEns$boot_orders = oneNetEns$boot_orders[nIdx]
    ##edges, coef objects remains from the whole ensemble, this is a bug for now
    oneNetEns$fitmodels = oneNetEns$fitmodels[nIdx]
    oneNetEns$allnet = oneNetEns$allnet[nidx]
    return(oneNetEns)
}

##create cross-validation folds
folds <- createFolds(trainSmall$CAG, k = 5, list = TRUE, returnTrain = FALSE)
latPred = list()
##cross-validate learning of CAG
for (fi in 1:length(folds)) {
    ##Make a data frame with CAG latent
    trainlv = trainSmall[-folds[[fi]], ] %>% dplyr::select(-c(CAG))
    testlv = trainSmall[folds[[fi]], ]

    ##make a new blacklist
    fixedVars = c('Gender', 'Age') #latent - CAG
    expAndProt = setdiff(colnames(trainlv), c(fixedVars, 'CAG'))
    blacklistlv = rbind(data.frame(from = rep(expAndProt, each=length(fixedVars)),
                                   to = fixedVars,
                                   stringsAsFactors=T))


    resMissing = getEnsemble2(trainlv, blacklist = blacklistlv,
                               restart = 100, Nboot = NBOOT,
                               prior = "vsp",
                               score = "bge",
                               algorithm = 'tabu',
                               parallel = TRUE
                               )

    ##truecoef=getCoef(resMissing, 'Weight')

    discovered = latentDiscovery(
        resMissing,
        nItera=NITER,
        data = trainlv,
        seed=RAND_SEED,
        workpath="latentDiscoveryRep",
        freqCutoff = 0.01,
        maxpath = 1,
        alpha = 0.05,
        scale. = TRUE,
        method = "linear",
        latent_iterations = NLATENT,
        truecoef = NULL, #truecoef,
        truelatent = NULL,
        include_downstream = TRUE,
        multiple_comparison_correction = T,
        debug = F,
        include_output=FALSE,  
        output=NULL,
        parallel = TRUE,
        testdata=testlv %>% dplyr::select(-c(CAG)),
        isOrdinal=T,
        batch_size=50,
        nSamples=30,
        recompute_vars=F #using all vars anyway, since output is NULL
    )

    ##perform variable selection to find the relevant PCs on training data
    predicted = predict(discovered, trainlv)
    
    ##predicted = predict(discovered, testlv)
    print(paste(ncol(predicted$confounders), 'latent variables found'))
    dframe = as.data.frame(cbind(predicted$confounders,
                                 trainSmall[-folds[[fi]], ] %>% dplyr::select(CAG)))

    mCag = lm(CAG ~ ., data = dframe)
    mCag2 = stepAIC(mCag, k = 2) #use AIC as better for predicting

    ##now predict out of sample
    predicted = predict(discovered, testlv)
    
    dframe = as.data.frame(cbind(predicted$confounders, testlv %>% dplyr::select(CAG)))
    print(cor.test(testlv$CAG, predict(mCag2, newdata = dframe)))

    
    latPred[[fi]] = list(
        CAG = testlv$CAG,
        CAGpred = predict(mCag2, newdata = dframe),
        prediction = predicted,
        ##bestModelCAG = mCag2,
        testData = testlv)

    print('clearing GPU')
    clear_gpu()
    ##reticulate::use_condaenv("TFv1")

    ##source_python("../aupy.py")

}

allCAG = unlist(sapply(latPred, function(x) { x$CAG }))
allCAGpred = unlist(sapply(latPred, function(x) { x$CAGpred }))

write.csv(data.frame(allCAG, allCAGpred), 'linear.csv')

library(visreg)

m = lm(allCAGpred ~ allCAG)
r2 = signif(cor.test(allCAGpred, allCAG, method='pearson')$estimate^2, 2)
pdf('CAG_linear_boot3_iter10.pdf', width=11)
visreg(m, xlab = 'CAG', ylab = 'Predicted CAG', main=paste('CAG predicted R2:', r2))
dev.off()
