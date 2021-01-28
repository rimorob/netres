## Simple Example
##   :PROPERTIES:
##   :ID:       A2741453-06EB-414B-9679-496BFDBCF0DF
##   :END:
## First we load the codes and dataset.

## [[id:A2741453-06EB-414B-9679-496BFDBCF0DF][Simple Example:1]]
## Load libraries
source("LatentConfounderBNlearnv2.R")
## load test data
suppressWarnings(load("final_model_nolvp_novp.RData", verbose = T))

figtrue = igraph::graph_from_data_frame(datalist$coef)
## Simple Example:1 ends here



## For the first example the true causal structure is shown below:

## [[id:A2741453-06EB-414B-9679-496BFDBCF0DF][Simple Example:2]]
plotIgraph(figtrue, layout = 'neato', sep = 0.0001,
                       fill = list("U.*" = 'darksalmon',
                                   "Z" = 'yellow',
                                   "^V1$|V2$" = 'skyblue'))
## Simple Example:2 ends here

## Structure learning with no missing  variables
##    :PROPERTIES:
##    :ID:       4D522327-83A0-40C8-BB3A-550ECC577C96
##    :END:
## If we observe all the variables we can get a good estimate the graph structure:

## [[id:4D522327-83A0-40C8-BB3A-550ECC577C96][Structure learning with no missing  variables:1]]
train = datalist$data_noisy
library(bnlearn)

blacklist = rbind(data.frame(from = "Z", to = colnames(train)),
		      data.frame(from = colnames(train), to = "U1.out"),
		      data.frame(from = colnames(train), to = "U2.out")
                  )

library(doParallel)
cl <- makeCluster(10) ## for multi-threading
registerDoParallel(cl)
## Structure learning with no missing  variables:1 ends here

## [[id:4D522327-83A0-40C8-BB3A-550ECC577C96][Structure learning with no missing  variables:2]]
res_alldata = getEnsemble2(train, blacklist = blacklist,
                           restart = 100, Nboot = 10,
                           prior = "vsp",
                           score = "bge",
                           algorithm = 'hc',
                           parallel = TRUE
                           )
## Structure learning with no missing  variables:2 ends here



## #+RESULTS:
## : [1] "Distributing ensemble learning"
## : Bootstrapping bnlearn 10 times.
## : calculate coefficients


## [[id:4D522327-83A0-40C8-BB3A-550ECC577C96][Structure learning with no missing  variables:3]]
plot(res_alldata,
     "Z",
     cutoff = 0.3,
     freqth = 0.3, 
     maxpath=3,
     nodesep=0.01,
     sep = 0.01,
     layout = 'neato',
     edge_labels = 'frequency',
     edgeweights=T,
     edgelabelsFilter = 0.5,
     edge_color=tribble(~input,~output,~color,"V1","Z","red","V2","Z","red",
                        "U1.out", "Z", "red",
                        "U2.out", "Z", "red"),
     fill = list("U.*" = 'darksalmon',
                 "Z" = 'yellow',
                 "^V1$|V2$" = 'skyblue'))
## Structure learning with no missing  variables:3 ends here

## Reconstruction with missing confounders.
##    :PROPERTIES:
##    :ID:       8D9528DA-7AD2-4F73-B959-CB3704DBDB38
##    :END:
## Suppose that now we try to reconstruct the structure without observing U1.out and U2.out.
## The resulting graph will be more complicated since many variables will try to compensate for the 
## missing U1.out and U2.out.

## [[id:8D9528DA-7AD2-4F73-B959-CB3704DBDB38][Reconstruction with missing confounders.:1]]
trainlv = train %>% dplyr::select(-U1.out, -U2.out)

blacklistlv = rbind(data.frame(from = "Z", to = colnames(trainlv)))
seed = 123
set.seed(seed)
res_missing = getEnsemble2(trainlv, blacklist = blacklistlv,
			  restart = 100, Nboot = 10,
			  prior = "vsp",
			  score = "bge",
			  algorithm = 'tabu',
			  parallel = TRUE
			  )
## Reconstruction with missing confounders.:1 ends here

## [[id:8D9528DA-7AD2-4F73-B959-CB3704DBDB38][Reconstruction with missing confounders.:2]]
plot(res_missing,"Z",
     cutoff = 0.3,
     freqth = 0.3, 
     maxpath=3,
     nodesep=0.01,
     sep = 0.01,
     layout = 'dot',
     edge_labels = 'frequency',
     edgeweights=T,
     edgelabelsFilter = 0.5,
     edge_color=tribble(~input,~output,~color,"V1","Z","red","V2","Z","red",
                        "U1.out", "Z", "red",
                        "U2.out", "Z", "red"),
     fill = list("U.*" = 'darksalmon',
                 "Z" = 'yellow',
                 "^V1$|V2$" = 'skyblue'))
## Reconstruction with missing confounders.:2 ends here

## Reconstruction estimating missing confounders
##    :PROPERTIES:
##    :ID:       71A1CEAA-2F6B-4EDC-BB69-C011BC8D306E
##    :END:
## Here we used the propose algorithm to learn the missing variables and improve the structure search.

## [[id:71A1CEAA-2F6B-4EDC-BB69-C011BC8D306E][Reconstruction estimating missing confounders:1]]
seed = 123
set.seed(seed)
graphics.off()
latvar_simple =  latentDiscovery(
    res_missing,
    nItera=5,
    data = trainlv,
    "Z",
    workpath="pca_simple",
    method = "linear",
    truecoef = datalist$coef %>% filter(output=="Z"),
    truelatent=datalist$data %>% dplyr::select("U1.out","U2.out"),
    parallel = TRUE
)
## Reconstruction estimating missing confounders:1 ends here

## [[id:71A1CEAA-2F6B-4EDC-BB69-C011BC8D306E][Reconstruction estimating missing confounders:2]]
plot(latvar_simple$details$final_ensemble,
     "Z",
     cutoff = 0.3,
     freqth = 0.3, 
     maxpath=3,
     nodesep=0.01,
     sep = 0.01,
     layout = 'dot',
     edge_labels = 'frequency',
     edgeweights=T,
     edgelabelsFilter = 0.5,
     edge_color=tribble(~input,~output,~color,"V1","Z","red","V2","Z","red",
                        "U1.out", "Z", "red",
                        "U2.out", "Z", "red"),
     fill = list("U.*" = 'darksalmon',
                 "Z" = 'yellow',
                 "^V1$|V2$" = 'skyblue'))
## Reconstruction estimating missing confounders:2 ends here

## Cleanup
##    :PROPERTIES:
##    :ID:       981E1400-F763-4255-808C-69D290E0DBDA
##    :END:
## Close the parallelization cluster:

## [[id:981E1400-F763-4255-808C-69D290E0DBDA][Cleanup:1]]
stopCluster(cl)
## Cleanup:1 ends here

## True network
##    :PROPERTIES:
##    :ID:       333CAB53-F7CF-4880-ABBC-1155942E84F0
##    :END:

## [[id:333CAB53-F7CF-4880-ABBC-1155942E84F0][True network:1]]
source("LatentConfounderBNlearnv2.R")
suppressWarnings(load("final_model_nolvp_withvp.RData", verbose = T))
datalist_med = datalist

figtrue_med = igraph::graph_from_data_frame(
			  datalist_med$coef %>%
			  filter(input %in% c(paste0("P", 1:10),
					      paste0("V", 45:50),
					      paste0("V", 1:5),
					      "U1.out", "U2.out", 'Z'),
				 output %in% c("Z","U1.out", "U2.out",
					       paste0("V", 40:50),
					      paste0("V", 1:10)))
)
## True network:1 ends here

## [[id:333CAB53-F7CF-4880-ABBC-1155942E84F0][True network:2]]
plotIgraph(figtrue_med, layout = 'dot', nodesep = 0.00001,
		       fill = list("U.*" = 'darksalmon',
				   "Z" = 'yellow',
				   "^V1$|V2$" = 'skyblue')
           )
## True network:2 ends here

## Structure search with missing variables
##    :PROPERTIES:
##    :ID:       E415ABAC-C801-411D-9BAA-BC28C885AD65
##    :END:


## [[id:E415ABAC-C801-411D-9BAA-BC28C885AD65][Structure search with missing variables:1]]
train_med = datalist_med$data_noisy

trainlv_med = train_med %>% dplyr::select(-U1.out, -U2.out)

blacklistlv_med = rbind(data.frame(from = "Z", to = colnames(trainlv_med)))

seed = 123
set.seed(seed)
res_missing_med_small = getEnsemble2(trainlv_med, blacklist = blacklistlv_med,
			    restart = 100, Nboot = 10,
			    prior = "vsp",
			    score = "bge",
			    algorithm = 'tabu',
			    parallel = TRUE
			    )
## Structure search with missing variables:1 ends here

## [[id:E415ABAC-C801-411D-9BAA-BC28C885AD65][Structure search with missing variables:2]]
plot(res_missing_med_small,
     "Z",
     cutoff = 0.5,
     freqth = 0.5, 
     maxpath=2,
     nodesep=0.01,
     sep = 0.01,
     layout = 'dot',
     edge_labels = 'frequency',
     edgeweights=T,
     edgelabelsFilter = 0.5,
     edge_color=tribble(~input,~output,~color,"V1","Z","red","V2","Z","red",
                        "U1.out", "Z", "red",
                        "U2.out", "Z", "red"),
     fill = list("U.*" = 'darksalmon',
                 "Z" = 'yellow',
                 "^V1$|V2$" = 'skyblue'))
## Structure search with missing variables:2 ends here

## Structure search estimating latent variables
##    :PROPERTIES:
##    :ID:       DC5D2AE4-2C27-4DA6-B0F8-89FE7A3E8172
##    :END:


## [[id:DC5D2AE4-2C27-4DA6-B0F8-89FE7A3E8172][Structure search estimating latent variables:1]]
seed = 123
set.seed(seed)
graphics.off()
medium_evo = latentDiscovery(
	res_missing_med_small,
	nItera=5,
	data = trainlv_med,
	"Z",
	seed=seed,
	workpath="latentDiscovery_med_linear",
	method = "linear",
	truecoef = datalist_med$coef %>% filter(output=="Z"),
	truelatent=datalist_med$data %>% dplyr::select("U1.out","U2.out"),
	parallel = TRUE
    )
## Structure search estimating latent variables:1 ends here

## [[id:DC5D2AE4-2C27-4DA6-B0F8-89FE7A3E8172][Structure search estimating latent variables:2]]
plot(medium_evo$details$final_ensemble,
     "Z",
     cutoff = 0.5,
     freqth = 0.5, 
     maxpath=2,
     nodesep=0.01,
     sep = 0.01,
     layout = 'dot',
     edge_labels = 'frequency',
     edgeweights=T,
     edgelabelsFilter = 0.5,
     edge_color=tribble(~input,~output,~color,"V1","Z","red","V2","Z","red",
                        "U1.out", "Z", "red",
                        "U2.out", "Z", "red"),
     fill = list("U.*" = 'darksalmon',
                 "Z" = 'yellow',
                 "^V1$|V2$" = 'skyblue'))
## Structure search estimating latent variables:2 ends here

## Final Example
##   :PROPERTIES:
##   :ID:       8F9AF2D0-3530-461B-8322-A183E15D07DA
##   :END:

## [[id:8F9AF2D0-3530-461B-8322-A183E15D07DA][Final Example:1]]
suppressWarnings(load("final_model_withlvp_withvp.RData", verbose = T))
datalist_com = datalist
figtrue_com = igraph::graph_from_data_frame(
			  datalist_com$coef %>%
			  filter(input %in% c(paste0("P", 1:10),
					      paste0("Up", 1:10),
					      paste0("V", 45:50),
					      paste0("V", 1:5),
					      "U1.out", "U2.out", 'Z'),
				 output %in% c("Z","U1.out", "U2.out",
					       paste0("V", 40:50),
					      paste0("V", 1:10)))
)
## Final Example:1 ends here

## [[id:8F9AF2D0-3530-461B-8322-A183E15D07DA][Final Example:2]]
plotIgraph(figtrue_com, layout = 'dot', nodesep = 0.00001,
		       fill = list("U\\d.*" = 'darksalmon',
				   "Z" = 'yellow',
				   "^V1$|V2$" = 'skyblue'),
	   saveToFile=F,
	   filename="complicated_model.pdf")
## Final Example:2 ends here

## Structure search with missing confounders.
##    :PROPERTIES:
##    :ID:       6F9DEBD2-700A-402B-8D1B-F3E2B60E2032
##    :END:

## [[id:6F9DEBD2-700A-402B-8D1B-F3E2B60E2032][Structure search with missing confounders.:1]]
train_com = datalist_com$data_noisy



trainlv_com = train_com %>% dplyr::select(-U1.out, -U2.out)

blacklistlv_com = rbind(data.frame(from = "Z",
				   to = colnames(trainlv_com)))

seed = 123
set.seed(seed)
res_missing_com = getEnsemble2(trainlv_com, blacklist = blacklistlv_com,
			    restart = 100, Nboot = 10,
			    prior = "vsp",
			    score = "bge",
			    algorithm = 'tabu',
			    parallel = TRUE
			    )
## Structure search with missing confounders.:1 ends here

## [[id:6F9DEBD2-700A-402B-8D1B-F3E2B60E2032][Structure search with missing confounders.:2]]
plot(res_missing_com,
     "Z",
     cutoff = 0.5,
     freqth = 0.5, 
     maxpath=2,
     nodesep=0.01,
     sep = 0.01,
     layout = 'dot',
     edge_labels = 'frequency',
     edgeweights=T,
     edgelabelsFilter = 0.5,
     edge_color=tribble(~input,~output,~color,"V1","Z","red","V2","Z","red",
                        "U1.out", "Z", "red",
                        "U2.out", "Z", "red"),
     fill = list("U1.out|U2.out" = 'darksalmon',
                 "Z" = 'yellow',
                 "^V1$|V2$" = 'skyblue'))
## Structure search with missing confounders.:2 ends here

## Reconstruction with estimated latent variables
##    :PROPERTIES:
##    :ID:       84079C63-F315-4062-BC0E-8000C025142D
##    :END:

## [[id:84079C63-F315-4062-BC0E-8000C025142D][Reconstruction with estimated latent variables:1]]
seed = 123
set.seed(seed)
graphics.off()
complicated_evo = latentDiscovery(
	res_missing_com,
	nItera=5,
	data = trainlv_com,
	"Z",
	seed=seed,
	workpath="latentDiscovery_com",
	method = "robustLinear",
	latent_iterations = 10, ## reduced for speed reasons
	truecoef = datalist_com$coef %>% filter(output=="Z"),
	truelatent=datalist_com$data %>% dplyr::select("U1.out","U2.out"),
	parallel = TRUE
    )
## Reconstruction with estimated latent variables:1 ends here

## [[id:84079C63-F315-4062-BC0E-8000C025142D][Reconstruction with estimated latent variables:2]]
plot(complicated_evo$details$final_ensemble, 
     "Z",
     cutoff = 0.5,
     freqth = 0.5, 
     maxpath=2,
     nodesep=0.01,
     sep = 0.01,
     layout = 'dot',
     edge_labels = 'frequency',
     edgeweights=T,
     edgelabelsFilter = 0.5,
     edge_color=tribble(~input,~output,~color,"V1","Z","red","V2","Z","red",
                        "U1.out", "Z", "red",
                        "U2.out", "Z", "red"),
     fill = list("U1.out|U2.out" = 'darksalmon',
                 "Z" = 'yellow',
                 "^V1$|V2$" = 'skyblue'))
## Reconstruction with estimated latent variables:2 ends here
