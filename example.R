## Load libraries
source("LatentConfounderBNlearnv2.R")
## load test data
suppressWarnings(load("final_model_nolvp_novp.RData", verbose = T))

figtrue = igraph::graph_from_data_frame(datalist$coef)

plotIgraph(figtrue, layout = 'neato', sep = 0.0001,
                       fill = list("U.*" = 'darksalmon',
                                   "Z" = 'yellow',
                                   "^V1$|V2$" = 'skyblue'))

train = datalist$data_noisy
library(bnlearn)

blacklist = rbind(data.frame(from = "Z", to = colnames(train)),
		      data.frame(from = colnames(train), to = "U1.out"),
		      data.frame(from = colnames(train), to = "U2.out")
                  )

library(doParallel)
cl <- makeCluster(10) ## for multi-threading
registerDoParallel(cl)

res_alldata = getEnsemble2(train, blacklist = blacklist,
                           restart = 100, Nboot = 10,
                           prior = "vsp",
                           score = "bge",
                           algorithm = 'hc',
                           parallel = TRUE
                           )

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

stopCluster(cl)

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

plotIgraph(figtrue_med, layout = 'dot', nodesep = 0.00001,
		       fill = list("U.*" = 'darksalmon',
				   "Z" = 'yellow',
				   "^V1$|V2$" = 'skyblue')
           )

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

plotIgraph(figtrue_com, layout = 'dot', nodesep = 0.00001,
		       fill = list("U\\d.*" = 'darksalmon',
				   "Z" = 'yellow',
				   "^V1$|V2$" = 'skyblue'),
	   saveToFile=F,
	   filename="complicated_model.pdf")

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
