## Simple example

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Simple%20example][Simple example:1]]
##setwd("~/Documents/projects//latentconfounder")
##setwd("~/Documents/gitRepos/latentconfounder")
  source("LatentConfounderBNlearn.R")

    ##source("/users/fred/Documents/gitRepos/latentconfounder/LatentConfounder.R")
  ##  devtools::load_all("~/Documents/gitRepos/gnsutils/gnsutils")
  ##devtools::load_all("/Users/fred/Documents/projects/gnsutils/gnsutils")
  load("final_model_nolvp_novp.RData", verbose = T)
figtrue = igraph::graph_from_data_frame(datalist$coef)

igtruelv = igraph::graph_from_data_frame(
                       filter(
                           datalist$coef,
                           !(input == "U1.out" & output == "V1"),
                           !(input == "U2.out" & output == "V1"),
                           !(input == "U1.out" & output == "V2"),
                           !(input == "U2.out" & output == "V2")))

igtruered = igraph::graph_from_data_frame(
                        filter(
                            datalist$coef,
                            output == 'Z'
                        )
                    )

truecoef = filter(
    datalist$coef,
    output == 'Z'
)
## Simple example:1 ends here

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Simple%20example][Simple example:2]]
plotIgraph(figtrue, layout = 'neato', sep = 0.0001,
                       fill = list("U.*" = 'darksalmon',
                                   "Z" = 'yellow',
                                   "^V1$|V2$" = 'skyblue'))
## Simple example:2 ends here

## With all data
## Here learn the structure using bnlearn and bootstrapping

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*With%20all%20data][With all data:1]]
train = datalist$data_noisy
library(bnlearn)



    blacklist = rbind(data.frame(from = "Z", to = colnames(train)),
		      data.frame(from = colnames(train), to = "U1.out"),
		      data.frame(from = colnames(train), to = "U2.out")
		      )

##REFSfs:::registerDoSGE() ## this is how you would do it in GNS to use nodes, Should remove this line before we share this publicly
library(doParallel)
  cl <- makeCluster(5) ## for multi-threading
  registerDoParallel(cl)

library(pracma)
  tic()
  res_alldata = getEnsemble2(train, blacklist = blacklist,
			      restart = 100, Nboot = 50,
			      prior = "vsp",
			      score = "bge",
			      algorithm = 'hc',
			      parallel = TRUE
			      )
  toc() ## about 50 seconds

  save(res_alldata, file = "~/latentconfounderotf/latent_Discovery/res_alldata.RData")

pdf("Paper/images/res_alldata.pdf", width = 5, height = 5)
##  plot(res_alldata, 'Z', layout = 'neato', sep = 0.01)

##  res_alldata = getEnsemble(train, blacklist = blacklist, restart = 200, Nboot = 20)



    ## method = hc
    ## dag = method(train, blacklist = blacklist, restart = 200)
## With all data:1 ends here

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*With%20all%20data][With all data:2]]
plot(res_alldata,"Z",0, cutoff = 0.5,maxpath=2,nodesep=0.01,sep = 0.01,
     layout = 'neato', edge_labels = 'frequency',
     edgeweights=T,
     edgelabelsFilter = 0.5,
     edge_color=tribble(~input,~output,~color,"V1","Z","red","V2","Z","red",
                        "U1.out", "Z", "red",
                        "U2.out", "Z", "red"),
     fill = list("U.*" = 'darksalmon',
                 "Z" = 'yellow',
                 "^V1$|V2$" = 'skyblue'))
dev.off()
## With all data:2 ends here

## Missing U1 and U2
## THis is learning structure with U1 and U2 unobserved,

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Missing%20U1%20and%20U2][Missing U1 and U2:1]]
trainlv = train %>% dplyr::select(-U1.out, -U2.out)

blacklistlv = rbind(data.frame(from = "Z", to = colnames(trainlv)))



##REFSfs:::registerDoSGE()
tic()
res_missing = getEnsemble2(trainlv, blacklist = blacklistlv,
			  restart = 100, Nboot = 50,
			  prior = "vsp",
			  score = "bge",
			  algorithm = 'hc',
			  parallel = TRUE
			  )
toc() ## about 50 seconds

save(res_missing, file = "~/latentconfounderotf/latent_Discovery/res_missing.RData")

load("res_missing.RData", verbose = T)

test = latentDiscovery(
    res_missing,
    nItera=100,
    data = trainlv,
    "Z",
    workpath="aetestnew2",
    freqCutoff = 0.01,
    maxpath = 1,
    alpha = 0.05,
    scale. = TRUE,
    method = "autoencoder",
    latent_iterations = 100,
    truecoef = datalist$coef %>% filter(output=="Z"),
    truelatent=datalist$data %>% dplyr::select("U1.out","U2.out"),
    include_downstream = TRUE,
    include_output = TRUE,
    multiple_comparison_correction = T,
    debug = FALSE,
    parallel = TRUE,
    drRate = 0.2,
    batch_size = 16,
    optimizer = "RMSprop",
    activation = 'sigmoid',
    activation_coding = "sigmoid",
    activation_output = "linear",
    metrics = 'mse'
)

## Missing U1 and U2:1 ends here

## Run Latent discovery with repetitions

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Run%20Latent%20discovery%20with%20repetitions][Run Latent discovery with repetitions:1]]
Nreps = 30
niter = 10
allseeds = sample(seq(1, 10000,1),size = Nreps)
res = NULL
for(iis in 1:Nreps){
    message("Repeat ",iis, " of ", Nreps)
    seed = allseeds[iis]
    set.seed(seed)
    simple_evo_repeat = latentDiscovery(
        res_missing,
        nItera=niter,
        data = trainlv,
        "Z",
	seed=seed,
        workpath="latentDiscoveryRep",
        freqCutoff = 0.01,
        maxpath = 1,
        alpha = 0.05,
        scale. = TRUE,
        method = "linear",
        latent_iterations = 100,
        truecoef = datalist$coef %>% filter(output=="Z"),
        truelatent=datalist$data %>% dplyr::select("U1.out","U2.out"),
        include_downstream = TRUE,
        multiple_comparison_correction = T,
        debug = F,include_output=FALSE,
        parallel = TRUE
    )
    cres=simple_evo_repeat$details$Diagnostics %>%
         mutate(Repeat=iis)
    if(iis==1)
       res=cres
    else
       res=rbind(res,cres)
    saveRDS(res,file="/home/fred/Documents/gitRepos/latentconfounder/latentDiscoveryRep/allres_rep.RDS")
}
## Run Latent discovery with repetitions:1 ends here

## Images from evolutions

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Images%20from%20evolutions][Images from evolutions:1]]
library(gnsutils)
devtools::load_all("/Users/fred/Documents/projects/gnsutils/gnsutils/")
source("~/Documents/projects//latentconfounder/bnlearnLatent.R")
setwd("/Users/fred/Documents/projects/latentconfounder/")

load("latent_Discovery/res_alldata.RData", verbose = T)
load("latent_Discovery/res_missing.RData", verbose = T)
latvar = readRDS("latVars.RDS")
load("final_model_nolvp_novp.RData", verbose = T)

figtrue = igraph::graph_from_data_frame(datalist$coef)
## Images from evolutions:1 ends here

## True Reconstruction

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*True%20Reconstruction][True Reconstruction:1]]
ppPlotIgraph(figtrue, layout = 'neato', sep = 0.001,
                       fill = list("U.*" = 'red',
                                   "Z" = 'green',
                                   "^V1$|^V2$" = 'skyblue',
                                   "V.*" = "gray"
                                   ),
             node_labels = list("U1.out" = 'U1',
                        "U2.out" = "U2"),
             edgelabels = TRUE,
             edgelabelsFilter=0.49,
             fontsize=14,
             filename="Paper/images/true_network.pdf",
             width=10,
             height=10,
             start = 3421,
             edgelabelsFontSize=20,
             saveToFile=TRUE)

ppPlotIgraph(figtrue, layout = 'neato',nodesep = 0.001, sep = 0.001,
                       fill = list("U.*" = 'red',
                                   "Z" = 'green',
                                   "^V1$|^V2$" = 'skyblue',
                                   "V.*" = "gray"
                                   ),
             edgelabels=TRUE,
             start = 3421,
             node_labels = list("U1.out" = 'U1',
                                "U2.out" = "U2"),
             edgelabelsFilter=0.49,edgelabelsFontSize=14,
             fontsize=10)
## True Reconstruction:1 ends here

## Reconstruction with all data

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Reconstruction%20with%20all%20data][Reconstruction with all data:1]]
plot(res_alldata, "Z", layout = 'neato',
     cutoff = 0.4,
     sep = 0.00001,
     pack = TRUE,
     maxpath = 1,
     fill = list("U.*" = 'red',
                 "Z" = 'green',
                 "^V1$|^V2$" = 'skyblue',
                 "V.*" = "gray"
                 ),
     node_labels = list("U1.out" = 'U1',
                        "U2.out" = "U2"),
     edge_labels = "coefficients",
     labels_regex = "Z|^V1$|^V2$",
     label_pad = 0,
     start = 1235,
     edgeweights=F,
     lwdmin=0,
     edgelabelsFontSize=20,
     title_size = 2,
     fontsize=12,
     filename="Paper/images/estimated_network_alldata.pdf",
     saveToFile=TRUE,
     width=10,
     height=10)


plot(res_alldata, "Z", layout = 'neato',
     cutoff = 0.4,
     sep = 0.00001,
     fill = list("U.*" = 'red',
                 "Z" = 'green',
                 "^V1$|^V2$" = 'skyblue',
                 "V.*" = "gray"
                 ),
     edge_labels = "coefficients",
     labels_regex = "Z|^V1$|^V2$",
     node_labels = list("U1.out" = 'U1',
                        "U2.out" = "U2"),
     label_pad = 0.1,
     start = 1324,
     edgeweights=F,
     lwdmin=0,
     edgelabelsFontSize=15,
     title="Estimated Network with All Data",
     title_size = 2,
     fontsize=12)
## Reconstruction with all data:1 ends here

## Reconstruction with missing variables

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Reconstruction%20with%20missing%20variables][Reconstruction with missing variables:1]]
plot(res_missing, "Z",
     layout = 'dot',
     cutoff = 0.4,
     maxpath = 2,
     freqth = 0.1,
     nodesep = 0.01,
     pack = FALSE,
     fill = list("U.*" = 'red',
                 "Z" = 'green',
                 "^V1$|^V2$" = 'skyblue',
                 "V.*" = "gray"
                 ),
     edge_labels="coef",
     edge_color = tribble(
         ~ inp, ~ out, ~ color,
         "V1", "Z", 'red',
         "V2", "Z", "red"
     ),
     labels_regex="Z|^V1$|^V2$",
     start = 1334324,
     edgeweights=F,
     lwdmin=0,
     edgelabelsFontSize=32,
     fontsize=12,
     filename ="Paper/images/estimated_network_missingdata.pdf",
     saveToFile=TRUE,
     width=10,
     height=10)


plot(res_missing, output = "Z",
     layout = 'dot',
     cutoff = 0.4,
     maxpath = 2,
     freqth = 0.1,
     nodesep = 0.01,
     sep = 0.001,
     fill = list("U.*" = 'red',
                 "Z" = 'green',
                 "^V1$|^V2$" = 'skyblue',
                 "V.*" = "gray"
                 ),
     edge_labels="coef",
     labels_regex="Z",
     label_pad = 0,
     start = 13324,
     edgeweight=T,
     lwdmin=0,
     edgelabelsFontSize=15,
     fontsize=14
     )
## Reconstruction with missing variables:1 ends here

## Reconstruction with last iteration

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Reconstruction%20with%20last%20iteration][Reconstruction with last iteration:1]]
res_lastitera = latvar$details$Diagnostics$Ensemble[[11]]

plot(res_lastitera, "Z",
     layout = 'dot',
     cutoff = 0.4,
     freqth = 0.1,
     maxpath = 2,
     nodesep = 0.01,
     pack = F,
     sep = 0.01,
     fill = list("LV.*" = 'red',
                 "Z" = 'green',
                 "^V1$|^V2$" = 'skyblue',
                 "V.*" = "gray"
                 ),
     edge_labels="coef",
     edge_color = tribble(
         ~ inp, ~ out, ~ color,
         "V1", "Z", 'red',
         "V2", "Z", "red",
         "LV.*", "Z", "red"
     ),
     labels_regex="Z|^V1$|^V2$",
     label_pad = 0.4,
     start = 133434,
     edgeweights=F,
     lwdmin=0,
     edgelabelsFontSize=30,
     fontsize=12,
     filename ="Paper/images/estimated_network_infered.pdf",
     saveToFile=TRUE,
     width=10,
     height=8)


plot(res_lastitera, "Z",
     layout = 'dot',
     cutoff = 0.4,
     freqth = 0.2,
     maxpath = 2,
     nodesep = 0.01,
     sep = 0.001,
     fill = list("LV.*" = 'red',
                 "Z" = 'green',
                 "^V1$|^V2$" = 'skyblue',
                 "V.*" = "gray"
                 ),
     edge_labels="coef",
     labels_regex="Z|^V1$|^V2$",
     label_pad = 0.01,
     pack = TRUE,
     start = 134,
     edgeweights = T,
     lwdmin=0,
     edgelabelsFontSize=15,
     fontsize=14
     )
## Reconstruction with last iteration:1 ends here

## R2 in Latent variables estimation

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*R2%20in%20Latent%20variables%20estimation][R2 in Latent variables estimation:1]]
library("latex2exp")

label_latex <- function(vs, am) {
    browser()
    tst = map(vs$var, TeX)
    names(tst) = vs$var
    tst
}

## mylabels = c(
##     "R2a_U1.out" = ("$R^2$ for $U_1$"),
##     "R2a_U2.out" = ("$R^2$ for $U_2$")
## )

mylabels = c(
     "R2a_U1.out" = TeX("$U_1$"),
    "R2a_U2.out" = TeX("$U_2$")
)

tst = latvar$details$Diagnostics %>%
    dplyr::select(Iteration,
                  R2a_U1.out,
                  R2a_U2.out) %>%
    gather(var, val, -Iteration) %>%
    mutate(var = factor(var))

## levels(tst$var) = c(
##     TeX("$U_1$"),
##     TeX("$U_2$")
## )

if(0){
ggp = tst %>%
    ggplot(aes(x = Iteration, y = val)) +
    geom_point()  + geom_smooth(span = 0.8) +
    facet_wrap( ~ var, nrow = 1, labeller = label_parsed) +
    ggpubr::theme_pubclean() +    ylab(TeX("R^2"))  +
##    ggtitle(TeX("Latent Variables Prediction $R^2$")) +
    theme(strip.text=element_text(size = rel(1.5)),
          axis.title=element_text(size = rel(1.5)),
          plot.title = element_text(size = rel(1.5), hjust=0.5),
          )
}

ggp = tst %>%
    ggplot(aes(x = Iteration, y = val, colour = var, fill = var)) +
    geom_point()  + geom_smooth(span = 0.8, alpha = 0.1) +
    ggpubr::theme_pubclean() +    ylab(TeX("R^2"))  +
    ##ggtitle(TeX("Latent Variables Prediction $R^2$")) +
    scale_colour_manual(values = c("skyblue", "orange"),
                        labels = expression(U[1], U[2])) +
    scale_fill_manual(values = c("skyblue", "orange"),
                      labels = expression(U[1], U[2])) +
    theme(strip.text=element_text(size = rel(1.5)),
          legend.text=element_text(size = rel(1.1)),
          legend.position = c(0.5, 0.2),
          legend.direction ="horizontal",
          axis.title=element_text(size = rel(1)),
          plot.title = element_text(size = rel(1.5), hjust=0.5),
          legend.title=element_blank( )
          )




pdf("Paper/images/fig_paper_r2lat.pdf", width = 5, height = 5)
print(ggp)
dev.off()

print(ggp)
## R2 in Latent variables estimation:1 ends here

## Error rmse

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Error%20rmse][Error rmse:1]]
library("latex2exp")


mylabels = c(
    "R2a_U1.out" = ("$R^2$ for $U_1$"),
    "R2a_U2.out" = ("$R^2$ for $U_2$")
)

errall = latvar$details$Diagnostics %>%
    dplyr::select(Repeat, Iteration, Coef) %>%
    unnest(Coef)

err = errall %>%
    filter(True_coef >= 0.02)
	  ##            input != colnames(true)) %>%
	  ##     dplyr::select(-coef, -Estimate) %>%
	  ##     mutate(input=otfname2name(varmap,input))

err = err %>%
    mutate(Error = abs(error)) %>%
    dplyr::select(
               Iteration, Variable = input, Error
           )

rmsedf = latvar$details$Diagnostics %>%
    dplyr::select(Iteration, Error_rmse) %>%
    mutate(Variable = "RMSE") %>%
    dplyr::select(Iteration, Variable, Error = Error_rmse)

allerror = rbind(err, rmsedf)

## get performance of model with all data
score_alldata = getScores(ens=res_alldata,output='Z',truecoef=datalist$coef %>% filter(output=="Z"),lat_estimates=NULL)
error_alldata = score_alldata$Coef[[1]] %>% filter(input %in% c("V1", "V2"))
rmse_alldata = score_alldata$Error_rmse
allerror_alldata = rbind(
    error_alldata %>%
    mutate(error = abs(error)) %>%
    dplyr::select(Variable = input, Error = error),
    tribble(
        ~ Variable, ~ Error,
        'RMSE', rmse_alldata
    )
)

## plot
ggp = allerror %>%
    ggplot(aes(x = Iteration, y = Error, colour = Variable, fill = Variable)) +
    geom_point() +
    geom_line() +
    geom_hline(data = allerror_alldata,
               aes(yintercept = Error, colour = Variable), linetype = 'dashed') +
    ## annotate("text", x = 0, y = max(allerror_alldata$Error), label = "Obs. All Variables",
    ##          vjust = -1) +
    ylab("Error in Coefficients") +
    ggpubr::theme_pubclean() +
    theme(strip.text=element_text(size = rel(1)),
          axis.title=element_text(size = rel(1)),
          legend.position = c(0.5, 0.7),
          legend.text = element_text(size = rel(1)),
          legend.direction='horizontal',
          plot.title = element_text(size = rel(1.5), hjust=0.5))

## ggpe = err %>%
## 	      ggplot(aes(x = Iteration, y = abs(error), colour = input)) +
## 	      geom_point() +
## 	      geom_line() +
## 	      ##ylim(0, max(err$error)) +
##     ggtitle("Error In Drivers of Outcome") + ##+ geom_smooth()+
##     ylab("Absolute Error") +
##     ggpubr::theme_pubclean() +
##     theme(strip.text=element_text(size = rel(1.5)),
##           axis.title=element_text(size = rel(2)),
##           plot.title = element_text(size = rel(1.5), hjust=0.5))

## ggpet = latvar$details$Diagnostics %>%
## 	      ggplot(aes(x = Iteration, y = Error_rmse)) +
## 	      geom_point() +
## 	      geom_line() +
##     ggtitle("Total Error") +
##     ylab("RMSE") +
##     ggpubr::theme_pubclean() +
##     theme(strip.text=element_text(size = rel(1.5)),
##           axis.title=element_text(size = rel(2)),
##           plot.title = element_text(size = rel(1.5), hjust=0.5))


## ggpall = cowplot::plot_grid(ggpe, ggpet, ncol = 2)


pdf("Paper/images/fig_paper_errors.pdf", width = 5, height = 5)
print(ggp)
dev.off()

print(ggp)
## Error rmse:1 ends here

## Images with repetition

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Images%20with%20repetition][Images with repetition:1]]
latvars_resp = readRDS("/Users/fred/Documents/projects/latentconfounder/Paper/images/allres_rep.RDS")
## Images with repetition:1 ends here

## R2 in Latent variables estimation

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*R2%20in%20Latent%20variables%20estimation][R2 in Latent variables estimation:1]]
library("latex2exp")

## mylabels = c(
##     "R2a_U1.out" = ("$R^2$ for $U_1$"),
##     "R2a_U2.out" = ("$R^2$ for $U_2$")
## )

mylabels = c(
     "R2a_U1.out" = TeX("$U_1$"),
    "R2a_U2.out" = TeX("$U_2$")
)

tstrep = latvars_resp %>%
    dplyr::select(Iteration,
                  R2a_U1.out,
                  R2a_U2.out) %>%
    gather(var, val, -Iteration) %>%
    mutate(var = factor(var))

## levels(tst$var) = c(
##     TeX("$U_1$"),
##     TeX("$U_2$")
## )


ggp = tstrep %>%
    ggplot(aes(x = Iteration, y = val, colour = var, fill = var)) +
    ggbeeswarm::geom_beeswarm()  + geom_smooth(span = 0.9, alpha = 0.1) +
    ggpubr::theme_pubclean() +    ylab(TeX("R^2"))  +
    ##ggtitle(TeX("Latent Variables Prediction $R^2$")) +
    scale_colour_manual(values = c("skyblue", "orange"),
                        labels = expression(U[1], U[2])) +
    scale_fill_manual(values = c("skyblue", "orange"),
                      labels = expression(U[1], U[2])) +
    theme(strip.text=element_text(size = rel(1.5)),
          legend.text=element_text(size = rel(1.1)),
          legend.position = c(0.5, 0.2),
          legend.direction ="horizontal",
          axis.title=element_text(size = rel(1)),
          plot.title = element_text(size = rel(1.5), hjust=0.5),
          legend.title=element_blank( )
          )




pdf("Paper/images/fig_paper_r2lat_rep.pdf", width = 5, height = 5)
print(ggp)
dev.off()

print(ggp)
## R2 in Latent variables estimation:1 ends here

## Error rmse

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Error%20rmse][Error rmse:1]]
library("latex2exp")


mylabels = c(
    "R2a_U1.out" = ("$R^2$ for $U_1$"),
    "R2a_U2.out" = ("$R^2$ for $U_2$")
)

errallresp = latvars_resp %>%
    dplyr::select(Repeat, Iteration, Coef) %>%
    unnest(Coef)

errresp = errallresp %>%
    filter(True_coef >= 0.02)
	  ##            input != colnames(true)) %>%
	  ##     dplyr::select(-coef, -Estimate) %>%
	  ##     mutate(input=otfname2name(varmap,input))

errresp = errresp %>%
    mutate(Error = abs(error)) %>%
    dplyr::select(
               Iteration, Variable = input, Error
           )

rmsedfresp = latvars_resp %>%
    dplyr::select(Iteration, Error_rmse) %>%
    mutate(Variable = "RMSE") %>%
    dplyr::select(Iteration, Variable, Error = Error_rmse)

allerrorresp = rbind(errresp, rmsedfresp)

## get performance of model with all data
score_alldata = getScores(ens=res_alldata,output='Z',truecoef=datalist$coef %>% filter(output=="Z"),lat_estimates=NULL)
error_alldata = score_alldata$Coef[[1]] %>% filter(input %in% c("V1", "V2"))
rmse_alldata = score_alldata$Error_rmse
allerror_alldata = rbind(
    error_alldata %>%
    mutate(error = abs(error)) %>%
    dplyr::select(Variable = input, Error = error),
    tribble(
        ~ Variable, ~ Error,
        'RMSE', rmse_alldata
    )
)

## plot
(ggp = allerrorresp %>%
    ggplot(aes(x = Iteration, y = Error, colour = Variable, fill = Variable)) +
    ggbeeswarm::geom_beeswarm() +
    geom_smooth(se = TRUE) +
    geom_hline(data = allerror_alldata,
               aes(yintercept = Error, colour = Variable), linetype = 'dashed') +
    ## annotate("text", x = 0, y = max(allerror_alldata$Error), label = "Obs. All Variables",
    ##          vjust = -1) +
    ylab("Error in Coefficients") +
    ggpubr::theme_pubclean() +
    theme(strip.text=element_text(size = rel(1)),
          axis.title=element_text(size = rel(1)),
          legend.position = c(0.5, 0.7),
          legend.text = element_text(size = rel(1)),
          legend.direction='horizontal',
          plot.title = element_text(size = rel(1.5), hjust=0.5))
)
## ggpe = err %>%
## 	      ggplot(aes(x = Iteration, y = abs(error), colour = input)) +
## 	      geom_point() +
## 	      geom_line() +
## 	      ##ylim(0, max(err$error)) +
##     ggtitle("Error In Drivers of Outcome") + ##+ geom_smooth()+
##     ylab("Absolute Error") +
##     ggpubr::theme_pubclean() +
##     theme(strip.text=element_text(size = rel(1.5)),
##           axis.title=element_text(size = rel(2)),
##           plot.title = element_text(size = rel(1.5), hjust=0.5))

## ggpet = latvar$details$Diagnostics %>%
## 	      ggplot(aes(x = Iteration, y = Error_rmse)) +
## 	      geom_point() +
## 	      geom_line() +
##     ggtitle("Total Error") +
##     ylab("RMSE") +
##     ggpubr::theme_pubclean() +
##     theme(strip.text=element_text(size = rel(1.5)),
##           axis.title=element_text(size = rel(2)),
##           plot.title = element_text(size = rel(1.5), hjust=0.5))


## ggpall = cowplot::plot_grid(ggpe, ggpet, ncol = 2)


pdf("Paper/images/fig_paper_errors_rep.pdf", width = 5, height = 5)
print(ggp)
dev.off()

print(ggp)
## Error rmse:1 ends here

## other stuff

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*other%20stuff][other stuff:1]]
plot.bnlearn_ens(res_missing,"Z",0, cutoff = 0.5,maxpath=2,nodesep=0.01,sep = 0.01,
     layout = 'neato', edgelabels = T,
     edgeweights=T,
     edgelabelsFilter = 0.5,
     edge_color=tribble(~inp,~out,~color,"V1","Z","red","V2","Z","red"
			),
     fill = list("U.*" = 'darksalmon',
		 "Z" = 'yellow',
		 "^V1$|V2$" = 'skyblue'))
## other stuff:1 ends here

## Medium Example

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Medium%20Example][Medium Example:1]]
setwd("~/Documents/gitRepos/latentconfounder")
source("bnlearnLatent.R")


load("final_model_nolvp_withvp.RData", verbose = T)

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

igtruelv_med = igraph::graph_from_data_frame(
		       filter(
			   datalist_med$coef,
			   !(input == "U1.out" & output == "V1"),
			   !(input == "U2.out" & output == "V1"),
			   !(input == "U1.out" & output == "V2"),
			   !(input == "U2.out" & output == "V2")))

igtruered_med = igraph::graph_from_data_frame(
			filter(
			    datalist_med$coef,
			    output == 'Z'
			)
		    )

truecoef_med = filter(
    datalist_med$coef,
    output == 'Z'
)
## Medium Example:1 ends here

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Medium%20Example][Medium Example:2]]
plotIgraph(figtrue_med, layout = 'dot', nodesep = 0.00001,
		       fill = list("U.*" = 'darksalmon',
				   "Z" = 'yellow',
				   "^V1$|V2$" = 'skyblue'),
	   saveToFile=F,
	   filename="medium_model.pdf")
## Medium Example:2 ends here

## WIth all data

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*WIth%20all%20data][WIth all data:1]]
train_med = datalist_med$data_noisy

blacklist = rbind(data.frame(from = "Z", to = colnames(train_med)),
		  data.frame(from = colnames(train_med), to = "U1.out"),
		  data.frame(from = colnames(train_med), to = "U2.out")
		  )

##blacklist_med = rbind(data.frame(from = "Z", to = colnames(train_med)))

##REFSfs:::registerDoSGE()

res_alldata_med = getEnsemble2(train_med, blacklist = blacklist,
			   restart = 100, Nboot = 50,
			   prior = "vsp",
			   score = "bge",
			   algorithm = 'hc',
			   parallel = T
			   )

save(res_alldata_med, file = "/home/fred/Documents/gitRepos/latentconfounder/latent_Discovery/res_alldata_med.RData")

plot(res_alldata_med, 'Z', layout = 'dot', sep = 0.01, edge_labels = "coef")


##res_alldata = getEnsemble(train, blacklist = blacklist, restart = 200, Nboot = 20)



    ## method = hc
    ## dag = method(train, blacklist = blacklist, restart = 200)
## WIth all data:1 ends here

## Missing latent variables

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Missing%20latent%20variables][Missing latent variables:1]]
trainlv_med = train_med %>% dplyr::select(-U1.out, -U2.out)

blacklistlv_med = rbind(data.frame(from = "Z", to = colnames(trainlv_med)))



library(tictoc)


REFSfs:::registerDoSGE()

tic()
res_missing_med = getEnsemble2(trainlv_med, blacklist = blacklistlv_med,
			    restart = 100, Nboot = 50,
			    prior = "vsp",
			    score = "bge",
			    algorithm = 'hc',
			    parallel = TRUE
			    )
toc() ## about

save(res_missing_med, file = "/home/fred/Documents/gitRepos/latentconfounder/latent_Discovery/res_missing_med.RData")
## Missing latent variables:1 ends here

## Causal Discovery

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Causal%20Discovery][Causal Discovery:1]]
niter = 10
seed = 123
set.seed(seed)
medium_evo = latentDiscovery(
	res_missing_med,
	nItera=niter,
	data = trainlv_med,
	"Z",
	seed=seed,
	workpath="latentDiscovery_med",
	freqCutoff = 0.01,
	maxpath = 1,
	alpha = 0.05,
	scale. = TRUE,
	method = "robustLinear",
	latent_iterations = 100,
	truecoef = datalist_med$coef %>% filter(output=="Z"),
	truelatent=datalist_med$data %>% dplyr::select("U1.out","U2.out"),
	include_downstream = TRUE,
	multiple_comparison_correction = TRUE,
	debug = F,
	parallel = TRUE
    )


save(medium_evo, file = "/home/fred/Documents/gitRepos/latentconfounder/latent_Discovery/evo_discovery_med.RData")
## Causal Discovery:1 ends here

## Complicated Example

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Complicated%20Example][Complicated Example:1]]
setwd("~/Documents/gitRepos/latentconfounder")
source("bnlearnLatent.R")


load("final_model_withlvp_withvp.RData", verbose = T)
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

igtruelv_com = igraph::graph_from_data_frame(
		       filter(
			   datalist_com$coef,
			   !(input == "U1.out" & output == "V1"),
			   !(input == "U2.out" & output == "V1"),
			   !(input == "U1.out" & output == "V2"),
			   !(input == "U2.out" & output == "V2")))

igtruered_com = igraph::graph_from_data_frame(
			filter(
			    datalist_com$coef,
			    output == 'Z'
			)
		    )

truecoef_com = filter(
    datalist_com$coef,
    output == 'Z'
)
## Complicated Example:1 ends here

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Complicated%20Example][Complicated Example:2]]
plotIgraph(figtrue_com, layout = 'dot', nodesep = 0.00001,
		       fill = list("U\\d.*" = 'darksalmon',
				   "Z" = 'yellow',
				   "^V1$|V2$" = 'skyblue'),
	   saveToFile=F,
	   filename="complicated_model.pdf")
## Complicated Example:2 ends here

## WIth all data

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*WIth%20all%20data][WIth all data:1]]
train_com = datalist_com$data_noisy


blacklist_com = rbind(data.frame(from = "Z", to = colnames(train_com)))


##blacklist_med = rbind(data.frame(from = "Z", to = colnames(train_med)))

REFSfs:::registerDoSGE()
res_alldata_com = getEnsemble2(train_com, blacklist = blacklist_com,
			   restart = 100, Nboot = 50,
			   prior = "vsp",
			   score = "bge",
			   algorithm = 'hc',
			   parallel = T
			   )

save(res_alldata_com, file = "/home/fred/Documents/gitRepos/latentconfounder/latent_Discovery/res_alldata_com.RData")

plot(res_alldata_com, 'Z', layout = 'neato', sep = 0.01, edge_labels = "coef")


##res_alldata = getEnsemble(train, blacklist = blacklist, restart = 200, Nboot = 20)



    ## method = hc
    ## dag = method(train, blacklist = blacklist, restart = 200)
## WIth all data:1 ends here

## Missing latent variables

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Missing%20latent%20variables][Missing latent variables:1]]
trainlv_com = train_com %>% dplyr::select(-U1.out, -U2.out)

blacklistlv_com = rbind(data.frame(from = "Z",
				   to = colnames(trainlv_com)))



library(tictoc)


REFSfs:::registerDoSGE()

tic()
res_missing_com = getEnsemble2(trainlv_com, blacklist = blacklistlv_com,
			    restart = 100, Nboot = 50,
			    prior = "vsp",
			    score = "bge",
			    algorithm = 'hc',
			    parallel = TRUE
			    )
toc() ## about

save(res_missing_com, file = "/home/fred/Documents/gitRepos/latentconfounder/latent_Discovery/res_missing_com.RData")
## Missing latent variables:1 ends here

## Causal Discovery

## [[file:~/Documents/gitRepos/latentconfounder/LatentConfounderBNlearn.org::*Causal%20Discovery][Causal Discovery:1]]
niter = 10
seed = 123
set.seed(seed)

complicated_evo = latentDiscovery(
	res_missing_com,
	nItera=niter,
	data = trainlv_com,
	"Z",
	seed=seed,
	workpath="latentDiscovery_com",
	freqCutoff = 0.01,
	maxpath = 1,
	alpha = 0.05,
	scale. = TRUE,
	method = "robustLinear",
	latent_iterations = 100,
	truecoef = datalist_com$coef %>% filter(output=="Z"),
	truelatent=datalist_com$data %>% dplyr::select("U1.out","U2.out"),
	include_downstream = TRUE,
	multiple_comparison_correction = TRUE,
	debug = F,
	parallel = TRUE
    )


save(complicated_evo, file = "/home/fred/Documents/gitRepos/latentconfounder/latent_Discovery/evo_discovery_com.RData")
## Causal Discovery:1 ends here

