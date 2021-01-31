

genData = function(N,
		     V = 20, ## number of nodes downstream of latent variables
		     Npx = 3, ## number of parents of drivers beside latent variables
		     latent_type = "Gaussian", ## type of latent variables
		     Np = 10, ## total number of parents                 
		     parentcoef = 0.3, ## max coef of parents
		     Npu = 0, ## number of parents of latent variables
		     uparentcoef = 0.3, ## max coef of parent of latent variables
		     latsd = 1, 
		     latdriver_coef1 = 1, ## coeficient value for latent variable on driver.
		     latdriver_coef2 = -1, 
		     latout_coef = 0.5, ## define the coeficient of latent variables on output
		     latother_max = 0.2, ## maximum coefficient of latent varible on other children
		     num_out_drivers = 2, ## number of drivers of output
		     coef_driver1 = 1,
		     coef_driver2 = 1, 
		     addnoisesd = 0.5, ## additional noise sd
		     outbias = 0, 
		     binary_endpoint = FALSE 
		     ){
      Ns =  N
      N = 2 * N
      allcoef = tribble( ~ input, ~ output, ~ coef)
      if(Npu == 0){
	  ## define input only latent variables
	  if(latent_type == 'Gaussian'){
	      U2.out = rnorm(N, 1, latsd)
	      U1.out = rnorm(N, 1, latsd)
	  }else if(latent_type == 'Beta'){
	      U2.out = rbeta(N, 10, 2)
	      U1.out = rbeta(N, 2,10)        
	  }else if(latent_type == "Lognormal"){
	      U1.out = rlnorm(N,1,0.5)
	      U2.out = rlnorm(N,3,0.1)
	  }else if(latent_type == "Multimodal"){
	      U1.out = map_dbl(sample(1:3, N, TRUE),
			       function(ii){
				   if(ii == 1)
				       rnorm(1, -1, 1)
				   else if(ii == 2)
				       rnorm(1, 1, 1)
				   else
				       rnorm(1, 3, 1)
			       })
	      U2.out = map_dbl(sample(1:2, N, TRUE),
			       function(ii){
				   if(ii == 1)
				       rnorm(1, 0, 1)
				   else if(ii == 2)
				       rnorm(1, 3, 1)
			       })
	      U1.out = scale(U1.out)[, 1]
	      U2.out = scale(U2.out)[, 1]
	  }
      }else{
	  ## create parents of latent variables
	  Npuu = Npu
	  Npu = 2 * Npu
	  Up = array(0, dim = c(N, Npu))
	  colnames(Up) = paste0("Up", 1:Npu)
	  for(pp in 1:Npu){
	      Up[, pp] = rnorm(N, runif(Npu, -2, 2), 1)
	  }
	  U1.out = rnorm(N, 0, addnoisesd)
	  U2.out = rnorm(N, 0, addnoisesd)
	  selparent = sample(1:Npu, Npu)
	  for(iipp in 1:Npuu){
	      selpp = selparent[iipp]
	      coefsel = runif(1, -uparentcoef, uparentcoef)
	      U1.out = U1.out + coefsel * Up[, selpp]
	      allcoef = rbind(allcoef,
			      tribble( ~ input, ~ output, ~ coef,
				      colnames(Up)[selpp], "U1.out", coefsel
				      ))
	      selpp = selparent[(Npuu) + iipp]
	      coefsel = runif(1, -uparentcoef, uparentcoef)
	      U2.out = U2.out + coefsel * Up[, selpp]
	      allcoef = rbind(allcoef,
			      tribble( ~ input, ~ output, ~ coef,
				      colnames(Up)[selpp], "U2.out", coefsel
				      ))
	  }
      }
      ## Generate parents
      if(Np > 0){
	    P = array(0, dim = c(N, Np))
	    colnames(P) = paste0("P", 1:Np)
	    for(pp in 1:Np){
		P[, pp] = rnorm(N, runif(Np, -2, 2), 1)
	    }
	}
	## Generate variables downstream of latent variables
	X = array(0, dim=c(N, V))
	for (ii in 1:V){
	    ## coefficient of latent variables on variable
	    if(ii <= num_out_drivers){
		## u1xcoef = rnorm(1, latdriver_coef1, 0.2)
		u1xcoef = latdriver_coef1
		##u2xcoef = rnorm(1, latdriver_coef2, 0.2)
		u2xcoef = latdriver_coef2
	    }else{
		u1xcoef = runif(1, -latother_max, latother_max)
		u2xcoef = runif(1, -latother_max, latother_max)
	    }
	    if(ii > num_out_drivers){
		if(runif(1,0, 1) <= 0.5){
		    X[, ii] =  u1xcoef * U1.out  + ## latent
			rnorm(nrow(X), 0, addnoisesd) ## additional unexplained drivers
		    allcoef = rbind(allcoef,
				    tribble( ~ input, ~ output, ~ coef,
					    "U1.out", paste0("V", ii), u1xcoef
					    ))
		}else{
		    X[, ii] = u2xcoef * U2.out + rnorm(nrow(X), 0, addnoisesd) ## additional unexplained drivers
		    allcoef = rbind(allcoef,
				    tribble( ~ input, ~ output, ~ coef,
					    "U2.out", paste0("V", ii), u2xcoef
					    ))
		}
	    }else{
		X[, ii] =  u1xcoef * U1.out  +
		    u2xcoef * U2.out + ## latent
		    rnorm(nrow(X), 0, addnoisesd) ## additional unexplained drivers
		allcoef = rbind(allcoef,
				tribble( ~ input, ~ output, ~ coef,
					"U1.out", paste0("V", ii), u1xcoef,
					"U2.out", paste0("V", ii), u2xcoef                                 
					))
	    }
	    ## add parents to each variable
	    ## decide which parents
	    if(Np > 0){
		parents = sample(1:Np, min(Np, Npx))
		for(pp in parents){
		    pxcoef = runif(1, -parentcoef, parentcoef)
		    X[, ii] = X[, ii] + pxcoef * P[, pp]
		    allcoef = rbind(allcoef,
				    tribble( ~ input, ~ output, ~ coef,
					    paste0("P",pp ), paste0("V", ii), pxcoef
					    )
				    )
		}
	    }
	}
	## define output
	Z = latout_coef * U1.out + latout_coef * U2.out
	allcoef = rbind(allcoef,
			tribble( ~ input, ~ output, ~ coef,
				"U1.out", "Z", latout_coef,
				"U2.out", "Z", latout_coef)
			)
      for (ii in 1:num_out_drivers){
	  if(ii == 1)
	      val = coef_driver1
	  else
	      val = coef_driver2
	    allcoef = rbind(allcoef,
			    tribble( ~ input, ~ output, ~ coef,
				    paste0("V", ii), "Z",val
				    )
			    )
	    Z = Z + X[, ii] * val
	}
      Z = Z + rnorm(length(Z), outbias, addnoisesd)
      if(binary_endpoint){
	  Zlogi = 1 / (1 + exp(-Z))
	  Z = rbinom(length(Z), 1, Zlogi)
      }
	##Z = Z + rnorm(length(Z), 0, addnoisesd) ## add additional unexplained variation
	colnames(X) = paste0("V", 1:ncol(X))
	if(Np > 0){
	    colnames(P) = paste0("P", 1:ncol(P))
	    dframe = as.data.frame(cbind(X, P, U1.out, U2.out, Z))
	}else{
	    dframe = as.data.frame(cbind(X, U1.out, U2.out, Z))
	}
      if(Npu > 0){
	  dframe = cbind(dframe, Up)
      }
      attributes(dframe)$allcoef = allcoef
      return(dframe)
  }






dframe = genData(
    N = 1000
    V = 25, 
    num_out_drivers = 2
    Npx = 3
    latent_type = "Gaussian"
    Npu = 0, 
    uparentcoef=0.3, 
    Np = 10
    parentcoef = 0.3
    latout_coef = 0.5
    latdriver_coef1 = 1
    latdriver_coef2 = 1
    latother_max = allconval$latother_max,
    addnoisesd = allconval$addnoisesd, 
    coef_driver1 = 1,
    coef_driver2 = 1,
    outbias = 0
)


## see load("final_model_nolvp_withvp.RData",verbose=T) for coeficcient values

