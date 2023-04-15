#' @include LatSpace.R


#' 
## The Variational Autoencoder implementation of latent space encoding
## Based on the example implementation found at:
## https://github.com/schmons/torch_R_examples/tree/main/

# @description A utility function to calculate a smooth fan-in and fan-out towards the latent space
# @param n_features number of "outward-facing" features to encode/decode
# @param n_latent_space number of latent space terms to encode/decode
# @param direction in/out, determines the order of operations
NetRes$set("private", "VAE_fan_inout", function(n_features, n_latent_space, direction, n_steps=3) {
  ##n_steps = max(floor(log(n_features) - log(n_latent_space)), 2)
  print('fan_inout')
  if (direction == 'in') {
    n_features_per_layer = floor(seq(sqrt(n_features), sqrt(n_latent_space), length.out=n_steps)^2)
    ##now ensure that the "edges" are equal to exactly the expected numbers because seq may not respect that    
    n_features_per_layer[1] = n_features
    n_features_per_layer[length(n_features_per_layer)] = n_latent_space
  } else {
    n_features_per_layer = floor(seq(sqrt(n_latent_space), sqrt(n_features), length.out=n_steps)^2)
    ##now ensure that the "edges" are equal to exactly the expected numbers because seq may not respect that
    n_features_per_layer[1] = n_latent_space
    n_features_per_layer[length(n_features_per_layer)] = n_features
  }
  print(paste('Direction of coding:', direction, 'layers:', n_features_per_layer, '\n', sep=' ', collapse=''))
  return(n_features_per_layer)
}
)

# @description
# The encoder constructor
# @param in_features 
NetRes$set("private", "encoder", function(n_features_per_layer) {
  print('encoder')
  encoder = nn_module(
    "encoder",
    initialize = function(n_features_per_layer){
      ##nn_sequential allows adding to "NULL" transparently, so set to NULL and build up from there
      self$modelFit <- NULL
      ##the last automatically-generated index *precedes* the latent layer
      for (layerSizeIndex in 1:(length(n_features_per_layer) - 2)) {
        if (layerSizeIndex == 1) { #first layer is linear to support continuous observables
          self$modelFit <- nn_sequential(
            self$modelFit,
            nn_linear(n_features_per_layer[layerSizeIndex],
                      n_features_per_layer[layerSizeIndex + 1]))
        } else { #the rest are relus
          self$modelFit <- nn_sequential(
            self$modelFit,
            nn_linear(n_features_per_layer[layerSizeIndex],
                      n_features_per_layer[layerSizeIndex + 1]),
            nn_relu())
        }
      }
      
      ##I think this is mu - CONFIRM
      self$linear1 <- nn_linear(n_features_per_layer[length(n_features_per_layer) - 1], 
                                n_features_per_layer[length(n_features_per_layer)])
      ##I think this is sigma - CONFIRM
      self$linear2 <- nn_linear(n_features_per_layer[length(n_features_per_layer) - 1], 
                                n_features_per_layer[length(n_features_per_layer)])
      
    },
    forward = function(x){
      hidden = self$modelFit(x)
      mu = self$linear1(hidden)
      log_var = self$linear2(hidden)
      
      list(mu, log_var)
    }
  )
  return(encoder(n_features_per_layer))
}
)

# @description
# The decoder module - note that there are more layers here than in the example this was extended from -
# the decoder now mirrors the encoder in geometry.  Is this a good thing for VAEs?
NetRes$set("private", "decoder", function(n_features_per_layer) {
  print('decoder')
  decoder = nn_module(
    "decoder",
    initialize = function(n_features_per_layer){
      ##nn_sequential allows adding to "NULL" transparently, so set to NULL and build up from there
      self$modelFit = NULL
      ##the last automatically-generated index *precedes* the latent layer
      for (layerSizeIndex in 1:(length(n_features_per_layer) - 1)) {
        if(layerSizeIndex == length(n_features_per_layer) - 1) { #last layer is linear to support the continuous observables
          self$modelFit <- nn_sequential(
            self$modelFit,
            nn_linear(n_features_per_layer[layerSizeIndex],
                      n_features_per_layer[layerSizeIndex + 1]))
        } else { #inner layers should be non-linear
          self$modelFit <- nn_sequential(
            self$modelFit,
            nn_linear(n_features_per_layer[layerSizeIndex],
                      n_features_per_layer[layerSizeIndex + 1]),
            nn_relu())
        }
      }
    },
    forward = function(x){
      self$modelFit(x)
    }
  )
  return(decoder(n_features_per_layer))
}
)

# @description
# The overall VAE constructor
# Define VAE model using encoder and decoder from above
NetRes$set("private", "vae_module", function(n_features, latent_dim, n_steps=3) {
  vae_template = nn_module(
    encoder_geometry = NULL,
    decoder_geometry = NULL,
    latent_dim = NULL,
    encoder = NULL,
    decoder = NULL,
    initialize = function(n_features, latent_dim, n_steps=3) {
      self$encoder_geometry = ns_env()$NetRes$private_methods$VAE_fan_inout(n_features, latent_dim, 'in', n_steps=n_steps)
      self$decoder_geometry = ns_env()$NetRes$private_methods$VAE_fan_inout(n_features, latent_dim, 'out', n_steps=n_steps)
      self$latent_dim = latent_dim
      self$encoder = ns_env()$NetRes$private_methods$encoder(self$encoder_geometry)
      self$decoder = ns_env()$NetRes$private_methods$decoder(self$decoder_geometry)
    },
    
    forward = function(x) {
      f <- self$encoder(x)
      mu <- f[[1]]
      log_var <- f[[2]]
      ##NEED TO CONSIDER THE MEANING OF 0.5 AND WHETHER IT CAN BE SET 'RATIONALLY' INSTEAD
      z <- mu + torch_exp(log_var$mul(0.5))*torch_randn(c(dim(x)[1], self$latent_dim))
      reconst_x <- self$decoder(z)
      
      list(reconst_x, mu, log_var)
    }
  )
  return(vae_template(n_features, latent_dim, n_steps))
}
)

NetRes$set("private", "trainVAE", function(vae, train_data, epochs = 1000) {
  ##prepare the data as an iterator
  dataset <- torch::dataset(
    
    name = "train_dataset",
    data = NULL,
    initialize = function(data) {
      self$data <- torch_tensor(data)
    },
    
    .getitem = function(index) {
      
      x <- self$data[index, ]
      
      x
    },
    
    .length = function() {
      self$data$size()[[1]]
    }
  )

  dl <- torch::dataloader(dataset(as.matrix(train_data)), batch_size = floor(nrow(train_data)/5), shuffle = TRUE, drop_last=TRUE)
  ## Optimize. Note that a scheduler and/or a different learning rate could improve performance
  optimizer <- optim_adam(vae$parameters, lr = 0.001)
  ## Optimization loop
  losses = c()
  for(epoch in 1:epochs) {
    l = 0
    
    coro::loop(for (b in dl) {  # loop over all minibatches for one epoch
      forward = vae(b)
      
      loss = torch::nn_l1_loss(reduction = "mean")
      mu = forward[[2]]
      log_var = forward[[3]]
      
      output <- loss(forward[[1]], b)
      l = l + as.numeric(output)
      
      optimizer$zero_grad()
      output$backward()
      optimizer$step()
    })
    print(paste('epoch:', epoch))
    print(paste('loss:', l))
    losses[epoch] = l
  }
  epoch = 1:epochs
  visreg::visreg(loess(losses ~ epoch, span=0.1), main='Loss by epoch')
  return(vae)
}
)

NetRes$set("private", "calculateLatVarsVAE", function(residuals, resPpca, algorithm.args=NULL) {
  ##check that the latent space variables are excluded
  if (is.null(algorithm.args) || is.null(algorithm.args$autoencoder.depth)) {
    n_steps = 3
  } else {
    n_steps = algorithm.args$autoencoder.depth
  }
  print(paste('Autoencoder depth set to', n_steps))
  vae = private$vae_module(n_features = ncol(residuals), latent_dim = resPpca$n, n_steps = n_steps)
  
  ##training occurs by side effect since torch is all implemented in C
  private$trainVAE(vae, residuals)
  ##encoding = vae(torch_tensor(as.matrix(residuals))) #check correct dimensionality  

  ##create an R6 class that will encapsulate the predictive ability using the common R syntax
  predictor = R6Class('VAEPredictor',
                         public=list(
                           fitted.vae = NULL,
                           data = NULL,
                           initialize = function(fitted.vae, data) {
                             self$fitted.vae = fitted.vae
                             self$data = data
                           },
                           predict = function(newdata = NULL) {
                             if (is.null(newdata)) {
                               ##res = self$fitted.vae$encoder(torch_tensor(self$data))
                               res = self$fitted.vae(torch_tensor(self$data))[[2]]
                             } else {
                               res = self$fitted.vae(torch_tensor(newdata))[[2]]
                             }
                             res = as.matrix(res)
                             colnames(res) =  paste('PC', 1:ncol(res), sep='')
                             return(res)
                           }
                         )
  )
  predictorObj = predictor$new(vae, as.matrix(residuals))

  return(list(nLatVars = resPpca$n, latVars = predictorObj$predict(), lvPredictor = predictorObj))
})
