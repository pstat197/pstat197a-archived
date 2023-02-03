library(keras)
library(tensorflow)
library(mvtnorm)
library(tidyverse)

set.seed(110222)
x <- runif(n = 500, min = -5, max = 5)
beta <- c(1, -2)
eta <- cbind(1, x) %*% beta
p <- 1/(1 + exp(-eta))
y <- rbinom(n = 500, size = 1, prob = p)

sim_data <- as_tibble(x) %>% bind_cols(y = y)

fit <- glm(y ~ ., family = 'binomial', data = sim_data)
coef(fit)

library(keras) 

# input layer 
inputs <- layer_input(shape = 1) 

# outputs compose input + dense layers 
outputs <- inputs %>% 
  layer_dense(units = 1, activation = 'sigmoid') 

# create and compile model 
model <- keras_model(inputs = inputs, outputs = outputs) 

model %>%
  compile(optimizer = optimizer_adam(learning_rate = 0.1),
          loss = 'binary_crossentropy',
          metric = 'accuracy')

model %>%
  fit(x = x,
      y = y,
      epochs = 10)

betahat_nn <- get_weights(model) %>% unlist() %>% rev()
betahat_glm <- coef(fit)

rbind(glm = betahat_glm,
      nn = betahat_nn)
