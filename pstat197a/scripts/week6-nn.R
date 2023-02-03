library(tidyverse)
library(tidymodels)
library(keras)
library(tensorflow)

setwd("~/pstat197/pstat197a/materials/slides")
load('data/claims-clean.RData')



set.seed(102722)
partitions <- clean %>%
  initial_split(prop = 0.8)

train_text <- training(partitions) %>%
  pull(text_clean)
train_labels <- training(partitions) %>%
  pull(bclass) %>%
  as.numeric() - 1

## preprocessing layer

preprocess_layer <- layer_text_vectorization(
  standardize = NULL,
  split = 'whitespace',
  ngrams = NULL,
  max_tokens = NULL,
  output_mode = 'tf_idf'
)

preprocess_layer %>% adapt(train_text)

train_processed <- preprocess_layer(train_text)
inshape = train_processed$shape[2] %>% as.numeric()

# architecture

model <- keras_model_sequential(input_shape = inshape) %>%
  layer_dropout(0.2) %>%
  layer_dense(units = 25) %>%
  layer_dropout(0.2) %>%
  layer_dense(units = 20) %>%
  layer_dropout(0.2) %>%
  layer_dense(1)

summary(model)

model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = 'adam',
  metrics = 'binary_accuracy'
)

history <- model %>%
  fit(train_processed, 
      train_labels,
      validation_split = 0.2,
      epochs = 15)

test_processed <- testing(partitions) %>%
  pull(text_clean) %>%
  preprocess_layer()
test_labels <- testing(partitions) %>%
  pull(bclass) %>%
  as.numeric() - 1

model %>% evaluate(test_processed, test_labels)

probability_model <- keras_model_sequential() %>%
  model() %>%
  layer_activation(activation = 'sigmoid')

preds <- probability_model(test_processed) %>%
  as.numeric() 

panel <- metric_set(sensitivity, specificity, accuracy, roc_auc)
testing(partitions) %>%
  select(bclass) %>%
  bind_cols(pred = preds) %>%
  mutate(pred.bclass = factor(pred > 0.5, labels = levels(bclass))) %>%
  panel(truth = bclass,
        estimate = pred.bclass,
        pred,
        event_level = 'second')
            
