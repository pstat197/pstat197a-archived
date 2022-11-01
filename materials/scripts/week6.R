library(tidyverse)
library(tidymodels)
library(modelr)
library(Matrix)
library(sparsesvd)
library(glmnet)
setwd("~/pstat197/pstat197a/materials/slides")
load('data/claims-processed.RData')

# rename for merging
claims_words <- claims_multi %>%
  rename_with(~ paste('word.', .x, sep = ''))

claims_bigrams <- claims_multi_bigram %>%
  rename_with(~ paste('bigram.', .x, sep = ''))
  
# partition
set.seed(103022)
partitions <- bind_cols(claims_words, claims_bigrams) %>%
  initial_split(prop = 0.8)

test_dtm_words <- testing(partitions) %>%
  select(contains('word.')) %>%
  rename_with(~ gsub('word.', '', .x)) %>%
  select(-.id, -bclass, -mclass)

test_dtm_bigrams <- testing(partitions) %>%
  select(contains('bigram.')) %>%
  rename_with(~ gsub('bigram.', '', .x)) %>%
  select(-.id, -bclass, -mclass)

test_labels <- testing(partitions) %>%
  select(contains('word.')) %>%
  rename_with(~ gsub('word.', '', .x)) %>%
  select(.id, bclass, mclass)

train_dtm_words <- training(partitions) %>%
  select(contains('word.')) %>%
  rename_with(~ gsub('word.', '', .x)) %>%
  select(-.id, -bclass, -mclass)

train_dtm_bigrams <- training(partitions) %>%
  select(contains('bigram.')) %>%
  rename_with(~ gsub('bigram.', '', .x)) %>%
  select(-.id, -bclass, -mclass)

train_labels <- training(partitions) %>%
  select(contains('word.')) %>%
  rename_with(~ gsub('word.', '', .x)) %>%
  select(.id, bclass, mclass)

# some helper functions
projection_fn <- function(.dtm, .prop){
  # coerce feature matrix to sparse
  dtm_mx <- .dtm %>%
    as.matrix() %>%
    as('sparseMatrix')
  
  # compute svd
  svd_out <- sparsesvd(dtm_mx)
  
  # select number of projections
  var_df <- tibble(var = svd_out$d^2) %>%
    mutate(pc = row_number(),
           cumulative = cumsum(var)/sum(var))
  
  n_pc <- which.min(var_df$cumulative < .prop)
  
  # extract loadings
  loadings <- svd_out$v[, 1:n_pc] %>% as.matrix()
  
  # extract scores
  scores <- (dtm_mx %*% svd_out$v[, 1:n_pc]) %>% as.matrix()
  
  # adjust names
  colnames(loadings) <- colnames(scores) <- paste('pc', 1:n_pc, sep = '')
  
  # output
  out <- list(n_pc = n_pc,
              var = var_df,
              projection = loadings,
              data = as_tibble(scores))
  
  return(out)
}

reproject_fn <- function(.dtm, .projection_fn_out){
  as_tibble(as.matrix(.dtm) %*% .projection_fn_out$projection)
}

# project the data onto subspaces
proj_out_words <- projection_fn(.dtm = train_dtm_words,
                                .prop = 0.6)

proj_out_bigrams <- projection_fn(.dtm = train_dtm_bigrams,
                                .prop = 0.6)


# separate for model fitting with glmnet
x_train_words <- proj_out_words$data %>% as.matrix()
y_train <- train_labels %>% pull(bclass)

# fit first model (words)
alpha_words <- 0.3
fit_words <- glmnet(x = x_train_words, 
                  y = y_train, 
                  family = 'binomial',
                  alpha = alpha_words)

set.seed(103022)
cvout_words <- cv.glmnet(x = x_train_words, 
                   y = y_train, 
                   family = 'binomial',
                   alpha = alpha_words)

# evaluate accuracy of first model (words)
panel <- metric_set(sensitivity, specificity, accuracy, roc_auc)

x_test_words <- reproject_fn(test_dtm_words, 
                             proj_out_words) %>% 
  as.matrix()

accuracy_onestage <- test_labels %>%
  mutate(pred = predict(fit_words, 
                        s = cvout_words$lambda.min, 
                        newx = x_test_words,
                        type = 'response')[, 1]) %>%
  mutate(bclass.pred = factor(pred > 0.5, 
                              labels = levels(bclass))) %>%
  panel(truth = bclass,
        estimate = bclass.pred,
        pred,
        event_level = 'second')

# fit second model
pred.logodds <- predict(fit_words, 
        s = cvout_words$lambda.min, 
        newx = x_train_words,
        type = 'link')[, 1]

x_train_bigrams <- proj_out_bigrams$data %>% 
  as.matrix()

alpha_bigrams <- 0
fit_bigrams <- glmnet(x = x_train_bigrams, 
                    y = y_train, 
                    family = 'binomial',
                    alpha = alpha_bigrams,
                    intercept = F,
                    offset = pred.logodds)

set.seed(103022)
cvout_bigrams <- cv.glmnet(x = x_train_bigrams, 
                   y = y_train, 
                   family = 'binomial',
                   alpha = alpha_bigrams,
                   intercept = F,
                   offset = pred.logodds)


x_test_bigrams <- reproject_fn(test_dtm_bigrams, 
                               proj_out_bigrams) %>% 
  as.matrix()

accuracy_twostage <- test_labels %>%
  mutate(pred.offset = predict(fit_words, 
                        s = cvout_words$lambda.min, 
                        newx = x_test_words,
                        type = 'link')[, 1],
         pred = predict(fit_bigrams,
                        s = cvout_bigrams$lambda.min,
                        newx = x_test_bigrams,
                        newoffset = pred.offset, 
                        type = 'response')[, 1]) %>%
  mutate(bclass.pred = factor(pred > 0.5, 
                              labels = levels(bclass))) %>%
  panel(truth = bclass,
        estimate = bclass.pred,
        pred,
        event_level = 'second')
