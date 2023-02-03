library(tidyverse)
library(tidytext)
library(textstem)
library(rvest)
library(stopwords)
library(tokenizers)
library(qdapRegex)
library(tidymodels)
library(Matrix)
library(sparsesvd)
setwd("~/pstat197/pstat197a/materials/slides")
load('data/carpe-raw-subsample.RData')

## PREPROCESSING/NLP
####################

# read, lump and relabel
rawdata_relabeled <- rawdata %>%
  mutate(bclass = fct_lump(internal_feedback, 
                           prop = 0.5, 
                           other_level = 'relevant'),
         bclass = fct_recode(bclass, 
                             irrelevant = 'N/A: No relevant content.'),
         bclass = fct_infreq(bclass))  %>%
  mutate(.id = paste('url', row_number(), sep = '')) %>%
  select(.id, bclass, text_tmp)

# function to parse html
parse_fn <- function(.html){
  read_html(.html) %>%
    html_elements('p') %>%
    html_text2() %>%
    str_c(collapse = ' ') %>%
    rm_url() %>%
    rm_email() %>%
    str_remove_all('\'') %>%
    str_replace_all(paste(c('\n', 
                            '[[:punct:]]', 
                            'nbsp', 
                            '[[:digit:]]', 
                            '[[:symbol:]]'),
                          collapse = '|'), ' ') %>%
    str_replace_all("([a-z])([A-Z])", "\\1 \\2") %>%
    tolower() %>%
    str_replace_all("\\s+", " ")
}

# clean text
clean <- rawdata_relabeled %>%
  filter(str_detect(text_tmp, '<!')) %>%
  rowwise() %>%
  mutate(text_clean = parse_fn(text_tmp)) %>%
  select(-text_tmp) %>%
  unnest(text_clean)

# generate tfidf document term matrix
claims <- clean %>% 
  unnest_tokens(output = token, 
                input = text_clean, 
                token = 'words',
                stopwords = str_remove_all(stop_words$word, 
                                           '[[:punct:]]')) %>%
  mutate(token.lem = lemmatize_words(token)) %>%
  filter(str_length(token.lem) > 2) %>%
  count(.id, bclass, token.lem, name = 'n') %>%
  bind_tf_idf(term = token.lem, 
              document = .id,
              n = n) %>%
  pivot_wider(id_cols = c('.id', 'bclass'),
              names_from = 'token.lem',
              values_from = 'tf_idf',
              values_fill = 0)

# save(claims, file = 'data/claims-tfidf.RData')

# tokenize w/ bigrams
claims_bigrams <- clean %>% 
  unnest_tokens(output = token, 
                input = text_clean, 
                token = 'ngrams',
                n = 2,
                stopwords = str_remove_all(stop_words$word, 
                                           '[[:punct:]]')) %>%
  filter(str_length(token) > 6) %>%
  count(.id, bclass, token, name = 'n') %>%
  bind_tf_idf(term = token, 
              document = .id,
              n = n) %>%
  pivot_wider(id_cols = c('.id', 'bclass'),
              names_from = 'token',
              values_from = 'tf_idf',
              values_fill = 0)

# generate multiclass labels
multi_labels <- rawdata %>%
  transmute(mclass = fct_lump(internal_feedback,
                             prop = 0.1),
            mclass = fct_recode(mclass,
                               irrelevant = 'N/A: No relevant content.',
                               physical = 'Physical Activity',
                               fatality = 'Possible Fatality',
                               unlawful = 'Potentially unlawful activity',
                               other = 'Other'),
            .id = paste('url', row_number(), sep = '')) 

# merge with word token tfidf
claims_multi <- multi_labels %>% 
  right_join(claims, by = '.id') 

# merge with bigram tfidf
claims_multi_bigram <- multi_labels %>% 
  right_join(claims_bigrams, by = '.id') 

# save
save(list = c('claims_multi', 'claims_multi_bigram'), 
     file = 'data/claims-processed.RData')


claims_long <- clean %>% 
  unnest_tokens(output = token, 
                input = text_clean, 
                token = 'words',
                stopwords = str_remove_all(stop_words$word, 
                                           '[[:punct:]]')) %>%
  mutate(token.lem = lemmatize_words(token)) %>%
  filter(str_length(token.lem) > 2) %>%
  count(.id, bclass, token.lem, name = 'n') %>%
  bind_tf_idf(term = token.lem, 
              document = .id,
              n = n) 


save(list = c('claims_multi', 'claims_long', 'clean'),
     file = 'data/claims-clean.RData')

## PROJECTION
##############
library(tidyverse)
library(tidymodels)
library(modelr)
library(Matrix)
library(sparsesvd)
library(glmnet)
setwd("~/pstat197/pstat197a/materials/slides")
load('data/claims-tfidf.RData')

# partition
set.seed(102722)
partitions <- claims %>% initial_split(prop = 0.8)

test_dtm <- testing(partitions) %>%
  select(-.id, -bclass)
test_labels <- testing(partitions) %>%
  select(.id, bclass)

train_dtm <- training(partitions) %>%
  select(-.id, -bclass)
train_labels <- training(partitions) %>%
  select(.id, bclass)

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

# project the data onto a subspace
proj_out <- projection_fn(.dtm = train_dtm,
                          .prop = 0.7)

# merge with labels
train <- bind_cols(train_labels, proj_out$data)

# fit logistic regression
fit <- glm(bclass ~ ., data = train, family = 'binomial')  

# add a penalty to reduce overfitting
x_train <- proj_out$data %>% as.matrix()
y_train <- train_labels %>% pull(bclass)
fit_reg <- glmnet(x = x_train, 
                  y = y_train, 
                  data = train, 
                  family = 'binomial',
                  alpha = 0.3)

set.seed(102722)
cvout <- cv.glmnet(x = x_train, 
                y = y_train, 
                data = train, 
                family = 'binomial',
                alpha = 0.3)

# evaluate accuracy
panel <- metric_set(sensitivity, specificity, accuracy, roc_auc)

x_test <- reproject_fn(test_dtm, proj_out) %>% as.matrix()
test %>%
  mutate(pred = predict(fit_reg, 
                         s = cvout$lambda.1se, 
                         newx = x_test,
                        type = 'response')[, 1]) %>%
  mutate(bclass.pred = factor(pred > 0.5, 
                              labels = levels(bclass))) %>%
  panel(truth = bclass,
        estimate = bclass.pred,
        pred,
        event_level = 'second')

claims_multi %>%
  select(.id, mclass) %>%
  right_join(train_labels, by = '.id') 

test_mlabels <- claims_multi %>%
  select(.id, mclass) %>%
  right_join(testing(partitions), by = '.id') %>% 
  select(mclass)

mfit_reg <- glmnet(x = x_train, 
                  y = train_mlabels,
                  family = 'multinomial',
                  grouped = T,
                  alpha = 0.7)

set.seed(102722)
cvout <- cv.glmnet(x = x_train, 
                   y = train_mlabels,  
                   family = 'multinomial',
                   grouped = T,
                   alpha = 0.7)

pred_df <- as_tibble(predict(mfit_reg, 
                    s = cvout$lambda.1se, 
                    newx = x_train,
                    type = 'response')[, , 1]) 
pred_df

