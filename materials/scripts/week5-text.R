library(tidyverse)
library(tidytext)
library(textstem)
library(rvest)
setwd("~/pstat197/pstat197a/materials/datasets")

# grab a few rows of raw html
raw_sample <- read_csv('web-fraud/fraud-raw.csv', n_max = 50)

## PREPROCESSING
##################

# grab one page
page <- raw_sample %>% slice(1) %>% pull(text_tmp) 
read_csv('web-fraud/fraud-raw.csv', n_max = 2)


# view raw html
page
# write_file(page, file = 'page.html')

# parse raw html
read_html(page) %>% 
  html_text2()

# generate a regex for string removal
remove <- c('\n', 
            '[[:punct:]]', 
            'nbsp', 
            '[[:digit:]]', 
            '[[:symbol:]]') %>%
  paste(collapse = '|')

# process the page
read_html(page) %>%
  # select paragraph contents
  html_elements('p') %>%
  # parse html
  html_text2() %>%
  # collapse into one string
  str_c(collapse = ' ') %>%
  # drop links
  qdapRegex::rm_url() %>%
  # remove punctuation, numbers, etc.
  str_replace_all(remove, ' ') %>%
  # add whitespace between captial letters
  str_replace_all("([a-z])([A-Z])", "\\1 \\2") %>%
  # convert to lowercase letters
  tolower() %>%
  # remove extra whitespace
  str_replace_all("\\s+", " ")

# wrap in a function
parse_fn <- function(.html){
  read_html(.html) %>%
    # select paragraph elements
    html_elements('p') %>%
    # parse html
    html_text2() %>%
    # collapse into one string
    str_c(collapse = ' ') %>%
    # drop links
    qdapRegex::rm_url() %>%
    # remove punctuation, numbers, etc.
    str_replace_all(remove, ' ') %>%
    # add whitespace between captial letters
    str_replace_all("([a-z])([A-Z])", "\\1 \\2") %>%
    # convert to lowercase letters
    tolower() %>%
    # remove extra whitespace
    str_replace_all("\\s+", " ")
}

# try it on a new page
raw_sample %>% 
  slice(38:40) %>%
  pull(text_tmp) %>% 
  parse_fn()

## NLP
########

clean_sample <- raw_sample %>%
  filter(!str_detect(text_tmp, 'Error:')) %>%
  filter(str_detect(text_tmp, '<!')) %>%
  mutate(page.id = paste('url', row_number(), sep = ''),
         label = factor(internal_feedback)) %>%
  select(page.id, label, text_tmp) %>%
  rowwise() %>%
  mutate(text_clean = parse_fn(text_tmp)) %>%
  select(-text_tmp) %>%
  unnest(text_clean)

clean_sample_long <- clean_sample %>%
  # tokenize
  unnest_tokens(word, text_clean, token = 'words') %>%
  # remove stopwords
  anti_join(get_stopwords(), by = 'word') %>%
  # lemmatize
  mutate(lemma = lemmatize_words(word)) %>%
  # tally lemma counts
  count(page.id, lemma) %>%
  # compute tf-idf
  bind_tf_idf(lemma, page.id, n) 

# document term matrix
clean_sample_long %>%
  cast_dtm(page.id, lemma, n, weighting = tm::weightTfIdf) %>%
  as.matrix() %>% 
  as_tibble()

## EXPLORATORY ANALYSIS
########################

# read in more pages
raw_full <- read_csv('web-fraud/fraud-raw.csv')

# check label counts
raw_full %>%
  count(internal_feedback) %>%
  mutate(freq = n/sum(n))

# lump together infrequent labels
raw_full_labeled <- raw_full %>%
  mutate(.id = paste('url', row_number(), sep = ''),
         label = factor(internal_feedback),
         class = fct_lump_prop(label, 
                               prop = 0.05, 
                               other_level = 'Other claim content'),
         class_binary = fct_lump_n(label, 
                                   n = 1, 
                                   other_level = 'Relevant claim content')) 

# downsample
set.seed(102222)
raw_subsample <- raw_full_labeled %>%
  slice_sample(prop = 0.2)

# preprocess
clean_subsample <- raw_subsample %>%
  filter(str_detect(text_tmp, '<!')) %>%
  select(.id, label, class, class_binary, text_tmp) %>%
  rowwise() %>%
  mutate(text_clean = parse_fn(text_tmp)) %>%
  select(-text_tmp) %>%
  unnest(text_clean)

subsample_long <- clean_subsample %>%
  # tokenize
  unnest_tokens(word, text_clean, token = 'words') %>%
  # remove stopwords
  anti_join(get_stopwords(), by = 'word') %>%
  # lemmatize
  mutate(lemma = lemmatize_words(word)) %>%
  # tally lemma counts
  count(.id, class, lemma) %>%
  # compute tf-idf
  bind_tf_idf(lemma, .id, n) 

# document term matrix
subsample_dtm <- subsample_long %>%
  cast_dtm(document = .id, 
           term = lemma, 
           value = n, 
           weighting = tm::weightTfIdf) 

dtm_df <- bind_cols(.id = rownames(subsample_dtm), 
          as.matrix(subsample_dtm)) 

# bind to labels
subsample_df <- clean_subsample %>%
  select(.id, class_binary) %>%
  rename(class_fct = class_binary) %>%
  right_join(dtm_df, by = '.id')

library(randomForest)
library(tidymodels)

set.seed(102222)
partitions <- subsample_df %>%
  initial_split()

x_train <- training(partitions) %>%
  select(-.id, -class_fct) 

y_train <- training(partitions) %>% pull(class_fct)

x_test <- testing(partitions) %>%
  select(-.id, -class_fct) 

y_test <- testing(partitions) %>% pull(class_fct)

set.seed(102222)
rf_out <- randomForest(x = x_train, 
                       y = y_train,
                       xtest = x_test,
                       ytest = y_test,
                       mtry = 100,
                       maxnodes = 20,
                       ntree = 500)

(rf_out$confusion %>% diag() %>% sum())/nrow(x_train)

training(partitions) %>% 
  count(class_fct) %>%
  mutate(prop = n/sum(n))

summary(rf_out)
rf_out$importance %>%
  as_tibble() %>%
  bind_cols(lemma = rownames(rf_out$importance)) %>%
  arrange(desc(MeanDecreaseGini))
