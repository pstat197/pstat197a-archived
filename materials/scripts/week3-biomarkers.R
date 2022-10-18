library(tidyverse)
library(infer)
setwd("~/pstat197/pstat197a/materials/slides/data")

# # preprocessing
# read_csv('biomarker-raw.csv', skip = 1) %>%
#   rename(group = `...1`,
#          ados = `...1320`) %>%
#   select(-Target) %>%
#   mutate(ados = na_if(ados, '-')) %>%
#   mutate(ados = as.numeric(ados),
#          group = factor(group)) %>%
#   select(group, ados, everything()) %>%
#   filter(!is.na(group)) %>%
#   write_csv('biomarker-clean.csv')

# read data
asd <- read_csv('biomarker-clean.csv')

# check that matches data description
asd %>% count(group)
asd %>% ncol()

# preprocessing (following paper)

trim_fn <- function(x){
  x[x > 3] <- 3
  x[x < -3] <- -3
  
  return(x)
}

# trim_fn(-5:5)

asd_clean <- asd %>% 
  select(-ados) %>%
  # log transform
  mutate(across(.cols = -group, log10)) %>%
  # center and scale
  mutate(across(.cols = -group, ~ scale(.x)[, 1])) %>%
  # trim outliers (affects results??)
  mutate(across(.cols = -group, trim_fn))

## MULTIPLE TESTING
####################

# one t test
asd_clean %>%
  t_test(formula = CHIP ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)

# nesting for many t tests
asd_nested <- asd_clean %>%
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group))

# look at one group
asd_nested %>% 
  slice(1L) %>% 
  unnest(cols = data)

# test on one group
asd_nested %>% 
  slice(1L) %>% 
  unnest(cols = data) %>% 
  t_test(formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)

# check variance ratios
asd_clean %>% 
  pivot_longer(-group, 
               names_to = "protein", 
               values_to = "level") %>%
  group_by(protein, group) %>%
  summarize(level.var = var(level), .groups = 'drop') %>%
  pivot_wider(id_cols = protein, 
              names_from = 'group', 
              values_from = 'level.var') %>%
  mutate(var.ratio = ASD/TD) %>%
  ggplot(aes(x = var.ratio)) +
  geom_histogram(aes(y = ..density..), bins = 50) +
  scale_x_log10()

# compute for several groups
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

tt_out <- asd_nested %>% 
  # slice(1L:100L) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value)

# multiple testing corrections
m <- nrow(tt_out)
hm <- log(m) + 1/(2*m) - digamma(1)
  
tt_corrected <- tt_out %>%
  select(data, protein, p_value) %>%
  mutate(rank = row_number()) %>%
  mutate(p_bh = p_value*m/rank,
         p_by = p_value*m*hm/rank,
         p_bonf = p_value*m)

## comparing methods
tt_corrected %>%
  rename(bh.adj = p_bh,
         by.adj = p_by,
         bonf.adj = p_bonf,
         p.value = p_value) %>%
  mutate(p.raw = p.value) %>%
  pivot_longer(cols = c(contains('.adj'), p.raw), 
               names_to = 'correction',
               values_to = 'p.adj') %>%
  ggplot(aes(x = p.value, 
             y = p.adj, 
             color = correction,
             linetype = correction)) +
  geom_path() +
  scale_y_log10() +
  scale_x_sqrt() +
  labs(x = 'p value', y = 'adjusted p value')

## 'volcano' plot
ratio_fn <- function(.df){
  .df %>%
    group_by(group) %>%
    summarize(diff = mean(level_raw)) %>%
    spread(group, diff) %>%
    mutate(ratio = ASD/TD) %>%
    pull(ratio)
}

ratios <- asd %>%
  select(-ados) %>%
  pivot_longer(cols = -group,
               names_to = 'protein',
               values_to = 'level_raw') %>%
  nest(data = c(group, level_raw)) %>%
  mutate(ratio = map(data, ratio_fn)) %>%
  select(protein, ratio) %>%
  unnest(cols = ratio)

plot_df <- tt_corrected %>%
  select(p_by, protein) %>%
  left_join(ratios, by = 'protein')

plot_df %>%
  ggplot(aes(x = log2(ratio), y = -log10(p_by))) +
  geom_point(aes(color = factor(sign(log2(ratio))),
                 alpha = p_by < 0.05)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dotdash') +
  geom_vline(xintercept = 0, linetype = 'dotdash') +
  guides(color = guide_none(), alpha = guide_none())

# select significant proteins -- close but not exact match to paper
tt_corrected %>%
  select(protein, p_by) %>%
  slice_min(p_by, n = 10)

## CORRELATIONS
################

asd_clean <- asd %>% 
  # log transform
  mutate(across(.cols = -c(group, ados), log10)) %>%
  # center and scale
  mutate(across(.cols = -c(group, ados), ~ scale(.x)[, 1])) %>%
  # trim outliers (affects results??)
  mutate(across(.cols = -c(group, ados), trim_fn)) %>%
  # ados only for ASD group
  filter(group == 'ASD') %>%
  select(-group)

# look for linear association
asd_clean %>%
  select(1:17) %>%
  pivot_longer(cols = -ados,
               names_to = 'protein',
               values_to = 'level') %>%
  ggplot(aes(x = level, y = ados)) +
  geom_point() +
  geom_smooth(formula = 'y ~ x', 
              method = 'lm', 
              se = F) +
  facet_wrap(~ protein, nrow = 4, ncol = 4)

# compute correlations
ados_cors <- asd_clean %>%
  pivot_longer(cols = -ados,
               names_to = 'protein',
               values_to = 'level') %>%
  group_by(protein) %>%
  summarize(correlation = cor(ados, level)) %>%
  arrange(correlation) %>%
  mutate(abs.corr = abs(correlation),
         rank = row_number()) 

cor_test <- function(x, y){
  cor_out <- cor.test(x, y, method = 'pearson')
  tibble(estimate = cor_out$estimate,
         p.value = cor_out$p.value)
}

asd_clean %>%
  pivot_longer(cols = -ados,
               names_to = 'protein',
               values_to = 'level') %>%
  group_by(protein) %>%
  summarize(correlation = cor_test(ados, level)) %>%
  unnest(correlation) %>%
  arrange(p.value)

# plot correlations
ados_cors %>%
  ggplot(aes(x = rank,
             y = correlation)) +
  geom_path() +
  geom_point(data = slice_max(ados_cors, abs.corr, n = 10)) +
  geom_text(data = slice_max(ados_cors, abs.corr, n = 10),
            aes(label = protein, hjust = 'inward'),
            check_overlap = T,
            size = 3)

# top 10
ados_cors %>% slice_max(abs.corr, n = 10)

# regression approach (equivalent)
fit_fn <- function(.df){
  lm(ados ~ level, data = .df)
}

asd_clean %>%
  pivot_longer(cols = -ados,
               names_to = 'protein',
               values_to = 'level') %>%
  nest(data = c(ados, level)) %>%
  mutate(fit = map(data, fit_fn),
         fit_tidy = map(fit, broom::tidy)) %>%
  unnest(fit_tidy) %>%
  filter(term == 'level') %>%
  select(protein, estimate, p.value) %>%
  slice_min(p.value, n = 10)

# do these pass multiple testing?
m <- ncol(asd_clean) - 1
hm <- log(m) + 1/(2*m) - digamma(1)

asd_clean %>%
  pivot_longer(cols = -ados,
               names_to = 'protein',
               values_to = 'level') %>%
  nest(data = c(ados, level)) %>%
  mutate(fit = map(data, fit_fn),
         fit_tidy = map(fit, broom::tidy)) %>%
  unnest(fit_tidy) %>%
  filter(term == 'level') %>%
  select(protein, estimate, p.value) %>%
  arrange(p.value) %>%
  mutate(rank = row_number(),
         p.adj = p.value*m*hm/rank) %>%
  slice_min(p.value, n = 10) %>%
  select(protein, p.adj, p.value)

## visualize
selected_proteins <- ados_cors %>% 
  slice_max(abs.corr, n = 10) %>%
  pull(protein)

asd_clean %>%
  select(c(ados, any_of(selected_proteins))) %>%
  pivot_longer(cols = -ados,
               names_to = 'protein',
               values_to = 'level') %>%
  ggplot(aes(x = level, y = ados)) +
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = 'lm', se = F) +
  facet_wrap(~ protein, nrow = 2)

library(randomForest)
asd_clean <- asd %>% 
  select(-ados) %>%
  # log transform
  mutate(across(.cols = -group, log10)) %>%
  # center and scale
  mutate(across(.cols = -group, ~ scale(.x)[, 1])) %>%
  # trim outliers (affects results??)
  mutate(across(.cols = -group, trim_fn))

asd_preds <- asd_clean %>% select(-group)
asd_resp <- asd_clean %>% pull(group) %>% factor()

set.seed(101222)
rf_out <- randomForest(x = asd_preds, y = asd_resp,
                       mtry = 100, ntree = 1000, 
                       importance = T)

rf_out$confusion
rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  select(protein)
