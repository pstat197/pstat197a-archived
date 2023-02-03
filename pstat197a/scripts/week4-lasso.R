library(tidyverse)
library(glmnet)
library(tidymodels)
library(yardstick)
setwd("~/pstat197/pstat197a/materials/slides/data")
load('biomarker-clean.RData')

biomarker_alt <- read_csv('biomarker-clean.csv') %>%
  select(-ados) %>%
  mutate(across(-group, ~scale(log(.x))[,1]),
         class = as.numeric(group == 'ASD'))

biomarker_clean <- biomarker_clean %>%
  select(-ados) %>%
  mutate(class = as.numeric(group == 'ASD'))

set.seed(101622)
partitions <- biomarker_clean %>%
  initial_split(prop = 0.8)

x_train <- training(partitions) %>%
  select(-group, -class) %>%
  as.matrix()
y_train <- training(partitions) %>%
  pull(class)

fit <- glmnet(x_train, y_train, family = 'binomial')

cv_out <- cv.glmnet(x_train, y_train, family = 'binomial', nfolds = 5, type.measure = 'auc')

cv_summary <- tibble(lambda = cv_out$lambda,
       auc.mean = cv_out$cvm,
       auc.sd = cv_out$cvsd,
       df = cv_out$nzero)


ggplot(cv_summary,
       aes(x = -log(lambda), y = auc.mean)) +
  geom_path() +
  geom_ribbon(aes(ymin = auc.mean - auc.sd, 
                  ymax = auc.mean + auc.sd),
              alpha = 0.2) +
  geom_point(data = filter(cv_summary, lambda == cv_out$lambda.min),
             color = 'red') +
  geom_point(data = filter(cv_summary, lambda == cv_out$lambda.1se))

ggplot(cv_summary,
       aes(x = -log(lambda), y = df)) +
  geom_path() +
  geom_point(data = filter(cv_summary, lambda == cv_out$lambda.min),
             color = 'red') +
  geom_point(data = filter(cv_summary, lambda == cv_out$lambda.1se)) +
  scale_y_log10()



fit_coef <- coef(fit, s = cv_out$lambda.1se) %>% as.matrix()

tibble(term = rownames(fit_coef),
       estimate = fit_coef) %>%
  filter(estimate != 0)
        
x_test <- testing(partitions) %>%
  select(-group, -class) %>%
  as.matrix()

panel <- metric_set(sensitivity, specificity, accuracy, roc_auc)

testing(partitions) %>%
  select(group, class) %>%
  bind_cols(pred = predict(fit, 
                        newx = x_test, 
                        s = cv_out$lambda.1se, 
                        type = 'response')[,1]) %>%
  mutate(class.pred = as.numeric(pred > 0.5)) %>%
  panel(truth = factor(class, labels = c('TD', 'ASD')),
        estimate = factor(class.pred, labels = c('TD', 'ASD')),
        pred,
        event_level = 'second')

terms <- tibble(term = rownames(fit_coef),
       estimate = fit_coef) %>%
  filter(estimate != 0) %>%
  pull(term)


biomarker_sub <- biomarker_alt %>%
  select(class, any_of(terms))

fit_logit <- glm(class ~ ., data = biomarker_sub, family = 'binomial')

testing(partitions) %>%
  modelr::add_predictions(fit_logit, type = 'response') %>%
  mutate(class.pred = as.numeric(pred > 0.5)) %>%
  panel(truth = factor(class, labels = c('TD', 'ASD')),
        estimate = factor(class.pred, labels = c('TD', 'ASD')),
        pred,
        event_level = 'second')

fit_null <- glm(class ~ 1, data = biomarker_alt, family = 'binomial')

fit_null %>% MASS::stepAIC(scope = 'class ~ DERM + CHIP', direction = 'forward')
