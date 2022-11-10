library(splines)
library(tidyverse)
library(modelr)
setwd("~/pstat197/pstat197a/materials/slides")


knot_pts <- seq_range(1:365, n = 6)
x = seq(from = 1, 
        to = 365, 
        length = 200)
spl <- bs(x, 
            knots = knot_pts[-c(1, 6)], 
            degree = 4) 
  
mx <- spl[, ]
colnames(mx) <- paste('f', 1:ncol(mx), sep = '')

spl_out <- as_tibble(mx) %>%
  bind_cols(x = x)
  

tibble(data = spl_out) %>%
  unnest(data) %>%
  pivot_longer(cols = -x, 
               names_to = 'basis',
               values_to = 'value') %>%
  ggplot(aes(x = x, 
             y = value, 
             color = basis)) +
  geom_path() +
  scale_y_log10() +
  labs(x = 'x', y = expression(paste(f[j], '(x)'))) +
  geom_vline(xintercept = knot_pts, 
             linetype = 'dashed')

ggsave(filename = 'img/even-knots-bspline-log.png',
       width = 5, height = 3, units = 'in')


set.seed(110622)
knot_pts <- sample(1:365, size = 4, replace = F)
x = seq(from = 1, 
        to = 365, 
        length = 200)

spl <- bs(x, 
          knots = knot_pts, 
          degree = 4) 

mx <- spl[, ]
colnames(mx) <- paste('f', 1:ncol(mx), sep = '')

spl_out <- as_tibble(mx) %>%
  bind_cols(x = x)


tibble(data = spl_out) %>%
  unnest(data) %>%
  pivot_longer(cols = -x, 
               names_to = 'basis',
               values_to = 'value') %>%
  ggplot(aes(x = x, 
             y = value, 
             color = basis)) +
  geom_path() +
  labs(x = 'x', y = expression(paste(f[j], '(x)'))) +
  scale_y_log10() +
  geom_vline(xintercept = c(0, 365, knot_pts), 
             linetype = 'dashed')

ggsave(filename = 'img/uneven-knots-bspline-log.png',
       width = 5, height = 3, units = 'in')
