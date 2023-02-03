library(tidyverse)
setwd("~/pstat197/pstat197a/materials/slides/data")

background <- read_csv('background-clean.csv')
interest <- read_csv('interest-clean.csv')

## individual variable summaries
###############################

# one variable, one summary
background %>%
  select(math.comf) %>%
  summarize_all(mean)

# many variables, same summary
background %>%
  select(contains('comf')) %>%
  summarise_all(.funs = mean)

background %>%
  select(contains('comf')) %>%
  summarise_all(.funs = median)

# many variables, many summaries
background %>%
  select(contains('comf')) %>%
  summarise_all(.funs = list(mean = mean, 
                             median = median,
                             min = min, 
                             max = max)) %>%
  gather(stat, val) %>%
  separate(stat, into = c('variable', 'stat'), sep = '_') %>%
  spread(stat, val)

# proficiency responses are factors
background %>%
  pull(math.prof) %>%
  factor() %>%
  fct_count()

# but levels are ordered
background %>%
  pull(math.prof) %>%
  factor(levels = c('Beg', 'Int', 'Adv')) %>%
  fct_count()

# same summary, each variable
background %>%
  select(contains('.prof')) %>%
  mutate_all(~factor(.x, levels = c('Beg', 'Int', 'Adv'))) %>%
  summarize_all(fct_count)

# clean up names a little
background %>%
  select(contains('.prof')) %>%
  mutate_all(~factor(.x, levels = c('Beg', 'Int', 'Adv'))) %>%
  rename_with(~gsub('.prof', '', .x)) %>%
  summarize_all(fct_count)

# what about converting to numeric? is this meaningful??
background %>%
  select(contains('.prof')) %>%
  mutate_all(~factor(.x, levels = c('Beg', 'Int', 'Adv'))) %>%
  rename_with(~gsub('.prof', '', .x)) %>%
  mutate_all(as.numeric) %>%
  summarize_all(.funs = list(mean = mean, 
                             median = median)) %>%
  gather(stat, val) %>%
  separate(stat, into = c('variable', 'stat'), sep = '_') %>%
  spread(stat, val)

## multivariable thinking
###########################

# unique combinations of proficiency ratings
proficiency <- background %>%
  select(contains('.prof')) %>%
  mutate_all(~factor(.x, levels = c('Beg', 'Int', 'Adv'))) %>%
  mutate_all(as.numeric)

proficiency %>% 
  rename_with(~gsub('.prof', '', .x)) %>%
  group_by(prog, math, stat) %>%
  count()

# unique combinations of comfort ratings
comfort <- background %>%
  select(contains('comf'))

comfort %>% 
  rename_with(~gsub('.comf', '', .x)) %>%
  group_by(prog, math, stat) %>%
  count()

# cluster responses into three groups
set.seed(92922)
clust <- bind_cols(proficiency, comfort) %>%
  kmeans(centers = 3)

clust %>% broom::tidy()

# plot clusters 
bind_cols(proficiency, comfort) %>%
  svd() %>%
  broom::tidy(matrix = 'u') %>%
  mutate(PC = paste('pc', PC, sep = '')) %>%
  pivot_wider(names_from = PC, values_from = value) %>%
  select(pc1, pc2) %>%
  mutate(clust = clust$cluster) %>%
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = factor(clust))) +
  labs(x = 'projection 1', 
       y = 'projection 2',
       color = 'cluster')
