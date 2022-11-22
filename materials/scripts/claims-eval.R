library(tidyverse)

github_url <- "https://github.com/pstat197/pstat197a/raw/main/materials/slides/data/f22-claims-evals.RData"
load(url(github_url))

evals %>%
  rowwise() %>%
  filter(!is_tibble(eval))

accuracies <- evals %>%
  rowwise() %>%
  filter(is_tibble(eval)) %>%
  ungroup() %>%
  unnest(eval)

accuracies %>%
  ggplot(aes(x = factor(group), y = .estimate)) +
  geom_point(aes(alpha = n)) +
  facet_wrap(~class*.metric, nrow = 2)
