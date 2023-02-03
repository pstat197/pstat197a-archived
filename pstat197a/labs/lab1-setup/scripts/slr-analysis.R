library(tidyverse)

# load data
pollution <- read_csv('data/pollution.csv')

# examine scatterplot with SLR fit
ggplot(pollution,
       aes(x = log(SO2), y = Mort)) +
  geom_point() +
  geom_smooth(method = 'lm')

# compute SLR fit
fit <- lm(Mort ~ log(SO2), data = pollution)

# confidence interval
fit_ci <- confint(fit, parm = 'log(SO2)')*log(1.2)

# interpretation
paste('With 95% confidence, every 20% increase in sulfur dioxide pollution is associated with an increase in two-year mortality rate between', round(fit_ci[1], 2), 'and', round(fit_ci[2], 2), 'per 100k', sep = ' ') %>% print()
