library(tidyverse)
library(lubridate)
library(sf)
library(ggspatial)
library(ggmap)
setwd("~/pstat197/pstat197a/materials/slides")

sites <- read_csv(unz("data/soil-temp-data.zip", 
                      "data/USArray_Sites.csv"))
         
sites_sf <- sites %>%
  select(site, latitude, longitude, elevation, 
         start_time, end_time) %>%
  rename(lon = longitude, 
         lat = latitude) %>%
  st_as_sf(coords = c('lon', 'lat'))

box <- st_bbox(sites_sf) + 0.5*c(-1, -1, 1, 1)
names(box) <- c('left', 'bottom', 'right', 'top')

map <- get_map(location = box, 
        maptype = 'terrain',
        source = 'stamen')

fig_sitemap <- ggmap(map) +
  geom_point(data = sites, 
             aes(x = longitude,
                 y = latitude))

fig_sitemap


files <- sites_sf %>%
  mutate(filename = paste(site, 
                          start_time, 
                          end_time, 
                          sep = '_')) %>%
  pull(filename) %>%
  paste('.csv', sep = '')

data_list <- lapply(1:64, function(i){
  tryCatch({
      read_csv(unz("data/soil-temp-data.zip", 
               paste("data/", 
                     files[i], 
                     sep = '')),
              col_types = cols(.default = 'c')) %>%
        pivot_longer(contains('tsoil'),
                     names_to = 'depth',
                     values_to = 'temp') %>%
        separate(depth, 
                 into = c('discard', 'depth'),
                 sep = '_') %>%
        select(-timezone, -discard) %>%
        mutate(date_time = ymd_hm(date_time),
               across(c(depth, temp), as.numeric))
      },
    error = function(cond){return(NA)}
  )
})

merged_data <- sites %>%
  select(site, latitude, longitude, elevation) %>%
  bind_cols(tibble(data = data_list))


read_failures <- merged_data %>%
  filter(is.na(data))

data_clean <- merged_data %>%
  filter(!is.na(data)) %>%
  unnest(data) %>%
  na.omit() %>%
  filter(temp > -100)

data_clean %>% count(year(date_time))

obs_counts <- data_clean %>%
  filter(year(date_time) > 2016,
         year(date_time) < 2020) %>%
  group_by(site, depth) %>%
  count()

obs_filter <- obs_counts %>%
  filter(n > 365*4) %>%
  mutate(grp = interaction(site, depth))

obs_filter %>% pull(site) %>% unique()

soil <- data_clean %>%
  filter(interaction(site, depth) %in% 
           obs_filter$grp,
         year(date_time) > 2016,
         year(date_time) < 2020)

save(soil, file = 'data/soiltemp-clean.RData')

soil %>%
  filter(site == unique(soil$site)[20]) %>%
  ggplot(aes(x = date_time, y = temp)) +
  geom_path(aes(group = depth, color = depth), 
            alpha = 0.6) +
  labs(x = '')

soil %>% 
  distinct(site, depth) %>% 
  count(depth) %>% 
  filter(n > 20)

soil_200cm <- soil %>% 
  filter(depth == 0.2) %>%
  mutate(year = year(date_time),
         day = yday(date_time),
         date = date(date_time)) %>%
  group_by(site, year, day, date) %>%
  summarize(elev = mean(elevation),
            temp = mean(temp)) %>%
  ungroup()

exclude <- soil_200cm %>%
  arrange(site, date) %>%
  group_by(site) %>%
  summarize(maxdiff = max(diff(day))) %>%
  filter(maxdiff > 1) %>%
  select(site) %>%
  ungroup()

soil_200cm_clean <- anti_join(soil_200cm, exclude) %>% 
  filter(elev > -100)

fig_annual <- soil_200cm_clean %>%
  ggplot(aes(x = day, y = temp)) +
  geom_path(aes(group = interaction(site, year)),
            alpha = 0.1) +
  labs(x = 'day of year',
       y = 'average temp')

fig_annual

knotpt <- quantile(soil_200cm_clean$day, 
                   c(0.2, 0.4, 0.6, 0.8))

fig_annual +
  geom_vline(xintercept = c(0, 365, knotpt))

fit <- lm(temp ~ bs(day, 
                    knots = knotpt, 
                    degree = 3),
          data = soil_200cm_clean)

pred_grid <- soil_200cm_clean %>%
  ungroup() %>%
  data_grid(day = seq_range(day, 200)) 

pred_df <- predict(fit, 
                   pred_grid, 
                   interval = 'prediction',
                   level = 0.95) %>%
  bind_cols(pred_grid)

fig_annual +
  geom_path(data = pred_df, 
            aes(y = fit, 
                x = day),
            inherit.aes = F,
            color = 'blue') +
  geom_ribbon(data = pred_df,
              aes(ymin = lwr,
                  ymax = upr, 
                  x = day),
              alpha = 0.1,
              inherit.aes = F)

pred_grid_multiyear <- soil_200cm_clean %>%
  data_grid(date = seq_range(date, n = 500))

pred_df_multiyear <- pred_grid_multiyear %>%
  mutate(day = yday(date)) %>%
  add_predictions(fit)

soil_200cm_clean %>%
  ungroup() %>%
  ggplot(aes(x = date, y = temp)) +
  geom_path(aes(group = site),
            alpha = 0.1) +
  geom_path(data = pred_df_multiyear,
            aes(y = pred, 
                group = year(date)),
            color = 'blue')

fit_f <- lm(temp ~ fourier(day, 
                           nbasis = 4, 
                           period = 365) - 1,
            data = soil_200cm_clean)

fpred_df <- predict(fit_f, 
                    pred_grid, 
                    interval = 'prediction',
                    level = 0.95) %>%
  bind_cols(pred_grid)

fig_annual +
  geom_path(data = fpred_df, 
            aes(y = fit, 
                x = day),
            inherit.aes = F,
            color = 'blue') +
  geom_ribbon(data = fpred_df,
              aes(ymin = lwr,
                  ymax = upr, 
                  x = day),
              alpha = 0.1,
              inherit.aes = F)

soil_200cm_clean %>%
  ungroup() %>%
  ggplot(aes(x = date, y = temp)) +
  geom_path(aes(group = site),
            alpha = 0.1) +
  geom_path(data = fpred_df_multiyear,
            aes(y = pred),
            color = 'blue')

fit_df <- soil_200cm_clean %>%
  add_residuals(fit_f)

fit_df %>%
  ggplot(aes(x = date, y = resid)) +
  geom_path(aes(group = site), 
            alpha = 0.1)

fit_elev <- lm(resid ~ elev, data = fit_df)

fit_df %>%
  add_residuals(fit_elev, var = 'resid2') %>%
  ggplot(aes(x = date, y = resid2)) +
  geom_path(aes(group = site), 
            alpha = 0.1)

fit_df %>%
  add_residuals(fit_elev, var = 'resid2') %>%
  arrange(site, date) %>%
  mutate(resid2.lag = lag(resid2, n = 1)) %>%
  mutate(diff = resid2 - resid2.lag) %>%
  na.omit() %>%
  ggplot(aes(x = date, y = diff)) +
  geom_path(aes(group = site),
            alpha = 0.1)

fit_df %>%
  add_residuals(fit_elev, var = 'resid2') %>%
  arrange(site, date) %>%
  mutate(resid2.lag = lag(resid2, n = 1)) %>%
  mutate(diff = resid2 - resid2.lag) %>%
  na.omit() %>%
  ggplot(aes(x = lag(diff, 1), y = diff)) +
  geom_point(alpha = 0.1)

fit_df %>%
  add_residuals(fit_elev, var = 'resid2') %>%
  arrange(site, date) %>%
  mutate(resid2.lag = lag(resid2, n = 1)) %>%
  ggplot(aes(y = resid2, x = resid2.lag)) +
  geom_point(alpha = 0.2)
