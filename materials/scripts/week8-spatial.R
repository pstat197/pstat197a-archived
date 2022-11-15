library(tidyverse)
library(lubridate)
library(forecast)
library(sf)
library(sp)
library(gstat)
library(fda)
library(ggspatial)
library(ggmap)
setwd("~/pstat197/pstat197a/materials/slides")
soil <- read_csv('data/soiltemp-200cm.csv')
sites <- read_csv(unz("data/soil-temp-data.zip", 
                      "data/USArray_Sites.csv"))

site_df <- sites %>%
  dplyr::select(site, longitude, latitude, elevation)


soil %>% 
  group_by(site) %>% 
  summarize(start = min(date), 
            stop = max(date),
            step = max(diff(date))) %>%
  arrange(stop)

soil %>%
  filter(date < ymd('2018-05-01')) %>%
  group_by(site) %>%
  count() %>%
  arrange(n)

fit_fn <- function(.x, .y){
  out <- forecast::Arima(y = .y, 
               order = c(2, 0, 0), 
               xreg = .x, 
               include.mean = F, 
               method = 'ML')
  return(out)
}

pred_fn <- function(.fit, .reg){
  out <- forecast::forecast(.fit, h = nrow(.reg), xreg = .reg)
  return(out)
}

fit_df <- soil %>% 
  dplyr::select(-year, -elev) %>%
  filter(!str_starts(site, 'SHA')) %>%
  arrange(date) %>%
  nest(data = c(day, date, temp)) %>%
  mutate(train = map(data, ~filter(.x, date < ymd('2018-05-01'))),
         test = map(data, ~filter(.x, date >= ymd('2018-05-01'))),
         x = map(data, ~fourier(.x$day, nbasis = 4, period = 365)),
         y = map(data, ~pull(.x, temp)),
         fit = map2(x, y, fit_fn),
         xtest = map(test, ~fourier(.x$day, nbasis = 4, period = 365)),
         pred = map2(fit, xtest, pred_fn))
         

pred_df <- fit_df %>%
  mutate(y.pred = map(pred, ~.x$mean)) %>%
  dplyr::select(site, y.pred, test) %>%
  unnest(everything()) %>%
  left_join(site_df, by = 'site')

pred_df %>%
  ggplot(aes(x = date, y = temp, group = site)) +
  geom_path() +
  geom_path(aes(y = y.pred), color = 'blue') +
  facet_wrap(~site)



preds_sf <- st_as_sf(pred_df, 
                       coords = c('longitude', 'latitude'),
                       crs = st_crs(4326)) 

grid_fn <- function(sf, w, h){
  
  # determine boundary (convex hull of sampled locations)
  boundary <- sf %>% 
    distinct(site, .keep_all = T) %>%
    st_combine() %>% 
    st_convex_hull()
  
  # partition region within boundary into boxes
  grid_geo <- boundary %>%
    st_make_grid(n = c(w, h),
                 what = 'polygons',
                 crs = 4326) %>%
    st_intersection(boundary)
  
  # add box id (modelr::map fails otherwise due to length attribute)
  out <- st_sf(tibble(.id = 1:length(grid_geo)), grid_geo)
  
  return(out)
}

preds_1step <- preds_sf %>% filter(date == ymd('2018-05-01'))
grid <- grid_fn(preds_sf, 25, 25)

spatial_preds <- lapply(c(0.5, 0.8, 1),
                        function(pwr){
                          idw(y.pred ~ 1, preds_1step, grid, idp = pwr) %>%
                            mutate(idw.power = paste('power ', pwr, sep = ''))
})

spatial_pred_df <- Reduce(bind_rows, spatial_preds)

box <- st_bbox(preds_sf) + 0.5*c(-1, -1, 1, 1)
names(box) <- c('left', 'bottom', 'right', 'top')
map <- get_map(location = box, 
               maptype = 'terrain',
               source = 'stamen')

fig_sitemap <- ggmap(map) +
  geom_point(data = pred_df, 
             aes(x = longitude,
                 y = latitude))

fig_sitemap + 
  layer_spatial(spatial_pred_df, 
                aes(fill = var1.pred),
                alpha = 0.8) +
  facet_wrap(~idw.power, nrow = 1) +
  scale_fill_distiller(palette = 'YlGnBu') +
  guides(fill = guide_colorbar('soil temp'))
