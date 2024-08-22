## code to prepare the topical cyclone maximum intensity series

library(dplyr)
cycle_all <- read.csv("ibtracs.ALL.list.v04r00.csv")

cyclone <- as_tibble(cycle_all) %>%
  select(
    SEASON,
    USA_WIND,
    BASIN,
    SID
  ) %>%
  mutate(
    SEASON = as.numeric(SEASON),
    USA_WIND = as.numeric(USA_WIND),
    BASIN = ifelse(is.na(BASIN), "NHA", BASIN)
  ) %>%
  filter(
    !is.na(USA_WIND),
    NATURE == "TS",
  ) %>%
  group_by(SID) %>%
  summarise(
    max_wind = max(USA_WIND, na.rm = T),
    BASIN = unique(BASIN),
    SEASON = unique(SEASON)
  ) %>%
  ungroup()

usethis::use_data(cyclone, overwrite = TRUE)
