library(tidyverse)
library(here)

## Fonts and themes
library(systemfonts)
clear_registry()

register_variant(
  name = "Myriad Pro SemiCondensed",
  family = "Myriad Pro",
  width = "semicondensed",
  weight = c("normal", "semibold"),
)

library(showtext)
showtext_opts(dpi = 300)
showtext_auto()

library(myriad)
import_myriad_semi()
import_myriad_condensed()
theme_set(theme_myriad_semi())


colors <- ggokabeito::palette_okabe_ito()
scales::show_col(colors)


## Setup

season_lab <-  tibble(yrday = yday(as.Date(c("2019-03-01",
                                             "2019-06-01",
                                             "2019-09-01",
                                             "2019-12-01"))),
                      lab = c("Spring", "Summer", "Autumn", "Winter"))

main_oceans <- c(
  "North Atlantic Ocean",
  "South Atlantic Ocean",
  "North Pacific Ocean",
  "South Pacific Ocean",
  "Indian Ocean")

### Read the dfs

# North Atlantic
natl_df <- read_csv(here("data", "natl_means.csv"))

# All the seas
seameans_df <- read_csv(here("data", "oceans_means_zonal.csv"))

# Global average 60S/60N
world_df <- read_csv(here("data", "global_means_60S60N.csv"))

### Graphs
month_labs <- seameans_df |>
  filter(sea == "North Atlantic Ocean",
         year == 2023,
         day == 15) |>
  select(date, year, yrday, month, day) |>
  mutate(month_lab = month(date, label = TRUE, abbr = TRUE))



out <- seameans_df |>
  filter(sea %in% main_oceans) |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years"
  ),
  sea_f = factor(sea, levels = main_oceans, ordered = TRUE)) |>
  filter(sea != "Indian Ocean") |>
  ggplot(aes(x = yrday, y = sst, group = year, color = year_flag)) +
  geom_line(linewidth = rel(0.5)) +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_color_manual(values = colors[c(1,6,2)]) +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 1.4))
  ) +
  facet_wrap(~ sea_f, axes = "all_x", axis.labels = "all_y") +
  labs(x = "Month of the Year", y = "Mean Temperature (Celsius)",
       color = "Year",
       title = "Mean Daily Sea Surface Temperatures, 1981-2024",
       subtitle = "Area-weighted 0.25° grid estimates; NOAA OISST v2.1; IHO Sea Boundaries",
       caption = "Data processed with R; Figure made with ggplot by Kieran Healy / @kjhealy") +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)),
        strip.text = element_text(face = "bold", size = rel(1.4)),
        plot.title = element_text(size = rel(1.525)),
        plot.subtitle = element_text(size = rel(1.1)))

ggsave(here("figures", "four_oceans.pdf"), out, width = 10, height = 10)
ggsave(here("figures", "four_oceans.png"), out, width = 10, height = 10, dpi = 300)


## All the world's oceans and seas
out <- seameans_df |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years")) |>
  ggplot(aes(x = yrday, y = sst, group = year, color = year_flag)) +
  geom_line(linewidth = rel(0.5)) +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_color_manual(values = colors[c(1,6,2)]) +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 1.4))
  ) +
  facet_wrap(~ reorder(sea, sst), axes = "all_x", axis.labels = "all_y") +
  labs(x = "Month of the Year", y = "Mean Temperature (Celsius)",
       color = "Year",
       title = "Mean Daily Sea Surface Temperatures, 1981-2024",
       subtitle = "Area-weighted 0.25° grid estimates; NOAA OISST v2.1; IHO Sea Boundaries",
       caption = "Data processed with R; Figure made with ggplot by Kieran Healy / @kjhealy") +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)),
        strip.text = element_text(face = "bold", size = rel(1.4)),
        plot.title = element_text(size = rel(1.525)),
        plot.subtitle = element_text(size = rel(1.1)))

ggsave(here("figures", "all_seas.pdf"), out, width = 40, height = 40)
ggsave(here("figures", "all_seas.png"), out, width = 40, height = 40, dpi = 300)



## North Atlantic only. Could also use natl_df here for cruder N Atl bounds.
out_atlantic <- seameans_df |>
  filter(sea == "North Atlantic Ocean") |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years"
  )) |>
  ggplot(aes(x = yrday, y = sst, group = year, color = year_flag)) +
  geom_line(linewidth = rel(1.1)) +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_color_manual(values = colors[c(1,6,2)]) +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 2))
  ) +
  labs(x = "Month", y = "Mean Temperature (Celsius)",
       color = "Year",
       title = "Mean Daily Sea Surface Temperature, North Atlantic Ocean, 1981-2024",
       subtitle = "Gridded and weighted NOAA OISST v2.1 estimates",
       caption = "Kieran Healy / @kjhealy") +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)))

ggsave(here("figures", "north_atlantic.png"), out_atlantic, height = 7, width = 10, dpi = 300)


## World Graph

world_avg <- world_df |>
  filter(year > 1981 & year < 2012) |>
  group_by(yrday) |>
  filter(yrday != 366) |>
  summarize(mean_8211 = mean(sst, na.rm = TRUE),
            sd_8211 = sd(sst, na.rm = TRUE)) |>
  mutate(fill = colors[2],
         color = colors[2])

out_world <- world_df |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years"))


out_world_plot <- ggplot() +
  geom_ribbon(data = world_avg,
              mapping = aes(x = yrday,
                            ymin = mean_8211 - 2*sd_8211,
                            ymax = mean_8211 + 2*sd_8211,
                            fill = fill),
              alpha = 0.3,
              inherit.aes = FALSE) +
  geom_line(data = world_avg,
            mapping = aes(x = yrday,
                          y = mean_8211,
                          color = color),
            linewidth = 2,
            inherit.aes = FALSE) +
  scale_color_identity(name = "Mean Temp. 1982-2011, ±2SD", guide = "legend",
                       breaks = unique(world_avg$color), labels = "") +
  scale_fill_identity(name = "Mean Temp. 1982-2011, ±2SD", guide = "legend",
                      breaks = unique(world_avg$fill), labels = "") +
  ggnewscale::new_scale_color() +
  geom_line(data = out_world,
            mapping = aes(x = yrday, y = sst, group = year, color = year_flag),
            inherit.aes = FALSE) +
  #  scale_color_manual(values = c("orange", "firebrick1", "grey50")) +
  scale_color_manual(values = colors[c(1,6,8)]) +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_y_continuous(breaks = seq(19.5, 21.5, 0.5),
                     limits = c(19.5, 21.5),
                     expand = expansion(mult = c(-0.05, 0.05))) +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 2))
  ) +
  labs(x = "Month", y = "Mean Temperature (°Celsius)",
       color = "Year",
       title = "Mean Daily Global Sea Surface Temperature, 1981-2024",
       subtitle = "Latitudes 60°N to 60°S; Area-weighted NOAA OISST v2.1 estimates",
       caption = "Kieran Healy / @kjhealy") +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)),
        plot.title = element_text(size = rel(1.9)))

ggsave(here("figures", "global_mean.png"), out_world_plot, height = 7, width = 10, dpi = 300)




out_world_plot_kelvin <- ggplot() +
  geom_line(data = world_avg,
            mapping = aes(x = yrday,
                          y = mean_8211 + 273.15,
                          color = color),
            linewidth = 2,
            inherit.aes = FALSE) +
  scale_color_identity(name = "Mean Temp. 1982-2011, ±2SD", guide = "legend",
                       breaks = unique(world_avg$color), labels = "") +
  scale_fill_identity(name = "Mean Temp. 1982-2011, ±2SD", guide = "legend",
                      breaks = unique(world_avg$fill), labels = "") +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_y_continuous(breaks = c(0, 75, 175, 275, 375),
                     labels = c("0", "75", "175", "275", "375"),
                     limits = c(0, 375),
                     expand = c(0,0)) +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE),
    color = guide_legend(override.aes = list(linewidth = 2))
  ) +
  labs(x = "Month", y = "Mean Temperature (°Kelvin)",
       color = "Year",
       title = "Mean Daily Global Sea Surface Temperature, 1981-2024",
       subtitle = "Latitudes 60°N to 60°S; Area-weighted NOAA OISST v2.1 estimates",
       caption = "Kieran Healy / @kjhealy") +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)),
        plot.title = element_text(size = rel(1.9)))


ggsave(here("figures", "global_mean_kelvin_zero.png"),
       out_world_plot_kelvin, height = 7, width = 10, dpi = 300)


