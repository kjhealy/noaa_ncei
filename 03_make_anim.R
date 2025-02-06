library(tidyverse)
library(here)
library(gganimate)
library(transformr)

library(kjhmisc)
kjhmisc::setup_socviz()

colors <- ggokabeito::palette_okabe_ito()
scales::show_col(colors)

year_seq_colors <- colors[c(1,6,9,2)]


theme_figs <- function(){
  theme_socviz_semi() +
    theme(
      plot.background = element_rect(color = "white"),
      plot.title = element_text(size = rel(2)),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      axis.text.x = element_text(size = rel(1.5)),
      axis.text.y = element_text(size = rel(1.5)))
}

theme_set(theme_figs())

## Setup

season_lab <-  tibble(yrday = yday(as.Date(c("2019-03-01",
                                             "2019-06-01",
                                             "2019-09-01",
                                             "2019-12-01"))),
                      lab = c("Spring", "Summer", "Autumn", "Winter"))

### Read the dfs

# All the seas
seameans_df <- read_csv(here("data", "oceans_means_zonal.csv"))

# Global average 60S/60N
world_df <- read_csv(here("data", "global_means_60S60N.csv"))

### For labels
month_labs <- seameans_df |>
  filter(sea == "North Atlantic Ocean",
         year == 2023,
         day == 15) |>
  select(date, year, yrday, month, day) |>
  mutate(month_lab = month(date, label = TRUE, abbr = TRUE))

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
    year == 2025 ~ "2025",
    .default = "All other years"),
    days_elapsed = as.integer(date - first(date)))


out_world_mov <- out_world |>
  ggplot() +
  geom_line(mapping = aes(x = yrday, y = sst, group = year, color = year_flag),
            linewidth = 1.2) +
  ggrepel::geom_text_repel(aes(x = 320,
                      y = 21.5,
                      label = year),
                  size = 25,
                  max.overlaps = 50,
                  max.time = 0.001,
                  max.iter = 1,
                  seed = 1234,
                  force = 0,
                  force_pull = 0,
                  segment.color = "transparent",
                  family = "Socviz Condensed",
                  color = colors[3],
                  alpha = 0.6,
                  show.legend = FALSE) +
  scale_color_manual(values = year_seq_colors) +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_y_continuous(breaks = seq(19.5, 21.5, 0.5),
                     limits = c(19.5, 21.5),
                     expand = expansion(mult = c(-0.05, 0.05))) +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 4))
  ) +
  labs(x = "Month", y = "Mean Temperature (°Celsius)",
       color = "Year",
       title = "Mean Daily Global Sea Surface Temperature, {out_world$year[which.min(abs(out_world$days_elapsed-frame_along))]}",
       subtitle = "Latitudes 60°N to 60°S; Area-weighted NOAA OISST v2.1 estimates",
       caption = "Kieran Healy / kieranhealy.org") +
  transition_reveal(days_elapsed) +
  theme(axis.line = element_line(color = "gray30", linewidth = rel(1)),
        legend.text = element_text(size = rel(1.75)),
        legend.title = element_text(size = rel(1.75)),
        plot.title = element_text(size = rel(2.5), face = "bold"))

animate(out_world_mov, fps = 30,
        duration = 60, width = 3840/3, height = 2160/3,
        renderer = av_renderer("figures/sst_anim_1min_1280x720.mp4"))


