library(ggplot2)
library(ggrepel)
library(ggtext)
library(dplyr)
library(forcats)
library(rcartocolor)
library(sf)
library(arcgis)
library(rnaturalearth)
library(rnaturalearthhires)

theme_set(theme_minimal(base_size = 10,
  base_family = "Arial"))

theme_update(
plot.caption.position = "plot",
# change color of caption
plot.caption = element_markdown(color = "grey30"),
panel.grid.minor = element_blank() # no minor grid lines
)

akhs_stock_url <- "https://services2.arcgis.com/C8EMgrsFcRFL6LrL/arcgis/rest/services/pv_dist/FeatureServer"
akhs_stocks <- arc_open(akhs_stock_url) |> 
  get_layer(id=0) |> 
  arc_select(fields = c("stock_num","stockname"))

akhs_stocks <- akhs_stocks |> 
  sf::st_transform(3338)

akhs_stocks <- akhs_stocks |> 
  mutate(stockname = forcats::as_factor(stockname),
         stockname = forcats::fct_reorder(stockname,stock_num)
        )

alaska_base <-
  rnaturalearth::ne_states("United States of America", return = "sf") |>
  dplyr::filter(name == "Alaska") |>
  sf::st_transform(3338)
russia_base <-
  rnaturalearth::ne_states("Russia", return = "sf")  |>
  sf::st_transform(3338)
canada_base <- rnaturalearth::ne_states("Canada", return = "sf") |>
  sf::st_transform(3338)

akhs_stock_colors <- c(rcartocolor::carto_pal(12, 'Bold')[1:11],
                    rcartocolor::carto_pal(2, 'Pastel')[1])

map_labels <- tibble(
  label = c("Bering \nSea","Gulf \nof Alaska"),
  longitude = c(-178,-145),
  latitude = c(57,54)
) |> 
  st_as_sf(coords = c('longitude','latitude'),
           crs = 4326) |> 
  st_transform(3338)

ggplot() +
  geom_sf(
    data = akhs_stocks,
    aes(fill = stockname),
    linewidth = 0
  ) +
  scale_fill_manual(values = akhs_stock_colors) +
  geom_sf(
    data = alaska_base,
    fill = "gray70",
    linewidth = 0.1
  ) +
  geom_sf_text(data = alaska_base, aes(label = name)) +
  geom_sf(
    data = russia_base,
    fill = "gray85",
    linewidth = 0
  ) +
  geom_sf_text(
    data = russia_base |> filter(gn_id == 2126099),
    aes(label = toupper(admin)),
    color = "gray50"
  ) +
  geom_sf(
    data = canada_base,
    fill = "gray85",
    linewidth = 0
  ) +
  geom_sf_text(
    data = canada_base |> filter(gn_id == 6185811),
    aes(label = toupper(admin)),
    color = "gray50"
  ) +
  geom_sf_text(
    data = map_labels,
    aes(label = label),
    color = "gray70",
    fontface = "italic"
  ) +
    geom_label_repel(
      data = akhs_stocks,
      aes(label = stockname, geometry = geometry, fill = stockname),
      color = "white",
      stat = "sf_coordinates",
      max.overlaps = 100
    ) +
  coord_sf(
    xlim = c(-2.25e+06, 1.75e+06),
    ylim = c(0.1e+06, 2.6e+06),
    expand = FALSE
  ) +
  scale_x_continuous(breaks = c(180, -170, -160, -150, -140)) +
  guides(fill = guide_legend(ncol = 4, nrow = 3,
          override.aes = aes(color = NA))) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_blank()
  )

ggsave(here::here("figures/png/akhs_stock_map.png"),
  device = agg_png,
  width = 6.5, height = 4.0,
  units = "in",scale =2,
  bg="white",
  create.dir = TRUE
)
