library(here)
load(here::here('data/akpv_datacube.rda'))
load(here::here('data/dGlac.rda'))
load(here::here('data/dTerr.rda'))

dTerr = dTerr[dTerr$yr > 1995,]
dTerr$polyid = as.factor(as.character(dTerr$polyid))
dTerr$stockid = as.factor(as.character(dTerr$stockid))
dGlac = dGlac[dGlac$yr > 1995,]
dGlac$polyid = as.factor(as.character(dGlac$polyid))
dGlac$stockid = as.factor(as.character(dGlac$stockid))
# rename column so consistent with terrestrial datasets
names(dGlac)[5] = 'count'
# remove records with missing count values
dGlac = dGlac[!is.na(dGlac$count),]
dstk = rbind(dTerr[,c('polyid', 'stockid', 'yr')],
	dGlac[,c('polyid', 'stockid', 'yr')])

library(tidyverse)

get_total_abundance_tidy <- function(datacube) {
    years <- colnames(datacube[[1]])
    stocks <- unique(attr(datacube[[1]], 'stockid'))
    mcmc_samples <- length(datacube)
    
    result <- map_dfr(1:mcmc_samples, function(i) {
        map_dfr(years, function(year) {
            map_dfr(stocks, function(stock) {
                tibble(
                    mcmc_sample = i,
                    year = year,
                    stock = stock,
                    total_abundance = sum(datacube[[i]][attr(datacube[[i]], 'stockid') == stock, year])
                )
            })
        })
    })
    
    return(result)
}

total_abundance_df <- get_total_abundance_tidy(akpv_datacube)
library(ggplot2)

# Calculate median and credible intervals
summary_df <- total_abundance_df %>%
  group_by(year, stock) %>%
  summarise(
    mean_abundance = mean(total_abundance),
    lower_ci = quantile(total_abundance, 0.025),
    upper_ci = quantile(total_abundance, 0.975)
  )

# Create the plot
ggplot(summary_df, aes(x = as.numeric(year), y = mean_abundance, group = stock)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
  facet_wrap(~ stock, scales = "free_y") +
  labs(
    title = "Mean Abundance Estimates with 95% Credible Intervals",
    x = "Year",
    y = "Mean Abundance"
  ) +
geom_point(color = "black")

# Calculate the fraction of seals surveyed each year
calculate_fraction_sampled <- function(datacube, dstk) {
    years <- colnames(datacube[[1]])
    stocks <- unique(attr(datacube[[1]], 'stockid'))
    
    result <- map_dfr(stocks, function(stock_id) {
        fraci <- map_dbl(1:length(years), function(i) {
            ind <- rownames(datacube[[1]]) %in% 
                levels(as.factor(as.character(dstk[dstk$stockid == stock_id & dstk$yr == 1995 + i, 'polyid'])))
            sum(datacube[[1]][ind, i]) / sum(datacube[[1]][attr(datacube[[1]], 'stockid') == stock_id, i])
        })
        
        tibble(
            year = years,
            stock = stock_id,
            fraction_sampled = fraci
        )
    })
    
    return(result)
}

fraction_sampled_df <- calculate_fraction_sampled(akpv_datacube, dstk)

# Create the bar plot
library(patchwork)

# Define a mapping of stock IDs to stock names
stock_names <- c(
  "1" = "Aleutian Islands",
  "2" = "Pribilof Islands",
  "3" = "Bristol Bay",
  "4" = "North Kodiak",
  "5" = "South Kodiak",
  "6" = "Prince William Sound",
  "7" = "Cook Inlet/Shelikof Strait",
  "8" = "Glacier Bay/Icy Strait",
  "9" = "Lynn Canal/Stephens Passage",
  "10" = "Sitka/Chatham Strait",
  "11" = "Dixon/Cape Decision",
  "12" = "Clarence Strait"
  # Add more mappings as needed
)

plots <- lapply(unique(summary_df$stock), function(stock_id) {
  stock_name <- stock_names[as.character(stock_id)]
  
  abundance_plot <- ggplot(subset(summary_df, stock == stock_id), aes(x = as.numeric(year), y = mean_abundance)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
    geom_vline(xintercept = c(2000, 2005, 2010, 2015, 2020), linetype = "dashed", color = "grey") +
    labs(
      title = stock_name,
      y = "Abundance"
    ) +
    geom_point(color = "black") +
    scale_x_continuous(limits = c(1996, NA), breaks = c(2000, 2005, 2010, 2015, 2020)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.title.x = element_blank()
    )
  
  fraction_sampled_plot <- ggplot(subset(fraction_sampled_df, stock == stock_id), aes(x = as.numeric(year), y = fraction_sampled)) +
    geom_bar(stat = "identity", fill = "grey", color = "black", linewidth = 0.2) +
    geom_vline(xintercept = c(2000, 2005, 2010, 2015, 2020), linetype = "dashed", color = "grey") +
    labs(
      x = "Year",
      y = "Effort"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(limits = c(1996, NA), breaks = c(2000, 2005, 2010, 2015, 2020)) +
    scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.position = "none"
    )
  
  combined_plot <- abundance_plot / fraction_sampled_plot + plot_layout(heights = c(2, 1))
  return(combined_plot)
})

# Combine all individual plots into a single plot
final_plot <- wrap_plots(plots, ncol = 3)

# Save the final combined plot with specified dimensions
ggsave("figures/png/final_combined_plot.png", plot = final_plot, 
       width = 6.35, height = 7, dpi = 300,
       units = "in", scale = 2)

# Display the final combined plot
print(final_plot)
