# Change point analysis

```r
library(tidyverse)
library(bcp)

# The neo-W chromosome has poorly mapping regions. Consequently, when computing neo-W neo-Z divergence we removed the 20% of windows with the lowest neo-W coverage, as the neo-W was ∼ 80% covered.

# Find the coverage threshold to filter the 80% of windows with highest coverage. 

meta <- read.csv("metadata.csv") %>% 
  select(id, class)

median_coverage <- read_csv("neosex.coverage.50kb.csv", show_col_types = FALSE) %>%
  left_join(meta, by = "id") %>%
  filter(class == "Female")
  group_by(mid) %>%
  summarise(median = median(median))

threshold_df <- tibble(threshold = numeric(), filtered = numeric())

n <- 1000

for (i in seq(n)) {
  
  threshold <- i / n
  filtered <- median_coverage %>% filter(median <= threshold) %>% nrow() / nrow(d)
  threshold_df <- threshold_df %>% rbind(tibble(threshold = threshold, filtered = filtered))
  
}

# Threshold = 0.84

covered_wins <- median_coverage %>%
  filter(median > 0.84) %>% pull(mid)

df <- read_csv("neosex.iberian.female_het.50kb.csv", show_col_types = FALSE)

windowed_divergence <- df %>%
  filter(mid %in% covered_wins)

# change point analysis with bcp
res <- bcp(y = windowed_divergence$divergence, burnin = 500, mcmc = 5000)

# plotting the results (Figure 4)

windowed_divergence %>%
    cbind(res$posterior.mean, res$posterior.prob) %>%
    rename(mean = X1, post = `res$posterior.prob`) %>%
    ggplot() +
    geom_point(aes(mid / 1e6, divergence), size = 1.2) +
    geom_line(aes(mid / 1e6, mean, col='red')) +
    labs(x = "Position (Mb)", y = "Divergence", title = "Neo-sex chromosome") +
    scale_x_continuous(breaks = seq(0, 12, 2)) +
    ylim(c(0, 0.035)) +
    annotate("segment", x = 0, xend = 10.45, y = 0.0025, yend = 0.0025, colour = "black", linewidth = 1) +
    annotate("text", x = 5.5, y = 0.004, label = "Plateau 1", size = 4) +
    annotate("segment", x = 10.55, xend = 12, y = 0.0025, yend = 0.0025, colour = "black", linewidth = 1) +
    annotate("text", x = 11.25, y = 0.004, label = "Plateau 2", size = 4) +
    annotate("segment", x = 12.0925, xend = 13, y = 0.0025, yend = 0.0025, colour = "black", linewidth = 1) +
    annotate("text", x = 12.6, y = 0.0012, label = "Plateau 3", size = 4) +
    theme_light() + 
    theme(legend.position = "None", plot.title = element_text(hjust = 0.5), 
          axis.title = element_text(size = 14))
```