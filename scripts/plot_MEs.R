library(tidyverse)

# Generate test data

Mfactor <- factor(c(paste0("M", 1:31), "MZ"), levels =c("MZ", paste0("M", 31:1)))
Gfactor <- factor(paste0("group", 1:31))
n_samples <- 100
randM <- sample(Mfactor, n_samples, replace = TRUE)
randG <- sample(Gfactor, n_samples, replace = TRUE)
data <- tibble(ME=randM, chr=randG)

# Alternatively, read in table from bash companion script
#data <- read_tsv("O_bidentata_merian_elements.tsv")
#data <- read_tsv("P_pilosaria_merian_elements.tsv")
#data <- read_tsv("E_epiphron_merian_elements.tsv")


# Calculate fraction of each ME comprising each chromosome

df_summary <- data %>%
  group_by(chr, ME) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(chr) %>%
  mutate(fraction = count / sum(count)) %>%
  ungroup()

# Identify chromosomes that are >50% of a single ME and sort them starting from M1

chr_order <- df_summary %>%
  filter(fraction > 0.5) %>%
  arrange(factor(ME, levels = c(paste0("M", 1:31), "MZ"))) %>%
  pull(chr) # Return just the chr column rather than the entire df

# Add in any missing chromosomes that were <50% of any ME
chr_order <- unique(c(chr_order, df_summary$chr))
# Reorder the chromosomes and MEs to start with the highest fraction of M1 in the top left
df_summary$chr <- factor(df_summary$chr, levels = chr_order)
df_summary$ME <- factor(df_summary$ME, levels = c("MZ", paste0("M", 31:1)))

# Create a heatmap
ggplot(df_summary, aes(x = chr, y = ME, fill = fraction)) +
  geom_tile(height = 0.8, width = 0.8) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(y = "Merian element", x = "Chromosome", fill = "Fraction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  

