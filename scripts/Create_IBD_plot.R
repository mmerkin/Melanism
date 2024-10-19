library(tidyverse)

Pair <- read.table("SPECIESID_pairwise_comparisons.txt")$V1
Mean_FST <- read.table("SPECIESID_mean_fsts.txt")$V1
Weighted_FST <- read.table("SPECIESID_weighted_fsts.txt")$V1
Bb_dist <- read.csv("SPECIESID_filtered_distances.csv")
Bb_dist$extra <- 1


dist <- dist %>% full_join(dist,c("extra"="extra"), relationship = "many-to-many") %>% 
  filter(Site.x != Site.y) %>% 
  mutate(distance=sqrt((E.x-E.y)^2 + (N.x-N.y)^2))

Distance <- dist$distance

df <- tibble(Pair, Distance, Mean_FST, Weighted_FST) %>% 
  distinct(Distance, .keep_all= TRUE)

df$Mean_FST[df$Mean_FST<0] <- 0

ggplot(data=df, aes(x=log(Distance, base=exp(1)), y=Mean_FST/(1-Mean_FST))) + geom_point() +
  geom_smooth(se=T)
