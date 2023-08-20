## Load Packages
library(plyr)
library(tidyverse) 
library(lubridate)
library(vegan)
library(ggrepel)
library(ggpubr)
library(nlme)
library(MuMIn)
library(sjPlot)
library(multcomp)
library(agricolae)
library(car)
library(Rmisc)

#set working directory
path = "C:\\Users\\Owner\\Documents\\Surf Zone\\Koval_etal_2023_SurfZoneFish"
setwd(path)

###Load Data####################################################################
#read in species data
bruv_month <- read_csv(".\\data\\maxn_period_seasonal_month.csv")
bruv_month$Date <- as.Date(bruv_month$Date, tryFormats = c("%m/%d/%Y"))

#read in interpolated data
bruv_environment_inter <- read_csv(".\\data\\bruv_meta_seasonal_interpolated_filled_indices.csv")
bruv_environment_inter$Date <- as.Date(bruv_environment_inter$Date, tryFormats = c("%m/%d/%Y"))

#create vector for not in fuction
`%notin%` <- Negate(`%in%`)

#filter out mammals and unidentified species
bruv_month_fish <- bruv_month %>% 
  dplyr::filter(Type == "Fish" & 
                  Sci_Name_Abbr %notin% c("NA. spp") & 
                  Common_Name %notin% c("No Species"))

#count species and filter fish species observed less than 10 times
species_10_fish <- plyr::count(bruv_month_fish, vars = "Sci_Name_Full")
species_10_fish <- dplyr::filter(species_10_fish, freq >= 10)
species_10_list_fish <- c(species_10_fish$Sci_Name_Full)

bruv_month_fish <- subset(bruv_month_fish, Sci_Name_Full %in% species_10_list_fish)

#calculate average max N by species per BRUV
bruv_month_subset_fish <- bruv_month_fish %>% 
  dplyr::group_by(Site, Date, BRUV, Analysis_Month, Sci_Name_Abbr) %>%
  dplyr::summarise(MaxN = max(MaxN))

#calculate maxN per site
bruv_month_subset_fish <- bruv_month_subset_fish %>% 
  dplyr::group_by(Site, Date, Analysis_Month, Sci_Name_Abbr) %>% 
  dplyr::summarise(MaxN = sum(MaxN)/6) %>% 
  print()

#convert to wide for spp table
bruv_month_wide_fish <- pivot_wider(bruv_month_subset_fish, names_from = Sci_Name_Abbr, values_from = MaxN)

#add in missing sampling days
bruv_month_wide_fish <- ungroup(bruv_month_wide_fish)
bruv_month_wide_fish <- bruv_month_wide_fish %>% 
  add_row(Site = "Carmel", Date = as.Date("2021-05-23"), Analysis_Month = "May") %>% 
  add_row(Site = "Spanish Bay", Date = as.Date("2020-08-15"), Analysis_Month = "Aug") %>% 
  add_row(Site = "Spanish Bay", Date = as.Date("2021-06-11"), Analysis_Month = "Jun")

#Fill in NA with 0
bruv_month_wide_fish[is.na(bruv_month_wide_fish)] <- 0

#summarize environmental data
bruv_environment_subset <- bruv_environment_inter %>% 
  dplyr::group_by(Site, Date, Analysis_Month, Season) %>% 
  dplyr::summarize(SkyConditions = unique(SkyConditions),
                   TideHeight = mean(TideHeight),
                   WindAve = mean(WindAve),
                   Salinity = mean(Salinity),
                   Depth = mean(Depth),
                   Visibility = mean(Visibility),
                   Temp = mean(Temp),
                   DistShore = mean(DistShore),
                   WaveHeight = mean(WaveHeight),
                   WavePeriod = mean(WavePeriod),
                   Percent_Cover = mean(Percent_Cover),
                   Red_Percent_Cover = mean(Red_Percent_Cover),
                   Brown_Percent_Cover = mean(Brown_Percent_Cover),
                   Green_Percent_Cover = mean(Green_Percent_Cover),
                   Seagrass_Percent_Cover = mean(Seagrass_Percent_Cover),
                   WDIR_ave = mean(WDIR_ave),
                   MWD_ave = mean(MWD_ave),
                   CUTI_ave = mean(CUTI),
                   BEUTI_ave = mean(BEUTI))

#join species data with environmental
bruv_month_wide_fish <- ungroup(bruv_month_wide_fish)
bruv_month_wide_fish <- dplyr::select(bruv_month_wide_fish, -c(Analysis_Month))

bruv_full_fish <- full_join(bruv_environment_subset, bruv_month_wide_fish, by = c("Site","Date"))

#add abundance column 
bruv_full_fish <- ungroup(bruv_full_fish)
bruv_full_fish <- bruv_full_fish %>% 
  dplyr::mutate(Abundance = rowSums(dplyr::select(bruv_full_fish,`C. stigmaeus`:`S. chrysomelas`)),
                .before = `C. stigmaeus`)

#calculate relative abundance
min_ab = min(bruv_full_fish$Abundance)
max_ab = max(bruv_full_fish$Abundance)
bruv_full_fish <- bruv_full_fish %>% 
  dplyr::mutate(Rel_Abundance = (Abundance-min_ab)/(max_ab-min_ab),
                .after = Abundance)

#check for NAs in data
colSums(is.na(bruv_full_fish[1:23]))

#extract groups
bruv_full_group_fish <- bruv_full_fish[1:25]
bruv_full_spp_fish <- bruv_full_fish[-c(1:25)]

#remove all non-environmental variables
bruv_full_group_fish <- ungroup(bruv_full_group_fish)
bruv_full_env_fish <- dplyr::select(bruv_full_group_fish, c("TideHeight","WindAve","Salinity",
                                                            "Depth", "Visibility","Temp",
                                                            "WaveHeight","WavePeriod","Red_Percent_Cover",
                                                            "Green_Percent_Cover", "Brown_Percent_Cover",
                                                            "Seagrass_Percent_Cover","WDIR_ave","MWD_ave"))


#add very small number to species data
bruv_full_spp_001_fish <- (bruv_full_spp_fish + 0.0001)

###Distanced Based Redundancy Analysis##########################################
#run distance based redundancy analysis by square rooting dissimilarities
dbRDA_sqrt_fish = capscale(bruv_full_spp_001_fish ~ TideHeight + WindAve + Salinity + Depth + Visibility +
                             Temp + WaveHeight + WavePeriod + Red_Percent_Cover + Brown_Percent_Cover + Green_Percent_Cover + 
                             Seagrass_Percent_Cover + WDIR_ave + MWD_ave, 
                           bruv_full_env_fish, dist = "bray", sqrt.dist = TRUE)

anova(dbRDA_sqrt_fish)
summary(dbRDA_sqrt_fish)

#test environmental variables significance
anova(dbRDA_sqrt_fish, by = "terms", permu = 200)

#get scores
scores_dbRDA = scores(dbRDA_sqrt_fish)
scores_dbRDA

#calculate loading and environmental correlations with axes
site_scores_environment = cbind(scores_dbRDA$sites, bruv_full_env_fish)
correlations = cor(site_scores_environment)
correlations2 = as.data.frame(correlations[3:16,1:2])
correlations2$environmental = rownames(correlations2)

#filter significant environmentals 
correlations_sig <- correlations2 %>% 
  dplyr::filter(environmental %in% c("Depth","Visibility","Temp","WaveHeight"))

#replace wave height with breaker height
correlations_sig$environmental <- str_replace(correlations_sig$environmental, "WaveHeight", "Breaker Height")
correlations_sig$environmental <- str_replace(correlations_sig$environmental, "Temp", "Water Temperature")

#convert the scores to a dataframe
species_scores = as.data.frame(scores_dbRDA$species)
species_scores$species = rownames(species_scores)

#calculate the species vectors
vec.sp <- envfit(dbRDA_sqrt_fish, bruv_full_spp_001_fish, perm = 1000, p.max = 0.001)
vec.sp.df <- as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species <- rownames(vec.sp.df)

#add the other data back to the scores
site_scores = as.data.frame(scores_dbRDA$sites)
site_scores$site = bruv_full_group_fish$Site
site_scores$season = bruv_full_group_fish$Season
site_scores$rel_ab = bruv_full_group_fish$Rel_Abundance

site_scores$season <- factor(site_scores$season, levels = c("Winter","Spring","Summer","Fall"))

#graph with environmental variables
ggplot(site_scores) +
  geom_point(aes(x = CAP1, y = CAP2, color = season, size = rel_ab)) +
  geom_segment(data = correlations_sig, aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.25, "cm")), color = "grey50", inherit.aes = FALSE, size = 0.65) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  #geom_text(data = correlations_sig, aes(x = CAP1*1.25, y = CAP2*1.25, label = environmental), size = 6) +
  geom_text_repel(data = correlations_sig, aes(x = CAP1, y = CAP2, label = environmental), size = 7,
                  segment.color = "white") +
  scale_color_manual(name = "Season", values = c("dodgerblue","orchid","forestgreen","darkorange")) +
  scale_size_continuous(name = "Relative \nAbundance", range = c(4,10), guide = 'none') +
  # labs(title = "Distance Based Redundancy Analysis of Fish Species by Season",
  #      subtitle = "Square Root Transformed with Significant Environmentals") +
  theme(plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 22),
        plot.caption = element_text(size = 12)) +   
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 18),
        legend.position = c(0.9, 0.15)) +
  guides(colour = guide_legend(override.aes = list(size=4)))


#graph with species variables
ggplot(site_scores) +
  geom_point(aes(x = CAP1, y = CAP2, color = season, size = rel_ab)) +
  geom_segment(data = vec.sp.df, aes(x = 0, xend = CAP1*2.5, y = 0, yend = CAP2*2.5),
               arrow = arrow(length = unit(0.25, "cm")), color = "grey50", inherit.aes = FALSE, size = 0.65) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  geom_text_repel(data = vec.sp.df, aes(x = CAP1*2.5, y = CAP2*2.5, label = species), size = 6) +
  #geom_text(data = vec.sp.df, aes(x = CAP1, y = CAP2, label = species), size = 6) + 
  scale_color_manual(name = "Season", values = c("dodgerblue","orchid","forestgreen","darkorange")) +
  scale_size_continuous(name = "Relative \nAbundance", range = c(4,10), guide = 'none') +
  # labs(title = "Distance Based Redundancy Analysis of Fish Species by Season",
  #      subtitle = "Square Root Transformed with Species") +
  theme(plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 22),
        plot.caption = element_text(size = 12)) +   
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 18),
        legend.position = c(0.9, 0.15)) +
  guides(colour = guide_legend(override.aes = list(size=4)))



