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
library(scales)

#set working directory
path = "C:\\Users\\Owner\\Documents\\Surf Zone\\Koval_etal_2023_SurfZoneFish"
setwd(path)

###Load Data####################################################################
#import updated data
bruv_month <- read_csv("data\\maxn_period_seasonal_month_site.csv")
family <- read_csv("data\\FamilyList.csv", na = ".")
bruv_month$Date <- as.Date(bruv_month$Date, tryFormats = c("%m/%d/%Y"))

#add new field for unique bruvs
bruv_month_bruv <- bruv_month %>% 
  mutate(Unique_bruv = paste(bruv_month$Site, "_", bruv_month$Date, "_", bruv_month$BRUV, sep = ""))

#group by type and species and summarize the number of bruvs each was seen on
unique_species <- bruv_month_bruv %>% 
  group_by(Type, Sci_Name, Common_Name) %>% 
  summarize(num_bruvs = length(unique(Unique_bruv)))

#summarize the values to get total number of species per type
unique_species %>% 
  group_by(Type) %>% 
  summarize(count = n())

#create list of BRUVS without species
bruv_missing <- dplyr::filter(bruv_month, Type == "None") 
bruv_missing <- unique(bruv_month[c("Site","Date")])

#Filter out unidentified species and fish 
bruv_fish <- dplyr::filter(bruv_month, (Species != "spp") & (Type == 'Fish'))

#count species and filter fish species observed less than 10 times
fish_10 <- plyr::count(bruv_fish, vars = "Sci_Name_Full")
fish_10 <- dplyr::filter(fish_10, freq >= 10)
fish_species_list <- c(fish_10$Sci_Name_Full)

bruv_fish_subset <- subset(bruv_fish, Sci_Name_Full %in% fish_species_list)

#calculate average max N by species per BRUV
bruv_fish_subset <- bruv_fish_subset %>% 
  dplyr::group_by(Site, Date, BRUV, Analysis_Month, Sci_Name_Abbr) %>%
  dplyr::summarise(MaxN_ave = (sum(MaxN))/max(Periods),
                   MaxN = max(MaxN)) %>% 
  print()

#calculate average max N by species per sampling
bruv_fish_subset <- bruv_fish_subset %>% 
  dplyr::group_by(Site, Date, Analysis_Month, Sci_Name_Abbr) %>% 
  dplyr::summarise(MaxN_ave = (sum(MaxN_ave))/6,
                   MaxN = sum(MaxN)/6) %>% 
  print()


bruv_fish_subset_cn <- subset(bruv_fish, Sci_Name_Full %in% fish_species_list)

#calculate average max N by species per BRUV
bruv_fish_subset_cn <- bruv_fish_subset_cn %>% 
  dplyr::group_by(Site, Date, BRUV, Analysis_Month, Common_Name) %>%
  dplyr::summarise(MaxN_ave = (sum(MaxN))/max(Periods),
                   MaxN = max(MaxN)) %>% 
  print()

#calculate average max N by species per sampling
bruv_fish_subset_cn <- bruv_fish_subset_cn %>% 
  dplyr::group_by(Site, Date, Analysis_Month, Common_Name) %>% 
  dplyr::summarise(MaxN_ave = (sum(MaxN_ave))/6,
                   MaxN = sum(MaxN)/6) %>% 
  print()

###Non-Metric Multidimensional Scaling#########################################
#add numeric date and season to the data
bruv_fish_nmds <-  bruv_fish_subset %>% 
  add_column(Numeric_Date = NA,
             Season = NA) 

bruv_fish_nmds <- dplyr::select(bruv_fish_nmds, -c(MaxN_ave))

for (i in 1:length(bruv_fish_nmds$Analysis_Month)) {
  if (bruv_fish_nmds$Analysis_Month[i] == "Jan") {
    bruv_fish_nmds$Numeric_Date[i] <- 1
    bruv_fish_nmds$Season[i] <- "Winter"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "Feb") {
    bruv_fish_nmds$Numeric_Date[i] <- 2
    bruv_fish_nmds$Season[i] <- "Winter"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "Mar") {
    bruv_fish_nmds$Numeric_Date[i] <- 3
    bruv_fish_nmds$Season[i] <- "Spring"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "Apr") {
    bruv_fish_nmds$Numeric_Date[i] <- 4
    bruv_fish_nmds$Season[i] <- "Spring"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "May") {
    bruv_fish_nmds$Numeric_Date[i] <- 5
    bruv_fish_nmds$Season[i] <- "Spring"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "Jun") {
    bruv_fish_nmds$Numeric_Date[i] <- 6
    bruv_fish_nmds$Season[i] <- "Summer"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "Jul") {
    bruv_fish_nmds$Numeric_Date[i] <- 7
    bruv_fish_nmds$Season[i] <- "Summer"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "Aug") {
    bruv_fish_nmds$Numeric_Date[i] <- 8
    bruv_fish_nmds$Season[i] <- "Summer"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "Sep") {
    bruv_fish_nmds$Numeric_Date[i] <- 9
    bruv_fish_nmds$Season[i] <- "Fall"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "Oct") {
    bruv_fish_nmds$Numeric_Date[i] <- 10
    bruv_fish_nmds$Season[i] <- "Fall"
  } else if (bruv_fish_nmds$Analysis_Month[i] == "Nov") {
    bruv_fish_nmds$Numeric_Date[i] <- 11
    bruv_fish_nmds$Season[i] <- "Fall"
  } else { 
    bruv_fish_nmds$Numeric_Date[i] <- 12
    bruv_fish_nmds$Season[i] <- "Winter"
  }
}

#convert to wide for spp table
maxn_std_seasonal_fish_sqrt <- spread(bruv_fish_nmds, Sci_Name_Abbr, MaxN)

#Fill in NA with 0
maxn_std_seasonal_fish_sqrt[is.na(maxn_std_seasonal_fish_sqrt)] <- 0

#extract groups
maxn_std_seasonal_fish_sqrt_group <- maxn_std_seasonal_fish_sqrt[1:5]
maxn_std_seasonal_fish_sqrt_spp <- maxn_std_seasonal_fish_sqrt[-c(1:5)]

#calculate square root
maxn_std_seasonal_fish_sqrt_spp <- sqrt(maxn_std_seasonal_fish_sqrt_spp)

#calculate relative abundance
maxn_std_seasonal_fish_sqrt_rel <- decostand(maxn_std_seasonal_fish_sqrt_spp, method = 'total')

#calculate distance matrix
maxn_std_seasonal_fish_sqrt_distmat <- vegdist(maxn_std_seasonal_fish_sqrt_rel, method = "bray")
maxn_std_seasonal_fish_sqrt_distmat <- as.matrix(maxn_std_seasonal_fish_sqrt_distmat, labels = TRUE)
#write.csv(maxn_ave_seasonal_distmat, "..\\Results\\maxn_ave_distmat.csv")

#run nMDS
maxn_std_seasonal_fish_sqrt_nmds <- metaMDS(maxn_std_seasonal_fish_sqrt_distmat,
                                            distance = "bray",
                                            k = 3,
                                            maxit = 999,
                                            trymax = 250,
                                            wascores = TRUE)

maxn_std_seasonal_fish_sqrt_nmds       

#Shepards test
goodness(maxn_std_seasonal_fish_sqrt_nmds)

#stress plot
stressplot(maxn_std_seasonal_fish_sqrt_nmds)

#extract stores from nMDS
data.scores <- as.data.frame(scores(maxn_std_seasonal_fish_sqrt_nmds, display = c('sites')))
data.scores$site <- rownames(data.scores)
data.scores$grp <- maxn_std_seasonal_fish_sqrt_group$Season
data.scores$grp3 <- maxn_std_seasonal_fish_sqrt_group$Site
data.scores$grp1 <- maxn_std_seasonal_fish_sqrt_group$Analysis_Month
data.scores$grp2 <- maxn_std_seasonal_fish_sqrt_group$Numeric_Date

#convert data to factors
data.scores$grp <- factor(data.scores$grp, levels = c("Winter","Spring","Summer","Fall"))

maxn_std_seasonal_fish_sqrt$Season <- factor(maxn_std_seasonal_fish_sqrt$Season, 
                                             levels = c("Winter","Spring","Summer","Fall"))
data.scores$grp3 <- factor(data.scores$grp3, 
                           levels = c("Spanish Bay", "Carmel", "Whalers Cove", "Stillwater Cove"))

maxn_std_seasonal_fish_sqrt$Site <- factor(maxn_std_seasonal_fish_sqrt$Site, 
                                           levels = c("Spanish Bay", "Carmel", "Whalers Cove", "Stillwater Cove"))

#add species vectors 
vec.sp <- envfit(maxn_std_seasonal_fish_sqrt_nmds$points, maxn_std_seasonal_fish_sqrt_rel, perm = 1000)
vec.sp.df <- as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species <- rownames(vec.sp.df)

#nMDS by season with species vectors
ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = grp), level = 0.5, size = 1) +
  geom_point(aes(x = NMDS1, y = NMDS2,color = grp), size = 5) + 
  annotate(geom = "text", x = 0.4, y = -0.68, label = "Stress = 0.066", size = 7) +
  geom_segment(data = vec.sp.df, aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.5, "cm")), color = "grey50", inherit.aes = FALSE, size = 0.65) +
  scale_color_manual(name = "Season", values = c("dodgerblue","orchid","forestgreen","darkorange")) +
  #geom_text(data = vec.sp.df, aes(x = MDS1, y = MDS2, label = Common_Name), size = 6) + 
  geom_text_repel(data = vec.sp.df, aes(x = MDS1, y = MDS2, label = species), size = 6,
                  segment.color = "white") +
  theme_bw() + 
  ylim(-0.8, 0.8) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  # labs(title = "nMDS of Fish Species Abundances by Season \nwith Species Vectors",
  #      subtitle = "Square Root Transformed") +
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
        legend.position = c(0.89, 0.85)) 


#nMDS by site with species vectors
ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = maxn_std_seasonal_fish_sqrt_group$Site), 
               level = 0.5, size = 1) +
  geom_point(aes(x = NMDS1, y = NMDS2,color = maxn_std_seasonal_fish_sqrt_group$Site), size = 5) + 
  annotate(geom = "text", x = 0.4, y = -0.68, label = "Stress = 0.066", size = 7) +
  geom_segment(data = vec.sp.df, aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.5, "cm")), color = "grey50", inherit.aes = FALSE, size = 0.65) +
  scale_color_discrete(name = "Site") +
  #geom_text(data = vec.sp.df, aes(x = MDS1, y = MDS2, label = Common_Name), size = 6) + 
  geom_text_repel(data = vec.sp.df, aes(x = MDS1, y = MDS2, label = species), size = 6, 
                  segment.color = "white") +
  theme_bw() + 
  ylim(-0.8, 0.8) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  # labs(title = "nMDS of Fish Species Abundances by Site \nwith Species Vectors",
  #      subtitle = "Square Root Transformed") +
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
        legend.position = c(0.89, 0.85))

#run PERMANOVA by season * site
adon.results <- adonis2(maxn_std_seasonal_fish_sqrt_rel ~ maxn_std_seasonal_fish_sqrt_group$Season *
                          maxn_std_seasonal_fish_sqrt_group$Site,method = "bray", perm = 999)
print(adon.results) #p = 0.085 

#calculate multivariate dispersion
maxn_ave_fish_distmat_sqrt <- vegdist(maxn_std_seasonal_fish_sqrt_rel, method = "bray")
mod <- betadisper(maxn_ave_fish_distmat_sqrt, maxn_std_seasonal_fish_sqrt_group$Season)
mod
