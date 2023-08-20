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

###Shannon-Wiener Diversity#####################################################
#create a new variable for diversity
bruv_fish_subset_diver <- dplyr::select(bruv_fish_subset, -c(MaxN_ave))

#Convert to wide
bruv_fish_wide_std <- spread(bruv_fish_subset_diver, Sci_Name_Abbr, MaxN)

#Fill in NA with 0
bruv_fish_wide_std[is.na(bruv_fish_wide_std)] <- 0

#Calculate Shannon Wiener Diversity
H <- diversity(bruv_fish_wide_std[-c(1:3)], index="shannon")

#add diversity back to data
bruv_fish_wide_std <- bruv_fish_wide_std %>% add_column(H) 

#add back in the additional data and the missing dates
shannon_fish <- bruv_fish_wide_std[c("Site", "Date", "Analysis_Month", "H")]
shannon_fish <- ungroup(shannon_fish) 
shannon_fish <- shannon_fish %>% 
  add_row(Site = "Carmel", Date = as.Date("2021-05-01"), Analysis_Month = "May", H = 0) %>% 
  add_row(Site =  "Spanish Bay", Date = as.Date("2021-06-01"), Analysis_Month = "Jun", H = 0)

head(shannon_fish)

shannon_fish %>% 
  group_by(Site) %>% 
  summarize(count = n())

#add season into the data
shannon_fish <-  shannon_fish %>% 
  add_column(Numeric_Date = NA,
             Season = NA) 

for (i in 1:length(shannon_fish$Analysis_Month)) {
  if (shannon_fish$Analysis_Month[i] == "Jan") {
    shannon_fish$Numeric_Date[i] <- 1
    shannon_fish$Season[i] <- "Winter"
  } else if (shannon_fish$Analysis_Month[i] == "Feb") {
    shannon_fish$Numeric_Date[i] <- 2
    shannon_fish$Season[i] <- "Winter"
  } else if (shannon_fish$Analysis_Month[i] == "Mar") {
    shannon_fish$Numeric_Date[i] <- 3
    shannon_fish$Season[i] <- "Spring"
  } else if (shannon_fish$Analysis_Month[i] == "Apr") {
    shannon_fish$Numeric_Date[i] <- 4
    shannon_fish$Season[i] <- "Spring"
  } else if (shannon_fish$Analysis_Month[i] == "May") {
    shannon_fish$Numeric_Date[i] <- 5
    shannon_fish$Season[i] <- "Spring"
  } else if (shannon_fish$Analysis_Month[i] == "Jun") {
    shannon_fish$Numeric_Date[i] <- 6
    shannon_fish$Season[i] <- "Summer"
  } else if (shannon_fish$Analysis_Month[i] == "Jul") {
    shannon_fish$Numeric_Date[i] <- 7
    shannon_fish$Season[i] <- "Summer"
  } else if (shannon_fish$Analysis_Month[i] == "Aug") {
    shannon_fish$Numeric_Date[i] <- 8
    shannon_fish$Season[i] <- "Summer"
  } else if (shannon_fish$Analysis_Month[i] == "Sep") {
    shannon_fish$Numeric_Date[i] <- 9
    shannon_fish$Season[i] <- "Fall"
  } else if (shannon_fish$Analysis_Month[i] == "Oct") {
    shannon_fish$Numeric_Date[i] <- 10
    shannon_fish$Season[i] <- "Fall"
  } else if (shannon_fish$Analysis_Month[i] == "Nov") {
    shannon_fish$Numeric_Date[i] <- 11
    shannon_fish$Season[i] <- "Fall"
  } else { 
    shannon_fish$Numeric_Date[i] <- 12
    shannon_fish$Season[i] <- "Winter"
  }
}

#Relabel Factors
shannon_fish$Season <- factor(shannon_fish$Season, levels = c("Winter","Spring","Summer","Fall"))


shannon_fish$Analysis_Month <- factor(shannon_fish$Analysis_Month, levels = c("Dec","Jan","Feb","Mar","Apr","May",
                                                                              "Jun","Jul","Aug","Sep","Oct","Nov"))
#Diversity bar by Site and Season
shannon_fish_summary_site_season <- summarySE(shannon_fish, measurevar="H", groupvars=c("Site", "Season"))
shannon_fish_summary_site_season$Site <- factor(shannon_fish_summary_site_season$Site, 
                                                levels = c("Spanish Bay","Carmel","Whalers Cove","Stillwater Cove"))
shannon_fish_summary_site_season$Season <- factor(shannon_fish_summary_site_season$Season, 
                                                  levels = c("Winter","Spring","Summer","Fall"))


ggplot(shannon_fish_summary_site_season, aes(x = Season, y = H, fill = Season)) +
  facet_wrap(.~Site, nrow = 1) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=H-se, ymax=H+se), width=.2, position=position_dodge(.9)) +
  theme_bw() +
  ylab("Diversity") +
  xlab("Season") + 
  labs(title = "Shannon-Wiener Diversity Index of Fish by Site and Season") +
  scale_y_continuous(limits=c(0, 2), oob = rescale_none) +
  scale_fill_manual(values = c("dodgerblue","orchid","forestgreen","darkorange")) +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        plot.caption = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16, angle = 330, vjust = 0.5, hjust = 0.1),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 20)) 

#two-way ANOVA for Diversity
#two way anova
m3 <- aov(H ~ Site * Season, data = shannon_fish)
summary(m3) 

#get CLDs
#TukeyHSD(m3)
post1=HSD.test(m3, "Season", group = TRUE)
post1$groups

##Plot Interaction
with(shannon_fish, {
  interaction.plot(Site, Season, H)
})

###Abundance####################################################################
#variable for average abundance
bruv_fish_ab <- dplyr::select(bruv_fish_subset, -c(MaxN_ave))
bruv_fish_ab <- dplyr::rename(bruv_fish_ab, Abundance = MaxN)

#add in missing dates
bruv_fish_ab <- ungroup(bruv_fish_ab) 
bruv_fish_ab <- bruv_fish_ab %>% 
  add_row(Site = "Carmel", Date = as.Date("2021-05-01"), Analysis_Month = "May", Abundance = 0) %>% 
  add_row(Site = "Spanish Bay", Date = as.Date("2020-08-01"), Analysis_Month = "Aug", Abundance = 0) %>% 
  add_row(Site =  "Spanish Bay", Date = as.Date("2021-06-01"), Analysis_Month = "Jun", Abundance = 0)

#add numeric date and season to the data
bruv_fish_ab <-  bruv_fish_ab %>% 
  add_column(Numeric_Date = NA,
             Season = NA) 

for (i in 1:length(bruv_fish_ab$Analysis_Month)) {
  if (bruv_fish_ab$Analysis_Month[i] == "Jan") {
    bruv_fish_ab$Numeric_Date[i] <- 1
    bruv_fish_ab$Season[i] <- "Winter"
  } else if (bruv_fish_ab$Analysis_Month[i] == "Feb") {
    bruv_fish_ab$Numeric_Date[i] <- 2
    bruv_fish_ab$Season[i] <- "Winter"
  } else if (bruv_fish_ab$Analysis_Month[i] == "Mar") {
    bruv_fish_ab$Numeric_Date[i] <- 3
    bruv_fish_ab$Season[i] <- "Spring"
  } else if (bruv_fish_ab$Analysis_Month[i] == "Apr") {
    bruv_fish_ab$Numeric_Date[i] <- 4
    bruv_fish_ab$Season[i] <- "Spring"
  } else if (bruv_fish_ab$Analysis_Month[i] == "May") {
    bruv_fish_ab$Numeric_Date[i] <- 5
    bruv_fish_ab$Season[i] <- "Spring"
  } else if (bruv_fish_ab$Analysis_Month[i] == "Jun") {
    bruv_fish_ab$Numeric_Date[i] <- 6
    bruv_fish_ab$Season[i] <- "Summer"
  } else if (bruv_fish_ab$Analysis_Month[i] == "Jul") {
    bruv_fish_ab$Numeric_Date[i] <- 7
    bruv_fish_ab$Season[i] <- "Summer"
  } else if (bruv_fish_ab$Analysis_Month[i] == "Aug") {
    bruv_fish_ab$Numeric_Date[i] <- 8
    bruv_fish_ab$Season[i] <- "Summer"
  } else if (bruv_fish_ab$Analysis_Month[i] == "Sep") {
    bruv_fish_ab$Numeric_Date[i] <- 9
    bruv_fish_ab$Season[i] <- "Fall"
  } else if (bruv_fish_ab$Analysis_Month[i] == "Oct") {
    bruv_fish_ab$Numeric_Date[i] <- 10
    bruv_fish_ab$Season[i] <- "Fall"
  } else if (bruv_fish_ab$Analysis_Month[i] == "Nov") {
    bruv_fish_ab$Numeric_Date[i] <- 11
    bruv_fish_ab$Season[i] <- "Fall"
  } else { 
    bruv_fish_ab$Numeric_Date[i] <- 12
    bruv_fish_ab$Season[i] <- "Winter"
  }
}

#Relabel Factors
bruv_fish_ab$Season <- factor(bruv_fish_ab$Season, levels = c("Winter","Spring","Summer","Fall"))

bruv_fish_ab$Analysis_Month <- factor(bruv_fish_ab$Analysis_Month, levels = c("Dec","Jan","Feb","Mar","Apr","May",
                                                                              "Jun","Jul","Aug","Sep","Oct","Nov"))
#Abundance bar by Site and Season
abundance_fish_summary_site_season <- summarySE(bruv_fish_ab, measurevar="Abundance", groupvars=c("Site", "Season"))
abundance_fish_summary_site_season$Site <- factor(abundance_fish_summary_site_season$Site, 
                                                  levels = c("Spanish Bay","Carmel","Whalers Cove","Stillwater Cove"))
abundance_fish_summary_site_season$Season <- factor(abundance_fish_summary_site_season$Season, 
                                                    levels = c("Winter","Spring","Summer","Fall"))

#graph of site and season for fish abundance
ggplot(abundance_fish_summary_site_season, aes(x = Season, y = Abundance, fill = Season)) +
  facet_wrap(.~Site, nrow = 1) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Abundance-se, ymax=Abundance+se), width=.2, position=position_dodge(.9)) +
  theme_bw() +
  ylab("MaxN") +
  xlab("Season") + 
  labs(title = "Abundance of Fish by Site and Season") +
  scale_fill_manual(values = c("dodgerblue","orchid","forestgreen","darkorange")) +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        plot.caption = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16, angle = 330, vjust = 0.5, hjust = 0.1),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 20))

#calculate a monthly abundance value for the anova test
bruv_fish_anova <- bruv_fish_ab %>% 
  group_by(Site, Season, Analysis_Month) %>% 
  summarize(Abundance = mean(Abundance)) %>% 
  ungroup()

#two way anova
m3 <- aov(Abundance ~ Site * Season, data = bruv_fish_anova)
summary(m3) #if interaction is not significant I would remove

#get CLDs
TukeyHSD(m3)
post1=HSD.test(m3, "Site")
post1$groups

##Plot Interaction
with(shannon_fish, {
  interaction.plot(Site, Season, H)
})
