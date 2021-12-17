# originaly started 2021 09 30
#last update 2021 12 02
# analysis of 041GFP and FV3-AF488 particle uptake by BV2 cells

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()) 




#---------------------- libraries and extras -----------------------------

library(readr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(doBy)
my_comparisons <- list( c("control", "25ugmL"),c("control", "50ugmL"),c("control", "100ugmL"))


#-----------------------------------041-GFP-------------------------------------

rm(list = ls()) 
my_comparisons <- list( c("control", "25ugmL"),c("control", "50ugmL"),c("control", "100ugmL"))
#load data
perinuclear <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/041/BV2_041_GFP_perinuclear.csv")
perinuclear$TYPE <- c('perinuclear')
cytoplasm <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/041/BV2_041_GFP_Cytoplasm.csv")
cytoplasm$TYPE <- c('cytoplasm')
FM <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/041/BV2_041_GFP_FM.csv")
FM$TYPE <- c('FM')
names(perinuclear)[names(perinuclear) == "Metadata_Concentration...6"] <- "Metadata_Concentration"
names(cytoplasm)[names(cytoplasm) == "Metadata_Concentration...6"] <- "Metadata_Concentration"
names(FM)[names(FM) == "Metadata_Concentration...6"] <- "Metadata_Concentration"
perinuclear$Metadata_Concentration <- factor(perinuclear$Metadata_Concentration , levels=c("control", "25ugmL", "50ugmL", "100ugmL"))
cytoplasm$Metadata_Concentration <- factor(cytoplasm$Metadata_Concentration , levels=c("control", "25ugmL", "50ugmL", "100ugmL"))
FM$Metadata_Concentration <- factor(FM$Metadata_Concentration , levels=c("control", "25ugmL", "50ugmL", "100ugmL"))
# Mergedata into a single frame
DF041<-rbind(FM,cytoplasm,perinuclear)
DF041$Metadata_Experiment_Number <- factor(DF041$Metadata_Experiment_Number)
# Note: the FM signal is in the whole cell. including the nuceleus.



# calculate the mean integrated intensity per every image
m1<-summaryBy(Intensity_IntegratedIntensity_GFP ~ ImageNumber+Metadata_Concentration+TYPE+Metadata_Experiment_Number+Metadata_LPS+Metadata_Well, data = DF041, FUN = mean)
m1.5<-summaryBy(Intensity_IntegratedIntensity_GFP.mean ~ Metadata_Well+Metadata_Concentration+TYPE+Metadata_Experiment_Number+Metadata_LPS, data = m1, FUN = mean)
m2<-summaryBy(Intensity_IntegratedIntensity_GFP.mean.mean ~ Metadata_Concentration+Metadata_Experiment_Number+TYPE+Metadata_LPS, data = m1.5, FUN = mean)

m2_p <- filter(m2, TYPE == 'perinuclear')
ggplot(m2_p, aes(x=Metadata_Concentration,y=Intensity_IntegratedIntensity_GFP.mean.mean.mean, fill=Metadata_LPS))+
  geom_boxplot()+stat_compare_means(comparisons = my_comparisons, method='t.test') +   geom_jitter() + ggtitle('041 perinuclear (statistics from 3 experiments)')

m2_c <- filter(m2, TYPE == 'cytoplasm')
ggplot(m2_c, aes(x=Metadata_Concentration,y=Intensity_IntegratedIntensity_GFP.mean.mean.mean, fill=Metadata_LPS))+
  geom_boxplot()+stat_compare_means(comparisons = my_comparisons, method='t.test') +   geom_jitter() + ggtitle('041 cytoplasm (statistics from 3 experiments)')

m2_FM <- filter(m2, TYPE == 'FM')
ggplot(m2_FM, aes(x=Metadata_Concentration,y=Intensity_IntegratedIntensity_GFP.mean.mean.mean, fill=Metadata_LPS))+
  geom_boxplot()+stat_compare_means(comparisons = my_comparisons, method='t.test') +   geom_jitter() + ggtitle('041 FM (statistics from 3 experiments)' )




#---------------------FV3-SiR--------------------------------------------------
rm(list = ls()) 
my_comparisons <- list( c("control", "25ugmL"),c("control", "50ugmL"),c("control", "100ugmL"))


# FV3-SiR 
# Load data
perinuclear <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/fv3/MyExpt_perinuclear.csv")
perinuclear$TYPE <- c('perinuclear')
cytoplasm <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/fv3/MyExpt_Cytoplasm.csv")
cytoplasm$TYPE <- c('cytoplasm')
FM <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/fv3/MyExpt_FM.csv")
FM$TYPE <- c('FM')
names(perinuclear)[names(perinuclear) == "Metadata_Concentration...6"] <- "Metadata_Concentration"
names(cytoplasm)[names(cytoplasm) == "Metadata_Concentration...6"] <- "Metadata_Concentration"
names(FM)[names(FM) == "Metadata_Concentration...6"] <- "Metadata_Concentration"
perinuclear$Metadata_Concentration <- factor(perinuclear$Metadata_Concentration , levels=c("control", "25ugmL", "50ugmL", "100ugmL"))
cytoplasm$Metadata_Concentration <- factor(cytoplasm$Metadata_Concentration , levels=c("control", "25ugmL", "50ugmL", "100ugmL"))
FM$Metadata_Concentration <- factor(FM$Metadata_Concentration , levels=c("control", "25ugmL", "50ugmL", "100ugmL"))
# Mergedata into a single frame
DF041<-rbind(FM,cytoplasm,perinuclear)
DF041$Metadata_Experiment_Number <- factor(DF041$Metadata_Experiment_Number)
# Note: the FM signal is in the whole cell. including the nuceleus.


# calculate the mean integrated intensity per every image
m1<-summaryBy(Intensity_IntegratedIntensity_GFP ~ ImageNumber+Metadata_Concentration+TYPE+Metadata_Experiment_Number+Metadata_LPS+Metadata_Well, data = DF041, FUN = mean)
m1.5<-summaryBy(Intensity_IntegratedIntensity_GFP.mean ~ Metadata_Well+Metadata_Concentration+TYPE+Metadata_Experiment_Number+Metadata_LPS, data = m1, FUN = mean)
m2<-summaryBy(Intensity_IntegratedIntensity_GFP.mean.mean ~ Metadata_Concentration+Metadata_Experiment_Number+TYPE+Metadata_LPS, data = m1.5, FUN = mean)

m2_p <- filter(m2, TYPE == 'perinuclear')
ggplot(m2_p, aes(x=Metadata_Concentration,y=Intensity_IntegratedIntensity_GFP.mean.mean.mean, fill=Metadata_LPS))+
  geom_boxplot()+stat_compare_means(comparisons = my_comparisons, method='t.test') +   geom_jitter() + ggtitle('FV3 perinuclear (statistics from 3 experiments)')

m2_c <- filter(m2, TYPE == 'cytoplasm')
ggplot(m2_c, aes(x=Metadata_Concentration,y=Intensity_IntegratedIntensity_GFP.mean.mean.mean, fill=Metadata_LPS))+
  geom_boxplot()+stat_compare_means(comparisons = my_comparisons, method='t.test') +   geom_jitter() + ggtitle('FV3 cytoplasm (statistics from 3 experiments)')

m2_FM <- filter(m2, TYPE == 'FM')
ggplot(m2_FM, aes(x=Metadata_Concentration,y=Intensity_IntegratedIntensity_GFP.mean.mean.mean, fill=Metadata_LPS))+
  geom_boxplot()+stat_compare_means(comparisons = my_comparisons, method='t.test') +   geom_jitter() + ggtitle('FV3 FM (statistics from 3 experiments)')





