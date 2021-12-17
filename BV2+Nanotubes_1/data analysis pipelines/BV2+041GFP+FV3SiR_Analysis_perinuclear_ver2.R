setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()) 

# 2021 12 10

# data analysis of 


#---------------------- libraries and extras -----------------------------

library(readr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(doBy)
my_comparisons <- list( c("control", "25ugmL"),c("control", "50ugmL"),c("control", "100ugmL"))


#---------- Functions-----------


#this one dworks only on lists
remove_outlayer<-function(x){
  x1<-x[!x %in% boxplot.stats(x)$out]
  return(x1)
}






# Analysing the signal from the cytoplasm:

#---------------------------  041 GFP-------------------------------------


cytoplasm <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/041/BV2_041_GFP_perinuclear.csv")
cytoplasm$TYPE <- c('cytoplasm')


names(cytoplasm)[names(cytoplasm) == "Metadata_Concentration...6"] <- "Metadata_Concentration"
cytoplasm$Metadata_Concentration <- factor(cytoplasm$Metadata_Concentration , levels=c("control", "25ugmL", "50ugmL", "100ugmL"))


data_woLPS<- filter(cytoplasm, Metadata_LPS == '041GFP_wo')
data_wiLPS<- filter(cytoplasm, Metadata_LPS == '041GFP_wi')

# wi LPS

cntrl11<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
cntrl12<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
cntrl21<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
cntrl22<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
cntrl31<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
cntrl32<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_25_11 <- filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_25_12<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_25_21<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_25_22<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_25_31<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_25_32<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_50_11 <- filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_50_12<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_50_21<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_50_22<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_50_31<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_50_32<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_100_11 <- filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_100_12<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_100_21<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_100_22<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_100_31<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_100_32<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)


# here we analyse every cell in a well (all images are accumulated)
# clean outlayers from each well

cntrl11_clean<-mean(remove_outlayer(cntrl11$Intensity_IntegratedIntensity_GFP))
cntrl12_clean<-mean(remove_outlayer(cntrl12$Intensity_IntegratedIntensity_GFP))
cntrl21_clean<-mean(remove_outlayer(cntrl21$Intensity_IntegratedIntensity_GFP))
cntrl22_clean<-mean(remove_outlayer(cntrl22$Intensity_IntegratedIntensity_GFP))
cntrl31_clean<-mean(remove_outlayer(cntrl31$Intensity_IntegratedIntensity_GFP))
cntrl32_clean<-mean(remove_outlayer(cntrl32$Intensity_IntegratedIntensity_GFP))

gfp_25_11_clean<-mean(remove_outlayer(gfp_25_11$Intensity_IntegratedIntensity_GFP))
gfp_25_12_clean<-mean(remove_outlayer(gfp_25_12$Intensity_IntegratedIntensity_GFP))
gfp_25_21_clean<-mean(remove_outlayer(gfp_25_21$Intensity_IntegratedIntensity_GFP))
gfp_25_22_clean<-mean(remove_outlayer(gfp_25_22$Intensity_IntegratedIntensity_GFP))
gfp_25_31_clean<-mean(remove_outlayer(gfp_25_31$Intensity_IntegratedIntensity_GFP))
gfp_25_32_clean<-mean(remove_outlayer(gfp_25_32$Intensity_IntegratedIntensity_GFP))

gfp_50_11_clean<-mean(remove_outlayer(gfp_50_11$Intensity_IntegratedIntensity_GFP))
gfp_50_12_clean<-mean(remove_outlayer(gfp_50_12$Intensity_IntegratedIntensity_GFP))
gfp_50_21_clean<-mean(remove_outlayer(gfp_50_21$Intensity_IntegratedIntensity_GFP))
gfp_50_22_clean<-mean(remove_outlayer(gfp_50_22$Intensity_IntegratedIntensity_GFP))
gfp_50_31_clean<-mean(remove_outlayer(gfp_50_31$Intensity_IntegratedIntensity_GFP))
gfp_50_32_clean<-mean(remove_outlayer(gfp_50_32$Intensity_IntegratedIntensity_GFP))

gfp_100_11_clean<-mean(remove_outlayer(gfp_100_11$Intensity_IntegratedIntensity_GFP))
gfp_100_12_clean<-mean(remove_outlayer(gfp_100_12$Intensity_IntegratedIntensity_GFP))
gfp_100_21_clean<-mean(remove_outlayer(gfp_100_21$Intensity_IntegratedIntensity_GFP))
gfp_100_22_clean<-mean(remove_outlayer(gfp_100_22$Intensity_IntegratedIntensity_GFP))
gfp_100_31_clean<-mean(remove_outlayer(gfp_100_31$Intensity_IntegratedIntensity_GFP))
gfp_100_32_clean<-mean(remove_outlayer(gfp_100_32$Intensity_IntegratedIntensity_GFP))

# normalise by experiment

normalising_factor1 <- mean(c(cntrl11_clean,cntrl12_clean))
normalising_factor2 <- mean(c(cntrl21_clean,cntrl22_clean))
normalising_factor3 <- mean(c(cntrl31_clean,cntrl32_clean))

cntrl_wilps1 <- mean(c(cntrl11_clean,cntrl12_clean))/normalising_factor1
cntrl_wilps2 <- mean(c(cntrl21_clean,cntrl22_clean))/normalising_factor2
cntrl_wilps3 <- mean(c(cntrl31_clean,cntrl32_clean))/normalising_factor3

gfp_25_wilps1 <- mean(c(gfp_25_11_clean,gfp_25_12_clean))/normalising_factor1
gfp_25_wilps2 <- mean(c(gfp_25_21_clean,gfp_25_22_clean))/normalising_factor2
gfp_25_wilps3 <- mean(c(gfp_25_31_clean,gfp_25_32_clean))/normalising_factor3

gfp_50_wilps1 <- mean(c(gfp_50_11_clean,gfp_50_12_clean))/normalising_factor1
gfp_50_wilps2 <- mean(c(gfp_50_21_clean,gfp_50_22_clean))/normalising_factor2
gfp_50_wilps3 <- mean(c(gfp_50_31_clean,gfp_50_32_clean))/normalising_factor3

gfp_100_wilps1 <- mean(c(gfp_100_11_clean,gfp_100_12_clean))/normalising_factor1
gfp_100_wilps2 <- mean(c(gfp_100_21_clean,gfp_100_22_clean))/normalising_factor2
gfp_100_wilps3 <- mean(c(gfp_100_31_clean,gfp_100_32_clean))/normalising_factor3

LPS <- c('+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS')
Concentration <- c('0 ug/mL','0 ug/mL','0 ug/mL','25 ug/mL','25 ug/mL','25 ug/mL','50 ug/mL','50 ug/mL','50 ug/mL','100 ug/mL','100 ug/mL','100 ug/mL')
Integrated_Intensity <- c(cntrl_wilps1,cntrl_wilps2,cntrl_wilps3, 
                          gfp_25_wilps1,gfp_25_wilps2,gfp_25_wilps3,
                          gfp_50_wilps1,gfp_50_wilps2,gfp_50_wilps3,
                          gfp_100_wilps1,gfp_100_wilps2,gfp_100_wilps3)

WILPS_041<-data.frame(Concentration,Integrated_Intensity,LPS)

ggplot(WILPS_041, aes(x=Concentration,y=Integrated_Intensity))+
  geom_boxplot()+stat_compare_means(comparisons = my_comparisons, method='t.test') +   geom_jitter() 

Integrated_Intensity_nonorm<- c(mean(c(cntrl11_clean,cntrl12_clean)), mean(c(cntrl21_clean,cntrl22_clean)), mean(c(cntrl31_clean,cntrl32_clean)), 
                                mean(c(gfp_25_11_clean,gfp_25_12_clean)), mean(c(gfp_25_21_clean,gfp_25_22_clean)), mean(c(gfp_25_31_clean,gfp_25_32_clean)),
                                mean(c(gfp_50_11_clean,gfp_50_12_clean)), mean(c(gfp_50_21_clean,gfp_50_22_clean)),mean(c(gfp_50_31_clean,gfp_50_32_clean)),
                                mean(c(gfp_100_11_clean,gfp_100_12_clean)), mean(c(gfp_100_21_clean,gfp_100_22_clean)),mean(c(gfp_100_31_clean,gfp_100_22_clean)))
WILPS_041_nonorm<-data.frame(Concentration,Integrated_Intensity_nonorm,LPS)

# wo LPS

cntrl11<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
cntrl12<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
cntrl21<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
cntrl22<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
cntrl31<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
cntrl32<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_25_11 <- filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_25_12<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_25_21<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_25_22<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_25_31<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_25_32<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_50_11 <- filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_50_12<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_50_21<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_50_22<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_50_31<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_50_32<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_100_11 <- filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_100_12<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_100_21<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_100_22<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_100_31<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_100_32<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)


# here we analyse every cell in a well (all images are accumulated)
# clean outlayers from each well

cntrl11_clean<-mean(remove_outlayer(cntrl11$Intensity_IntegratedIntensity_GFP))
cntrl12_clean<-mean(remove_outlayer(cntrl12$Intensity_IntegratedIntensity_GFP))
cntrl21_clean<-mean(remove_outlayer(cntrl21$Intensity_IntegratedIntensity_GFP))
cntrl22_clean<-mean(remove_outlayer(cntrl22$Intensity_IntegratedIntensity_GFP))
cntrl31_clean<-mean(remove_outlayer(cntrl31$Intensity_IntegratedIntensity_GFP))
cntrl32_clean<-mean(remove_outlayer(cntrl32$Intensity_IntegratedIntensity_GFP))

gfp_25_11_clean<-mean(remove_outlayer(gfp_25_11$Intensity_IntegratedIntensity_GFP))
gfp_25_12_clean<-mean(remove_outlayer(gfp_25_12$Intensity_IntegratedIntensity_GFP))
gfp_25_21_clean<-mean(remove_outlayer(gfp_25_21$Intensity_IntegratedIntensity_GFP))
gfp_25_22_clean<-mean(remove_outlayer(gfp_25_22$Intensity_IntegratedIntensity_GFP))
gfp_25_31_clean<-mean(remove_outlayer(gfp_25_31$Intensity_IntegratedIntensity_GFP))
gfp_25_32_clean<-mean(remove_outlayer(gfp_25_32$Intensity_IntegratedIntensity_GFP))

gfp_50_11_clean<-mean(remove_outlayer(gfp_50_11$Intensity_IntegratedIntensity_GFP))
gfp_50_12_clean<-mean(remove_outlayer(gfp_50_12$Intensity_IntegratedIntensity_GFP))
gfp_50_21_clean<-mean(remove_outlayer(gfp_50_21$Intensity_IntegratedIntensity_GFP))
gfp_50_22_clean<-mean(remove_outlayer(gfp_50_22$Intensity_IntegratedIntensity_GFP))
gfp_50_31_clean<-mean(remove_outlayer(gfp_50_31$Intensity_IntegratedIntensity_GFP))
gfp_50_32_clean<-mean(remove_outlayer(gfp_50_32$Intensity_IntegratedIntensity_GFP))

gfp_100_11_clean<-mean(remove_outlayer(gfp_100_11$Intensity_IntegratedIntensity_GFP))
gfp_100_12_clean<-mean(remove_outlayer(gfp_100_12$Intensity_IntegratedIntensity_GFP))
gfp_100_21_clean<-mean(remove_outlayer(gfp_100_21$Intensity_IntegratedIntensity_GFP))
gfp_100_22_clean<-mean(remove_outlayer(gfp_100_22$Intensity_IntegratedIntensity_GFP))
gfp_100_31_clean<-mean(remove_outlayer(gfp_100_31$Intensity_IntegratedIntensity_GFP))
gfp_100_32_clean<-mean(remove_outlayer(gfp_100_32$Intensity_IntegratedIntensity_GFP))

# normalise by experiment

normalising_factor1 <- mean(c(cntrl11_clean,cntrl12_clean))
normalising_factor2 <- mean(c(cntrl21_clean,cntrl22_clean))
normalising_factor3 <- mean(c(cntrl31_clean,cntrl32_clean))

cntrl_wilps1 <- mean(c(cntrl11_clean,cntrl12_clean))/normalising_factor1
cntrl_wilps2 <- mean(c(cntrl21_clean,cntrl22_clean))/normalising_factor2
cntrl_wilps3 <- mean(c(cntrl31_clean,cntrl32_clean))/normalising_factor3

gfp_25_wilps1 <- mean(c(gfp_25_11_clean,gfp_25_12_clean))/normalising_factor1
gfp_25_wilps2 <- mean(c(gfp_25_21_clean,gfp_25_22_clean))/normalising_factor2
gfp_25_wilps3 <- mean(c(gfp_25_31_clean,gfp_25_32_clean))/normalising_factor3

gfp_50_wilps1 <- mean(c(gfp_50_11_clean,gfp_50_12_clean))/normalising_factor1
gfp_50_wilps2 <- mean(c(gfp_50_21_clean,gfp_50_22_clean))/normalising_factor2
gfp_50_wilps3 <- mean(c(gfp_50_31_clean,gfp_50_32_clean))/normalising_factor3

gfp_100_wilps1 <- mean(c(gfp_100_11_clean,gfp_100_12_clean))/normalising_factor1
gfp_100_wilps2 <- mean(c(gfp_100_21_clean,gfp_100_22_clean))/normalising_factor2
gfp_100_wilps3 <- mean(c(gfp_100_31_clean,gfp_100_32_clean))/normalising_factor3


LPS <- c('-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS')
Concentration <- c('0 ug/mL','0 ug/mL','0 ug/mL','25 ug/mL','25 ug/mL','25 ug/mL','50 ug/mL','50 ug/mL','50 ug/mL','100 ug/mL','100 ug/mL','100 ug/mL')
Integrated_Intensity <- c(cntrl_wilps1,cntrl_wilps2,cntrl_wilps3, 
                          gfp_25_wilps1,gfp_25_wilps2,gfp_25_wilps3,
                          gfp_50_wilps1,gfp_50_wilps2,gfp_50_wilps3,
                          gfp_100_wilps1,gfp_100_wilps2,gfp_100_wilps3)

WOLPS_041<-data.frame(Concentration,Integrated_Intensity,LPS)


Integrated_Intensity_nonorm<- c(mean(c(cntrl11_clean,cntrl12_clean)), mean(c(cntrl21_clean,cntrl22_clean)), mean(c(cntrl31_clean,cntrl32_clean)), 
                                mean(c(gfp_25_11_clean,gfp_25_12_clean)), mean(c(gfp_25_21_clean,gfp_25_22_clean)), mean(c(gfp_25_31_clean,gfp_25_32_clean)),
                                mean(c(gfp_50_11_clean,gfp_50_12_clean)), mean(c(gfp_50_21_clean,gfp_50_22_clean)),mean(c(gfp_50_31_clean,gfp_50_32_clean)),
                                mean(c(gfp_100_11_clean,gfp_100_12_clean)), mean(c(gfp_100_21_clean,gfp_100_22_clean)),mean(c(gfp_100_31_clean,gfp_100_22_clean)))
WOLPS_041_nonorm<-data.frame(Concentration,Integrated_Intensity_nonorm,LPS)

#------------------------------- FV3-AF488------
cytoplasm <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/fv3/MyExpt_perinuclear.csv")
cytoplasm$TYPE <- c('cytoplasm')


names(cytoplasm)[names(cytoplasm) == "Metadata_Concentration...6"] <- "Metadata_Concentration"
cytoplasm$Metadata_Concentration <- factor(cytoplasm$Metadata_Concentration , levels=c("control", "25ugmL", "50ugmL", "100ugmL"))


data_woLPS<- filter(cytoplasm, Metadata_LPS == 'FV3SiR_wo')
data_wiLPS<- filter(cytoplasm, Metadata_LPS == 'FV3SiR_wi')

# wi LPS

cntrl11<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
cntrl12<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
cntrl21<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
cntrl22<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
cntrl31<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
cntrl32<-filter(data_wiLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_25_11 <- filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_25_12<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_25_21<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_25_22<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_25_31<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_25_32<-filter(data_wiLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_50_11 <- filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_50_12<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_50_21<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_50_22<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_50_31<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_50_32<-filter(data_wiLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_100_11 <- filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_100_12<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_100_21<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_100_22<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_100_31<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_100_32<-filter(data_wiLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)



# here we analyse every cell in a well (all images are accumulated)
# clean outlayers from each well

cntrl11_clean<-mean(remove_outlayer(cntrl11$Intensity_IntegratedIntensity_GFP))
cntrl12_clean<-mean(remove_outlayer(cntrl12$Intensity_IntegratedIntensity_GFP))
cntrl21_clean<-mean(remove_outlayer(cntrl21$Intensity_IntegratedIntensity_GFP))
cntrl22_clean<-mean(remove_outlayer(cntrl22$Intensity_IntegratedIntensity_GFP))
cntrl31_clean<-mean(remove_outlayer(cntrl31$Intensity_IntegratedIntensity_GFP))
cntrl32_clean<-mean(remove_outlayer(cntrl32$Intensity_IntegratedIntensity_GFP))

gfp_25_11_clean<-mean(remove_outlayer(gfp_25_11$Intensity_IntegratedIntensity_GFP))
gfp_25_12_clean<-mean(remove_outlayer(gfp_25_12$Intensity_IntegratedIntensity_GFP))
gfp_25_21_clean<-mean(remove_outlayer(gfp_25_21$Intensity_IntegratedIntensity_GFP))
gfp_25_22_clean<-mean(remove_outlayer(gfp_25_22$Intensity_IntegratedIntensity_GFP))
gfp_25_31_clean<-mean(remove_outlayer(gfp_25_31$Intensity_IntegratedIntensity_GFP))
gfp_25_32_clean<-mean(remove_outlayer(gfp_25_32$Intensity_IntegratedIntensity_GFP))

gfp_50_11_clean<-mean(remove_outlayer(gfp_50_11$Intensity_IntegratedIntensity_GFP))
gfp_50_12_clean<-mean(remove_outlayer(gfp_50_12$Intensity_IntegratedIntensity_GFP))
gfp_50_21_clean<-mean(remove_outlayer(gfp_50_21$Intensity_IntegratedIntensity_GFP))
gfp_50_22_clean<-mean(remove_outlayer(gfp_50_22$Intensity_IntegratedIntensity_GFP))
gfp_50_31_clean<-mean(remove_outlayer(gfp_50_31$Intensity_IntegratedIntensity_GFP))
gfp_50_32_clean<-mean(remove_outlayer(gfp_50_32$Intensity_IntegratedIntensity_GFP))

gfp_100_11_clean<-mean(remove_outlayer(gfp_100_11$Intensity_IntegratedIntensity_GFP))
gfp_100_12_clean<-mean(remove_outlayer(gfp_100_12$Intensity_IntegratedIntensity_GFP))
gfp_100_21_clean<-mean(remove_outlayer(gfp_100_21$Intensity_IntegratedIntensity_GFP))
gfp_100_22_clean<-mean(remove_outlayer(gfp_100_22$Intensity_IntegratedIntensity_GFP))
gfp_100_31_clean<-mean(remove_outlayer(gfp_100_31$Intensity_IntegratedIntensity_GFP))
gfp_100_32_clean<-mean(remove_outlayer(gfp_100_32$Intensity_IntegratedIntensity_GFP))

# normalise by experiment

normalising_factor1 <- mean(c(cntrl11_clean,cntrl12_clean))
normalising_factor2 <- mean(c(cntrl21_clean,cntrl22_clean))
normalising_factor3 <- mean(c(cntrl31_clean,cntrl32_clean))

cntrl_wilps1 <- mean(c(cntrl11_clean,cntrl12_clean))/normalising_factor1
cntrl_wilps2 <- mean(c(cntrl21_clean,cntrl22_clean))/normalising_factor2
cntrl_wilps3 <- mean(c(cntrl31_clean,cntrl32_clean))/normalising_factor3

gfp_25_wilps1 <- mean(c(gfp_25_11_clean,gfp_25_12_clean))/normalising_factor1
gfp_25_wilps2 <- mean(c(gfp_25_21_clean,gfp_25_22_clean))/normalising_factor2
gfp_25_wilps3 <- mean(c(gfp_25_31_clean,gfp_25_32_clean))/normalising_factor3

gfp_50_wilps1 <- mean(c(gfp_50_11_clean,gfp_50_12_clean))/normalising_factor1
gfp_50_wilps2 <- mean(c(gfp_50_21_clean,gfp_50_22_clean))/normalising_factor2
gfp_50_wilps3 <- mean(c(gfp_50_31_clean,gfp_50_32_clean))/normalising_factor3

gfp_100_wilps1 <- mean(c(gfp_100_11_clean,gfp_100_12_clean))/normalising_factor1
gfp_100_wilps2 <- mean(c(gfp_100_21_clean,gfp_100_22_clean))/normalising_factor2
gfp_100_wilps3 <- mean(c(gfp_100_31_clean,gfp_100_32_clean))/normalising_factor3

LPS <- c('+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS','+LPS')
Concentration <- c('0 ug/mL','0 ug/mL','0 ug/mL','25 ug/mL','25 ug/mL','25 ug/mL','50 ug/mL','50 ug/mL','50 ug/mL','100 ug/mL','100 ug/mL','100 ug/mL')
Integrated_Intensity <- c(cntrl_wilps1,cntrl_wilps2,cntrl_wilps3, 
                          gfp_25_wilps1,gfp_25_wilps2,gfp_25_wilps3,
                          gfp_50_wilps1,gfp_50_wilps2,gfp_50_wilps3,
                          gfp_100_wilps1,gfp_100_wilps2,gfp_100_wilps3)

WILPS_FV3<-data.frame(Concentration,Integrated_Intensity,LPS)

Integrated_Intensity_nonorm<- c(mean(c(cntrl11_clean,cntrl12_clean)), mean(c(cntrl21_clean,cntrl22_clean)), mean(c(cntrl31_clean,cntrl32_clean)), 
                                mean(c(gfp_25_11_clean,gfp_25_12_clean)), mean(c(gfp_25_21_clean,gfp_25_22_clean)), mean(c(gfp_25_31_clean,gfp_25_32_clean)),
                                mean(c(gfp_50_11_clean,gfp_50_12_clean)), mean(c(gfp_50_21_clean,gfp_50_22_clean)),mean(c(gfp_50_31_clean,gfp_50_32_clean)),
                                mean(c(gfp_100_11_clean,gfp_100_12_clean)), mean(c(gfp_100_21_clean,gfp_100_22_clean)),mean(c(gfp_100_31_clean,gfp_100_22_clean)))
WILPS_FV3_nonorm<-data.frame(Concentration,Integrated_Intensity_nonorm,LPS)



# wo LPS

cntrl11<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
cntrl12<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
cntrl21<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
cntrl22<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
cntrl31<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
cntrl32<-filter(data_woLPS, Metadata_Concentration == 'control' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_25_11 <- filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_25_12<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_25_21<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_25_22<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_25_31<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_25_32<-filter(data_woLPS, Metadata_Concentration == '25ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_50_11 <- filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_50_12<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_50_21<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_50_22<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_50_31<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_50_32<-filter(data_woLPS, Metadata_Concentration == '50ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)

gfp_100_11 <- filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 1)
gfp_100_12<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 1 & Metadata_Well == 2)
gfp_100_21<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 1)
gfp_100_22<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 2 & Metadata_Well == 2)
gfp_100_31<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 1)
gfp_100_32<-filter(data_woLPS, Metadata_Concentration == '100ugmL' &   Metadata_Experiment_Number == 3 & Metadata_Well == 2)


# here we analyse every cell in a well (all images are accumulated)
# clean outlayers from each well

cntrl11_clean<-mean(remove_outlayer(cntrl11$Intensity_IntegratedIntensity_GFP))
cntrl12_clean<-mean(remove_outlayer(cntrl12$Intensity_IntegratedIntensity_GFP))
cntrl21_clean<-mean(remove_outlayer(cntrl21$Intensity_IntegratedIntensity_GFP))
cntrl22_clean<-mean(remove_outlayer(cntrl22$Intensity_IntegratedIntensity_GFP))
cntrl31_clean<-mean(remove_outlayer(cntrl31$Intensity_IntegratedIntensity_GFP))
cntrl32_clean<-mean(remove_outlayer(cntrl32$Intensity_IntegratedIntensity_GFP))

gfp_25_11_clean<-mean(remove_outlayer(gfp_25_11$Intensity_IntegratedIntensity_GFP))
gfp_25_12_clean<-mean(remove_outlayer(gfp_25_12$Intensity_IntegratedIntensity_GFP))
gfp_25_21_clean<-mean(remove_outlayer(gfp_25_21$Intensity_IntegratedIntensity_GFP))
gfp_25_22_clean<-mean(remove_outlayer(gfp_25_22$Intensity_IntegratedIntensity_GFP))
gfp_25_31_clean<-mean(remove_outlayer(gfp_25_31$Intensity_IntegratedIntensity_GFP))
gfp_25_32_clean<-mean(remove_outlayer(gfp_25_32$Intensity_IntegratedIntensity_GFP))

gfp_50_11_clean<-mean(remove_outlayer(gfp_50_11$Intensity_IntegratedIntensity_GFP))
gfp_50_12_clean<-mean(remove_outlayer(gfp_50_12$Intensity_IntegratedIntensity_GFP))
gfp_50_21_clean<-mean(remove_outlayer(gfp_50_21$Intensity_IntegratedIntensity_GFP))
gfp_50_22_clean<-mean(remove_outlayer(gfp_50_22$Intensity_IntegratedIntensity_GFP))
gfp_50_31_clean<-mean(remove_outlayer(gfp_50_31$Intensity_IntegratedIntensity_GFP))
gfp_50_32_clean<-mean(remove_outlayer(gfp_50_32$Intensity_IntegratedIntensity_GFP))

gfp_100_11_clean<-mean(remove_outlayer(gfp_100_11$Intensity_IntegratedIntensity_GFP))
gfp_100_12_clean<-mean(remove_outlayer(gfp_100_12$Intensity_IntegratedIntensity_GFP))
gfp_100_21_clean<-mean(remove_outlayer(gfp_100_21$Intensity_IntegratedIntensity_GFP))
gfp_100_22_clean<-mean(remove_outlayer(gfp_100_22$Intensity_IntegratedIntensity_GFP))
gfp_100_31_clean<-mean(remove_outlayer(gfp_100_31$Intensity_IntegratedIntensity_GFP))
gfp_100_32_clean<-mean(remove_outlayer(gfp_100_32$Intensity_IntegratedIntensity_GFP))

# normalise by experiment

normalising_factor1 <- mean(c(cntrl11_clean,cntrl12_clean))
normalising_factor2 <- mean(c(cntrl21_clean,cntrl22_clean))
normalising_factor3 <- mean(c(cntrl31_clean,cntrl32_clean))

cntrl_wilps1 <- mean(c(cntrl11_clean,cntrl12_clean))/normalising_factor1
cntrl_wilps2 <- mean(c(cntrl21_clean,cntrl22_clean))/normalising_factor2
cntrl_wilps3 <- mean(c(cntrl31_clean,cntrl32_clean))/normalising_factor3

gfp_25_wilps1 <- mean(c(gfp_25_11_clean,gfp_25_12_clean))/normalising_factor1
gfp_25_wilps2 <- mean(c(gfp_25_21_clean,gfp_25_22_clean))/normalising_factor2
gfp_25_wilps3 <- mean(c(gfp_25_31_clean,gfp_25_32_clean))/normalising_factor3

gfp_50_wilps1 <- mean(c(gfp_50_11_clean,gfp_50_12_clean))/normalising_factor1
gfp_50_wilps2 <- mean(c(gfp_50_21_clean,gfp_50_22_clean))/normalising_factor2
gfp_50_wilps3 <- mean(c(gfp_50_31_clean,gfp_50_32_clean))/normalising_factor3

gfp_100_wilps1 <- mean(c(gfp_100_11_clean,gfp_100_12_clean))/normalising_factor1
gfp_100_wilps2 <- mean(c(gfp_100_21_clean,gfp_100_22_clean))/normalising_factor2
gfp_100_wilps3 <- mean(c(gfp_100_31_clean,gfp_100_32_clean))/normalising_factor3


LPS <- c('-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS','-LPS')
Concentration <- c('0 ug/mL','0 ug/mL','0 ug/mL','25 ug/mL','25 ug/mL','25 ug/mL','50 ug/mL','50 ug/mL','50 ug/mL','100 ug/mL','100 ug/mL','100 ug/mL')
Integrated_Intensity <- c(cntrl_wilps1,cntrl_wilps2,cntrl_wilps3, 
                          gfp_25_wilps1,gfp_25_wilps2,gfp_25_wilps3,
                          gfp_50_wilps1,gfp_50_wilps2,gfp_50_wilps3,
                          gfp_100_wilps1,gfp_100_wilps2,gfp_100_wilps3)

WOLPS_FV3<-data.frame(Concentration,Integrated_Intensity,LPS)

Integrated_Intensity_nonorm<- c(mean(c(cntrl11_clean,cntrl12_clean)), mean(c(cntrl21_clean,cntrl22_clean)), mean(c(cntrl31_clean,cntrl32_clean)), 
                                mean(c(gfp_25_11_clean,gfp_25_12_clean)), mean(c(gfp_25_21_clean,gfp_25_22_clean)), mean(c(gfp_25_31_clean,gfp_25_32_clean)),
                                mean(c(gfp_50_11_clean,gfp_50_12_clean)), mean(c(gfp_50_21_clean,gfp_50_22_clean)),mean(c(gfp_50_31_clean,gfp_50_32_clean)),
                                mean(c(gfp_100_11_clean,gfp_100_12_clean)), mean(c(gfp_100_21_clean,gfp_100_22_clean)),mean(c(gfp_100_31_clean,gfp_100_22_clean)))
WOLPS_FV3_nonorm<-data.frame(Concentration,Integrated_Intensity_nonorm,LPS)

#----------- FINALISE - CYTOPLASM----------------------------------------------------------

FInal.Frame.041<- rbind(WOLPS_041,WILPS_041)
FInal.Frame.041$Concentration <- factor(FInal.Frame.041$Concentration , levels=c("0 ug/mL", "25 ug/mL", "50 ug/mL", "100 ug/mL"))

ggplot(FInal.Frame.041, aes(x=Concentration,y=Integrated_Intensity, fill=LPS))+
  geom_boxplot()+  ggtitle('041-GFP (PERINUCLEAR) n=3 normalised' )+
  xlab(c('041-GFP concentration'))+
  ylab(c('Fluorescence intensity fold increase'))#+   geom_jitter()

FInal.Frame.fv3<- rbind(WOLPS_FV3,WILPS_FV3)
FInal.Frame.fv3$Concentration <- factor(FInal.Frame.fv3$Concentration , levels=c("0 ug/mL", "25 ug/mL", "50 ug/mL", "100 ug/mL"))

ggplot(FInal.Frame.fv3, aes(x=Concentration,y=Integrated_Intensity, fill=LPS))+
  geom_boxplot() + ggtitle('FV3-SiR (PERINUCLEAR) n=3 normalised') +
  xlab(c('041-GFP concentration'))+
  ylab(c('Fluorescence intensity fold increase'))




# Finalised data w/o normalisation 

FInal.Frame.041_nonorm<- rbind(WOLPS_041_nonorm,WILPS_041_nonorm)
FInal.Frame.041_nonorm$Concentration <- factor(FInal.Frame.041_nonorm$Concentration , levels=c("0 ug/mL", "25 ug/mL", "50 ug/mL", "100 ug/mL"))

ggplot(FInal.Frame.041_nonorm, aes(x=Concentration,y=Integrated_Intensity_nonorm, fill=LPS))+
  geom_boxplot()+  ggtitle('041-GFP (PERINUCLEAR) n=3 w/o normalisation')+
  xlab(c('041-GFP concentration'))+
  ylab(c('Mean integrated fluorescence intensity'))#+   geom_jitter()

FInal.Frame.fv3_nonorm<- rbind(WOLPS_FV3_nonorm,WILPS_FV3_nonorm)
FInal.Frame.fv3_nonorm$Concentration <- factor(FInal.Frame.fv3_nonorm$Concentration , levels=c("0 ug/mL", "25 ug/mL", "50 ug/mL", "100 ug/mL"))

ggplot(FInal.Frame.fv3_nonorm, aes(x=Concentration,y=Integrated_Intensity_nonorm, fill=LPS))+
  geom_boxplot() + ggtitle('FV3-SiR (PERINUCLEAR) n=3 w/o normalisation') +
  xlab(c('041-GFP concentration'))+
  ylab(c('Mean integrated fluorescence intensity'))


#---------------- stuff
#m1<-summaryBy(Intensity_IntegratedIntensity_GFP ~ ImageNumber+Metadata_Concentration+TYPE+Metadata_Experiment_Number+Metadata_LPS+Metadata_Well, data = cytoplasm, FUN = mean)
#m1.5<-summaryBy(Intensity_IntegratedIntensity_GFP.mean ~ Metadata_Well+Metadata_Concentration+TYPE+Metadata_Experiment_Number+Metadata_LPS, data = m1, FUN = mean)
#m2<-summaryBy(Intensity_IntegratedIntensity_GFP.mean.mean ~ Metadata_Concentration+Metadata_Experiment_Number+TYPE+Metadata_LPS, data = m1.5, FUN = mean)


#m2_c <- filter(m2, TYPE == 'cytoplasm')
#ggplot(m2_c, aes(x=Metadata_Concentration,y=Intensity_IntegratedIntensity_GFP.mean.mean.mean, fill=Metadata_LPS))+geom_boxplot()+stat_compare_means(comparisons = my_comparisons, method='t.test') +   geom_jitter() + ggtitle('041 cytoplasm (statistics from 3 experiments)')








