setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()) 

# 2021 12 16


# minor bugs fixed.
# the 041 data intensity is rescaled between 0.15-1.0
# the fv3 data intensity is rescaled betwqeen 0.1-1.0 (?)- check 

# data analysis of  041 and fv3 internalisation
#cytoplasm and membrane


#---------------------- libraries and extras -----------------------------

library(readr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(doBy)


my_comparisons <- list( c("25 ug/mL", "50 ug/mL"),c("25 ug/mL", "100 ug/mL"),c("50 ug/mL", "100 ug/mL"))
my_comparisons1 <- list( c("Without treatment", "100ng/mL LPS"))

remove_outlayer<-function(x){
  x1<-x[!x %in% boxplot.stats(x)$out]
  return(x1)
}



# Analysing the signal from the cytoplasm:


# cytoplasm 
#---------------------------  041 GFP-------------------------------------

cytoplasm <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/final/BV2_041_GFP_FINALCytoplasm.csv")

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

cntrl11_clean<-mean(remove_outlayer(cntrl11$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl12_clean<-mean(remove_outlayer(cntrl12$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl21_clean<-mean(remove_outlayer(cntrl21$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl22_clean<-mean(remove_outlayer(cntrl22$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl31_clean<-mean(remove_outlayer(cntrl31$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl32_clean<-mean(remove_outlayer(cntrl32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_25_11_clean<-mean(remove_outlayer(gfp_25_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_12_clean<-mean(remove_outlayer(gfp_25_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_21_clean<-mean(remove_outlayer(gfp_25_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_22_clean<-mean(remove_outlayer(gfp_25_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_31_clean<-mean(remove_outlayer(gfp_25_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_32_clean<-mean(remove_outlayer(gfp_25_32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_50_11_clean<-mean(remove_outlayer(gfp_50_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_12_clean<-mean(remove_outlayer(gfp_50_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_21_clean<-mean(remove_outlayer(gfp_50_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_22_clean<-mean(remove_outlayer(gfp_50_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_31_clean<-mean(remove_outlayer(gfp_50_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_32_clean<-mean(remove_outlayer(gfp_50_32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_100_11_clean<-mean(remove_outlayer(gfp_100_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_12_clean<-mean(remove_outlayer(gfp_100_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_21_clean<-mean(remove_outlayer(gfp_100_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_22_clean<-mean(remove_outlayer(gfp_100_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_31_clean<-mean(remove_outlayer(gfp_100_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_32_clean<-mean(remove_outlayer(gfp_100_32$Intensity_IntegratedIntensity_GFP_aftermath))

# normalise by experiment

normalising_factor1 <- mean(c(cntrl11_clean,cntrl12_clean))+0.0001
normalising_factor2 <- mean(c(cntrl21_clean,cntrl22_clean))+0.0001
normalising_factor3 <- mean(c(cntrl31_clean,cntrl32_clean))+0.0001

cntrl_wilps1 <- normalising_factor1/normalising_factor1
cntrl_wilps2 <- normalising_factor2/normalising_factor2
cntrl_wilps3 <- normalising_factor3/normalising_factor3

gfp_25_wilps1 <- mean(c(gfp_25_11_clean,gfp_25_12_clean))/normalising_factor1
gfp_25_wilps2 <- mean(c(gfp_25_21_clean,gfp_25_22_clean))/normalising_factor2
gfp_25_wilps3 <- mean(c(gfp_25_31_clean,gfp_25_32_clean))/normalising_factor3

gfp_50_wilps1 <- mean(c(gfp_50_11_clean,gfp_50_12_clean))/normalising_factor1
gfp_50_wilps2 <- mean(c(gfp_50_21_clean,gfp_50_22_clean))/normalising_factor2
gfp_50_wilps3 <- mean(c(gfp_50_31_clean,gfp_50_32_clean))/normalising_factor3

gfp_100_wilps1 <- mean(c(gfp_100_11_clean,gfp_100_12_clean))/normalising_factor1
gfp_100_wilps2 <- mean(c(gfp_100_21_clean,gfp_100_22_clean))/normalising_factor2
gfp_100_wilps3 <- mean(c(gfp_100_31_clean,gfp_100_32_clean))/normalising_factor3

LPS <- c('100ng/mL LPS','100ng/mL LPS','100ng/mL LPS',
         '100ng/mL LPS','100ng/mL LPS','100ng/mL LPS',
         '100ng/mL LPS','100ng/mL LPS','100ng/mL LPS',
         '100ng/mL LPS','100ng/mL LPS','100ng/mL LPS')
Concentration <- c('0 ug/mL','0 ug/mL','0 ug/mL','25 ug/mL','25 ug/mL','25 ug/mL','50 ug/mL','50 ug/mL','50 ug/mL','100 ug/mL','100 ug/mL','100 ug/mL')
Integrated_Intensity <- c(cntrl_wilps1,cntrl_wilps2,cntrl_wilps3, 
                          gfp_25_wilps1,gfp_25_wilps2,gfp_25_wilps3,
                          gfp_50_wilps1,gfp_50_wilps2,gfp_50_wilps3,
                          gfp_100_wilps1,gfp_100_wilps2,gfp_100_wilps3)

WILPS_041<-data.frame(Concentration,Integrated_Intensity,LPS)


Integrated_Intensity_nonorm<- c(mean(c(cntrl11_clean,cntrl12_clean)), mean(c(cntrl21_clean,cntrl22_clean)), mean(c(cntrl31_clean,cntrl32_clean)), 
                                mean(c(gfp_25_11_clean,gfp_25_12_clean)), mean(c(gfp_25_21_clean,gfp_25_22_clean)), mean(c(gfp_25_31_clean,gfp_25_32_clean)),
                                mean(c(gfp_50_11_clean,gfp_50_12_clean)), mean(c(gfp_50_21_clean,gfp_50_22_clean)),mean(c(gfp_50_31_clean,gfp_50_32_clean)),
                                mean(c(gfp_100_11_clean,gfp_100_12_clean)), mean(c(gfp_100_21_clean,gfp_100_22_clean)),mean(c(gfp_100_31_clean,gfp_100_32_clean)))
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

cntrl11_clean<-mean(remove_outlayer(cntrl11$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl12_clean<-mean(remove_outlayer(cntrl12$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl21_clean<-mean(remove_outlayer(cntrl21$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl22_clean<-mean(remove_outlayer(cntrl22$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl31_clean<-mean(remove_outlayer(cntrl31$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl32_clean<-mean(remove_outlayer(cntrl32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_25_11_clean<-mean(remove_outlayer(gfp_25_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_12_clean<-mean(remove_outlayer(gfp_25_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_21_clean<-mean(remove_outlayer(gfp_25_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_22_clean<-mean(remove_outlayer(gfp_25_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_31_clean<-mean(remove_outlayer(gfp_25_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_32_clean<-mean(remove_outlayer(gfp_25_32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_50_11_clean<-mean(remove_outlayer(gfp_50_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_12_clean<-mean(remove_outlayer(gfp_50_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_21_clean<-mean(remove_outlayer(gfp_50_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_22_clean<-mean(remove_outlayer(gfp_50_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_31_clean<-mean(remove_outlayer(gfp_50_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_32_clean<-mean(remove_outlayer(gfp_50_32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_100_11_clean<-mean(remove_outlayer(gfp_100_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_12_clean<-mean(remove_outlayer(gfp_100_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_21_clean<-mean(remove_outlayer(gfp_100_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_22_clean<-mean(remove_outlayer(gfp_100_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_31_clean<-mean(remove_outlayer(gfp_100_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_32_clean<-mean(remove_outlayer(gfp_100_32$Intensity_IntegratedIntensity_GFP_aftermath))

# normalise by experiment

normalising_factor1 <- mean(c(cntrl11_clean,cntrl12_clean))+0.0001
normalising_factor2 <- mean(c(cntrl21_clean,cntrl22_clean))+0.0001
normalising_factor3 <- mean(c(cntrl31_clean,cntrl32_clean))+0.0001

cntrl_wilps1 <- normalising_factor1/normalising_factor1
cntrl_wilps2 <- normalising_factor2/normalising_factor2
cntrl_wilps3 <- normalising_factor3/normalising_factor3

gfp_25_wilps1 <- mean(c(gfp_25_11_clean,gfp_25_12_clean))/normalising_factor1
gfp_25_wilps2 <- mean(c(gfp_25_21_clean,gfp_25_22_clean))/normalising_factor2
gfp_25_wilps3 <- mean(c(gfp_25_31_clean,gfp_25_32_clean))/normalising_factor3

gfp_50_wilps1 <- mean(c(gfp_50_11_clean,gfp_50_12_clean))/normalising_factor1
gfp_50_wilps2 <- mean(c(gfp_50_21_clean,gfp_50_22_clean))/normalising_factor2
gfp_50_wilps3 <- mean(c(gfp_50_31_clean,gfp_50_32_clean))/normalising_factor3

gfp_100_wilps1 <- mean(c(gfp_100_11_clean,gfp_100_12_clean))/normalising_factor1
gfp_100_wilps2 <- mean(c(gfp_100_21_clean,gfp_100_22_clean))/normalising_factor2
gfp_100_wilps3 <- mean(c(gfp_100_31_clean,gfp_100_32_clean))/normalising_factor3


LPS <- c('Without treatment','Without treatment','Without treatment',
         'Without treatment','Without treatment','Without treatment',
         'Without treatment','Without treatment','Without treatment',
         'Without treatment','Without treatment','Without treatment')
Concentration <- c('0 ug/mL','0 ug/mL','0 ug/mL','25 ug/mL','25 ug/mL','25 ug/mL','50 ug/mL','50 ug/mL','50 ug/mL','100 ug/mL','100 ug/mL','100 ug/mL')
Integrated_Intensity <- c(cntrl_wilps1,cntrl_wilps2,cntrl_wilps3, 
                          gfp_25_wilps1,gfp_25_wilps2,gfp_25_wilps3,
                          gfp_50_wilps1,gfp_50_wilps2,gfp_50_wilps3,
                          gfp_100_wilps1,gfp_100_wilps2,gfp_100_wilps3)

WOLPS_041<-data.frame(Concentration,Integrated_Intensity,LPS)


Integrated_Intensity_nonorm<- c(mean(c(cntrl11_clean,cntrl12_clean)), mean(c(cntrl21_clean,cntrl22_clean)), mean(c(cntrl31_clean,cntrl32_clean)), 
                                mean(c(gfp_25_11_clean,gfp_25_12_clean)), mean(c(gfp_25_21_clean,gfp_25_22_clean)), mean(c(gfp_25_31_clean,gfp_25_32_clean)),
                                mean(c(gfp_50_11_clean,gfp_50_12_clean)), mean(c(gfp_50_21_clean,gfp_50_22_clean)),mean(c(gfp_50_31_clean,gfp_50_32_clean)),
                                mean(c(gfp_100_11_clean,gfp_100_12_clean)), mean(c(gfp_100_21_clean,gfp_100_22_clean)),mean(c(gfp_100_31_clean,gfp_100_32_clean)))
WOLPS_041_nonorm<-data.frame(Concentration,Integrated_Intensity_nonorm,LPS)







#---------------------------  FV3-AF488-------------------------------------

cytoplasm <- read_csv("/Volumes/DG_Backup_E/SAPNAS/2021-09-29_BV2+041GFP_and_FV3/BV2+041GFP/data_output/final/BV2_FV3_GFP_FINALCytoplasm.csv")
cytoplasm$TYPE <- c('cytoplasm')
remove_outlayer<-function(x){
  x1<-x[!x %in% boxplot.stats(x)$out]
  return(x1)
}

names(cytoplasm)[names(cytoplasm) == "Metadata_Concentration...6"] <- "Metadata_Concentration"
cytoplasm$Metadata_Concentration <- factor(cytoplasm$Metadata_Concentration , levels=c("control", "25ugmL", "50ugmL", "100ugmL"))


data_woLPS<- filter(cytoplasm, Metadata_LPS == 'FV3SiR_wo')
data_wiLPS<- filter(cytoplasm, Metadata_LPS == 'FV3SiR_wi')



mean(cytoplasm$Intensity_MeanIntensity_GFP)
mean(cytoplasm$Intensity_MeanIntensity_GFP_aftermath)

mean(cytoplasm$Intensity_MedianIntensity_GFP)
mean(cytoplasm$Intensity_MedianIntensity_GFP_aftermath)

hist(cytoplasm$Intensity_MeanIntensity_GFP)
hist(cytoplasm$Intensity_MeanIntensity_GFP_aftermath)

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

cntrl11_clean<-mean(remove_outlayer(cntrl11$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl12_clean<-mean(remove_outlayer(cntrl12$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl21_clean<-mean(remove_outlayer(cntrl21$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl22_clean<-mean(remove_outlayer(cntrl22$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl31_clean<-mean(remove_outlayer(cntrl31$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl32_clean<-mean(remove_outlayer(cntrl32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_25_11_clean<-mean(remove_outlayer(gfp_25_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_12_clean<-mean(remove_outlayer(gfp_25_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_21_clean<-mean(remove_outlayer(gfp_25_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_22_clean<-mean(remove_outlayer(gfp_25_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_31_clean<-mean(remove_outlayer(gfp_25_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_32_clean<-mean(remove_outlayer(gfp_25_32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_50_11_clean<-mean(remove_outlayer(gfp_50_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_12_clean<-mean(remove_outlayer(gfp_50_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_21_clean<-mean(remove_outlayer(gfp_50_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_22_clean<-mean(remove_outlayer(gfp_50_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_31_clean<-mean(remove_outlayer(gfp_50_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_32_clean<-mean(remove_outlayer(gfp_50_32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_100_11_clean<-mean(remove_outlayer(gfp_100_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_12_clean<-mean(remove_outlayer(gfp_100_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_21_clean<-mean(remove_outlayer(gfp_100_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_22_clean<-mean(remove_outlayer(gfp_100_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_31_clean<-mean(remove_outlayer(gfp_100_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_32_clean<-mean(remove_outlayer(gfp_100_32$Intensity_IntegratedIntensity_GFP_aftermath))

# normalise by experiment

normalising_factor1 <- mean(c(cntrl11_clean,cntrl12_clean))+0.0001
normalising_factor2 <- mean(c(cntrl21_clean,cntrl22_clean))+0.0001
normalising_factor3 <- mean(c(cntrl31_clean,cntrl32_clean))+0.0001

cntrl_wilps1 <- normalising_factor1/normalising_factor1
cntrl_wilps2 <- normalising_factor2/normalising_factor2
cntrl_wilps3 <- normalising_factor3/normalising_factor3

gfp_25_wilps1 <- mean(c(gfp_25_11_clean,gfp_25_12_clean))/normalising_factor1
gfp_25_wilps2 <- mean(c(gfp_25_21_clean,gfp_25_22_clean))/normalising_factor2
gfp_25_wilps3 <- mean(c(gfp_25_31_clean,gfp_25_32_clean))/normalising_factor3

gfp_50_wilps1 <- mean(c(gfp_50_11_clean,gfp_50_12_clean))/normalising_factor1
gfp_50_wilps2 <- mean(c(gfp_50_21_clean,gfp_50_22_clean))/normalising_factor2
gfp_50_wilps3 <- mean(c(gfp_50_31_clean,gfp_50_32_clean))/normalising_factor3

gfp_100_wilps1 <- mean(c(gfp_100_11_clean,gfp_100_12_clean))/normalising_factor1
gfp_100_wilps2 <- mean(c(gfp_100_21_clean,gfp_100_22_clean))/normalising_factor2
gfp_100_wilps3 <- mean(c(gfp_100_31_clean,gfp_100_32_clean))/normalising_factor3

LPS <- c('100ng/mL LPS','100ng/mL LPS','100ng/mL LPS',
         '100ng/mL LPS','100ng/mL LPS','100ng/mL LPS',
         '100ng/mL LPS','100ng/mL LPS','100ng/mL LPS',
         '100ng/mL LPS','100ng/mL LPS','100ng/mL LPS')
Concentration <- c('0 ug/mL','0 ug/mL','0 ug/mL','25 ug/mL','25 ug/mL','25 ug/mL','50 ug/mL','50 ug/mL','50 ug/mL','100 ug/mL','100 ug/mL','100 ug/mL')
Integrated_Intensity <- c(cntrl_wilps1,cntrl_wilps2,cntrl_wilps3, 
                          gfp_25_wilps1,gfp_25_wilps2,gfp_25_wilps3,
                          gfp_50_wilps1,gfp_50_wilps2,gfp_50_wilps3,
                          gfp_100_wilps1,gfp_100_wilps2,gfp_100_wilps3)

WILPS_fv3<-data.frame(Concentration,Integrated_Intensity,LPS)



Integrated_Intensity_nonorm<- c(mean(c(cntrl11_clean,cntrl12_clean)), mean(c(cntrl21_clean,cntrl22_clean)), mean(c(cntrl31_clean,cntrl32_clean)), 
                                mean(c(gfp_25_11_clean,gfp_25_12_clean)), mean(c(gfp_25_21_clean,gfp_25_22_clean)), mean(c(gfp_25_31_clean,gfp_25_32_clean)),
                                mean(c(gfp_50_11_clean,gfp_50_12_clean)), mean(c(gfp_50_21_clean,gfp_50_22_clean)),mean(c(gfp_50_31_clean,gfp_50_32_clean)),
                                mean(c(gfp_100_11_clean,gfp_100_12_clean)), mean(c(gfp_100_21_clean,gfp_100_22_clean)),mean(c(gfp_100_31_clean,gfp_100_32_clean)))
WILPS_fv3_nonorm<-data.frame(Concentration,Integrated_Intensity_nonorm,LPS)

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

cntrl11_clean<-mean(remove_outlayer(cntrl11$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl12_clean<-mean(remove_outlayer(cntrl12$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl21_clean<-mean(remove_outlayer(cntrl21$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl22_clean<-mean(remove_outlayer(cntrl22$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl31_clean<-mean(remove_outlayer(cntrl31$Intensity_IntegratedIntensity_GFP_aftermath))
cntrl32_clean<-mean(remove_outlayer(cntrl32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_25_11_clean<-mean(remove_outlayer(gfp_25_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_12_clean<-mean(remove_outlayer(gfp_25_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_21_clean<-mean(remove_outlayer(gfp_25_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_22_clean<-mean(remove_outlayer(gfp_25_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_31_clean<-mean(remove_outlayer(gfp_25_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_25_32_clean<-mean(remove_outlayer(gfp_25_32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_50_11_clean<-mean(remove_outlayer(gfp_50_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_12_clean<-mean(remove_outlayer(gfp_50_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_21_clean<-mean(remove_outlayer(gfp_50_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_22_clean<-mean(remove_outlayer(gfp_50_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_31_clean<-mean(remove_outlayer(gfp_50_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_50_32_clean<-mean(remove_outlayer(gfp_50_32$Intensity_IntegratedIntensity_GFP_aftermath))

gfp_100_11_clean<-mean(remove_outlayer(gfp_100_11$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_12_clean<-mean(remove_outlayer(gfp_100_12$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_21_clean<-mean(remove_outlayer(gfp_100_21$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_22_clean<-mean(remove_outlayer(gfp_100_22$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_31_clean<-mean(remove_outlayer(gfp_100_31$Intensity_IntegratedIntensity_GFP_aftermath))
gfp_100_32_clean<-mean(remove_outlayer(gfp_100_32$Intensity_IntegratedIntensity_GFP_aftermath))

# normalise by experiment

normalising_factor1 <- mean(c(cntrl11_clean,cntrl12_clean))+0.0001
normalising_factor2 <- mean(c(cntrl21_clean,cntrl22_clean))+0.0001
normalising_factor3 <- mean(c(cntrl31_clean,cntrl32_clean))+0.0001

cntrl_wilps1 <- normalising_factor1/normalising_factor1
cntrl_wilps2 <- normalising_factor2/normalising_factor2
cntrl_wilps3 <- normalising_factor3/normalising_factor3

gfp_25_wilps1 <- mean(c(gfp_25_11_clean,gfp_25_12_clean))/normalising_factor1
gfp_25_wilps2 <- mean(c(gfp_25_21_clean,gfp_25_22_clean))/normalising_factor2
gfp_25_wilps3 <- mean(c(gfp_25_31_clean,gfp_25_32_clean))/normalising_factor3

gfp_50_wilps1 <- mean(c(gfp_50_11_clean,gfp_50_12_clean))/normalising_factor1
gfp_50_wilps2 <- mean(c(gfp_50_21_clean,gfp_50_22_clean))/normalising_factor2
gfp_50_wilps3 <- mean(c(gfp_50_31_clean,gfp_50_32_clean))/normalising_factor3

gfp_100_wilps1 <- mean(c(gfp_100_11_clean,gfp_100_12_clean))/normalising_factor1
gfp_100_wilps2 <- mean(c(gfp_100_21_clean,gfp_100_22_clean))/normalising_factor2
gfp_100_wilps3 <- mean(c(gfp_100_31_clean,gfp_100_32_clean))/normalising_factor3


LPS <- c('Without treatment','Without treatment','Without treatment',
         'Without treatment','Without treatment','Without treatment',
         'Without treatment','Without treatment','Without treatment',
         'Without treatment','Without treatment','Without treatment')
Concentration <- c('0 ug/mL','0 ug/mL','0 ug/mL','25 ug/mL','25 ug/mL','25 ug/mL','50 ug/mL','50 ug/mL','50 ug/mL','100 ug/mL','100 ug/mL','100 ug/mL')
Integrated_Intensity <- c(cntrl_wilps1,cntrl_wilps2,cntrl_wilps3, 
                          gfp_25_wilps1,gfp_25_wilps2,gfp_25_wilps3,
                          gfp_50_wilps1,gfp_50_wilps2,gfp_50_wilps3,
                          gfp_100_wilps1,gfp_100_wilps2,gfp_100_wilps3)

WOLPS_fv3<-data.frame(Concentration,Integrated_Intensity,LPS)


Integrated_Intensity_nonorm<- c(mean(c(cntrl11_clean,cntrl12_clean)), mean(c(cntrl21_clean,cntrl22_clean)), mean(c(cntrl31_clean,cntrl32_clean)), 
                                mean(c(gfp_25_11_clean,gfp_25_12_clean)), mean(c(gfp_25_21_clean,gfp_25_22_clean)), mean(c(gfp_25_31_clean,gfp_25_32_clean)),
                                mean(c(gfp_50_11_clean,gfp_50_12_clean)), mean(c(gfp_50_21_clean,gfp_50_22_clean)),mean(c(gfp_50_31_clean,gfp_50_32_clean)),
                                mean(c(gfp_100_11_clean,gfp_100_12_clean)), mean(c(gfp_100_21_clean,gfp_100_22_clean)),mean(c(gfp_100_31_clean,gfp_100_32_clean)))
WOLPS_fv3_nonorm<-data.frame(Concentration,Integrated_Intensity_nonorm,LPS)










#----------- FINALISE - CYTOPLASM----------------------------------------------------------

#041
#
FInal.Frame.041<- rbind(WOLPS_041,WILPS_041)
FInal.Frame.041$Concentration <- factor(FInal.Frame.041$Concentration , levels=c("0 ug/mL", "25 ug/mL", "50 ug/mL", "100 ug/mL"))
FInal.Frame.041 <- filter(FInal.Frame.041, Concentration != '0 ug/mL')
FInal.Frame.041.b<-filter(FInal.Frame.041, LPS == 'Without treatment')
FInal.Frame.041.c<-filter(FInal.Frame.041, LPS != 'Without treatment')
FInal.Frame.041_25 <- filter(FInal.Frame.041, Concentration == '25 ug/mL')
FInal.Frame.041_50 <- filter(FInal.Frame.041, Concentration == '50 ug/mL')
FInal.Frame.041_100 <- filter(FInal.Frame.041, Concentration == '100 ug/mL')



a<-ggplot(FInal.Frame.041, aes(x=Concentration,y=Integrated_Intensity, color=LPS))+
  geom_boxplot()+  ggtitle('041-GFP (cytoplasm) n=3 normalised' )+
  xlab(c('041-GFP concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point(position=position_jitterdodge())

b<-ggplot(FInal.Frame.041.b, aes(x=Concentration,y=Integrated_Intensity,color=LPS))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  ggtitle('Without treatment' )+
  xlab(c('041-GFP concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()

c<-ggplot(FInal.Frame.041.c, aes(x=Concentration,y=Integrated_Intensity,color=LPS))+
  geom_boxplot()+  ggtitle('100 ng/mL LPS' )+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  xlab(c('041-GFP concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()

figure1 <- ggarrange(a, b, c,
                    labels = c("A", "B", "C"),
                    ncol = 1, nrow = 3)
figure1

d<-ggplot(FInal.Frame.041_25, aes(x=LPS,y=Integrated_Intensity))+
  geom_boxplot()+  ggtitle('25 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()

e<-ggplot(FInal.Frame.041_50, aes(x=LPS,y=Integrated_Intensity))+
  geom_boxplot()+  ggtitle('25 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()
f<-ggplot(FInal.Frame.041_100, aes(x=LPS,y=Integrated_Intensity))+
  geom_boxplot()+  ggtitle('25 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()

figure1.5 <- ggarrange(d,e,f,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3)
figure1.5

# Finalised data w/o normalisation 

FInal.Frame.041_nonorm<- rbind(WOLPS_041_nonorm,WILPS_041_nonorm)
FInal.Frame.041_nonorm$Concentration <- factor(FInal.Frame.041_nonorm$Concentration , levels=c("0 ug/mL", "25 ug/mL", "50 ug/mL", "100 ug/mL"))
FInal.Frame.041_nonorm <- filter(FInal.Frame.041_nonorm, Concentration != '0 ug/mL')
FInal.Frame.041.b_nonorm<-filter(FInal.Frame.041_nonorm, LPS == 'Without treatment')
FInal.Frame.041.c_nonorm<-filter(FInal.Frame.041_nonorm, LPS != 'Without treatment')
FInal.Frame.041_25_nonorm <- filter(FInal.Frame.041_nonorm, Concentration == '25 ug/mL')
FInal.Frame.041_50_nonorm <- filter(FInal.Frame.041_nonorm, Concentration == '50 ug/mL')
FInal.Frame.041_100_nonorm <- filter(FInal.Frame.041_nonorm, Concentration == '100 ug/mL')


aaa <- ggplot(FInal.Frame.041_nonorm, aes(x=Concentration,y=Integrated_Intensity_nonorm, color=LPS))+
  geom_boxplot()+  ggtitle('041-GFP (cytoplasm) n=3' )+
  xlab(c('041-GFP concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point(position=position_jitterdodge())

bbb <- ggplot(FInal.Frame.041.b_nonorm, aes(x=Concentration,y=Integrated_Intensity_nonorm,color=LPS))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  ggtitle('Without treatment' )+
  xlab(c('041-GFP concentration'))+
  ylab(c('Fluorescence intensity '))+
  theme_bw()+geom_point()

ccc <- ggplot(FInal.Frame.041.c_nonorm, aes(x=Concentration,y=Integrated_Intensity_nonorm,color=LPS))+
  geom_boxplot()+  ggtitle('100 ng/mL LPS' )+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  xlab(c('041-GFP concentration'))+
  ylab(c('Fluorescence intensity '))+
  theme_bw()+geom_point()

figure2 <- ggarrange(aaa, bbb, ccc,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3)
figure2

dd<-ggplot(FInal.Frame.041_25_nonorm, aes(x=LPS,y=Integrated_Intensity_nonorm))+
  geom_boxplot()+  ggtitle('25 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity  '))+
  theme_bw()+geom_point()

ee<-ggplot(FInal.Frame.041_50_nonorm, aes(x=LPS,y=Integrated_Intensity_nonorm))+
  geom_boxplot()+  ggtitle('50 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity  '))+
  theme_bw()+geom_point()
ff<-ggplot(FInal.Frame.041_100_nonorm, aes(x=LPS,y=Integrated_Intensity_nonorm))+
  geom_boxplot()+  ggtitle('100 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity  '))+
  theme_bw()+geom_point()

figure2.5 <- ggarrange(dd,ee,ff,
                       labels = c("A", "B", "C"),
                       ncol = 1, nrow = 3)
figure2.5

# FV3
#

FInal.Frame.fv3<- rbind(WOLPS_fv3,WILPS_fv3)
FInal.Frame.fv3$Concentration <- factor(FInal.Frame.fv3$Concentration , levels=c("0 ug/mL", "25 ug/mL", "50 ug/mL", "100 ug/mL"))
FInal.Frame.fv3 <- filter(FInal.Frame.fv3, Concentration != '0 ug/mL')
FInal.Frame.fv3.b<-filter(FInal.Frame.fv3, LPS == 'Without treatment')
FInal.Frame.fv3.c<-filter(FInal.Frame.fv3, LPS != 'Without treatment')
FInal.Frame.fv3_25 <- filter(FInal.Frame.fv3, Concentration == '25 ug/mL')
FInal.Frame.fv3_50 <- filter(FInal.Frame.fv3, Concentration == '50 ug/mL')
FInal.Frame.fv3_100 <- filter(FInal.Frame.fv3, Concentration == '100 ug/mL')

aa <- ggplot(FInal.Frame.fv3, aes(x=Concentration,y=Integrated_Intensity, color=LPS))+
  geom_boxplot()+  ggtitle('FV3-AF488 (cytoplasm) n=3 normalised' )+
  xlab(c('FV3-AF488 concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point(position=position_jitterdodge())

bb <- ggplot(FInal.Frame.fv3.b, aes(x=Concentration,y=Integrated_Intensity,color=LPS))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, method = 't.test')+
  ggtitle('Without treatment' )+
  xlab(c('FV3-AF488 concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()

cc <- ggplot(FInal.Frame.fv3.c, aes(x=Concentration,y=Integrated_Intensity,color=LPS))+
  geom_boxplot()+  ggtitle('100 ng/mL LPS' )+
  stat_compare_means(comparisons = my_comparisons, method = 't.test')+
  xlab(c('FV3-AF488 concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()

figure3 <- ggarrange(aa, bb, cc,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3)
figure3

ddd<-ggplot(FInal.Frame.fv3_25, aes(x=LPS,y=Integrated_Intensity))+
  geom_boxplot()+  ggtitle('25 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()

eee<-ggplot(FInal.Frame.fv3_50, aes(x=LPS,y=Integrated_Intensity))+
  geom_boxplot()+  ggtitle('50 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()
fff<-ggplot(FInal.Frame.fv3_100, aes(x=LPS,y=Integrated_Intensity))+
  geom_boxplot()+  ggtitle('100 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point()

figure3.5 <- ggarrange(ddd,eee,fff,
                       labels = c("A", "B", "C"),
                       ncol = 1, nrow = 3)
figure3.5

# Finalised data w/o normalisation 

FInal.Frame.fv3_nonorm<- rbind(WOLPS_fv3_nonorm,WILPS_fv3_nonorm)
FInal.Frame.fv3_nonorm$Concentration <- factor(FInal.Frame.fv3_nonorm$Concentration , levels=c("0 ug/mL", "25 ug/mL", "50 ug/mL", "100 ug/mL"))
FInal.Frame.fv3_nonorm <- filter(FInal.Frame.fv3_nonorm, Concentration != '0 ug/mL')
FInal.Frame.fv3.b_nonorm<-filter(FInal.Frame.fv3_nonorm, LPS == 'Without treatment')
FInal.Frame.fv3.c_nonorm<-filter(FInal.Frame.fv3_nonorm, LPS != 'Without treatment')
FInal.Frame.fv3_25_nonorm <- filter(FInal.Frame.fv3_nonorm, Concentration == '25 ug/mL')
FInal.Frame.fv3_50_nonorm <- filter(FInal.Frame.fv3_nonorm, Concentration == '50 ug/mL')
FInal.Frame.fv3_100_nonorm <- filter(FInal.Frame.fv3_nonorm, Concentration == '100 ug/mL')

aaaa <- ggplot(FInal.Frame.fv3_nonorm, aes(x=Concentration,y=Integrated_Intensity_nonorm, color=LPS))+
  geom_boxplot()+  ggtitle('FV3-AF488 (cytoplasm) n=3 ' )+
  xlab(c('FV3-AF488 concentration'))+
  ylab(c('Fluorescence intensity fold increase'))+
  theme_bw()+geom_point(position=position_jitterdodge())

bbbb <- ggplot(FInal.Frame.fv3.b_nonorm, aes(x=Concentration,y=Integrated_Intensity_nonorm,color=LPS))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  ggtitle('Without treatment' )+
  xlab(c('FV3-AF488 concentration'))+
  ylab(c('Fluorescence intensity  '))+
  theme_bw()+geom_point()

cccc <- ggplot(FInal.Frame.fv3.c_nonorm, aes(x=Concentration,y=Integrated_Intensity_nonorm,color=LPS))+
  geom_boxplot()+  ggtitle('100 ng/mL LPS' )+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  xlab(c('FV3-AF488 concentration'))+
  ylab(c('Fluorescence intensity  '))+
  theme_bw()+geom_point()

figure4 <- ggarrange(aaaa, bbbb, cccc,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3)
figure4

dddd<-ggplot(FInal.Frame.fv3_25, aes(x=LPS,y=Integrated_Intensity_nonorm))+
  geom_boxplot()+  ggtitle('25 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity  '))+
  theme_bw()+geom_point()

eeee<-ggplot(FInal.Frame.fv3_50, aes(x=LPS,y=Integrated_Intensity_nonorm))+
  geom_boxplot()+  ggtitle('50 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity  '))+
  theme_bw()+geom_point()
ffff<-ggplot(FInal.Frame.fv3_100, aes(x=LPS,y=Integrated_Intensity_nonorm))+
  geom_boxplot()+  ggtitle('100 ug/mL 041-GFP')+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test")+
  xlab(c('LPS concentration'))+
  ylab(c('Fluorescence intensity  '))+
  theme_bw()+geom_point()

figure4.5 <- ggarrange(ddd,eee,fff,
                       labels = c("A", "B", "C"),
                       ncol = 1, nrow = 3)
figure4.5

# also show how LPS effect works on different concentrations
#


#-------------------------------Analysing the signal from the membrane:------------
