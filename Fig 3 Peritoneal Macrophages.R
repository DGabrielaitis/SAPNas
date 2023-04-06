# 2022 05 25

# SAPNAS
#Figure 3 
# Peritoneal macrophages incubated with 041-GFP and 053-SNAP-AF488
# Dovydas Gabrielaitis

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls()) 

#------------libraries and functinos --------
library(readr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(doBy)
library(svMisc)
library(data.table)
library(sciplot)



# quick and dirty normalisation function
normalize_values <- function(df) {
  # extract all the experiment numbers
  all_exp_numbers <- unique(df$Metadata_ExperimentNumber)
  
  new_df <- data.frame()
  #loop through the experiments and normalize values; normalized values are saved in a new data frame together with metadata
  for (i in all_exp_numbers) {
    temp_df <- df[df$Metadata_ExperimentNumber==i,1:4]
    temp_df$Normalized_values <- df[df$Metadata_ExperimentNumber==i,4]/
      as.double(df[df$Metadata_ExperimentNumber==i&df$Metadata_Treatment=='CNT',4])
    new_df <- rbind(new_df,temp_df)
    temp_df <- data.frame()
  }
  new_df
}


#----------------------------Peritoneal macrophages ------------------

# import data
#TOTAL CELL
#DAT <- read_csv("data_out/TotalCell_filtered.csv")
#Cytoplasm 
DAT <- read_csv("data_out/Cytoplasm_filtered.csv")
# Membrane
#DAT <- read_csv("data_out/Membrane_filtered.csv")
#



# if you want to select the data by age or sex, un-comment these lines:
#only old animals
#DAT_old <- filter(DAT, Metadata_Age == 'OLD' )
#only young animals
#DAT <- filter(DAT, Metadata_Age == 'OLD' )
#DAT <- filter(DAT, Metadata_Sex == 'FEM' )

#filter 041
PM_041<-filter(DAT, Metadata_Protein == '41' )
#Filter only fv3
PM_fv3<-filter(DAT, Metadata_Protein == 'fv3' )


# filter out bad experiments (in this case nly the CHD portion of the experiment):

#PM_fv3<-filter(PM_fv3, Metadata_Treatment != 'CHD' | Metadata_ExperimentNumber != '15' )
#PM_fv3<-filter(PM_fv3, Metadata_Treatment != 'CHD' | Metadata_ExperimentNumber != '13' )
#PM_fv3<-filter(PM_fv3, Metadata_Treatment != 'CHD' | Metadata_ExperimentNumber != '2' )
#PM_fv3<-filter(PM_fv3, Metadata_Treatment != 'CHD' | Metadata_ExperimentNumber != '3' )

# or remove the whole experiment:
PM_fv3<-filter(PM_fv3, Metadata_ExperimentNumber != '15' )
PM_fv3<-filter(PM_fv3,  Metadata_ExperimentNumber != '13' )
PM_fv3<-filter(PM_fv3,  Metadata_ExperimentNumber != '2' )
PM_fv3<-filter(PM_fv3,  Metadata_ExperimentNumber != '3' )




PM_041<-filter(PM_041, Metadata_ExperimentNumber != '8' )


# how many experiments are left?
table(PM_fv3$Metadata_ExperimentNumber) # 12 experiments are left (YM 3 YF 4 OM 3  OF 2)

# ----average 041 data and remove outlayers
int1<-setDT(PM_041)[, list(Intensity_IntegratedIntensity_Nanotube_Rescaled= mean(Intensity_IntegratedIntensity_Nanotube_Rescaled)) , .(Metadata_ImageNumber, Metadata_WellNumber,Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age)]
int2<-setDT(int1)[, list(Intensity_IntegratedIntensity_Nanotube_Rescaled= mean(Intensity_IntegratedIntensity_Nanotube_Rescaled)) , .( Metadata_WellNumber,Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age)]
int3<-setDT(int2)[, list(Intensity_IntegratedIntensity_Nanotube_Rescaled= mean(Intensity_IntegratedIntensity_Nanotube_Rescaled)) , .( Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age)]
# ----Remove outlayers for 041 data
# Split dataframe by treatment ant NT nocentration 
int1_YNG_CNT <- int1 %>%filter(Metadata_Treatment %in% c("CNT"), Metadata_Age == 'YNG' )
int1_YNG_LPS<- int1 %>%filter(Metadata_Treatment %in% c("LPS"), Metadata_Age == 'YNG' )
int1_YNF_CHD <- int1 %>%filter(Metadata_Treatment %in% c("CHD"), Metadata_Age == 'YNG' )  
int1_OLD_CNT <- int1 %>%filter(Metadata_Treatment %in% c("CNT"), Metadata_Age == 'OLD' )
int1_OLD_LPS <- int1 %>%filter(Metadata_Treatment %in% c("LPS"), Metadata_Age == 'OLD' )
int1_OLD_CHD <- int1 %>%filter(Metadata_Treatment %in% c("CHD"), Metadata_Age == 'OLD' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int1_YNG_CNT$is_outlier <- ifelse(int1_YNG_CNT$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int1_YNG_CNT$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int1_YNG_LPS$is_outlier <- ifelse(int1_YNG_LPS$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int1_YNG_LPS$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int1_YNF_CHD$is_outlier <- ifelse(int1_YNF_CHD$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int1_YNF_CHD$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int1_OLD_CNT$is_outlier <- ifelse(int1_OLD_CNT$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int1_OLD_CNT$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int1_OLD_LPS$is_outlier <- ifelse(int1_OLD_LPS$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int1_OLD_LPS$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int1_OLD_CHD$is_outlier <- ifelse(int1_OLD_CHD$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int1_OLD_CHD$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)

# Bind the DF back
DATA_PM_041 <- rbind(int1_YNG_CNT,
                     int1_YNG_LPS,
                     int1_YNF_CHD,
                     int1_OLD_CNT,
                     int1_OLD_LPS,
                     int1_OLD_CHD)

# Remove the values that have outlayers
DATA_PM_041<- DATA_PM_041 %>%filter(is_outlier %in% c("0"))
#-------End of remove outlayer


#----------------------------Peritoneal macrophages - 053-SNAP-AF488 DATA----


# average fv3 data and remove outlayers
int5<-setDT(PM_fv3)[, list(Intensity_IntegratedIntensity_Nanotube_Rescaled= mean(Intensity_IntegratedIntensity_Nanotube_Rescaled)) , .(Metadata_ImageNumber, Metadata_WellNumber,Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age)]
int6<-setDT(int5)[, list(Intensity_IntegratedIntensity_Nanotube_Rescaled= mean(Intensity_IntegratedIntensity_Nanotube_Rescaled)) , .( Metadata_WellNumber,Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age)]
int7<-setDT(int6)[, list(Intensity_IntegratedIntensity_Nanotube_Rescaled= mean(Intensity_IntegratedIntensity_Nanotube_Rescaled)) , .( Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age)]
# ----Remove outlayers for fv3 data
# Split dataframe by treatment ant NT nocentration 
int5_YNG_CNT <- int5 %>%filter(Metadata_Treatment %in% c("CNT"), Metadata_Age == 'YNG' )
int5_YNG_LPS<- int5 %>%filter(Metadata_Treatment %in% c("LPS"), Metadata_Age == 'YNG' )
int5_YNF_CHD <- int5 %>%filter(Metadata_Treatment %in% c("CHD"), Metadata_Age == 'YNG' )  
int5_OLD_CNT <- int5 %>%filter(Metadata_Treatment %in% c("CNT"), Metadata_Age == 'OLD' )
int5_OLD_LPS <- int5 %>%filter(Metadata_Treatment %in% c("LPS"), Metadata_Age == 'OLD' )
int5_OLD_CHD <- int5 %>%filter(Metadata_Treatment %in% c("CHD"), Metadata_Age == 'OLD' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int5_YNG_CNT$is_outlier <- ifelse(int5_YNG_CNT$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int5_YNG_CNT$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int5_YNG_LPS$is_outlier <- ifelse(int5_YNG_LPS$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int5_YNG_LPS$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int5_YNF_CHD$is_outlier <- ifelse(int5_YNF_CHD$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int5_YNF_CHD$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int5_OLD_CNT$is_outlier <- ifelse(int5_OLD_CNT$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int5_OLD_CNT$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int5_OLD_LPS$is_outlier <- ifelse(int5_OLD_LPS$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int5_OLD_LPS$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)
int5_OLD_CHD$is_outlier <- ifelse(int5_OLD_CHD$Intensity_IntegratedIntensity_Nanotube_Rescaled %in% boxplot.stats(int5_OLD_CHD$Intensity_IntegratedIntensity_Nanotube_Rescaled)$out, 1, 0)

# Bind the DF back
DATA_PM_fv3 <- rbind(int5_YNG_CNT,
                     int5_YNG_LPS,
                     int5_YNF_CHD,
                     int5_OLD_CNT,
                     int5_OLD_LPS,
                     int5_OLD_CHD)

# Remove the values that have outlayers
DATA_PM_fv3<- DATA_PM_fv3 %>%filter(is_outlier %in% c("0"))
#-------End of remove outlayer


#------------ PLOT DATA----------------


#normalise by experiment mean and plot 041 data
x1<- normalize_values(int3)
ggplot(x1, aes(x=Metadata_Treatment,y=Normalized_values, color= Metadata_Age))+
  geom_boxplot()+geom_point( pch = 4,position=position_jitterdodge()) +
  #stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  labs(title = "PERITONEAL MACROPHAGES + 041-GFP", x = "Treatment", y = "Fold integrated intensity")+
  theme_bw()+ labs(color='Mice Age')+ylim(c(0,2.5))+
  theme(axis.text = element_text(size = 17))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size = 20)) 





#normalise and plot fv3 data
x2<- normalize_values(int7)
ggplot(x2, aes(x=Metadata_Treatment,y=Normalized_values, color= Metadata_Age))+
  geom_boxplot()+geom_point(pch = 4, position=position_jitterdodge()) +
  #stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  labs(title = "PERITONEAL MACROPHAGES + FV3-AF488", x = "Treatment", y = "Fold integrated intensity")+
  theme_bw()+
  theme(axis.text = element_text(size = 17))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size = 20)) + ylim(c(0,4))


# Plot every experiment:



# 041:

# ----average 041 data and remove outlayers

ggplot(x1, aes(x=Metadata_Treatment,y=Normalized_values, color= Metadata_ExperimentNumber,group = Metadata_ExperimentNumber))+
  geom_point( size = 3) + geom_line()+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  labs(title = "PERITONEAL MACROPHAGES + 041-GFP", x = "Treatment", y = "Fold integrated intensity")+
  theme_bw() + labs(color='EXPERIMENT NUMBER') + geom_hline(yintercept=1, linetype="dashed", color = "red")+ylim(c(0,4.5))


ggplot(x2, aes(x=Metadata_Treatment,y=Normalized_values, color= Metadata_ExperimentNumber,group = Metadata_ExperimentNumber))+
  geom_point( size = 3) + geom_line()+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  labs(title = "PERITONEAL MACROPHAGES + FV3-AF488", x = "Treatment", y = "Fold integrated intensity")+
  theme_bw() + labs(color='EXPERIMENT NUMBER') + geom_hline(yintercept=1, linetype="dashed", color = "red")+ylim(c(0,3))




#------------ make log2 data table for 053SNAP
library(Rmisc)
x2 # this is the normalised 053 data by experiemnt
x2$log2_dat<- log2(x2$Normalized_values)
x2_stat_log2 <- summarySE(x2, measurevar="log2_dat", groupvars=c("Metadata_Age","Metadata_Treatment"))
x2_stat_log2

#------------ make log2 data table for 041GFP

x1 # this is the normalised 053 data by experiemnt
x1$log2_dat<- log2(x1$Normalized_values)
x1_stat_log2 <- summarySE(x1, measurevar="log2_dat", groupvars=c("Metadata_Age","Metadata_Treatment"))
x1_stat_log2
