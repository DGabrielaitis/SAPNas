# 2022 05 26

# SAPNAS
#Figure S1
# Microglia incubated with 041-GFP and 053-SNAP-AF488

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



#-------------Microglia----

#Import Total Cell data
#DAT_mg <- read_csv("data_out/MG_TotalCell_Filtered.csv")
#Import Cytoplasm
#DAT_mg <- read_csv("data_out/MG_Cytoplasm_Filtered.csv")
#Import Membrane
DAT_mg <- read_csv("data_out/MG_Membrane_Filtered.csv")


#filter by sex or age
#DAT_mg <- filter(DAT_mg, Metadata_Age...3 == 'OLD' )
#DAT_mg <- filter(DAT_mg, Metadata_Sex == 'FEM' )

#Filter the dataframe into two df's depending on the protein used
MG_041<-filter(DAT_mg, Metadata_Protein == '41' )
MG_fv3<-filter(DAT_mg, Metadata_Protein == 'fv3' )


# filter out bad experiments:

MG_fv3<-filter(MG_fv3, Metadata_ExperimentNumber != '6' )
MG_fv3<-filter(MG_fv3, Metadata_ExperimentNumber != '5' )
MG_fv3<-filter(MG_fv3, Metadata_ExperimentNumber != '11' )
MG_fv3<-filter(MG_fv3, Metadata_ExperimentNumber != '18' )

MG_041<-filter(MG_041, Metadata_ExperimentNumber != '6' )
MG_041<-filter(MG_041, Metadata_ExperimentNumber != '15' )
MG_041<-filter(MG_041, Metadata_ExperimentNumber != '18' )
MG_041<-filter(MG_041, Metadata_ExperimentNumber != '13' )
MG_041<-filter(MG_041, Metadata_ExperimentNumber != '8' )



#--------------- Microglia 041-GFP Data-----------------
# average 041 data
int8<-setDT(MG_041)[, list(Intensity_IntegratedIntensity_NanoTube_Rescaled= mean(Intensity_IntegratedIntensity_NanoTube_Rescaled)) , .(Metadata_ImageNumber, Metadata_WellNumber,Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age...3)]
int9<-setDT(int8)[, list(Intensity_IntegratedIntensity_NanoTube_Rescaled= mean(Intensity_IntegratedIntensity_NanoTube_Rescaled)) , .( Metadata_WellNumber,Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age...3)]
int10<-setDT(int9)[, list(Intensity_IntegratedIntensity_NanoTube_Rescaled= mean(Intensity_IntegratedIntensity_NanoTube_Rescaled)) , .( Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age...3)]

# ----Remove outlayers for 041 data
# Split dataframe by treatment ant NT nocentration 
int8_YNG_CNT <- int8 %>%filter(Metadata_Treatment %in% c("CNT"), Metadata_Age...3 == 'YNG' )
int8_YNG_LPS<- int8 %>%filter(Metadata_Treatment %in% c("LPS"), Metadata_Age...3 == 'YNG' )
int8_YNF_CHD <- int8 %>%filter(Metadata_Treatment %in% c("CHD"), Metadata_Age...3 == 'YNG' )  
int8_OLD_CNT <- int8%>%filter(Metadata_Treatment %in% c("CNT"), Metadata_Age...3 == 'OLD' )
int8_OLD_LPS <- int8 %>%filter(Metadata_Treatment %in% c("LPS"), Metadata_Age...3 == 'OLD' )
int8_OLD_CHD <- int8 %>%filter(Metadata_Treatment %in% c("CHD"), Metadata_Age...3 == 'OLD' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int8_YNG_CNT$is_outlier <- ifelse(int8_YNG_CNT$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int8_YNG_CNT$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int8_YNG_LPS$is_outlier <- ifelse(int8_YNG_LPS$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int8_YNG_LPS$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int8_YNF_CHD$is_outlier <- ifelse(int8_YNF_CHD$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int8_YNF_CHD$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int8_OLD_CNT$is_outlier <- ifelse(int8_OLD_CNT$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int8_OLD_CNT$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int8_OLD_LPS$is_outlier <- ifelse(int8_OLD_LPS$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int8_OLD_LPS$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int8_OLD_CHD$is_outlier <- ifelse(int8_OLD_CHD$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int8_OLD_CHD$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)

# Bind the DF back
DATA_MG_041 <- rbind(int8_YNG_CNT,
                     int8_YNG_LPS,
                     int8_YNF_CHD,
                     int8_OLD_CNT,
                     int8_OLD_LPS,
                     int8_OLD_CHD)

# Remove the values that have outlayers
DATA_MG_041<- DATA_MG_041 %>%filter(is_outlier %in% c("0"))
#-------End of remove outlayer

#--------------- Microglia 053-SNAP-AF488 Data-----------
# average fv3 data
int11<-setDT(MG_fv3)[, list(Intensity_IntegratedIntensity_NanoTube_Rescaled= mean(Intensity_IntegratedIntensity_NanoTube_Rescaled)) , .(Metadata_ImageNumber, Metadata_WellNumber,Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age...3)]
int12<-setDT(int11)[, list(Intensity_IntegratedIntensity_NanoTube_Rescaled= mean(Intensity_IntegratedIntensity_NanoTube_Rescaled)) , .( Metadata_WellNumber,Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age...3)]
int13<-setDT(int12)[, list(Intensity_IntegratedIntensity_NanoTube_Rescaled= mean(Intensity_IntegratedIntensity_NanoTube_Rescaled)) , .( Metadata_Treatment,Metadata_ExperimentNumber,Metadata_Age...3)]
# ----Remove outlayers for 041 data
# Split dataframe by treatment ant NT nocentration 
int11_YNG_CNT <- int11 %>%filter(Metadata_Treatment %in% c("CNT"), Metadata_Age...3 == 'YNG' )
int11_YNG_LPS<- int11 %>%filter(Metadata_Treatment %in% c("LPS"), Metadata_Age...3 == 'YNG' )
int11_YNF_CHD <- int11 %>%filter(Metadata_Treatment %in% c("CHD"), Metadata_Age...3 == 'YNG' )  
int11_OLD_CNT <- int11%>%filter(Metadata_Treatment %in% c("CNT"), Metadata_Age...3 == 'OLD' )
int11_OLD_LPS <- int11 %>%filter(Metadata_Treatment %in% c("LPS"), Metadata_Age...3 == 'OLD' )
int11_OLD_CHD <- int11 %>%filter(Metadata_Treatment %in% c("CHD"), Metadata_Age...3 == 'OLD' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int11_YNG_CNT$is_outlier <- ifelse(int11_YNG_CNT$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int11_YNG_CNT$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int11_YNG_LPS$is_outlier <- ifelse(int11_YNG_LPS$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int11_YNG_LPS$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int11_YNF_CHD$is_outlier <- ifelse(int11_YNF_CHD$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int11_YNF_CHD$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int11_OLD_CNT$is_outlier <- ifelse(int11_OLD_CNT$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int11_OLD_CNT$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int11_OLD_LPS$is_outlier <- ifelse(int11_OLD_LPS$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int11_OLD_LPS$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)
int11_OLD_CHD$is_outlier <- ifelse(int11_OLD_CHD$Intensity_IntegratedIntensity_NanoTube_Rescaled %in% boxplot.stats(int11_OLD_CHD$Intensity_IntegratedIntensity_NanoTube_Rescaled)$out, 1, 0)

# Bind the DF back
DATA_MG_fv3 <- rbind(int11_YNG_CNT,
                     int11_YNG_LPS,
                     int11_YNF_CHD,
                     int11_OLD_CNT,
                     int11_OLD_LPS,
                     int11_OLD_CHD)

# Remove the values that have outlayers
DATA_MG_fv3<- DATA_MG_fv3 %>%filter(is_outlier %in% c("0"))
#-------End of remove outlayer

# --> Go to other file for data plotting and statistics


#------------------- PLOT THE NORMALISED DATA--------------------------------------
x3<- normalize_values(int10)
x3$Metadata_Age...3 <- factor(x3$Metadata_Age...3, levels=c("YNG", "OLD"))
ggplot(x3, aes(x=Metadata_Treatment,y=Normalized_values, color= Metadata_Age...3))+
  geom_boxplot()+geom_point( pch = 4,size=8,position=position_jitterdodge()) +
  #stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  labs(title = "Microglia + 041-GFP", x = "Treatment", y = "Fold integrated intensity")+
  theme_bw()+ labs(color='Mice Age')+
  theme(axis.text = element_text(size = 17))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size = 20)) #+ ylim(c(0,3.5))



#normalise and plot fv3 data
x4<- normalize_values(int13)
x4$Metadata_Age...3 <- factor(x4$Metadata_Age...3, levels=c("YNG", "OLD"))
ggplot(x4, aes(x=Metadata_Treatment,y=Normalized_values, color= Metadata_Age...3))+
  geom_boxplot()+geom_point(pch = 4, size=8, position=position_jitterdodge()) +
  #stat_compare_means(comparisons = my_comparisons, method = "t.test")+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  labs(title = "Microglia + FV3-AF488", x = "Treatment", y = "Fold integrated intensity")+
  theme_bw()+
  theme(axis.text = element_text(size = 17))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size = 20)) #+ ylim(c(0,3.5))




# ---- Print the normalised values by experiment


ggplot(x3, aes(x=Metadata_Treatment,y=Normalized_values, color= Metadata_ExperimentNumber,group = Metadata_ExperimentNumber))+
  geom_point( size = 3) + geom_line()+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  labs(title = "Microglia + 041-GFP", x = "Treatment", y = "Fold integrated intensity")+
  theme_bw() + labs(color='EXPERIMENT NUMBER') + geom_hline(yintercept=1, linetype="dashed", color = "red")+ylim(c(0,5))
ggplot(x4, aes(x=Metadata_Treatment,y=Normalized_values, color= Metadata_ExperimentNumber,group = Metadata_ExperimentNumber))+
  geom_point( size = 3) + geom_line()+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  labs(title = "Microglia + FV3-AF488", x = "Treatment", y = "Fold integrated intensity")+
  theme_bw() + labs(color='EXPERIMENT NUMBER') + geom_hline(yintercept=1, linetype="dashed", color = "red")+ylim(c(0,10))

#-------
library(Rmisc)
x4 # 053SNAP
x4_stat <- summarySE(x4, measurevar="Normalized_values", groupvars=c("Metadata_Age...3","Metadata_Treatment"))
x4_stat

#041GFP
x3_stat <- summarySE(x3, measurevar="Normalized_values", groupvars=c("Metadata_Age...3","Metadata_Treatment"))
x3_stat

# log2 data - 053SNAP
x4$log2_dat<- log2(x4$Normalized_values)
x4_stat_log2 <- summarySE(x4, measurevar="log2_dat", groupvars=c("Metadata_Age...3","Metadata_Treatment"))
x4_stat_log2

# log2 data - 041GFP
x3$log2_dat<- log2(x3$Normalized_values)
x3_stat_log2 <- summarySE(x3, measurevar="log2_dat", groupvars=c("Metadata_Age...3","Metadata_Treatment"))
x3_stat_log2


