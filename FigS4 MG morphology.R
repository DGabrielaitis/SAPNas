# 2022 05 25

# SAPNAS
#Figure 3 
# Peritoneal macrophages incubated with 041-GFP and 053-SNAP-AF488
# Dovydas Gabrielaitis

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls()) 

#------------libraries and functinos --------
library(readr)
library(Rmisc)
library(dplyr)
library(doBy)
library(svMisc)
library(data.table)
library(sciplot)



# import data

DAT <- read_csv("data_out/SAPNAS_MG_FinalWholeCell.csv")
MICROGLIA<-DAT

#Filter only fv3
#MICROGLIA<-filter(DAT, Metadata_Protein != '041' )

#Filter only 041
#MICROGLIA<-filter(DAT, Metadata_Protein != 'fv3' )
# filter only CNT experiments
MICROGLIA<-filter(MICROGLIA, Metadata_Treatment == 'CNT' )




#MICROGLIA<-filter(MICROGLIA, Metadata_ExperimentNumber != '6' )
#MICROGLIA<-filter(MICROGLIA, Metadata_ExperimentNumber != '5' )
#MICROGLIA<-filter(MICROGLIA, Metadata_ExperimentNumber != '11' )
#MICROGLIA<-filter(MICROGLIA, Metadata_ExperimentNumber != '18' )


#filter out outlayer segments
MICROGLIA<-filter(MICROGLIA, AreaShape_Area > 450 )
MICROGLIA<-filter(MICROGLIA, AreaShape_Area < 10000 )



#---------------------------- Comparing the Area------


# average fv3 data and remove outlayers
int5<-setDT(MICROGLIA)[, list(AreaShape_Area= mean(AreaShape_Area)) , .(Metadata_Treatment,Metadata_ImageNumber, Metadata_WellNumber,Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3)]
int6<-setDT(int5)[, list(AreaShape_Area= mean(AreaShape_Area)) , .(Metadata_Treatment, Metadata_WellNumber,Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3,Metadata_Protein)]
int7<-setDT(int6)[, list(AreaShape_Area= mean(AreaShape_Area)) , .(Metadata_Treatment, Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3,Metadata_Protein)]
# ----Remove outlayers for fv3 data
# Split dataframe by treatment ant NT nocentration 
int5_YNG_fv3 <- int5 %>%filter(Metadata_Protein %in% c("fv3"), Metadata_Age...3 == 'YNG' )
int5_YNG_NEG<- int5 %>%filter(Metadata_Protein %in% c("NEG"), Metadata_Age...3 == 'YNG' )
int5_OLD_fv3 <- int5 %>%filter(Metadata_Protein %in% c("fv3"), Metadata_Age...3 == 'OLD' )
int5_OLD_NEG <- int5 %>%filter(Metadata_Protein %in% c("NEG"), Metadata_Age...3 == 'OLD' )
int5_YNG_041<- int5 %>%filter(Metadata_Protein %in% c("041"), Metadata_Age...3 == 'YNG' )
int5_OLD_041 <- int5 %>%filter(Metadata_Protein %in% c("041"), Metadata_Age...3 == 'OLD' )



# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int5_YNG_fv3$is_outlier <- ifelse(int5_YNG_fv3$AreaShape_Area %in% boxplot.stats(int5_YNG_fv3$AreaShape_Area)$out, 1, 0)
int5_YNG_NEG$is_outlier <- ifelse(int5_YNG_NEG$AreaShape_Area %in% boxplot.stats(int5_YNG_NEG$AreaShape_Area)$out, 1, 0)
int5_OLD_fv3$is_outlier <- ifelse(int5_OLD_fv3$AreaShape_Area %in% boxplot.stats(int5_OLD_fv3$AreaShape_Area)$out, 1, 0)
int5_OLD_NEG$is_outlier <- ifelse(int5_OLD_NEG$AreaShape_Area %in% boxplot.stats(int5_OLD_NEG$AreaShape_Area)$out, 1, 0)
int5_YNG_041$is_outlier <- ifelse(int5_YNG_041$AreaShape_Area %in% boxplot.stats(int5_YNG_041$AreaShape_Area)$out, 1, 0)
int5_OLD_041$is_outlier <- ifelse(int5_OLD_041$AreaShape_Area %in% boxplot.stats(int5_OLD_041$AreaShape_Area)$out, 1, 0)


# Bind the DF back
DATA_MICROGLIA <- rbind(int5_YNG_fv3,
                     int5_YNG_NEG,
                     int5_OLD_fv3,
                     int5_OLD_NEG,
                     int5_YNG_041,
                     int5_OLD_041
)

# Remove the values that have outlayers
PM_Area_2<- DATA_MICROGLIA %>%filter(is_outlier %in% c("0"))
PM_Area_2 <- summarySE(PM_Area_2, measurevar="AreaShape_Area", groupvars=c("Metadata_Age...3","Metadata_Protein",'Metadata_Treatment'))
PM_Area_2


#---------------------------- Comparing the Compactness------


# average fv3 data and remove outlayers
int5<-setDT(MICROGLIA)[, list(AreaShape_Compactness= mean(AreaShape_Compactness)) , .(Metadata_Treatment,Metadata_ImageNumber, Metadata_WellNumber,Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3)]
int6<-setDT(int5)[, list(AreaShape_Compactness= mean(AreaShape_Compactness)) , .(Metadata_Treatment, Metadata_WellNumber,Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3,Metadata_Protein)]
int7<-setDT(int6)[, list(AreaShape_Compactness= mean(AreaShape_Compactness)) , .(Metadata_Treatment, Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3,Metadata_Protein)]
# ----Remove outlayers for fv3 data
# Split dataframe by treatment ant NT nocentration 
int5_YNG_fv3 <- int5 %>%filter(Metadata_Protein %in% c("fv3"), Metadata_Age...3 == 'YNG' )
int5_YNG_NEG<- int5 %>%filter(Metadata_Protein %in% c("NEG"), Metadata_Age...3 == 'YNG' )
int5_OLD_fv3 <- int5 %>%filter(Metadata_Protein %in% c("fv3"), Metadata_Age...3 == 'OLD' )
int5_OLD_NEG <- int5 %>%filter(Metadata_Protein %in% c("NEG"), Metadata_Age...3 == 'OLD' )
int5_YNG_041<- int5 %>%filter(Metadata_Protein %in% c("041"), Metadata_Age...3 == 'YNG' )
int5_OLD_041 <- int5 %>%filter(Metadata_Protein %in% c("041"), Metadata_Age...3 == 'OLD' )



# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int5_YNG_fv3$is_outlier <- ifelse(int5_YNG_fv3$AreaShape_Compactness %in% boxplot.stats(int5_YNG_fv3$AreaShape_Compactness)$out, 1, 0)
int5_YNG_NEG$is_outlier <- ifelse(int5_YNG_NEG$AreaShape_Compactness %in% boxplot.stats(int5_YNG_NEG$AreaShape_Compactness)$out, 1, 0)
int5_OLD_fv3$is_outlier <- ifelse(int5_OLD_fv3$AreaShape_Compactness %in% boxplot.stats(int5_OLD_fv3$AreaShape_Compactness)$out, 1, 0)
int5_OLD_NEG$is_outlier <- ifelse(int5_OLD_NEG$AreaShape_Compactness %in% boxplot.stats(int5_OLD_NEG$AreaShape_Compactness)$out, 1, 0)
int5_YNG_041$is_outlier <- ifelse(int5_YNG_041$AreaShape_Compactness %in% boxplot.stats(int5_YNG_041$AreaShape_Compactness)$out, 1, 0)
int5_OLD_041$is_outlier <- ifelse(int5_OLD_041$AreaShape_Compactness %in% boxplot.stats(int5_OLD_041$AreaShape_Compactness)$out, 1, 0)


# Bind the DF back
DATA_MICROGLIA <- rbind(int5_YNG_fv3,
                        int5_YNG_NEG,
                        int5_OLD_fv3,
                        int5_OLD_NEG,
                        int5_YNG_041,
                        int5_OLD_041
)

# Remove the values that have outlayers
PM_Compactness_2<- DATA_MICROGLIA %>%filter(is_outlier %in% c("0"))
PM_Compactness_2 <- summarySE(PM_Compactness_2, measurevar="AreaShape_Compactness", groupvars=c("Metadata_Age...3","Metadata_Protein",'Metadata_Treatment'))
PM_Compactness_2


#---------------------------- Comparing the Eccentricity------


# average fv3 data and remove outlayers
int5<-setDT(MICROGLIA)[, list(AreaShape_Eccentricity= mean(AreaShape_Eccentricity)) , .(Metadata_Treatment,Metadata_ImageNumber, Metadata_WellNumber,Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3)]
int6<-setDT(int5)[, list(AreaShape_Eccentricity= mean(AreaShape_Eccentricity)) , .(Metadata_Treatment, Metadata_WellNumber,Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3,Metadata_Protein)]
int7<-setDT(int6)[, list(AreaShape_Eccentricity= mean(AreaShape_Eccentricity)) , .(Metadata_Treatment, Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3,Metadata_Protein)]
# ----Remove outlayers for fv3 data
# Split dataframe by treatment ant NT nocentration 
int5_YNG_fv3 <- int5 %>%filter(Metadata_Protein %in% c("fv3"), Metadata_Age...3 == 'YNG' )
int5_YNG_NEG<- int5 %>%filter(Metadata_Protein %in% c("NEG"), Metadata_Age...3 == 'YNG' )
int5_OLD_fv3 <- int5 %>%filter(Metadata_Protein %in% c("fv3"), Metadata_Age...3 == 'OLD' )
int5_OLD_NEG <- int5 %>%filter(Metadata_Protein %in% c("NEG"), Metadata_Age...3 == 'OLD' )
int5_YNG_041<- int5 %>%filter(Metadata_Protein %in% c("041"), Metadata_Age...3 == 'YNG' )
int5_OLD_041 <- int5 %>%filter(Metadata_Protein %in% c("041"), Metadata_Age...3 == 'OLD' )



# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int5_YNG_fv3$is_outlier <- ifelse(int5_YNG_fv3$AreaShape_Eccentricity %in% boxplot.stats(int5_YNG_fv3$AreaShape_Eccentricity)$out, 1, 0)
int5_YNG_NEG$is_outlier <- ifelse(int5_YNG_NEG$AreaShape_Eccentricity %in% boxplot.stats(int5_YNG_NEG$AreaShape_Eccentricity)$out, 1, 0)
int5_OLD_fv3$is_outlier <- ifelse(int5_OLD_fv3$AreaShape_Eccentricity %in% boxplot.stats(int5_OLD_fv3$AreaShape_Eccentricity)$out, 1, 0)
int5_OLD_NEG$is_outlier <- ifelse(int5_OLD_NEG$AreaShape_Eccentricity %in% boxplot.stats(int5_OLD_NEG$AreaShape_Eccentricity)$out, 1, 0)
int5_YNG_041$is_outlier <- ifelse(int5_YNG_041$AreaShape_Eccentricity %in% boxplot.stats(int5_YNG_041$AreaShape_Eccentricity)$out, 1, 0)
int5_OLD_041$is_outlier <- ifelse(int5_OLD_041$AreaShape_Eccentricity %in% boxplot.stats(int5_OLD_041$AreaShape_Eccentricity)$out, 1, 0)


# Bind the DF back
DATA_MICROGLIA <- rbind(int5_YNG_fv3,
                        int5_YNG_NEG,
                        int5_OLD_fv3,
                        int5_OLD_NEG,
                        int5_YNG_041,
                        int5_OLD_041
)

# Remove the values that have outlayers
PM_Eccentricity_2<- DATA_MICROGLIA %>%filter(is_outlier %in% c("0"))
PM_Eccentricity_2 <- summarySE(PM_Eccentricity_2, measurevar="AreaShape_Eccentricity", groupvars=c("Metadata_Age...3","Metadata_Protein",'Metadata_Treatment'))
PM_Eccentricity_2


#---------------------------- Comparing the Solidity------


# average fv3 data and remove outlayers
int5<-setDT(MICROGLIA)[, list(AreaShape_Solidity= mean(AreaShape_Solidity)) , .(Metadata_Treatment,Metadata_ImageNumber, Metadata_WellNumber,Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3)]
int6<-setDT(int5)[, list(AreaShape_Solidity= mean(AreaShape_Solidity)) , .(Metadata_Treatment, Metadata_WellNumber,Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3,Metadata_Protein)]
int7<-setDT(int6)[, list(AreaShape_Solidity= mean(AreaShape_Solidity)) , .(Metadata_Treatment, Metadata_Protein,Metadata_ExperimentNumber,Metadata_Age...3,Metadata_Protein)]
# ----Remove outlayers for fv3 data
# Split dataframe by treatment ant NT nocentration 
int5_YNG_fv3 <- int5 %>%filter(Metadata_Protein %in% c("fv3"), Metadata_Age...3 == 'YNG' )
int5_YNG_NEG<- int5 %>%filter(Metadata_Protein %in% c("NEG"), Metadata_Age...3 == 'YNG' )
int5_OLD_fv3 <- int5 %>%filter(Metadata_Protein %in% c("fv3"), Metadata_Age...3 == 'OLD' )
int5_OLD_NEG <- int5 %>%filter(Metadata_Protein %in% c("NEG"), Metadata_Age...3 == 'OLD' )
int5_YNG_041<- int5 %>%filter(Metadata_Protein %in% c("041"), Metadata_Age...3 == 'YNG' )
int5_OLD_041 <- int5 %>%filter(Metadata_Protein %in% c("041"), Metadata_Age...3 == 'OLD' )



# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int5_YNG_fv3$is_outlier <- ifelse(int5_YNG_fv3$AreaShape_Solidity %in% boxplot.stats(int5_YNG_fv3$AreaShape_Solidity)$out, 1, 0)
int5_YNG_NEG$is_outlier <- ifelse(int5_YNG_NEG$AreaShape_Solidity %in% boxplot.stats(int5_YNG_NEG$AreaShape_Solidity)$out, 1, 0)
int5_OLD_fv3$is_outlier <- ifelse(int5_OLD_fv3$AreaShape_Solidity %in% boxplot.stats(int5_OLD_fv3$AreaShape_Solidity)$out, 1, 0)
int5_OLD_NEG$is_outlier <- ifelse(int5_OLD_NEG$AreaShape_Solidity %in% boxplot.stats(int5_OLD_NEG$AreaShape_Solidity)$out, 1, 0)
int5_YNG_041$is_outlier <- ifelse(int5_YNG_041$AreaShape_Solidity %in% boxplot.stats(int5_YNG_041$AreaShape_Solidity)$out, 1, 0)
int5_OLD_041$is_outlier <- ifelse(int5_OLD_041$AreaShape_Solidity %in% boxplot.stats(int5_OLD_041$AreaShape_Solidity)$out, 1, 0)


# Bind the DF back
DATA_MICROGLIA <- rbind(int5_YNG_fv3,
                        int5_YNG_NEG,
                        int5_OLD_fv3,
                        int5_OLD_NEG,
                        int5_YNG_041,
                        int5_OLD_041
)

# Remove the values that have outlayers
PM_Solidity_2<- DATA_MICROGLIA %>%filter(is_outlier %in% c("0"))
PM_Solidity_2 <- summarySE(PM_Solidity_2, measurevar="AreaShape_Solidity", groupvars=c("Metadata_Age...3","Metadata_Protein",'Metadata_Treatment'))
PM_Solidity_2




#--------- data------------

PM_Solidity_2
PM_Eccentricity_2
PM_Compactness_2