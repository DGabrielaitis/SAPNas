setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls()) 

# 2022 05 12
# SAPNAS
# 041-GFP ir 053SNAP particle uptake into BV2 cells

#---------------------- libraries and extras -----------------------------
#library(Rcmdr)
library(Rmisc)
library(readr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(doBy)
library(svMisc)
library(data.table)
library(sciplot)




remove_outlayer<-function(x){
  x1<-x[!x %in% boxplot.stats(x)$out]
  return(x1)
}

# Analyse the Total cell data:

# ---- 041-GFP-Total-------------------
# measure the integrated intensity after math - which is the corrected intensity to the background
Total_cell_041 <- read_csv("data_output/final/BV2_041_GFP_FINAL_Total_cell.csv")


Total_cell_041[Total_cell_041 == "041GFP_wi" | Total_cell_041 == "FV3SiR_wi"] <- "100 ng/mL"
Total_cell_041[Total_cell_041 == "041GFP_wo" | Total_cell_041 == "FV3SiR_wo"] <- "No treatment"




# calculate means per image
int1<-setDT(Total_cell_041)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,ImageNumber,Metadata_Concentration...6,Metadata_LPS)]
write.csv(int1,"/Volumes/Dovydas HD2/SAPNAS-Paper/ppt and analysis files/2_BV2_Internalisation/data_output/After averaging for every image/total041.csv", row.names = FALSE)

# calculate means per well
int2<-setDT(int1)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
# calculate means per experiment
int3<-setDT(int2)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
# Fix the treatment order
int1$Metadata_LPS <- factor(int1$Metadata_LPS, levels=c("No treatment", "100 ng/mL"))

# Remove outlayers

# Split dataframe by treatment ant NT nocentration 
int1_NOLPS_25ugml <- int1 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '25ugmL' )
int1_NOLPS_50ugml<- int1 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '50ugmL' )
int1_NOLPS_100ugml<- int1 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '100ugmL' )
int1_LPS_25ugml <- int1 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '25ugmL' )
int1_LPS_50ugml<- int1 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '50ugmL' )
int1_LPS_100ugml<- int1 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '100ugmL' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int1_NOLPS_25ugml$is_outlier <- ifelse(int1_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_NOLPS_50ugml$is_outlier <- ifelse(int1_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_NOLPS_100ugml$is_outlier <- ifelse(int1_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_LPS_25ugml$is_outlier <- ifelse(int1_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_LPS_50ugml$is_outlier <- ifelse(int1_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_LPS_100ugml$is_outlier <- ifelse(int1_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
# Bind the DF back
DATA_041 <- rbind(int1_NOLPS_25ugml,
                  int1_NOLPS_50ugml,
                  int1_NOLPS_100ugml,
                  int1_LPS_25ugml,
                  int1_LPS_50ugml,
                  int1_LPS_100ugml
                  )
# Remove the values that have outlayers
DATA_041<- DATA_041 %>%filter(is_outlier %in% c("0"))

gg_041_Total <- ggplot(DATA_041, aes(x=Metadata_Concentration...6,y=log(Intensity_IntegratedIntensity_GFP_aftermath), color= Metadata_LPS))+
  geom_boxplot()+  theme_bw() +geom_point(pch = 4, position=position_jitterdodge()) +
  scale_x_discrete(limit = c( "25ugmL", "50ugmL", '100ugmL'),labels = c('25',"50", "100"))+labs(color='LPS Treatment')+
  labs(title = "BV2 cells incubated with 041-GFP - CYTOPLASM", x = "Nanotube Concentration, µg/mL", y = "log Mean integrated intensity, a.u.")+
  theme(text = element_text(size=18))+ ylim(c(-2,5))


DATA_041$Intensity_IntegratedIntensity_GFP_aftermath <- DATA_041$Intensity_IntegratedIntensity_GFP_aftermath*12.657
DATA_041$LOG_Intensity_IntegratedIntensity_GFP_aftermath <- log(DATA_041$Intensity_IntegratedIntensity_GFP_aftermath)
DATA_041$LOG_Intensity_IntegratedIntensity_GFP_aftermath  <- DATA_041$LOG_Intensity_IntegratedIntensity_GFP_aftermath *23.91032
DATA_041<-filter(DATA_041, LOG_Intensity_IntegratedIntensity_GFP_aftermath > 0.00001 )
stat_041_Total <- summarySE(DATA_041, measurevar="Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))
stat_041_Total_LOG <- summarySE(DATA_041, measurevar="LOG_Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))



# ---- 041-GFP-Cytoplasm-------------------
# measure the integrated intensity after math - which is the corrected intensity to the background
Cytoplasm_041 <- read_csv("data_output/final/BV2_041_GFP_FINAL_Cytoplasm.csv")


Cytoplasm_041[Cytoplasm_041 == "041GFP_wi" | Cytoplasm_041 == "FV3SiR_wi"] <- "100 ng/mL"
Cytoplasm_041[Cytoplasm_041 == "041GFP_wo" | Cytoplasm_041 == "FV3SiR_wo"] <- "No treatment"




# calculate means per image
int1<-setDT(Cytoplasm_041)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,ImageNumber,Metadata_Concentration...6,Metadata_LPS)]
write.csv(int1,"/Volumes/Dovydas HD2/SAPNAS-Paper/ppt and analysis files/2_BV2_Internalisation/data_output/After averaging for every image/cytoplasm041.csv", row.names = FALSE)

# calculate means per well
int2<-setDT(int1)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
# calculate means per experiment
int3<-setDT(int2)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
# Fix the treatment order
int1$Metadata_LPS <- factor(int1$Metadata_LPS, levels=c("No treatment", "100 ng/mL"))

# Remove outlayers

# Split dataframe by treatment ant NT nocentration 
int1_NOLPS_25ugml <- int1 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '25ugmL' )
int1_NOLPS_50ugml<- int1 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '50ugmL' )
int1_NOLPS_100ugml<- int1 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '100ugmL' )
int1_LPS_25ugml <- int1 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '25ugmL' )
int1_LPS_50ugml<- int1 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '50ugmL' )
int1_LPS_100ugml<- int1 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '100ugmL' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int1_NOLPS_25ugml$is_outlier <- ifelse(int1_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_NOLPS_50ugml$is_outlier <- ifelse(int1_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_NOLPS_100ugml$is_outlier <- ifelse(int1_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_LPS_25ugml$is_outlier <- ifelse(int1_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_LPS_50ugml$is_outlier <- ifelse(int1_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_LPS_100ugml$is_outlier <- ifelse(int1_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
# Bind the DF back
DATA_041 <- rbind(int1_NOLPS_25ugml,
                  int1_NOLPS_50ugml,
                  int1_NOLPS_100ugml,
                  int1_LPS_25ugml,
                  int1_LPS_50ugml,
                  int1_LPS_100ugml
)
# Remove the values that have outlayers
DATA_041<- DATA_041 %>%filter(is_outlier %in% c("0"))

gg_041_Cytoplasm <- ggplot(DATA_041, aes(x=Metadata_Concentration...6,y=log(Intensity_IntegratedIntensity_GFP_aftermath), color= Metadata_LPS))+
  geom_boxplot()+  theme_bw() +geom_point(pch = 4, position=position_jitterdodge()) +
  scale_x_discrete(limit = c( "25ugmL", "50ugmL", '100ugmL'),labels = c('25',"50", "100"))+labs(color='LPS Treatment')+
  labs(title = "BV2 cells incubated with 041-GFP - CYTOPLASM", x = "Nanotube Concentration, µg/mL", y = "log Mean integrated intensity, a.u.")+
  theme(text = element_text(size=18))+ ylim(c(-2,5))

DATA_041$Intensity_IntegratedIntensity_GFP_aftermath <- DATA_041$Intensity_IntegratedIntensity_GFP_aftermath*12.657
DATA_041$LOG_Intensity_IntegratedIntensity_GFP_aftermath <- log(DATA_041$Intensity_IntegratedIntensity_GFP_aftermath)
DATA_041$LOG_Intensity_IntegratedIntensity_GFP_aftermath  <- DATA_041$LOG_Intensity_IntegratedIntensity_GFP_aftermath *23.91032
DATA_041<-filter(DATA_041, LOG_Intensity_IntegratedIntensity_GFP_aftermath > 0.00001 )
stat_041_Cytoplasm <- summarySE(DATA_041, measurevar="Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))
stat_041_Cytoplasm_LOG <- summarySE(DATA_041, measurevar="LOG_Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))





# ---- 041-GFP-Membrane-------------------
# measure the integrated intensity after math - which is the corrected intensity to the background
Membrane_041 <- read_csv("data_output/final/BV2_041_GFP_FINAL_Membrane.csv")
Membrane_041[Membrane_041 == "041GFP_wi" | Membrane_041 == "FV3SiR_wi"] <- "100 ng/mL"
Membrane_041[Membrane_041 == "041GFP_wo" | Membrane_041 == "FV3SiR_wo"] <- "No treatment"




# calculate means per image
int1<-setDT(Membrane_041)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,ImageNumber,Metadata_Concentration...6,Metadata_LPS)]
write.csv(int1,"/Volumes/Dovydas HD2/SAPNAS-Paper/ppt and analysis files/2_BV2_Internalisation/data_output/After averaging for every image/membrane041.csv", row.names = FALSE)

# calculate means per well
int2<-setDT(int1)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
# calculate means per experiment
int3<-setDT(int2)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
# Fix the treatment order
int1$Metadata_LPS <- factor(int1$Metadata_LPS, levels=c("No treatment", "100 ng/mL"))

# Remove outlayers

# Split dataframe by treatment ant NT nocentration 
int1_NOLPS_25ugml <- int1 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '25ugmL' )
int1_NOLPS_50ugml<- int1 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '50ugmL' )
int1_NOLPS_100ugml<- int1 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '100ugmL' )
int1_LPS_25ugml <- int1 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '25ugmL' )
int1_LPS_50ugml<- int1 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '50ugmL' )
int1_LPS_100ugml<- int1 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '100ugmL' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int1_NOLPS_25ugml$is_outlier <- ifelse(int1_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_NOLPS_50ugml$is_outlier <- ifelse(int1_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_NOLPS_100ugml$is_outlier <- ifelse(int1_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_LPS_25ugml$is_outlier <- ifelse(int1_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_LPS_50ugml$is_outlier <- ifelse(int1_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int1_LPS_100ugml$is_outlier <- ifelse(int1_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int1_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
# Bind the DF back
DATA_041 <- rbind(int1_NOLPS_25ugml,
                  int1_NOLPS_50ugml,
                  int1_NOLPS_100ugml,
                  int1_LPS_25ugml,
                  int1_LPS_50ugml,
                  int1_LPS_100ugml
)
# Remove the values that have outlayers
DATA_041<- DATA_041 %>%filter(is_outlier %in% c("0"))

gg_041_Membrane <- ggplot(DATA_041, aes(x=Metadata_Concentration...6,y=log(Intensity_IntegratedIntensity_GFP_aftermath), color= Metadata_LPS))+
  geom_boxplot()+  theme_bw() +geom_point(pch = 4, position=position_jitterdodge()) +
  scale_x_discrete(limit = c( "25ugmL", "50ugmL", '100ugmL'),labels = c('25',"50", "100"))+labs(color='LPS Treatment')+
  labs(title = "BV2 cells incubated with 041-GFP - CYTOPLASM", x = "Nanotube Concentration, µg/mL", y = "log Mean integrated intensity, a.u.")+
  theme(text = element_text(size=18))+ ylim(c(-2,5))

DATA_041$Intensity_IntegratedIntensity_GFP_aftermath <- DATA_041$Intensity_IntegratedIntensity_GFP_aftermath*12.657
DATA_041$LOG_Intensity_IntegratedIntensity_GFP_aftermath <- log(DATA_041$Intensity_IntegratedIntensity_GFP_aftermath)
DATA_041$LOG_Intensity_IntegratedIntensity_GFP_aftermath  <- DATA_041$LOG_Intensity_IntegratedIntensity_GFP_aftermath *23.91032
DATA_041<-filter(DATA_041, LOG_Intensity_IntegratedIntensity_GFP_aftermath > 0.00001 )
stat_041_Membrane <- summarySE(DATA_041, measurevar="Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))
stat_041_Membrane_LOG <- summarySE(DATA_041, measurevar="LOG_Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))




#------fv3-AF488-Total--------------

Total_cell_053 <- read_csv("data_output/final/BV2_FV3_AF488_FINAL_Total_cell.csv")
Total_cell_053[Total_cell_053 == "041GFP_wi" | Total_cell_053 == "FV3SiR_wi"] <- "100 ng/mL"
Total_cell_053[Total_cell_053 == "041GFP_wo" | Total_cell_053 == "FV3SiR_wo"] <- "No treatment"

int4<-setDT(Total_cell_053)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,ImageNumber,Metadata_Concentration...6,Metadata_LPS)]
write.csv(int4,"/Volumes/Dovydas HD2/SAPNAS-Paper/ppt and analysis files/2_BV2_Internalisation/data_output/After averaging for every image/total053.csv", row.names = FALSE)

int5<-setDT(int4)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
int6<-setDT(int5)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
int6_new<-filter(int6, Metadata_Concentration...6 != 'control' )
int4$Metadata_LPS <- factor(int4$Metadata_LPS, levels=c("No treatment", "100 ng/mL"))

# Remove outlayers

# Split dataframe by treatment ant NT nocentration 
int4_NOLPS_25ugml <- int4 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '25ugmL' )
int4_NOLPS_50ugml<- int4 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '50ugmL' )
int4_NOLPS_100ugml<- int4 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '100ugmL' )
int4_LPS_25ugml <- int4 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '25ugmL' )
int4_LPS_50ugml<- int4 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '50ugmL' )
int4_LPS_100ugml<- int4 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '100ugmL' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int4_NOLPS_25ugml$is_outlier <- ifelse(int4_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_NOLPS_50ugml$is_outlier <- ifelse(int4_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_NOLPS_100ugml$is_outlier <- ifelse(int4_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_LPS_25ugml$is_outlier <- ifelse(int4_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_LPS_50ugml$is_outlier <- ifelse(int4_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_LPS_100ugml$is_outlier <- ifelse(int4_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)

DATA_FV3 <- rbind(int4_NOLPS_25ugml,
                  int4_NOLPS_50ugml,
                  int4_NOLPS_100ugml,
                  int4_LPS_25ugml,
                  int4_LPS_50ugml,
                  int4_LPS_100ugml)
DATA_FV3<- DATA_FV3 %>%filter(is_outlier %in% c("0"))



gg_053_Total <- ggplot(DATA_FV3, aes(x=Metadata_Concentration...6,y=log(Intensity_IntegratedIntensity_GFP_aftermath), color= Metadata_LPS))+
  geom_boxplot()+geom_point(pch = 4, position=position_jitterdodge()) + 
  scale_x_discrete(limit = c( "25ugmL", "50ugmL", '100ugmL'),labels = c("25","50", "100"))+
  theme_bw() + labs(color='LPS Treatment')+
  labs(title = "BV2 cells incubated with 053SNAP - Total Cell", x = "Nanotube Concentration, µg/mL", y = "log Mean integrated intensity, a.u.")+
  theme(text = element_text(size=18))+ ylim(c(0,5))

DATA_FV3$Intensity_IntegratedIntensity_GFP_aftermath <- DATA_FV3$Intensity_IntegratedIntensity_GFP_aftermath*9.869233
DATA_FV3$LOG_Intensity_IntegratedIntensity_GFP_aftermath <- log(DATA_FV3$Intensity_IntegratedIntensity_GFP_aftermath)
DATA_FV3$LOG_Intensity_IntegratedIntensity_GFP_aftermath  <- DATA_FV3$LOG_Intensity_IntegratedIntensity_GFP_aftermath *23.02941
DATA_FV3<-filter(DATA_FV3, LOG_Intensity_IntegratedIntensity_GFP_aftermath > 0.00001 )
stat_FV3_Total <- summarySE(DATA_FV3, measurevar="Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))
stat_FV3_Total_LOG <- summarySE(DATA_FV3, measurevar="LOG_Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))


#------fv3-AF488-Cytoplasm--------------

Cytoplasm_053 <- read_csv("data_output/final/BV2_FV3_AF488_FINAL_Cytoplasm.csv")
Cytoplasm_053[Cytoplasm_053 == "041GFP_wi" | Cytoplasm_053 == "FV3SiR_wi"] <- "100 ng/mL"
Cytoplasm_053[Cytoplasm_053 == "041GFP_wo" | Cytoplasm_053 == "FV3SiR_wo"] <- "No treatment"

int4<-setDT(Cytoplasm_053)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,ImageNumber,Metadata_Concentration...6,Metadata_LPS)]
write.csv(int4,"/Volumes/Dovydas HD2/SAPNAS-Paper/ppt and analysis files/2_BV2_Internalisation/data_output/After averaging for every image/cytoplasm053.csv", row.names = FALSE)

int5<-setDT(int4)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
int6<-setDT(int5)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
int6_new<-filter(int6, Metadata_Concentration...6 != 'control' )
int4$Metadata_LPS <- factor(int4$Metadata_LPS, levels=c("No treatment", "100 ng/mL"))

# Remove outlayers

# Split dataframe by treatment ant NT nocentration 
int4_NOLPS_25ugml <- int4 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '25ugmL' )
int4_NOLPS_50ugml<- int4 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '50ugmL' )
int4_NOLPS_100ugml<- int4 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '100ugmL' )
int4_LPS_25ugml <- int4 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '25ugmL' )
int4_LPS_50ugml<- int4 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '50ugmL' )
int4_LPS_100ugml<- int4 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '100ugmL' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int4_NOLPS_25ugml$is_outlier <- ifelse(int4_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_NOLPS_50ugml$is_outlier <- ifelse(int4_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_NOLPS_100ugml$is_outlier <- ifelse(int4_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_LPS_25ugml$is_outlier <- ifelse(int4_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_LPS_50ugml$is_outlier <- ifelse(int4_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_LPS_100ugml$is_outlier <- ifelse(int4_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)

DATA_FV3 <- rbind(int4_NOLPS_25ugml,
                  int4_NOLPS_50ugml,
                  int4_NOLPS_100ugml,
                  int4_LPS_25ugml,
                  int4_LPS_50ugml,
                  int4_LPS_100ugml)
DATA_FV3<- DATA_FV3 %>%filter(is_outlier %in% c("0"))



gg_053_Cytoplasm <- ggplot(DATA_FV3, aes(x=Metadata_Concentration...6,y=log(Intensity_IntegratedIntensity_GFP_aftermath), color= Metadata_LPS))+
  geom_boxplot()+geom_point(pch = 4, position=position_jitterdodge()) + 
  scale_x_discrete(limit = c( "25ugmL", "50ugmL", '100ugmL'),labels = c("25","50", "100"))+
  theme_bw() + labs(color='LPS Treatment')+
  labs(title = "BV2 cells incubated with053SNAP - Cytoplasm", x = "Nanotube Concentration, µg/mL", y = "log Mean integrated intensity, a.u.")+
  theme(text = element_text(size=18))+ ylim(c(0,5))

DATA_FV3$Intensity_IntegratedIntensity_GFP_aftermath <- DATA_FV3$Intensity_IntegratedIntensity_GFP_aftermath*9.869233
DATA_FV3$LOG_Intensity_IntegratedIntensity_GFP_aftermath <- log(DATA_FV3$Intensity_IntegratedIntensity_GFP_aftermath)
DATA_FV3$LOG_Intensity_IntegratedIntensity_GFP_aftermath  <- DATA_FV3$LOG_Intensity_IntegratedIntensity_GFP_aftermath *23.02941
DATA_FV3<-filter(DATA_FV3, LOG_Intensity_IntegratedIntensity_GFP_aftermath > 0.00001 )
stat_FV3_Cytoplasm <- summarySE(DATA_FV3, measurevar="Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))
stat_FV3_Cytoplasm_LOG <- summarySE(DATA_FV3, measurevar="LOG_Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))



#------fv3-AF488-Membrane--------------

Membrane_053 <- read_csv("data_output/final/BV2_FV3_AF488_FINAL_Membrane.csv")
Membrane_053[Membrane_053 == "041GFP_wi" | Membrane_053 == "FV3SiR_wi"] <- "100 ng/mL"
Membrane_053[Membrane_053 == "041GFP_wo" | Membrane_053 == "FV3SiR_wo"] <- "No treatment"

int4<-setDT(Membrane_053)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,ImageNumber,Metadata_Concentration...6,Metadata_LPS)]
write.csv(int4,"/Volumes/Dovydas HD2/SAPNAS-Paper/ppt and analysis files/2_BV2_Internalisation/data_output/After averaging for every image/membrane053.csv", row.names = FALSE)

int5<-setDT(int4)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Well,Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
int6<-setDT(int5)[, list(Intensity_IntegratedIntensity_GFP_aftermath= mean(Intensity_IntegratedIntensity_GFP_aftermath)) , .(Metadata_Experiment_Number,Metadata_Concentration...6,Metadata_LPS)]
int6_new<-filter(int6, Metadata_Concentration...6 != 'control' )
int4$Metadata_LPS <- factor(int4$Metadata_LPS, levels=c("No treatment", "100 ng/mL"))

# Remove outlayers

# Split dataframe by treatment ant NT nocentration 
int4_NOLPS_25ugml <- int4 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '25ugmL' )
int4_NOLPS_50ugml<- int4 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '50ugmL' )
int4_NOLPS_100ugml<- int4 %>%filter(Metadata_LPS %in% c("No treatment"), Metadata_Concentration...6 == '100ugmL' )
int4_LPS_25ugml <- int4 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '25ugmL' )
int4_LPS_50ugml<- int4 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '50ugmL' )
int4_LPS_100ugml<- int4 %>%filter(Metadata_LPS %in% c("100 ng/mL"), Metadata_Concentration...6 == '100ugmL' )


# Remove outlayers
# Here we make a new column, where we indicate if the value is an outlayer
int4_NOLPS_25ugml$is_outlier <- ifelse(int4_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_NOLPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_NOLPS_50ugml$is_outlier <- ifelse(int4_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_NOLPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_NOLPS_100ugml$is_outlier <- ifelse(int4_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_NOLPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_LPS_25ugml$is_outlier <- ifelse(int4_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_LPS_25ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_LPS_50ugml$is_outlier <- ifelse(int4_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_LPS_50ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)
int4_LPS_100ugml$is_outlier <- ifelse(int4_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath %in% boxplot.stats(int4_LPS_100ugml$Intensity_IntegratedIntensity_GFP_aftermath)$out, 1, 0)

DATA_FV3 <- rbind(int4_NOLPS_25ugml,
                  int4_NOLPS_50ugml,
                  int4_NOLPS_100ugml,
                  int4_LPS_25ugml,
                  int4_LPS_50ugml,
                  int4_LPS_100ugml)
DATA_FV3<- DATA_FV3 %>%filter(is_outlier %in% c("0"))



gg_053_Membrane <- ggplot(DATA_FV3, aes(x=Metadata_Concentration...6,y=log(Intensity_IntegratedIntensity_GFP_aftermath), color= Metadata_LPS))+
  geom_boxplot()+geom_point(pch = 4, position=position_jitterdodge()) + 
  scale_x_discrete(limit = c( "25ugmL", "50ugmL", '100ugmL'),labels = c("25","50", "100"))+
  theme_bw() + labs(color='LPS Treatment')+
  labs(title = "BV2 cells incubated with 053SNAP - Membrane", x = "Nanotube Concentration, µg/mL", y = "log Mean integrated intensity, a.u.")+
  theme(text = element_text(size=18))+ ylim(c(0,5))

DATA_FV3$Intensity_IntegratedIntensity_GFP_aftermath <- DATA_FV3$Intensity_IntegratedIntensity_GFP_aftermath*9.869233
DATA_FV3$LOG_Intensity_IntegratedIntensity_GFP_aftermath <- log(DATA_FV3$Intensity_IntegratedIntensity_GFP_aftermath)
DATA_FV3$LOG_Intensity_IntegratedIntensity_GFP_aftermath  <- DATA_FV3$LOG_Intensity_IntegratedIntensity_GFP_aftermath *23.02941
DATA_FV3<-filter(DATA_FV3, LOG_Intensity_IntegratedIntensity_GFP_aftermath > 0.00001 )
stat_FV3_Membrane <- summarySE(DATA_FV3, measurevar="Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))
stat_FV3_Membrane_LOG <- summarySE(DATA_FV3, measurevar="LOG_Intensity_IntegratedIntensity_GFP_aftermath", groupvars=c("Metadata_Concentration...6","Metadata_LPS"))



#-------------------- Plot Graphs--------------
'
gg_041_Total
gg_041_Cytoplasm
gg_041_Membrane

gg_053_Total
gg_053_Cytoplasm
gg_053_Membrane
'


#----Print stat data-------

#
stat_041_Total
stat_041_Total_LOG
#
stat_041_Cytoplasm
stat_041_Cytoplasm_LOG
#
stat_041_Membrane_LOG
stat_041_Membrane
#
stat_FV3_Total
stat_FV3_Total_LOG
#
#
stat_FV3_Cytoplasm
stat_FV3_Cytoplasm_LOG
#
stat_FV3_Membrane_LOG
stat_FV3_Membrane







