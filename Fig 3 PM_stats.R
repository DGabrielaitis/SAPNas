# 2022 05 26
# SAPNAS
# Figure 3 statistics
# Peritoneal macrophafes with 0411-GFP and 053-SNAP-AF488

library(Rmisc)
#----------------------------------- Peritoneal Macrophages -----------------------------------

#041
tgc1 <- summarySE(DATA_PM_041, measurevar="Intensity_IntegratedIntensity_Nanotube_Rescaled", groupvars=c("Metadata_Age","Metadata_Treatment"))
tgc1$Metadata_Age <- factor(tgc1$Metadata_Age, levels=c("YNG", "OLD"))
pd <- position_dodge(1) # move them .05 to the left and right
ggplot(tgc1, aes(x=Metadata_Treatment, y=Intensity_IntegratedIntensity_Nanotube_Rescaled, colour=Metadata_Age, group=Metadata_Age)) + 
  geom_errorbar(aes(ymin=Intensity_IntegratedIntensity_Nanotube_Rescaled-se, ymax=Intensity_IntegratedIntensity_Nanotube_Rescaled+se), colour="black", width=1,size=1,  position=pd) +
  geom_point(position=pd, size=5, colour=1)+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  ylim(c(0,300))+  theme_bw()  +
  ylab('')+xlab('') +
  theme(axis.text = element_text(size=20)) # export as a 700x700 square



#fv3
tgc2 <- summarySE(DATA_PM_fv3, measurevar="Intensity_IntegratedIntensity_Nanotube_Rescaled", groupvars=c("Metadata_Age","Metadata_Treatment"))
tgc2$Metadata_Age <- factor(tgc2$Metadata_Age, levels=c("YNG", "OLD"))
pd <- position_dodge(1) # move them .05 to the left and right
ggplot(tgc2, aes(x=Metadata_Treatment, y=Intensity_IntegratedIntensity_Nanotube_Rescaled, colour=Metadata_Age, group=Metadata_Age)) + 
  geom_errorbar(aes(ymin=Intensity_IntegratedIntensity_Nanotube_Rescaled-se, ymax=Intensity_IntegratedIntensity_Nanotube_Rescaled+se), colour="black", width=1,size=1,  position=pd) +
  geom_point(position=pd, size=5, colour=1)+
  scale_x_discrete(limit = c("CNT", "LPS", "CHD"),labels = c("Control","LPS","Cytochalasin D"))+
  ylim(c(0,300))+  theme_bw()  +
  ylab('')+xlab('') +
  theme(axis.text = element_text(size=20)) # export as a 700x700 square

# --------------------------------------------- STATISTICS ---------------------------------------



# 041-GFP and Peritoneal macrophages
# Filter data by condition - first by age, then by treatment
int1_old <- filter(DATA_PM_041, Metadata_Age == 'OLD' )
int1_old_cnt <- filter(int1_old,  Metadata_Treatment == 'CNT')
int1_old_cnt$Cathegory <- 'OLD_CNT'
int1_old_chd <- filter(int1_old,  Metadata_Treatment == 'CHD')
int1_old_chd$Cathegory <- 'OLD_CHD'
int1_old_lps <- filter(int1_old,  Metadata_Treatment == 'LPS')
int1_old_lps$Cathegory <- 'OLD_LPS'

int1_yng <- filter(DATA_PM_041, Metadata_Age == 'YNG' )
int1_yng_cnt <- filter(int1_yng,  Metadata_Treatment == 'CNT')
int1_yng_cnt$Cathegory <- 'YNG_CNT'
int1_yng_chd <- filter(int1_yng,  Metadata_Treatment == 'CHD')
int1_yng_chd$Cathegory <- 'YNG_CHD'
int1_yng_lps <- filter(int1_yng,  Metadata_Treatment == 'LPS')
int1_yng_lps$Cathegory <- 'YNG_LPS'

DATA_PM_041_stat <- rbind(int1_old_cnt,
                          int1_old_chd,
                          int1_old_lps,
                          int1_yng_cnt,
                          int1_yng_chd,
                          int1_yng_lps )

# two-way anova
summary(res.aov2 <- aov(Intensity_IntegratedIntensity_Nanotube_Rescaled ~ Metadata_Treatment + Metadata_Age, data = DATA_PM_041))
# Krushal Wallis
kruskal.test(Intensity_IntegratedIntensity_Nanotube_Rescaled ~ Cathegory, data = DATA_PM_041_stat)

pairwise.wilcox.test(DATA_PM_041_stat$Intensity_IntegratedIntensity_Nanotube_Rescaled, DATA_PM_041_stat$Cathegory,
                     p.adjust.method = "BH")


# 053-SNAP-AF488 and Peritoneal macrophages
# Filter data by condition - first by age, then by treatment
int5_old <- filter(DATA_PM_fv3, Metadata_Age == 'OLD' )
int5_old_cnt <- filter(int5_old,  Metadata_Treatment == 'CNT')
int5_old_cnt$Cathegory <- 'OLD_CNT'
int5_old_chd <- filter(int5_old,  Metadata_Treatment == 'CHD')
int5_old_chd$Cathegory <- 'OLD_CHD'
int5_old_lps <- filter(int5_old,  Metadata_Treatment == 'LPS')
int5_old_lps$Cathegory <- 'OLD_LPS'

int5_yng <- filter(DATA_PM_fv3, Metadata_Age == 'YNG' )
int5_yng_cnt <- filter(int5_yng,  Metadata_Treatment == 'CNT')
int5_yng_cnt$Cathegory <- 'YNG_CNT'
int5_yng_chd <- filter(int5_yng,  Metadata_Treatment == 'CHD')
int5_yng_chd$Cathegory <- 'YNG_CHD'
int5_yng_lps <- filter(int5_yng,  Metadata_Treatment == 'LPS')
int5_yng_lps$Cathegory <- 'YNG_LPS'

DATA_PM_fv3_stat <- rbind(int5_old_cnt,
                          int5_old_chd,
                          int5_old_lps,
                          int5_yng_cnt,
                          int5_yng_chd,
                          int5_yng_lps )

# two-way anova
summary(res.aov2 <- aov(Intensity_IntegratedIntensity_Nanotube_Rescaled ~ Metadata_Treatment + Metadata_Age, data = DATA_PM_fv3))
# Krushal Wallis
kruskal.test(Intensity_IntegratedIntensity_Nanotube_Rescaled ~ Cathegory, data = DATA_PM_fv3_stat)

pairwise.wilcox.test(DATA_PM_fv3_stat$Intensity_IntegratedIntensity_Nanotube_Rescaled, DATA_PM_fv3_stat$Cathegory,
                     p.adjust.method = "BH")


#--- testing normality
aaa<-DATA_PM_fv3[Metadata_Age=='OLD']
aaa<-aaa[Metadata_Treatment=='CHD']
aaa
shapiro.test(aaa$Intensity_IntegratedIntensity_Nanotube_Rescaled)


