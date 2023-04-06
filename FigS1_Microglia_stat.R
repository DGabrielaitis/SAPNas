# 2022 05 26
# SAPNAS
#Figure S1
# Microglia incubated with 041-GFP and 053-SNAP-AF488

# Dovydas Gabrielaitis
library(Rmisc)
# ---------------------------- Microglia -----------  

# mean +- Graphs
# 041-GFP
#remove CytD data
#DATA_MG_041<-filter(DATA_MG_041, Metadata_Treatment != 'CHD' )
tgc3 <- summarySE(DATA_MG_041, measurevar="Intensity_IntegratedIntensity_NanoTube_Rescaled", groupvars=c("Metadata_Age...3","Metadata_Treatment"))
tgc3$Metadata_Age...3 <- factor(tgc3$Metadata_Age...3, levels=c("YNG", "OLD"))
pd <- position_dodge(1) # move them .05 to the left and right
ggplot(tgc3, aes(x=Metadata_Treatment, y=Intensity_IntegratedIntensity_NanoTube_Rescaled, colour=Metadata_Age...3, group=Metadata_Age...3)) + 
  geom_errorbar(aes(ymin=Intensity_IntegratedIntensity_NanoTube_Rescaled-se, ymax=Intensity_IntegratedIntensity_NanoTube_Rescaled+se), colour="black", width=1,size=1,  position=pd) +
  geom_point(position=pd, size=5, colour=1)+
  scale_x_discrete(limit = c("CNT", "LPS"),labels = c("Control","LPS"))+
  ylim(c(0,750))+  theme_bw()  +
  ylab('')+xlab('') +
  theme(axis.text = element_text(size=20)) # export as a 700x700 square


# 053-SNAP-AF488
#DATA_MG_fv3<-filter(DATA_MG_fv3, Metadata_Treatment != 'CHD' )
tgc4 <- summarySE(DATA_MG_fv3, measurevar="Intensity_IntegratedIntensity_NanoTube_Rescaled", groupvars=c("Metadata_Age...3","Metadata_Treatment"))
tgc4$Metadata_Age...3 <- factor(tgc4$Metadata_Age...3, levels=c("YNG", "OLD"))
pd <- position_dodge(1) # move them .05 to the left and right
ggplot(tgc4, aes(x=Metadata_Treatment, y=Intensity_IntegratedIntensity_NanoTube_Rescaled, colour=Metadata_Age...3, group=Metadata_Age...3)) +
  geom_errorbar(aes(ymin=Intensity_IntegratedIntensity_NanoTube_Rescaled-se, ymax=Intensity_IntegratedIntensity_NanoTube_Rescaled+se), colour="black", width=1,size=1, position=pd)+
  geom_point(position=pd, size=5, colour=1)+scale_x_discrete(limit = c("CNT", "LPS"),labels = c("Control","LPS"))+
  ylim(c(0,500))  + theme_bw()+ylab('')+xlab('') + 
  theme(axis.text = element_text(size=20)) 



