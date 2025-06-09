########################## LOAD PACKAGES ##########################

library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyverse)
library(RMINC)
library(MRIcrotome)
library(reshape2)
library(glue)
library(grid)
library(gridExtra)

########################## LOAD DATA ##########################

# Use this statement below in all scripts to load everything in this environment
load("/projects/schoi/MTXcomparison/analysis/MTX17022024/med_MTX.RData")

df_dose <- med_df_MTX %>% filter(Dose == 5)
df_dose$Group <- paste(df_dose$Treatment, df_dose$Sex, sep = "_")


MTX24_IT0.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "M")
MTX24_IT0.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "F")
MTX24_IT1_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "M")
MTX24_IT1_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "F")
MTX24_IT2.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "M")
MTX24_IT2.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "F")
MTX24_IT5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "M")
MTX24_IT5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "F")

MTX42_IT0.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "M")
MTX42_IT0.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "F")
MTX42_IT1_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "M")
MTX42_IT1_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "F")
MTX42_IT2.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "M")
MTX42_IT2.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "F")
MTX42_IT5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "M")
MTX42_IT5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "F")

MTX63_IT0.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "M")
MTX63_IT0.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*(med_df_MTX$Sex == "F")
MTX63_IT1_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "M")
MTX63_IT1_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*(med_df_MTX$Sex == "F")
MTX63_IT2.5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "M")
MTX63_IT2.5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*(med_df_MTX$Sex == "F")
MTX63_IT5_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "M")
MTX63_IT5_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*(med_df_MTX$Sex == "F")


########################## SETTING UP IV COEFFICIENTS ##########################
################################################################################

MTX24_IV_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "M")
MTX24_IV_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "F")
MTX42_IV_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "M")
MTX42_IV_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "F")
MTX63_IV_M <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "M")
MTX63_IV_F <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IV")*(med_df_MTX$Sex == "F")


########################## GENERAL FLAGS ##########################

age_flag <- as.character(med_df_MTX$Age)
sex_flag <- med_df_MTX$Sex
dose_flag <- as.character(med_df_MTX$Dose)
study_id <- med_df_MTX$Study_ID
pregnancy_cage <- med_df_MTX$Pregnancy_cage
med_df_MTX$Group <- paste(med_df_MTX$Route, med_df_MTX$Dose, sep = "_")
group <- med_df_MTX$Group
route_flag <- med_df_MTX$Route



med_df_MTX_ranef <- med_df_MTX
for (j in 1:length(structure_list)){
  cstruct <- structure_list[j]
  print(paste(as.character(j), cstruct, sep=": "))
  first_form <- glue("`{cstruct}` ~ -1 + age_flag + age_flag:route_flag + age_flag:sex_flag +
                      MTX24_IT0.5_M + MTX24_IT0.5_F + 
                      MTX24_IT1_M + MTX24_IT1_F + 
                      MTX24_IT2.5_M + MTX24_IT2.5_F + 
                      MTX24_IT5_M + MTX24_IT5_F + 
                      MTX42_IT0.5_M + MTX42_IT0.5_F + 
                      MTX42_IT1_M + MTX42_IT1_F + 
                      MTX42_IT2.5_M + MTX42_IT2.5_F + 
                      MTX42_IT5_M + MTX42_IT5_F + 
                      MTX63_IT0.5_M + MTX63_IT0.5_F + 
                      MTX63_IT1_M + MTX63_IT1_F + 
                      MTX63_IT2.5_M + MTX63_IT2.5_F + 
                      MTX63_IT5_M + MTX63_IT5_F + 
                      MTX24_IV_M + MTX24_IV_F + 
                      MTX42_IV_M + MTX42_IV_F + 
                      MTX63_IV_M + MTX63_IV_F +
                      (1|study_id)")
  crlm <- lmer(as.formula(first_form), data = med_df_MTX_ranef, REML = FALSE)
  # Add the random effect of litter to the data frame
  med_df_MTX_ranef$ID_random_effect <- ranef(crlm)$study_id[as.character(med_df_MTX_ranef$Study_ID), ]
  # Subtract the random effect from the structures volumes
  med_df_MTX_ranef[cstruct] <- med_df_MTX_ranef[cstruct] - med_df_MTX_ranef$ID_random_effect
}



mean_ci <- function(indata){ mnse <- mean_se(indata); mnse$ymin <- mnse$y-1.96*(mnse$y-mnse$ymin); mnse$ymax <- mnse$y+1.96*(mnse$ymax-mnse$y); return(mnse) }
med_df_MTX_ranef$Treatment <- fct_rev(med_df_MTX_ranef$Treatment)
df_dose_ranef <- med_df_MTX_ranef %>% filter(Dose == 5)

########################## FIG 3A ##########################
# amygdala

plot_amygdala <- ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"amygdala"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Amygdala") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none", axis.title.x=element_blank())

########################## FIG 3B ##########################
# corpus_callosum

plot_CC <- ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"corpus_callosum"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Corpus Callosum") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())

########################## FIG 3C ##########################
# hypothalamus

plot_hypothalamus <- ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"hypothalamus"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Hypothalamus") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none", axis.title.x=element_blank())

########################## FIG 3D ##########################
# Pituitary_gland

plot_PG <- ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"Pituitary_gland"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Pituitary Gland") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())

########################## FIG 3E ##########################
# fourth_ventricle

plot_fourthventricle <- ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"fourth_ventricle"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Fourth Ventricle") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none")

########################## FIG 3F ##########################
# inferior_olivary_complex

plot_IOC <- ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"inferior_olivary_complex"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Inferior Olivary Complex") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none", axis.title.y=element_blank())


# RUN TO PRODUCE ALL PLOTS
plot_amygdala
plot_CC
plot_hypothalamus
plot_PG
plot_fourthventricle
plot_IOC


########################## PUT ALL PANELS TOGETHER ##########################

# png(filename="/projects/schoi/MTXcomparison/figures/fig3_v2.png", width = 12, height = 12, units = 'in', res = 600)
# grid.arrange(grobs = list(plot_amygdala, plot_CC,
#                           plot_hypothalamus, plot_PG,
#                           plot_fourthventricle, plot_IOC),
#              ncol=2,nrow=3)
# dev.off()

########################## TAKE LEGEND FROM HERE ##########################

# png(filename="/projects/schoi/MTXcomparison/figures/fig3_legend.png", width = 12, height = 12, units = 'in', res = 600)
# ggplot(df_dose_ranef,aes(x=as.character(Age),y=df_dose_ranef[,"amygdala"],colour=Treatment))+
#   geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
#   stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
#   stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
#   facet_wrap(~Route) +
#   scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
#   theme_bw()+
#   labs(x="Age", y="Volume (mm³)", title="Amygdala") +
#   theme(text=element_text(size=15,  family="Arial"))
# dev.off()



### Calculating percent change values

# Calculate mean for Saline and MTX groups, split by Route

# Calculate percent change for each structure
df_percent_change <- df_dose_ranef %>%
  group_by(Age, Route) %>%
  summarise(
    percent_change_amygdala = (mean(amygdala[Treatment == "MTX"], na.rm = TRUE) - mean(amygdala[Treatment == "Saline"], na.rm = TRUE)) / mean(amygdala[Treatment == "Saline"], na.rm = TRUE) * 100,
    percent_change_corpus_callosum = (mean(corpus_callosum[Treatment == "MTX"], na.rm = TRUE) - mean(corpus_callosum[Treatment == "Saline"], na.rm = TRUE)) / mean(corpus_callosum[Treatment == "Saline"], na.rm = TRUE) * 100,
    percent_change_hypothalamus = (mean(hypothalamus[Treatment == "MTX"], na.rm = TRUE) - mean(hypothalamus[Treatment == "Saline"], na.rm = TRUE)) / mean(hypothalamus[Treatment == "Saline"], na.rm = TRUE) * 100,
    percent_change_pituitary_gland = (mean(Pituitary_gland[Treatment == "MTX"], na.rm = TRUE) - mean(Pituitary_gland[Treatment == "Saline"], na.rm = TRUE)) / mean(Pituitary_gland[Treatment == "Saline"], na.rm = TRUE) * 100,
    percent_change_fourth_ventricle = (mean(fourth_ventricle[Treatment == "MTX"], na.rm = TRUE) - mean(fourth_ventricle[Treatment == "Saline"], na.rm = TRUE)) / mean(fourth_ventricle[Treatment == "Saline"], na.rm = TRUE) * 100,
    percent_change_inferior_olivary_complex = (mean(inferior_olivary_complex[Treatment == "MTX"], na.rm = TRUE) - mean(inferior_olivary_complex[Treatment == "Saline"], na.rm = TRUE)) / mean(inferior_olivary_complex[Treatment == "Saline"], na.rm = TRUE) * 100
  )

# View the result
print(df_percent_change)
