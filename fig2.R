# tissueTimecourse.R

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


########################## GENERAL FLAGS ##########################

age_flag <- as.character(med_df_MTX$Age)
sex_flag <- med_df_MTX$Sex
dose_flag <- as.character(med_df_MTX$Dose)
study_id <- med_df_MTX$Study_ID
pregnancy_cage <- med_df_MTX$Pregnancy_cage
med_df_MTX$Group <- paste(med_df_MTX$Route, med_df_MTX$Dose, sep = "_")
group <- med_df_MTX$Group
route_flag <- med_df_MTX$Route

########################## SETTING UP IT COEFFICIENTS ##########################
################################################################################

MTX24_IT0.5 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)
MTX24_IT0.5_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX24_IT1 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)
MTX24_IT1_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX24_IT2.5 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)
MTX24_IT2.5_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX24_IT5 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)
MTX24_IT5_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)

MTX42_IT0.5 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)
MTX42_IT0.5_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX42_IT1 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)
MTX42_IT1_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX42_IT2.5 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)
MTX42_IT2.5_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX42_IT5 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)
MTX42_IT5_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)

MTX63_IT0.5 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)
MTX63_IT0.5_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 0.5)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX63_IT1 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)
MTX63_IT1_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 1)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX63_IT2.5 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)
MTX63_IT2.5_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 2.5)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX63_IT5 <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)
MTX63_IT5_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IT")*(med_df_MTX$Dose == 5)*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)


########################## SETTING UP IV COEFFICIENTS ##########################
################################################################################

MTX24_IV <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IV")
MTX24_IV_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 24)*(med_df_MTX$Route == "IV")*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX42_IV <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IV")
MTX42_IV_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 42)*(med_df_MTX$Route == "IV")*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)
MTX63_IV <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IV")
MTX63_IV_sex <- (med_df_MTX$Treatment == "MTX")*(med_df_MTX$Age == 63)*(med_df_MTX$Route == "IV")*ifelse(med_df_MTX$Sex == "F",0.5,-0.5)


########################## TISSUE PLOTTING ##########################

df_dose <- med_df_MTX %>% filter(Dose == 5)

med_df_MTX_ranef <- med_df_MTX
for (j in 1:length(structure_list)){
  cstruct <- structure_list[j]
  print(paste(as.character(j), cstruct, sep=": "))
  first_form <- glue("`{cstruct}` ~ -1 + age_flag:Dose + age_flag:route_flag + age_flag:sex_flag + 
  MTX24_IT0.5 + MTX24_IT0.5_sex + 
  MTX24_IT1 + MTX24_IT1_sex + 
  MTX24_IT2.5 + MTX24_IT2.5_sex + 
  MTX24_IT5 + MTX24_IT5_sex + 
  MTX42_IT0.5 + MTX42_IT0.5_sex + 
  MTX42_IT1 + MTX42_IT1_sex + 
  MTX42_IT2.5 + MTX42_IT2.5_sex + 
  MTX42_IT5 + MTX42_IT5_sex + 
  MTX63_IT0.5 + MTX63_IT0.5_sex + 
  MTX63_IT1 + MTX63_IT1_sex + 
  MTX63_IT2.5 + MTX63_IT2.5_sex + 
  MTX63_IT5 + MTX63_IT5_sex + 
  MTX24_IV + MTX24_IV_sex + 
  MTX42_IV + MTX42_IV_sex + 
  MTX63_IV + MTX63_IV_sex + 
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

df_tissue <- data.frame(df_dose_ranef[,1:(which(colnames(df_dose_ranef) == "amygdala")-1)])
WM_list <- subset(dt.struc, tissue.type == "WM")[,1]
WM_list <- c(WM_list, "anterior_lobule_white_matter") #adding missing structure
df_tissue$WM <- rowSums(df_dose_ranef[,as.character(WM_list)], na.rm = TRUE)
GM_list <- subset(dt.struc, tissue.type == "GM")[,1]
GM_list <- c(GM_list, "midbrain", "medulla", "lateral_septum", "medial_septum") #adding missing structures
df_tissue$GM <- rowSums(df_dose_ranef[,GM_list], na.rm = TRUE)
HP_list <- subset(dt.struc, hierarchy == "Hippocampal region")[,1]
df_tissue$HP <- rowSums(df_dose_ranef[,HP_list], na.rm = TRUE)
OB_list <- subset(dt.struc, hierarchy == "Olfactory bulbs")[,1]
df_tissue$OB<- rowSums(df_dose_ranef[,OB_list], na.rm = TRUE)
AV_list <- subset(dt.struc, ABI == "arbor vitae")[,1]
df_tissue$AV<- rowSums(df_dose_ranef[,AV_list], na.rm = TRUE)

#Whole brain
WB_list <- dt.struc[, 1]
df_tissue$WB <- rowSums(df_dose_ranef[,WB_list], na.rm = TRUE)

########################## FIG 2A: WHOLE BRAIN ########################## 

plot_WB <- ggplot(df_tissue,aes(x=as.character(Age),y=df_tissue[,"WB"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Whole Brain") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none", axis.title.x=element_blank())


########################## FIG 2B: GREY MATTER ##########################

plot_GM <- ggplot(df_tissue,aes(x=as.character(Age),y=df_tissue[,"GM"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Grey Matter") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())

########################## FIG 2C: WHITE MATTER ##########################

plot_WM <- ggplot(df_tissue,aes(x=as.character(Age),y=df_tissue[,"WM"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="White Matter") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none")

########################## FIG 2D: HIPPOCAMPUS ##########################

plot_HP <- ggplot(df_tissue,aes(x=as.character(Age),y=df_tissue[,"HP"],colour=Treatment))+
  geom_point(alpha=0.4, aes(shape=Treatment, colour=Treatment), size=2, position=position_dodge(width=0.8))+
  stat_summary(fun=mean, geom="point", aes(shape=Treatment, colour=Treatment), size=3.5, position=position_dodge(width=0.8))+
  stat_summary(fun.data=mean_ci,geom="errorbar", width=0.8, size=0.4, alpha=0.9, position=position_dodge(width=0.8))+
  facet_wrap(~Route) +
  scale_color_manual(values = c("Saline" = "blue", "MTX" = "red")) +
  theme_bw()+
  labs(x="Age", y="Volume (mm³)", title="Hippocampus") +
  theme(text=element_text(size=15,  family="Arial")) +
  theme(legend.position = "none", axis.title.y=element_blank())


# RUN TO PRODUCE ALL PLOTS
plot_WB
plot_GM
plot_WM
plot_HP

########################## PUT ALL PANELS TOGETHER ##########################

#png(filename="/projects/schoi/MTXcomparison/figures/fig2_v1.png", width = 12, height = 8, units = 'in', res = 600)
#grid.arrange(grobs = list(plot_WB, plot_GM, 
#                          plot_WM, plot_HP),
#             ncol=2,nrow=2)
#dev.off()



########################## STATS CALCULATIONS FOR EACH PLOT ##########################

df_tissue_IV <- df_tissue %>% filter(Route == 'IV')
df_tissue_IT <- df_tissue %>% filter(Route == 'IT')

##### IV WB #####
# Perform t-test for each age time point
t_test_IV_WB <- by(df_tissue_IV, df_tissue_IV$Age, function(sub_df) {
  t_test <- t.test(WB ~ Treatment, data = sub_df)
  return(t_test)
})

# Print the results
print(t_test_IV_WB)

# Extract p-values
p_values_IV_WB <- sapply(t_test_IV_WB, function(t_test) t_test$p.value)

# Print p-values
print(p_values_IV_WB)

##### IV GM ##### 
# Perform t-test for each age time point
t_test_IV_GM <- by(df_tissue_IV, df_tissue_IV$Age, function(sub_df) {
  t_test <- t.test(GM ~ Treatment, data = sub_df)
  return(t_test)
})

# Print the results
print(t_test_IV_GM)

# Extract p-values
p_values_IV_GM <- sapply(t_test_IV_GM, function(t_test) t_test$p.value)

# Print p-values
print(p_values_IV_GM)

##### IV WM #####
# Perform t-test for each age time point
t_test_IV_WM <- by(df_tissue_IV, df_tissue_IV$Age, function(sub_df) {
  t_test <- t.test(WM ~ Treatment, data = sub_df)
  return(t_test)
})

# Print the results
print(t_test_IV_WM)

# Extract p-values
p_values_IV_WM <- sapply(t_test_IV_WM, function(t_test) t_test$p.value)

# Print p-values
print(p_values_IV_WM)


##### IV HP #####
# Perform t-test for each age time point
t_test_IV_HP <- by(df_tissue_IV, df_tissue_IV$Age, function(sub_df) {
  t_test <- t.test(HP ~ Treatment, data = sub_df)
  return(t_test)
})

# Print the results
print(t_test_IV_HP)

# Extract p-values
p_values_IV_HP <- sapply(t_test_IV_HP, function(t_test) t_test$p.value)

# Print p-values
print(p_values_IV_HP)


##### IT WB #####
# Perform t-test for each age time point
t_test_IT_WB <- by(df_tissue_IT, df_tissue_IT$Age, function(sub_df) {
  t_test <- t.test(WB ~ Treatment, data = sub_df)
  return(t_test)
})

# Print the results
print(t_test_IT_WB)

# Extract p-values
p_values_IT_WB <- sapply(t_test_IT_WB, function(t_test) t_test$p.value)

# Print p-values
print(p_values_IT_WB)

##### IT GM #####
# Perform t-test for each age time point
t_test_IT_GM <- by(df_tissue_IT, df_tissue_IT$Age, function(sub_df) {
  t_test <- t.test(GM ~ Treatment, data = sub_df)
  return(t_test)
})

# Print the results
print(t_test_IT_GM)

# Extract p-values
p_values_IT_GM <- sapply(t_test_IT_GM, function(t_test) t_test$p.value)

# Print p-values
print(p_values_IT_GM)

##### IT WM #####
# Perform t-test for each age time point
t_test_IT_WM <- by(df_tissue_IT, df_tissue_IT$Age, function(sub_df) {
  t_test <- t.test(WM ~ Treatment, data = sub_df)
  return(t_test)
})

# Print the results
print(t_test_IT_WM)

# Extract p-values
p_values_IT_WM <- sapply(t_test_IT_WM, function(t_test) t_test$p.value)

# Print p-values
print(p_values_IT_WM)

##### IT HP #####
# Perform t-test for each age time point
t_test_IT_HP <- by(df_tissue_IT, df_tissue_IT$Age, function(sub_df) {
  t_test <- t.test(HP ~ Treatment, data = sub_df)
  return(t_test)
})

# Print the results
print(t_test_IT_HP)

# Extract p-values
p_values_IT_HP <- sapply(t_test_IT_HP, function(t_test) t_test$p.value)

# Print p-values
print(p_values_IT_HP)


########################## TAKE LEGEND FROM HERE ##########################

#png(filename="/projects/schoi/MTXcomparison/figures/fig5_legend.png", width = 12, height = 12, units = 'in', res = 600)
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
mean_values <- df_tissue %>%
  group_by(Age, Treatment, Route) %>%
  summarise(
    mean_WB = mean(WB, na.rm = TRUE),
    mean_GM = mean(GM, na.rm = TRUE),
    mean_WM = mean(WM, na.rm = TRUE), 
    mean_HP = mean(HP, na.rm = TRUE), 
    .groups = "drop"
  )

# Pivot to get Saline and MTX means together, split by Route
mean_values_wide <- mean_values %>%
  pivot_wider(names_from = Treatment, values_from = c(mean_WB, mean_GM, mean_WM, mean_HP))

# Calculate percent change for each route
mean_values_wide <- mean_values_wide %>%
  mutate(
    percent_change_WB = (mean_WB_MTX - mean_WB_Saline) / mean_WB_Saline * 100,
    percent_change_GM = (mean_GM_MTX - mean_GM_Saline) / mean_GM_Saline * 100,
    percent_change_WM = (mean_WM_MTX - mean_WM_Saline) / mean_WM_Saline * 100,
    percent_change_HP = (mean_HP_MTX - mean_HP_Saline) / mean_HP_Saline * 100
  )

# View the result
print(mean_values_wide)
