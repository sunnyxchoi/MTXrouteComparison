### Weights analysis ###

########################## LOAD DATA ##########################

# Use this statement below in all scripts to load everything in this environment
load("/projects/schoi/MTXcomparison/analysis/MTX17022024/med_MTX.RData")

########################## LOAD PACKAGES ##########################

library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyverse)
library(RMINC)
library(MRIcrotome)
library(reshape2)
library(glue)
library(dplyr)
library(splines)


########################## FIG 1A: EXPERIMENTAL TIMELINE ########################## 
# Made timeline on powerpoint


########################## SET UP DATA ########################## 

# Making weights dfs
ITweights <- read.csv("/projects/schoi/MTXcomparison/analysis/sep2024analysis/MTXdata.csv")
# Use sub or gsub to remove the patterns
ITweights$P14.Dist.Corr <- gsub("/projects/schoi/IT_MTX/images/dist_corr/Gr../P../|/projects/schoi/IT_MTX/images/dist_corr/Gr.../P../", "", ITweights$P14.Dist.Corr)
# Keep only unique rows based on P14.Dist.Corr
ITweights <- ITweights[!duplicated(ITweights$P14.Dist.Corr), ]

IVweights <- read.csv("/projects/schoi/MTXcomparison/analysis/sep2024analysis/MTXdatasystemic.csv")
IVweights <- IVweights[!duplicated(IVweights$Dist.Corr.P14), ]

df_ITweights <- data.frame(ITweights$P14.Dist.Corr, ITweights$P13.Weight, ITweights$P17.Weight, ITweights$P18.Weight,
                           ITweights$P19.Weight, ITweights$P20.Weight, ITweights$P21.Weight, ITweights$P22.Weight, ITweights$P23.Weight,
                           ITweights$P25.Weight, ITweights$P26.Weight, ITweights$P27.Weight,
                           ITweights$P28.Weight, ITweights$P35.Weight, ITweights$P41.Weight, ITweights$P49.Weight,
                           ITweights$P56.Weight, ITweights$P62.Weight)

colnames(df_ITweights) <- c("Image", "13", "17", "18", "19", "20", "21", "22", "23", "25", "26", "27", "28", "35", "41", "49", "56", "62")

df_IVweights <- data.frame(IVweights$Dist.Corr.P14, IVweights$P13.Weight, IVweights$P17.Weight, IVweights$P19.Weight, 
                           IVweights$P20.Weight, IVweights$P21Weight, IVweights$P22Weight, IVweights$P23.Weight,
                           IVweights$P28.Weight, IVweights$P35.Weight, IVweights$P41.Weight, IVweights$P49.Weight,
                           IVweights$P56.Weight, IVweights$P62.Weight)

colnames(df_IVweights) <- c("Image", "13", "17", "19", "20", "21", "22", "23", "28", "35", "41", "49", "56", "62")

#Adding NAs instead of spaces
df_IVweights[, -1] <- lapply(df_IVweights[, -1], function(x) as.numeric(as.character(x)))


# Making mouse info df
mouseinfo_df <- med_df_MTX[, c(1:9, 11:13)]

mouseinfo_df <- mouseinfo_df[mouseinfo_df$Image %in% unique(subset(mouseinfo_df, Age==14)$Image),]

# Combining IT and IV weights data
combined_df <- bind_rows(df_ITweights, df_IVweights)

joined_df <- left_join(mouseinfo_df, combined_df, by="Image")


# Convert to long format
joinedlong_df <- joined_df %>% gather(Age, Weights, "13":"62")
joinedlong_df$Age <- as.numeric(joinedlong_df$Age)

joinedlong_df$Group <- paste(joinedlong_df$Route, joinedlong_df$Dose, joinedlong_df$Treatment, joinedlong_df$Sex, sep = "_")
joinedlong_df$Group <- factor(joinedlong_df$Group, levels = c("IT_0.5_Saline_M", "IT_0.5_Saline_F", "IT_0.5_MTX_M", "IT_0.5_MTX_F", 
                                                              "IT_1_Saline_M", "IT_1_Saline_F", "IT_1_MTX_M", "IT_1_MTX_F", 
                                                              "IT_2.5_Saline_M", "IT_2.5_Saline_F", "IT_2.5_MTX_M", "IT_2.5_MTX_F", 
                                                              "IT_5_Saline_M", "IT_5_Saline_F", "IT_5_MTX_M", "IT_5_MTX_F",
                                                              "IV_5_Saline_M", "IV_5_Saline_F", "IV_5_MTX_M", "IV_5_MTX_F"))


########################## FIG 1B: WEIGHTS CURVES ########################## 

# Filter for dose == 5 and create groups for plot
joinedlong_df5 <- joinedlong_df %>% filter(Dose ==5)
joinedlong_df5$GroupPlot <- paste(joinedlong_df5$Route, joinedlong_df5$Treatment, sep = "_")

# Plot
ggplot(joinedlong_df5, aes(x = Age, y = Weights, color = GroupPlot)) +
  geom_smooth(aes(group = GroupPlot), method = "loess", se = FALSE, size = 1) +
  geom_point(size = 1.5, alpha = 0.3) +  # Add points for the weights
  labs(x = "Age", y = "Weights", color = "GroupPlot") +  # Label the axes
  scale_color_discrete(labels=c("IT Saline", "IT MTX (5.0 mg/kg)", "IV Saline", "IV MTX (5.0 mg/kg)")) + 
  facet_grid(fct_rev(Sex) ~ Route) +
  theme_minimal() +  # Use a minimal theme
  theme(legend.position = "right")  # Position the legend



########################## FIG 1C: NORMALIZED WEIGHT PLOTS ########################## 


# Making weights dfs
ITweights <- read.csv("/projects/schoi/MTXcomparison/analysis/sep2024analysis/MTXdata.csv")
# Use sub or gsub to remove the patterns
ITweights$P14.Dist.Corr <- gsub("/projects/schoi/IT_MTX/images/dist_corr/Gr../P../|/projects/schoi/IT_MTX/images/dist_corr/Gr.../P../", "", ITweights$P14.Dist.Corr)
ITweights$Treatment <- gsub("Sal", "Saline", ITweights$Treatment)

IVweights <- read.csv("/projects/schoi/MTXcomparison/analysis/sep2024analysis/MTXdatasystemic.csv")
# Remove rows where the "ID" column is empty
IVweights <- IVweights[IVweights$ID != "", ]
IVweights$Sex <- gsub("Male", "M", IVweights$Sex)
IVweights$Sex <- gsub("Female", "F", IVweights$Sex)


df_ITweights <- data.frame(ITweights$ID, ITweights$Dose, ITweights$Route, ITweights$Treatment, ITweights$Sex, ITweights$P14.Dist.Corr, ITweights$P13.Weight, ITweights$P17.Weight, ITweights$P18.Weight,
                           ITweights$P19.Weight, ITweights$P20.Weight, ITweights$P21.Weight, ITweights$P22.Weight, ITweights$P23.Weight,
                           ITweights$P25.Weight, ITweights$P26.Weight, ITweights$P27.Weight,
                           ITweights$P28.Weight, ITweights$P35.Weight, ITweights$P41.Weight, ITweights$P49.Weight,
                           ITweights$P56.Weight, ITweights$P62.Weight)

colnames(df_ITweights) <- c("ID", "Dose", "Route", "Treatment", "Sex", "Image", "13", "17", "18", "19", "20", "21", "22", "23", "25", "26", "27", "28", "35", "41", "49", "56", "62")

df_IVweights <- data.frame(IVweights$ID, Dose=5.0, Route="IV", IVweights$Treatment, IVweights$Sex, IVweights$Dist.Corr.P14, IVweights$P13.Weight, IVweights$P17.Weight, ExtraCol = NA, IVweights$P19.Weight, 
                           IVweights$P20.Weight, IVweights$P21Weight, IVweights$P22Weight, IVweights$P23.Weight, ExtraCol = NA, ExtraCol = NA, ExtraCol = NA, 
                           IVweights$P28.Weight, IVweights$P35.Weight, IVweights$P41.Weight, IVweights$P49.Weight,
                           IVweights$P56.Weight, IVweights$P62.Weight)

colnames(df_IVweights) <- c("ID", "Dose", "Route", "Treatment", "Sex", "Image", "13", "17", "18", "19", "20", "21", "22", "23", "25", "26", "27", "28", "35", "41", "49", "56", "62")



# Making mouse info df
#mouseinfo_df <- df_MTX[, c(1:9, 11:13)]

#mouseinfo_df <- mouseinfo_df[mouseinfo_df$Image %in% unique(subset(mouseinfo_df, Age==14)$Image),]

# Combining IT and IV weights data
#combined_df <- bind_rows(df_ITweights, df_IVweights)
combined_df <- rbind(df_ITweights, df_IVweights)

#joined_df <- left_join(mouseinfo_df, combined_df, by="Image")
#joined_df <- cbind(mouseinfo_df, combined_df)


# Convert to long format
joinedlong_df <- combined_df %>% gather(Age, Weights, "13":"62")
joinedlong_df$Age <- as.numeric(joinedlong_df$Age)

joinedlong_df$Group <- paste(joinedlong_df$Route, joinedlong_df$Dose, joinedlong_df$Treatment, joinedlong_df$Sex, sep = "_")
joinedlong_df$Group <- factor(joinedlong_df$Group, levels = c("IT_0.5_Saline_M", "IT_0.5_Saline_F", "IT_0.5_MTX_M", "IT_0.5_MTX_F", 
                                                              "IT_1_Saline_M", "IT_1_Saline_F", "IT_1_MTX_M", "IT_1_MTX_F", 
                                                              "IT_2.5_Saline_M", "IT_2.5_Saline_F", "IT_2.5_MTX_M", "IT_2.5_MTX_F", 
                                                              "IT_5_Saline_M", "IT_5_Saline_F", "IT_5_MTX_M", "IT_5_MTX_F",
                                                              "IV_5_Saline_M", "IV_5_Saline_F", "IV_5_MTX_M", "IV_5_MTX_F"))

joinedlong_df$GroupNoSex <- paste(joinedlong_df$Route, joinedlong_df$Dose, joinedlong_df$Treatment, sep = "_")
joinedlong_df$GroupNoSex <- factor(joinedlong_df$GroupNoSex, levels = c("IT_0.5_Saline", "IT_0.5_MTX",
                                                                        "IT_1_Saline", "IT_1_MTX",
                                                                        "IT_2.5_Saline", "IT_2.5_MTX",
                                                                        "IT_5_Saline", "IT_5_MTX",
                                                                        "IV_5_Saline", "IV_5_MTX"))

joinedlong_df$GroupNoTreatment <- paste(joinedlong_df$Route, joinedlong_df$Dose, sep = "_")
joinedlong_df$GroupNoTreatment <- factor(joinedlong_df$GroupNoTreatment, levels = c("IT_0.5", "IT_1", "IT_2.5", "IT_5", "IV_5"))

weights_df <- joinedlong_df %>% filter(Age==13|Age==17|Age==19|Age==20|Age==21|Age==22|Age==23|Age==28|Age==35|Age==41|Age==49|Age==56|Age==62)
weights_df$Age <- as.factor(weights_df$Age)


joinedlong_df5 <- subset(joinedlong_df, Dose == 5)
joinedlong_df5$Group <- droplevels(joinedlong_df5$Group)


###


# Plotting at each timepoint

df_ages <- joinedlong_df %>% filter(Age %in% c(17, 23, 35, 62))

###
#Plots relative to control average

ctrl_weights <- df_ages %>%
  filter(Treatment == "Saline") %>%
  group_by(Age, Sex, GroupNoSex, GroupNoTreatment) %>%
  summarise(mean_weight = mean(Weights, na.rm = TRUE)) %>%
  ungroup()

df_ages <- df_ages %>%
  mutate(Weights = as.numeric(Weights))
ctrl_weights <- ctrl_weights %>%
  mutate(mean_weight = as.numeric(mean_weight))

ctrl_weights <- df_ages %>%
  filter(Treatment == "Saline") %>%
  group_by(Age, Sex, GroupNoSex, GroupNoTreatment) %>%
  summarise(mean_weight = mean(Weights, na.rm = TRUE)) %>%
  ungroup()

# Join df_ages with ctrl_weights to add the control (Saline) mean weight
df_ages_normalized <- df_ages %>%
  left_join(ctrl_weights, by = c("Age", "Sex", "GroupNoTreatment")) %>%
  mutate(normalized_weight = ifelse(Treatment == "MTX", Weights / mean_weight, NA))

df_summary_normalized <- df_ages_normalized %>%
  group_by(GroupNoTreatment, Age, Sex) %>%
  summarize(
    Mean = mean(normalized_weight, na.rm = TRUE),
    LowerCI = Mean - qt(0.975, df = n() - 1) * sd(normalized_weight, na.rm = TRUE) / sqrt(n()),
    UpperCI = Mean + qt(0.975, df = n() - 1) * sd(normalized_weight, na.rm = TRUE) / sqrt(n())
  )

x_labels <- c("IT (0.5 mg/kg)", "IT (1.0 mg/kg)", "IT (2.5 mg/kg)", "IT (5.0 mg/kg)", "IV (5.0 mg/kg)")

# Create the ggplot
ggplot(df_summary_normalized) +
  aes(x = GroupNoTreatment, y = Mean, color = GroupNoTreatment) + 
  # Add horizontal lines for means
  geom_crossbar(
    aes(ymin = Mean, ymax = Mean), 
    width = 0.5, position = position_dodge(width = 0.8), size = 0.7
  ) +
  # Add error bars for confidence intervals
  geom_errorbar(
    aes(ymin = LowerCI, ymax = UpperCI), 
    width = 0.3, position = position_dodge(width = 0.8), size = 1
  ) +
  # Add jittered points for normalized weights
  geom_jitter(
    data = df_ages_normalized, 
    aes(x = GroupNoTreatment, y = normalized_weight), 
    alpha = 0.4, position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8)
  ) +
  # Add a horizontal reference line at y = 1
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  # Add labels and custom x-axis labels
  labs(
    title = "",
    x = "Group",
    y = "Normalized Mass (Treatment/Control)",
    color = "Group"
  ) +
  scale_x_discrete(labels = x_labels) +
  # Facet by Age (rows) and reversed Sex (columns)
  facet_grid(rows = vars(Age), cols = vars(fct_rev(Sex))) +
  # Apply a clean theme
  theme_bw() +
  # Customize text and layout
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(color = "gray20", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14),  # Facet label size
    legend.position = "none"  # Remove the legend
  )


########################## STATS CALCULATIONS FOR EACH AGE ##########################

df_p17 <- df_ages %>% filter(Age == 17)
df_p23 <- df_ages %>% filter(Age == 23)
df_p35 <- df_ages %>% filter(Age == 35)
df_p62 <- df_ages %>% filter(Age == 62)

##### IV WB #####
# Perform t-test for each age time point

# Remove any rows where 'Weights' or 'Treatment' is missing
df_p17_clean <- df_p17 %>% filter(!is.na(Weights), !is.na(Treatment))
df_p23_clean <- df_p23 %>% filter(!is.na(Weights), !is.na(Treatment))
df_p35_clean <- df_p35 %>% filter(!is.na(Weights), !is.na(Treatment))
df_p62_clean <- df_p62 %>% filter(!is.na(Weights), !is.na(Treatment))

# PICK WHICH AGE DF TO USE HERE
df_selected <- df_p23_clean

# Run t-test
t_test_results <- df_selected %>%
  group_by(Sex, GroupNoTreatment) %>% 
  summarise(
    t_test = list(t.test(Weights ~ Treatment, data = cur_data())),
    .groups = "drop"
  )

# Extract p-values and other details if needed
t_test_summary <- t_test_results %>%
  mutate(
    p_value = map_dbl(t_test, ~ .x$p.value),
    statistic = map_dbl(t_test, ~ .x$statistic)
  )

print(t_test_summary)

# Calculate diff between treatment and control group
mean_diff = df_selected %>%
  group_by(Sex, GroupNoTreatment, Treatment) %>%
  summarise(mean_weight = mean(Weights, na.rm = TRUE), .groups = "drop") %>%
  spread(Treatment, mean_weight) %>%
  mutate(change = ((MTX - Saline)/Saline)*100)

print(mean_diff)


# Calculate diff between treatment and control group WITHOUT SEX GROUPS
mean_diff = df_selected %>%
  group_by(GroupNoTreatment, Treatment) %>%
  summarise(mean_weight = mean(Weights, na.rm = TRUE), .groups = "drop") %>%
  spread(Treatment, mean_weight) %>%
  mutate(change = ((MTX - Saline)/Saline)*100)

print(mean_diff)


#######
# Run t-test between sexes - replace first df_selected above by whichever age you want to test
# Split Saline and MTX groups

saline_data <- df_selected[df_selected$Treatment == "Saline", c("GroupNoTreatment", "Sex", "Weights")] %>%
  rename(Control_Weight = Weights)

mtx_data <- df_selected[df_selected$Treatment == "MTX", c("GroupNoTreatment", "Sex", "Weights")] %>%
  rename(Treatment_Weight = Weights)

# Create all possible pairs using a cross join
mtx_decrease <- saline_data %>%
  inner_join(mtx_data, by = c("GroupNoTreatment", "Sex")) %>%
  mutate(decrease = Control_Weight - Treatment_Weight)

# View the result
print(mtx_decrease)

# Run t-test
t_test_results <- mtx_decrease %>%
  group_by(GroupNoTreatment) %>% 
  summarise(
    t_test = list(t.test(decrease ~ Sex, data = cur_data())),
    .groups = "drop"
  )

# Extract p-values and other details if needed
t_test_summary <- t_test_results %>%
  mutate(
    p_value = map_dbl(t_test, ~ .x$p.value),
    statistic = map_dbl(t_test, ~ .x$statistic)
  )

print(t_test_summary)


### Comparing routes

# Run t-test between Route IV and Route IT for Weights, split by Sex
t_test_route_sex <- df_selected %>%
  filter(Route %in% c("IV", "IT")) %>%
  group_by(Sex) %>%
  summarise(
    t_test = list(t.test(Weights ~ Route, data = .)),
    .groups = "drop"
  )

# Extract p-values and other details if needed
t_test_summary_route_sex <- t_test_route_sex %>%
  mutate(
    p_value = map_dbl(t_test, ~ .x$p.value),
    statistic = map_dbl(t_test, ~ .x$statistic)
  )

# Print the results
print(t_test_summary_route_sex)
