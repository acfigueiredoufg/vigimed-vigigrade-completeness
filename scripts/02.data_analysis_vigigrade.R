library(tidyverse)
library(summarytools)

# DATA ANALYSIS ------------------------------------------------------------
## DESCRIPTIVE ANALYSIS ----------------------------------------------------

# Overview of the final ICSR dataset
glimpse(bd_icsr_final)

# Reporter type
freq(bd_icsr_final$Reporter_Type, order = "freq")

# State or region
freq(
  bd_icsr_final$State_Region,
  order = "freq",
  report.nas = FALSE
)

# Report quality classification
freq(bd_icsr_final$qualidade, order = "freq")

# Measures of central tendency for mean vigiGrade score
summary(bd_icsr_final$vigiGrade_medio)          # Base R
descr(bd_icsr_final$vigiGrade_medio)            # summarytools
describe(bd_icsr_final["vigiGrade_medio"])      # summarytools (data frame format)

# ATC classification (Level 1)
freq(bd_atc_final$ATC_Nivel1, order = "freq")

# -------------------------------------------------------------------------
# COMPLETENESS ANALYSIS
# -------------------------------------------------------------------------
# Completeness summary using renamed variables from bd_vigigrade

completeness_analysis <- bd_vigigrade %>%
  group_by(Report_ID) %>%
  summarise(across(starts_with("P_"), ~ first(.))) %>%
  ungroup() %>%
  summarise(
    `Time to Onset`  = mean(P_Time_to_Onset == 1.0) * 100,
    `Indication`     = mean(P_Indication == 1.0) * 100,
    `Outcome`        = mean(P_Outcome == 1.0) * 100,
    `Sex`            = mean(P_Sex == 1.0) * 100,
    `Age`            = mean(P_Age == 1.0) * 100,
    `Dose`           = mean(P_Dose == 1.0) * 100,
    `State`          = mean(P_Geographic_Region == 1.0) * 100,
    `Reporter`       = mean(P_Reporter_Qualif == 1.0) * 100,
    `Report Type`    = mean(P_Report_Type == 1.0) * 100,
    `Comments`       = mean(P_Narrative == 1.0) * 100
  ) %>%
  pivot_longer(
    everything(),
    names_to = "Dimension",
    values_to = "Completeness_Pct"
  ) %>%
  arrange(desc(Completeness_Pct))

# -------------------------------------------------------------------------
# FIGURE 1 – COMPLETENESS BAR PLOT
# -------------------------------------------------------------------------

fig_1 <- ggplot(
  completeness_analysis,
  aes(
    x = Completeness_Pct,
    y = reorder(Dimension, Completeness_Pct),
    fill = Dimension
  )
) +

  # Bars
  geom_col(width = 0.75, show.legend = FALSE) +

  # Color palette (ggsci – Nature Publishing Group)
  scale_fill_npg() +

  # Percentage labels at the end of bars
  geom_text(
    aes(label = sprintf("%.1f%%", Completeness_Pct)),
    hjust = -0.15,
    color = "black",
    fontface = "bold",
    size = 3.5
  ) +

  # Axis limits and scaling
  scale_x_continuous(limits = c(0, 115), expand = c(0, 0)) +
  theme_minimal(base_size = 12) +

  # Labels
  labs(
    title = "Completeness of Information in Individual Case Safety Reports (ICSR)",
    x = NULL,
    y = NULL
  ) +

  # Theme refinement
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      color = "black",
      hjust = 0.5
    ),
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_line(color = "black")
  )

# -------------------------------------------------------------------------
# EXPORT FIGURES
# -------------------------------------------------------------------------

# Save as PNG (high quality for presentations)
ggsave(
  filename = "Completeness_Analysis.png",
  plot = fig_1,
  device = "png",
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Save as TIFF (journal standard)
ggsave(
  filename = "Completeness_Analysis.tiff",
  plot = fig_1,
  device = "tiff",
  width = 10,
  height = 7,
  dpi = 300,
  compression = "lzw"
)

# Save as PDF (vector format)
ggsave(
  filename = "Completeness_Analysis.pdf",
  plot = fig_1,
  device = "pdf",
  width = 10,
  height = 7
)

# Proportion of well-documented reports
freq(bd_icsr_final$qualidade)

# Central tendency measures (redundancy kept for transparency)
summary(bd_icsr_final$vigiGrade_medio)
descr(bd_icsr_final$vigiGrade_medio)
describe(bd_icsr_final["vigiGrade_medio"])
psych::describe(bd_icsr_final$vigiGrade_medio)

# Distribution histogram
hist(bd_icsr_final$vigiGrade_medio)

## INFERENTIAL ANALYSIS -----------------------------------------------------

# Unit of analysis: ICSR

# Normality assessment for mean vigiGrade score
# Lilliefors test
lillie.test(bd_icsr_final$vigiGrade_medio)

# Shapiro-Wilk test (subsample due to large N)
shapiro.test(sample(bd_icsr_final$vigiGrade_medio, 5000))

# Histogram with density curve
ggplot(bd_icsr_final, aes(x = vigiGrade_medio)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightblue", color = "black") +
  geom_density(color = "red", size = 1) +
  labs(title = "Distribution of the Mean vigiGrade Score")

# Q-Q plot
qqnorm(bd_icsr_final$vigiGrade_medio)
qqline(bd_icsr_final$vigiGrade_medio, col = "red")

# -------------------------------------------------------------------------
# REPORTER TYPE × VIGIGRADE (KRUSKAL-WALLIS)
# -------------------------------------------------------------------------

# Kruskal-Wallis test
kruskal_result_reporter <- kruskal.test(vigiGrade_medio ~ Reporter_Type, data = bd_icsr_final)
print(kruskal_result_reporter)

# Dunn post-hoc test (Bonferroni correction)
dunn_result_reporter <- dunnTest(vigiGrade_medio ~ Reporter_Type, data = bd_icsr_final, method = "bonferroni")
print(dunn_result_reporter)

# Extract comparisons table
dunn_table_reporter <- dunn_result_reporter$res
print(dunn_table_reporter)

# Descriptive table with median and IQR
boot_median <- function(data, indices) {
  return(median(data[indices], na.rm = TRUE))
}

set.seed(123)

tabela_tend_central_reporter <- bd_icsr_final %>%
  group_by(Reporter_Type) %>%
  summarise(
    Median = round(median(vigiGrade_medio, na.rm = TRUE), 3),
    Q1 = round(quantile(vigiGrade_medio, 0.25, na.rm = TRUE), 3),
    Q3 = round(quantile(vigiGrade_medio, 0.75, na.rm = TRUE), 3)
  ) %>%
  mutate(IQR = paste0(Q1, " – ", Q3)) %>%
  select(Reporter_Type, Median, IQR)

# Format Dunn results into a readable table with adjusted p-values
tabela_dunn_reporter <- dunn_result_reporter$res %>%
  mutate(
    Group_A = str_trim(str_split_fixed(Comparison, " - ", 2)[, 1]),
    Group_B = str_trim(str_split_fixed(Comparison, " - ", 2)[, 2]),

    p_value = case_when(
      P.adj < 0.001 ~ "< 0.001",
      TRUE ~ formatC(P.adj, format = "f", digits = 3)
    ),

    Significant = ifelse(P.adj < 0.05, "Yes", "No")
  ) %>%
  select(Group_A, Group_B, Z, p_value, Significant)

# Effect size (Rosenthal's r)
N_total <- sum(!is.na(bd_icsr_final$vigiGrade_medio))

tabela_dunn_reporter <- tabela_dunn_reporter %>%
  mutate(
    Effect_size_r = round(abs(Z) / sqrt(N_total), 3),
    Effect_interpretation = case_when(
      Effect_size_r < 0.1 ~ "Very small",
      Effect_size_r < 0.3 ~ "Small",
      Effect_size_r < 0.5 ~ "Moderate",
      TRUE ~ "Large"
    )
  )

# Join median/IQR for Group A
tabela_final_reporter <- tabela_dunn_reporter %>%
  left_join(tabela_tend_central_reporter, by = c("Group_A" = "Reporter_Type")) %>%
  rename(Median_A = Median, IQR_A = IQR)

# Join median/IQR for Group B
tabela_final_reporter <- tabela_final_reporter %>%
  left_join(tabela_tend_central_reporter, by = c("Group_B" = "Reporter_Type")) %>%
  rename(Median_B = Median, IQR_B = IQR)

print(tabela_final_reporter)

# Final table for manuscript
tabela_reporter_multiple_comp <- tabela_final_reporter %>%
  select(
    `Group A` = Group_A,
    `Group B` = Group_B,
    `Median (A)` = Median_A,
    `IQR (A)` = IQR_A,
    `Median (B)` = Median_B,
    `IQR (B)` = IQR_B,
    `Z statistic` = Z,
    `p-value` = p_value,
    `Effect size (r)` = Effect_size_r,
    `Effect interpretation` = Effect_interpretation,
    `Significant` = Significant
  )

# Significant comparisons only (p < 0.05)
tabela_signif_reporter <- tabela_reporter_multiple_comp %>%
  filter(Significant == "Yes")

# -------------------------------------------------------------------------
# EXPORT TABLES (WORD)
# -------------------------------------------------------------------------

# Multiple comparisons table
multiple_comparisons_reporter <- flextable(tabela_reporter_multiple_comp)

multiple_comparisons_reporter_docx <- read_docx() %>%
  body_add_par("Multiple Comparisons Table - Reporter Type", style = "heading 1") %>%
  body_add_flextable(multiple_comparisons_reporter)

print(multiple_comparisons_reporter_docx, target = "Table_Multiple_Comparisons_Reporter.docx")

# Significant comparisons table
significant_comparisons_reporter <- flextable(tabela_signif_reporter)

significant_comparisons_reporter_docx <- read_docx() %>%
  body_add_par("Significant Comparisons Table - Reporter Type", style = "heading 1") %>%
  body_add_flextable(significant_comparisons_reporter)

print(significant_comparisons_reporter_docx, target = "Table_Significant_Comparisons_Reporter.docx")

# -------------------------------------------------------------------------
# VIOLIN PLOT WITH SIGNIFICANT DUNN COMPARISONS
# -------------------------------------------------------------------------

# Keep only significant pairs
significant_pairs_reporter <- tabela_signif_reporter %>%
  select(`Group A`, `Group B`, `p-value`) %>%
  rename(group1 = `Group A`, group2 = `Group B`, p_text = `p-value`) %>%
  mutate(
    p = case_when(
      p_text == "< 0.001" ~ 0.0009,
      TRUE ~ as.numeric(p_text)
    ),
    y.position = seq(1.02, 1.02 + 0.04 * (n() - 1), by = 0.04)
  )

# Filter dataset to only groups involved in significant comparisons
bd_violin_reporter <- bd_icsr_final %>%
  filter(Reporter_Type %in% unique(c(tabela_signif_reporter$`Group A`, tabela_signif_reporter$`Group B`))) %>%
  mutate(
    Reporter_Type = factor(
      Reporter_Type,
      levels = c(
        "Consumer or non-healthcare professional",
        "Lawyer",
        "Other healthcare professional",
        "Pharmacist",
        "Physician"
      )
    )
  )

# Median and IQR per group (to display on plot)
medians_iqr_reporter <- bd_violin_reporter %>%
  group_by(Reporter_Type) %>%
  summarise(
    Median = round(median(vigiGrade_medio, na.rm = TRUE), 2),
    Q1 = round(quantile(vigiGrade_medio, 0.25, na.rm = TRUE), 2),
    Q3 = round(quantile(vigiGrade_medio, 0.75, na.rm = TRUE), 2)
  ) %>%
  mutate(
    label = paste0("Median: ", Median, "\nIQR: ", Q1, "–", Q3),
    y = -0.05
  )

# Violin plot with p-values
violin_plot_kw_reporter <- ggplot(
  bd_violin_reporter,
  aes(x = Reporter_Type, y = vigiGrade_medio, fill = Reporter_Type)
) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  geom_text(
    data = medians_iqr_reporter,
    aes(x = Reporter_Type, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 1,
    fontface = "italic"
  ) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2, fill = "white") +
  stat_pvalue_manual(
    data = significant_pairs_reporter,
    label = "p_text",
    y.position = seq(1.05, 1.55, length.out = nrow(significant_pairs_reporter)),
    tip.length = 0.01,
    size = 4
  ) +
  scale_fill_npg() +
  scale_x_discrete(labels = c(
    "Consumer or non-healthcare professional" = "Consumer or\nnon-healthcare professional",
    "Other healthcare professional" = "Other\nhealthcare professional",
    "Lawyer" = "Lawyer",
    "Pharmacist" = "Pharmacist",
    "Physician" = "Physician"
  )) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Reporter Type",
    y = "vigiGrade Score",
    title = "Comparison of vigiGrade Scores Across Reporter Types",
    caption = "Only significant differences (Dunn-Bonferroni, p < 0.05) are shown"
  )

print(violin_plot_kw_reporter)

# Export plots
ggsave(
  filename = "violin_plot_reporter.png",
  plot = violin_plot_kw_reporter,
  width = 10,
  height = 6,
  dpi = 600,
  units = "in",
  bg = "white"
)

ggsave(
  filename = "violin_plot_reporter.tiff",
  plot = violin_plot_kw_reporter,
  width = 10,
  height = 6,
  dpi = 600,
  units = "in",
  bg = "white"
)

ggsave(
  filename = "violin_plot_reporter.pdf",
  plot = violin_plot_kw_reporter,
  width = 10,
  height = 6,
  dpi = 600,
  units = "in",
  bg = "white"
)
# -------------------------------------------------------------------------
# STATE × VIGIGRADE (KRUSKAL–WALLIS)
# -------------------------------------------------------------------------

# Remove invalid state codes
bd_state <- bd_icsr_final %>%
  mutate(
    State_Region = if_else(State_Region == "C", NA_character_, State_Region)
  ) %>%
  filter(!is.na(State_Region))

# Kruskal-Wallis test
kruskal_result_state <- kruskal.test(vigiGrade_medio ~ State_Region, data = bd_state)
print(kruskal_result_state)

# Dunn post-hoc test (Bonferroni correction)
dunn_result_state <- dunnTest(vigiGrade_medio ~ State_Region, data = bd_state, method = "bonferroni")
print(dunn_result_state)

# Extract comparisons table
dunn_table_state <- dunn_result_state$res
print(dunn_table_state)

# -------------------------------------------------------------------------
# DESCRIPTIVE STATISTICS (MEDIAN AND IQR)
# -------------------------------------------------------------------------

set.seed(123)

table_central_tendency_state <- bd_state %>%
  group_by(State_Region) %>%
  summarise(
    Median = round(median(vigiGrade_medio, na.rm = TRUE), 3),
    Q1 = round(quantile(vigiGrade_medio, 0.25, na.rm = TRUE), 3),
    Q3 = round(quantile(vigiGrade_medio, 0.75, na.rm = TRUE), 3)
  ) %>%
  mutate(IQR = paste0(Q1, " – ", Q3)) %>%
  select(State_Region, Median, IQR)

# Format Dunn table
table_dunn_state <- dunn_result_state$res %>%
  mutate(
    Group_A = str_trim(str_split_fixed(Comparison, " - ", 2)[, 1]),
    Group_B = str_trim(str_split_fixed(Comparison, " - ", 2)[, 2]),

    p_value = case_when(
      P.adj < 0.001 ~ "< 0.001",
      TRUE ~ formatC(P.adj, format = "f", digits = 3)
    ),

    Significant = ifelse(P.adj < 0.05, "Yes", "No")
  ) %>%
  select(Group_A, Group_B, Z, p_value, Significant)

# Effect size (Rosenthal's r)
N_total_state <- sum(!is.na(bd_state$vigiGrade_medio))

table_dunn_state <- table_dunn_state %>%
  mutate(
    Effect_size_r = round(abs(Z) / sqrt(N_total_state), 3),
    Effect_interpretation = case_when(
      Effect_size_r < 0.1 ~ "Very small",
      Effect_size_r < 0.3 ~ "Small",
      Effect_size_r < 0.5 ~ "Moderate",
      TRUE ~ "Large"
    )
  )

# Join descriptive statistics
final_table_state <- table_dunn_state %>%
  left_join(table_central_tendency_state, by = c("Group_A" = "State_Region")) %>%
  rename(Median_A = Median, IQR_A = IQR) %>%
  left_join(table_central_tendency_state, by = c("Group_B" = "State_Region")) %>%
  rename(Median_B = Median, IQR_B = IQR)

print(final_table_state)

# Final table for manuscript
table_state_multiple_comp <- final_table_state %>%
  select(
    `Group A` = Group_A,
    `Group B` = Group_B,
    `Median (A)` = Median_A,
    `IQR (A)` = IQR_A,
    `Median (B)` = Median_B,
    `IQR (B)` = IQR_B,
    `Z statistic` = Z,
    `p-value` = p_value,
    `Effect size (r)` = Effect_size_r,
    `Effect interpretation` = Effect_interpretation,
    `Significant` = Significant
  )

# Significant comparisons only
table_significant_state <- table_state_multiple_comp %>%
  filter(Significant == "Yes")
# -------------------------------------------------------------------------
# RIDGELINE PLOT – STATES WITH SIGNIFICANT DIFFERENCES
# -------------------------------------------------------------------------

# Minimum p-value per state
pvalues_state <- table_significant_state %>%
  pivot_longer(cols = c(`Group A`, `Group B`), names_to = "Group", values_to = "State_Region") %>%
  group_by(State_Region) %>%
  summarise(
    p_min = min(as.numeric(gsub("< ", "", `p-value`)), na.rm = TRUE)
  ) %>%
  mutate(
    p_text = ifelse(p_min < 0.001, "< 0.001", formatC(p_min, format = "f", digits = 3))
  )

# Descriptive statistics for plotting
table_state_plot <- bd_state %>%
  group_by(State_Region) %>%
  summarise(
    Median = median(vigiGrade_medio, na.rm = TRUE),
    Q1 = quantile(vigiGrade_medio, 0.25, na.rm = TRUE),
    Q3 = quantile(vigiGrade_medio, 0.75, na.rm = TRUE)
  ) %>%
  left_join(pvalues_state, by = "State_Region") %>%
  filter(State_Region %in% unique(c(table_significant_state$`Group A`,
                                   table_significant_state$`Group B`)))

# Dataset for plotting
bd_state_plot <- bd_state %>%
  filter(State_Region %in% table_state_plot$State_Region)

# Ridgeline plot
state_ridgeline_plot <- ggplot(
  bd_state_plot,
  aes(x = vigiGrade_medio, y = fct_rev(State_Region))
) +
  geom_density_ridges(
    scale = 1.2,
    rel_min_height = 0.01,
    fill = "#69b3a2",
    alpha = 0.6,
    color = "black"
  ) +
  geom_point(
    data = table_state_plot,
    aes(x = Median, y = State_Region),
    inherit.aes = FALSE,
    size = 2
  ) +
  geom_errorbarh(
    data = table_state_plot,
    aes(xmin = Q1, xmax = Q3, y = State_Region),
    inherit.aes = FALSE,
    height = 0.3,
    size = 0.8
  ) +
  geom_text(
    data = table_state_plot,
    aes(x = Median, y = State_Region, label = round(Median, 3)),
    inherit.aes = FALSE,
    vjust = -1,
    size = 3
  ) +
  geom_text(
    data = table_state_plot,
    aes(x = 0.05, y = State_Region, label = paste0("p < ", p_text)),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3
  ) +
  labs(
    x = "vigiGrade Score",
    y = "State",
    title = "Distribution of vigiGrade Scores by State"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(state_ridgeline_plot)
# -------------------------------------------------------------------------
# STATE × REGION MAPPING
# -------------------------------------------------------------------------
regions_map <- c(
  "Norte" = "North",
  "Nordeste" = "Northeast",
  "Centro-Oeste" = "Midwest",
  "Sudeste" = "Southeast",
  "Sul" = "South"
)

state_region_table <- tibble::tibble(
  State_Region = c("AC","AP","AM","PA","RO","RR","TO",
                   "AL","BA","CE","MA","PB","PE","PI","RN","SE",
                   "DF","GO","MT","MS",
                   "ES","MG","RJ","SP",
                   "PR","RS","SC"),
  Region = c(rep("Norte",7),
             rep("Nordeste",9),
             rep("Centro-Oeste",4),
             rep("Sudeste",4),
             rep("Sul",3))
) %>%
  mutate(Region_English = regions_map[Region])

bd_state_plot <- bd_state_plot %>%
  left_join(state_region_table, by = "State_Region")

table_state_plot <- table_state_plot %>%
  left_join(state_region_table, by = "State_Region")

# Faceted ridgeline plot by region
state_region_ridgeline <- ggplot(
  bd_state_plot,
  aes(x = vigiGrade_medio, y = fct_rev(State_Region), fill = Region_English)
) +
  geom_density_ridges(scale = 1.2, rel_min_height = 0.01, alpha = 0.6) +
  facet_wrap(~Region_English, scales = "free_y", ncol = 2) +
  scale_fill_npg() +
  labs(
    x = "vigiGrade Score",
    y = "State",
    title = "Distribution of vigiGrade Scores by State and Region"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(state_region_ridgeline)
### ATC × vigiGrade ---------------------------------------------------------
# Unit of analysis: drug–event pairs by ICSR

# Assess normality of vigiGrade score
# Lilliefors test
lillie.test(bd_atc_final$vigiGrade_med_especifico)

# Shapiro–Wilk test (subsample due to large N)
shapiro.test(sample(bd_atc_final$vigiGrade_med_especifico, 5000))

# Histogram with density curve
ggplot(bd_atc_final, aes(x = vigiGrade_med_especifico)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightblue", color = "black") +
  geom_density(color = "red", size = 1) +
  labs(title = "Distribution of mean vigiGrade score by drug")

# Q–Q plot
qqnorm(bd_atc_final$vigiGrade_med_especifico)
qqline(bd_atc_final$vigiGrade_med_especifico, col = "red")

# Inspect number of observations per ATC level
bd_atc_final %>%
  group_by(ATC_Nivel1) |>
  count(ATC_Nivel1, sort = TRUE)

# Kruskal–Wallis test
kruskal.test(vigiGrade_med_especifico ~ ATC_Nivel1, data = bd_atc_final)

# Dunn post hoc test with Bonferroni correction
dunn_result_atc <- dunnTest(
  vigiGrade_med_especifico ~ ATC_Nivel1,
  data = bd_atc_final,
  method = "bonferroni"
)

# Extract pairwise comparisons
dunn_table_atc <- dunn_result_atc$res

# Descriptive statistics (median and IQR) by ATC level
tabela_tend_central_atc <- bd_atc_final %>%
  group_by(ATC_Nivel1) %>%
  summarise(
    Median = round(median(vigiGrade_med_especifico, na.rm = TRUE), 3),
    Q1 = round(quantile(vigiGrade_med_especifico, 0.25, na.rm = TRUE), 3),
    Q3 = round(quantile(vigiGrade_med_especifico, 0.75, na.rm = TRUE), 3)
  ) %>%
  mutate(IQR = paste0(Q1, " – ", Q3))

# Format Dunn results
tabela_dunn_atc <- dunn_table_atc %>%
  mutate(
    Group_A = str_trim(str_split_fixed(Comparison, " - ", 2)[, 1]),
    Group_B = str_trim(str_split_fixed(Comparison, " - ", 2)[, 2]),
    p_value = case_when(
      P.adj < 0.001 ~ "< 0.001",
      TRUE ~ formatC(P.adj, format = "f", digits = 3)
    ),
    Significant = ifelse(P.adj < 0.05, "Yes", "No")
  ) %>%
  select(Group_A, Group_B, Z, p_value, Significant)

# Effect size (Rosenthal’s r)
N_total_atc <- sum(!is.na(bd_atc_final$vigiGrade_med_especifico))

tabela_dunn_atc <- tabela_dunn_atc %>%
  mutate(
    Effect_size_r = round(abs(Z) / sqrt(N_total_atc), 3),
    Effect_interpretation = case_when(
      Effect_size_r < 0.1 ~ "Very small",
      Effect_size_r < 0.3 ~ "Small",
      Effect_size_r < 0.5 ~ "Moderate",
      TRUE ~ "Large"
    )
  )

# Merge medians and IQRs
final_atc <- tabela_dunn_atc %>%
  left_join(tabela_tend_central_atc, by = c("Group_A" = "ATC_Nivel1")) %>%
  rename(Median_A = Median, IQR_A = IQR) %>%
  left_join(tabela_tend_central_atc, by = c("Group_B" = "ATC_Nivel1")) %>%
  rename(Median_B = Median, IQR_B = IQR)

# Ridgeline plot for significant ATC classes
grafico_atc_vigigrade <- ggplot(
  bd_atc_plot,
  aes(x = vigiGrade_med_especifico, y = fct_rev(ATC_Nivel1))
) +
  geom_density_ridges(
    scale = 1.2,
    rel_min_height = 0.01,
    fill = "#69b3a2",
    alpha = 0.6,
    color = "black"
  ) +
  labs(
    x = "vigiGrade score",
    y = "ATC class (level 1)",
    title = "Distribution of vigiGrade scores by ATC class (level 1)"
  ) +
  theme_minimal()

### Year × vigiGrade ---------------------------------------------------------
# Create year variable from system entry date
bd_icsr_final <- bd_icsr_final %>%
  mutate(
    Year = as.factor(substr(as.character(System_Entry_Date), 1, 4))
  )

# Kruskal–Wallis test
kruskal.test(vigiGrade_medio ~ Year, data = bd_icsr_final)

# Dunn post hoc test
dunn_result_year <- dunnTest(
  vigiGrade_medio ~ Year,
  data = bd_icsr_final,
  method = "bonferroni"
)

# Descriptive statistics by year
tabela_tend_central_year <- bd_icsr_final %>%
  group_by(Year) %>%
  summarise(
    Median = round(median(vigiGrade_medio, na.rm = TRUE), 3),
    Q1 = round(quantile(vigiGrade_medio, 0.25, na.rm = TRUE), 3),
    Q3 = round(quantile(vigiGrade_medio, 0.75, na.rm = TRUE), 3)
  ) %>%
  mutate(IQR = paste0(Q1, " – ", Q3))

# Violin plot by year
grafico_violin_kwallis_ano <- ggplot(
  bd_icsr_final,
  aes(x = Year, y = vigiGrade_medio, fill = Year)
) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  scale_fill_npg() +
  labs(
    x = "Year",
    y = "vigiGrade score",
    title = "Comparison of vigiGrade scores across years"
  ) +
  theme_minimal()

### Jonckheere–Terpstra trend analysis --------------------------------------
# Define ordered year factor
bd_icsr_final <- bd_icsr_final %>%
  mutate(
    Year = factor(Year, levels = sort(unique(Year)), ordered = TRUE)
  )

# Jonckheere–Terpstra test for temporal trend
teste_tendencia_jonckheere <- jonckheere.test(
  bd_icsr_final$vigiGrade_medio,
  bd_icsr_final$Year
)

### Flowchart: data processing workflow ------------------------------------
library(DiagrammeR)

flowchart_plot <- grViz("
digraph data_flow {
  graph [rankdir = TB, fontsize = 16]
  node  [shape = box, fontsize = 12]
  edge  [style = dashed, color = gray60]

  subgraph cluster_sources {
    label = '1. Data sources'
    style = rounded

    src_drug  [label = 'Medicines']
    src_icsr  [label = 'Notifications']
    src_reac  [label = 'Reactions']
  }

  subgraph cluster_clean {
    label = '2. Cleaning and deduplication'
    style = rounded

    dedup_drug [label = 'Medicines (deduplicated)']
    dedup_icsr [label = 'Notifications (deduplicated)']
    dedup_reac [label = 'Reactions (deduplicated)']
  }

  subgraph cluster_link {
    label = '3. Record linkage'
    style = rounded

    linked [label = 'Final linked dataset']
  }

  src_drug -> dedup_drug
  src_icsr -> dedup_icsr
  src_reac -> dedup_reac

  dedup_drug -> linked
  dedup_icsr -> linked
  dedup_reac -> linked
}
")
