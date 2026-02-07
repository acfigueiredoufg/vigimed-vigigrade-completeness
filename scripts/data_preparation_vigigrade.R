### Data Preparation and Integration ###
library(tidyverse)
library(readr)
library(readxl)
library(summarytools)
library(psych)
library(janitor)
library(writexl)
# -------------------------------------------------------------------------
# STUDY CONFIGURATION
# -------------------------------------------------------------------------
DATA_DOWNLOAD <- "2026-01-02"
DATA_SOURCE   <- "ANVISA Open Data / Catalogued: 2021-09-30 / Updated: 2025-12-26"
LOG_FILE_NAME <- paste0("processing_log_", format(Sys.Date(), "%d_%m_%Y"), ".csv")

# -------------------------------------------------------------------------
# 1. Load the three core datasets
# -------------------------------------------------------------------------
medications   <- VigiMed_Medicamentos
reactions     <- VigiMed_Reacoes
notifications <- VigiMed_Notificacoes

# -------------------------------------------------------------------------
# 2. Select variables of interest
# -------------------------------------------------------------------------
notifications <- notifications %>% 
  select(UF, RECEBIDO_DE, IDENTIFICACAO_NOTIFICACAO, DATA_ULTIMA_ATUALIZACAO,
         DATA_INCLUSAO_SISTEMA, DATA_NASCIMENTO, IDADE_MOMENTO_REACAO,
         SEXO, NOTIFICADOR, TIPO_NOTIFICACAO)

reactions <- reactions %>% 
  select(IDENTIFICACAO_NOTIFICACAO, PT, HLT, HLGT, SOC, 
         DATA_INICIO_HORA, DATA_FINAL_HORA, DESFECHO)

medications <- medications %>% 
  select(IDENTIFICACAO_NOTIFICACAO, NOME_MEDICAMENTO_WHODRUG,
         PRINCIPIOS_ATIVOS_WHODRUG, CODIGO_ATC, DOSE,
         INDICACAO_RELATADA_NOTIFICADOR_INICIAL,
         PROBLEMAS_ADICIONAIS_RELCIONADOS_MEDICAMENTO,
         INICIO_ADMINISTRACAO, FIM_ADMINISTRACAO)

# Inspect system entry date range
range(notifications$DATA_INCLUSAO_SISTEMA, na.rm = TRUE)

# -------------------------------------------------------------------------
# DATA NORMALIZATION AND DEDUPLICATION
# -------------------------------------------------------------------------

# Audit snapshot function: counts rows across datasets at each stage
generate_snapshot <- function(stage_label) {
  tibble(
    Stage          = stage_label,
    Notifications  = nrow(notifications),
    Medications    = nrow(medications),
    Reactions      = nrow(reactions)
  )
}

# Stage 1: Raw data (baseline)
log_stage1 <- generate_snapshot("1. Raw data")

# -------------------------------------------------------------------------
# 3. String normalization
# Trims whitespace and converts explicit "None" strings to NA
# -------------------------------------------------------------------------
clean_characters <- function(df) {
  df %>% mutate(across(
    where(is.character),
    ~ str_trim(.) |> na_if("None")
  ))
}

notifications <- clean_characters(notifications)
medications   <- clean_characters(medications)
reactions     <- clean_characters(reactions)

log_stage2 <- bind_rows(log_stage1,
                        generate_snapshot("2. After normalization (trim + None → NA)"))

# -------------------------------------------------------------------------
# 4. Final deduplication
# -------------------------------------------------------------------------
notifications <- notifications %>% distinct()
medications   <- medications   %>% distinct()
reactions     <- reactions     %>% distinct()

log_stage3 <- bind_rows(log_stage2,
                        generate_snapshot("3. Final deduplication"))

# -------------------------------------------------------------------------
# 5. Exclusion of immunobiological products (vaccines, sera, immunoglobulins)
# -------------------------------------------------------------------------
medications <- medications %>%
  filter(!str_starts(CODIGO_ATC, "J06|J07|L04AA"))

log_stage4 <- bind_rows(log_stage3,
                        generate_snapshot("4. After ATC exclusion (J06, J07, L04AA)"))

# -------------------------------------------------------------------------
# FINAL AUDIT LOG
# -------------------------------------------------------------------------
final_log <- log_stage4 %>%
  mutate(
    Processing_Date   = Sys.Date(),
    Source_Download   = DATA_DOWNLOAD,
    Notes             = paste("Data source:", DATA_SOURCE)
  )

write_csv(final_log, LOG_FILE_NAME)
cat("Audit log saved as:", LOG_FILE_NAME, "\n")

# -------------------------------------------------------------------------
# DATASET INTEGRATION (JOIN OPERATIONS)
# -------------------------------------------------------------------------

# Identify report IDs present in all three datasets
common_ids <- intersect(notifications$IDENTIFICACAO_NOTIFICACAO,
                         intersect(medications$IDENTIFICACAO_NOTIFICACAO,
                                   reactions$IDENTIFICACAO_NOTIFICACAO))

# Estimate expected row inflation due to many-to-many joins
join_forecast <- tibble(IDENTIFICACAO_NOTIFICACAO = common_ids) %>%
  left_join(medications %>% count(IDENTIFICACAO_NOTIFICACAO, name = "n_med"),
            by = "IDENTIFICACAO_NOTIFICACAO") %>%
  left_join(reactions %>% count(IDENTIFICACAO_NOTIFICACAO, name = "n_reac"),
            by = "IDENTIFICACAO_NOTIFICACAO") %>%
  mutate(expected_rows = n_med * n_reac)

expected_total_rows <- sum(join_forecast$expected_rows)

cat("Number of complete report IDs:", length(common_ids), "\n")
cat("Expected number of rows after join:", expected_total_rows, "\n")

# Final integrated dataset
bd <- notifications %>%
  filter(IDENTIFICACAO_NOTIFICACAO %in% common_ids) %>%
  inner_join(medications, by = "IDENTIFICACAO_NOTIFICACAO") %>%
  inner_join(reactions,   by = "IDENTIFICACAO_NOTIFICACAO")

cat("Rows in integrated dataset:", nrow(bd), "\n")

# Check for fully duplicated rows
duplicated_rows <- nrow(bd) - nrow(distinct(bd))
cat("Fully duplicated rows in final dataset:", duplicated_rows, "\n")

# -------------------------------------------------------------------------
# DATE HANDLING AND TIME-TO-ONSET DERIVATION
# -------------------------------------------------------------------------
bd <- bd %>%
  mutate(
    med_date_raw  = str_extract(as.character(INICIO_ADMINISTRACAO), "\\d{8}"),
    reac_date_raw = str_extract(as.character(DATA_INICIO_HORA), "\\d{8}"),

    med_granularity = case_when(
      nchar(med_date_raw) == 8 ~ "day",
      nchar(med_date_raw) == 6 ~ "month",
      nchar(med_date_raw) == 4 ~ "year",
      TRUE                     ~ "invalid"
    ),
    reac_granularity = case_when(
      nchar(reac_date_raw) == 8 ~ "day",
      nchar(reac_date_raw) == 6 ~ "month",
      nchar(reac_date_raw) == 4 ~ "year",
      TRUE                      ~ "invalid"
    ),

    med_date = case_when(
      med_granularity == "day"   ~ ymd(med_date_raw),
      med_granularity == "month" ~ ymd(paste0(med_date_raw, "01")),
      med_granularity == "year"  ~ ymd(paste0(med_date_raw, "0101")),
      TRUE                       ~ as.Date(NA)
    ),
    reac_date = case_when(
      reac_granularity == "day"   ~ ymd(reac_date_raw),
      reac_granularity == "month" ~ ymd(paste0(reac_date_raw, "01")),
      reac_granularity == "year"  ~ ymd(paste0(reac_date_raw, "0101")),
      TRUE                        ~ as.Date(NA)
    )
  ) %>%
  mutate(
    diff_days = as.numeric(reac_date - med_date),

    # Time-to-onset penalty according to vigiGrade principles
    P_Time_to_Onset = case_when(
      !is.na(diff_days) & diff_days >= 60 & diff_days <= 120 ~ 1.00,
      !is.na(diff_days) & diff_days >= -30 & diff_days <= 30 ~ 0.90,
      !is.na(diff_days)                                     ~ 0.70,
      TRUE                                                   ~ 0.50
    )
  )

# -------------------------------------------------------------------------
# vigiGrade SCORE CALCULATION (PAIR LEVEL)
# -------------------------------------------------------------------------
bd_vigigrade <- bd %>%
  mutate(
    P_Indication = if_else(!is.na(INDICACAO_RELATADA_NOTIFICADOR_INICIAL), 1.0, 0.7),
    P_Outcome    = if_else(!is.na(DESFECHO), 1.0, 0.7),
    P_Sex        = if_else(!is.na(SEXO) & SEXO != "Desconhecido", 1.0, 0.7),
    P_Age        = if_else(!is.na(IDADE_MOMENTO_REACAO), 1.0, 0.7),
    P_Dose       = if_else(!is.na(DOSE), 1.0, 0.9),
    P_Region     = if_else(!is.na(UF), 1.0, 0.9),
    P_Reporter   = case_when(
      NOTIFICADOR %in% c("Médico", "Farmacêutico") ~ 1.0,
      !is.na(NOTIFICADOR)                          ~ 0.9,
      TRUE                                         ~ 0.9
    ),
    P_ReportType = if_else(!is.na(TIPO_NOTIFICACAO), 1.0, 0.9),
    P_Narrative  = if_else(!is.na(PROBLEMAS_ADICIONAIS_RELCIONADOS_MEDICAMENTO), 1.0, 0.9),

    score_pair = P_Time_to_Onset * P_Indication * P_Outcome * P_Sex * P_Age *
                 P_Dose * P_Region * P_Reporter * P_ReportType * P_Narrative
  )

write_csv(bd_vigigrade, "bd_vigigrade_pair_level.csv")

cat("✅ Data preparation and vigiGrade calculation completed successfully.\n")
