#----------------------------------------------------------------------------------
# PREVALENCE MODELS
#----------------------------------------------------------------------------------

# GOAL: Model the probability that a bird is infected by haemosporidian
# parasites (binary outcome), as a function of bird order (Passeriformes vs.
# non-Passeriformes) plus anthropogenic, environmental and host-diversity
# covariates, with a random intercept for community cluster.
#
# MODEL TYPE: binomial GLMM (logit link) via glmmTMB, fitted to
# aggregated counts (cbind(successes, failures)).

# NAVIGATION:
#   A. Load libraries
#   B. Read and clean raw data  -> build inf_bin, group, q1_aves
#   C. Aggregate to (community_cluster x group) level
#   D. Merge in environmental covariates (one row per community)
#   E. Standardize (z-scale) continuous predictors
#   F. Define shared model formula (main effects + interactions + RE)
#   G. Fit model, extract parameters, run diagnostics, export tables

#----------------------------------------------------------------------------------
# PART 1 — ALL GENERA PREVALENCE MODEL (Haemoproteus + Plasmodium + Leucocytozoon)
#----------------------------------------------------------------------------------

# A. LIBRARIES --------------------------------------------------------------
# dplyr/tidyr/stringr : data manipulation
# readxl/writexl     : Excel I/O
# glmmTMB            : mixed-effects model fitting
# parameters         : tidy coefficient tables
# performance        : model diagnostics (singularity, dispersion, R²)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(glmmTMB)
library(parameters)
library(performance)

# B. READ AND CLEAN RAW DATA -----------------------------------------------
# Read the raw bird-level dataset and coerce community_cluster to character
# so it is treated as a grouping factor (random effect), not numeric.
data_raw <- read_excel("dataset.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

# Build analysis variables at the individual-bird level:
#   inf_bin : 1 if Diagnosis == "P" (positive for any parasite), else 0.
#   group   : Passeriformes vs. non_Passeriformes (reference = Passeriformes).
#   q1_aves : host diversity (q1) assigned from the order-specific column,
#   so each bird inherits the diversity value for its own group.
data_raw <- data_raw %>%
  mutate(
    inf_bin = ifelse(Diagnosis == "P", 1, 0),   # All genera: binary infection
    group   = ifelse(BirdTree.Order == "Passeriformes",
                     "Passeriformes", "non_Passeriformes"),
    group   = factor(group, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      group == "non_Passeriformes" ~ q1_hosts_nonpass,
      group == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    )
  ) %>%
  mutate(inf_bin = replace_na(inf_bin, 0))   # treat missing diagnosis as negative

# C. AGGREGATE TO COMMUNITY x GROUP LEVEL -----------------------------------------
# For binomial GLMM we need successes (n_inf_all) and trials (total) per
# community_cluster x group cell. q1_aves is taken as the first non-missing
# value (it is constant within a community x group cell).
agregated <- data_raw %>%
  group_by(community_cluster, group) %>%
  summarise(
    total     = n(),
    n_inf_all = sum(inf_bin, na.rm = TRUE),
    q1_aves   = first(na.omit(q1_aves)),
    .groups = "drop"
  )

# D. EXTRACT ENVIRONMENTAL COVARIATES (ONE ROW PER COMMUNITY) -------------------
# Environmental variables are community-level, so we keep one distinct row
# per community_cluster and merge it back to the aggregated table.
ambientais <- data_raw %>%
  distinct(community_cluster, .keep_all = TRUE) %>%
  dplyr::select(community_cluster,
                HII_por_ano, NDVI_por_ano,
                bio01, bio04, bio12)

data_agregated <- agregated %>%
  left_join(ambientais, by = "community_cluster")

# Export the final modelling table for inspection / reproducibility.
write_xlsx(data_agregated, "agregated_all_genera.xlsx")

# E. STANDARDIZE CONTINUOUS PREDICTORS -----------------------------------------
# z-scaling puts all predictors on a common scale, which (a) improves
# numerical stability and (b) makes interaction coefficients interpretable
# per 1 SD change.
data <- data_agregated %>%
  mutate(
    NDVI.scaled    = as.numeric(scale(NDVI_por_ano)),
    HII.scaled     = as.numeric(scale(HII_por_ano)),
    bio12.scaled   = as.numeric(scale(bio12)),
    bio04.scaled   = as.numeric(scale(bio04)),
    bio01.scaled   = as.numeric(scale(bio01)),
    q1_aves.scaled = as.numeric(scale(q1_aves))
  )

# F. DEFINE SHARED MODEL FORMULA -----------------------------------------------
# common_terms holds the RHS used across parasite models: main effects of
# group + anthropogenic + environment + host diversity, quadratic terms for HII/bio12/bio01,
# and group x covariate interactions (allows effects to differ by order).
# (1 | community_cluster) adds a random intercept for clustering.
common_terms <- "
  group +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  bio12.scaled + I(bio12.scaled^2) +
  bio04.scaled +
  bio01.scaled + I(bio01.scaled^2) +
  q1_aves.scaled +
  group:NDVI.scaled +
  group:HII.scaled + group:I(HII.scaled^2) +
  group:bio12.scaled + group:I(bio12.scaled^2) +
  group:bio04.scaled +
  group:bio01.scaled + group:I(bio01.scaled^2) +
  group:q1_aves.scaled +
  (1 | community_cluster)
"

# Combine with the binomial response cbind(successes, failures).
formula_all <- as.formula(paste("cbind(n_inf_all, total - n_inf_all) ~",
                                common_terms))

# G. FIT MODEL, DIAGNOSE, EXPORT -----------------------------------------------
# Fit binomial GLMM (logit link is the glmmTMB default for binomial).
model_all3.2 <- glmmTMB(formula_all, data = data, family = binomial)

# Tidy coefficient table (Estimate, SE, CI, z, p).
params_all <- parameters::parameters(model_all3.2)

# --- Diagnostics ---
# Singularity : flags over-parameterized random effects (variance ~ 0).
# Overdispersion : checks whether binomial variance assumption holds.
# Convergence : flags non-converged fits.
cat("\n===== SINGULARITY =====\n")
print(check_singularity(model_all3.2))

cat("\n===== OVERDISPERSION =====\n")
print(check_overdispersion(model_all3.2))

cat("\n===== CONVERGENCY =====\n")
print(check_convergence(model_all3.2))

# Marginal R² (fixed effects only) and conditional R² (fixed + random).
cat("\n===== R² =====\n")
r2_all <- performance::r2(model_all3.2)
cat("R² marginal   :", round(r2_all$R2_marginal, 3), "\n")
cat("R² condicional:", round(r2_all$R2_conditional, 3), "\n")

# Export model data + parameter table.
write_xlsx(
  list(
    data_model = data,
    params   = as.data.frame(params_all)
  ),
  "prevalence_all_genera_params.xlsx"
)

# --- Build a coefficient table ---------------------------
# Re-select and rename columns, round, and add significance stars.
COEFS <- as.data.frame(parameters::parameters(model_all3.2))

table_COEFS <- COEFS %>%
  select(
    Termo    = Parameter,
    Estimate = Coefficient,
    SE       = SE,
    IC_low   = CI_low,
    IC_high  = CI_high,
    z        = z,
    p        = p
  ) %>%
  mutate(
    Estimate = round(Estimate, 3),
    SE       = round(SE, 3),
    IC_low   = round(IC_low, 3),
    IC_high  = round(IC_high, 3),
    z        = round(z, 2),
    Signif   = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE      ~ ""
    ),
    p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

write_xlsx(table_COEFS, "COEFS_model_all_genera.xlsx")

#---------------------------------------------------------------------------------
# PART 2 — PLASMODIUM PREVALENCE MODEL
#---------------------------------------------------------------------------------
# GOAL: Same structure as PART 1, but the response is infection by
#       Plasmodium specifically (identified by "PL" in Genus.for.all).
#
# KEY DIFFERENCES vs PART 1:
#   - Response variable PL_bin (regex match "PL" in Genus.for.all).
#--------------------------------------------------------------------------------

# --- A. LIBRARIES --------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(glmmTMB)
library(parameters)

# B. READ AND CLEAN RAW DATA ---------------------------------------------------
# PL_bin = 1 if Genus.for.all contains "PL" (Plasmodium), else 0.
# group, q1_aves built identically to PART 1 for comparability.
data_raw <- read_excel("dataset.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    PL_bin = ifelse(str_detect(Genus.for.all, "PL"), 1, 0),
    group  = ifelse(BirdTree.Order == "Passeriformes",
                    "Passeriformes", "non_Passeriformes"),
    group  = factor(group, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      group == "non_Passeriformes" ~ q1_hosts_nonpass,
      group == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    )
  ) %>%
  mutate(PL_bin = replace_na(PL_bin, 0))

# C. AGGREGATE (REUSABLE HELPER) ---------------------------------------------
# Function aggregates any df that already has PL_bin/group/q1_aves into
# community x group counts. Reused here; could be reused in PART 3.
agregar <- function(df) {
  df %>%
    group_by(community_cluster, group) %>%
    summarise(
      total    = n(),
      n_inf_PL = sum(PL_bin, na.rm = TRUE),
      q1_aves  = first(na.omit(q1_aves)),
      .groups = "drop"
    )
}

agregated <- agregar(data_raw)

# D. MERGE ENVIRONMENTAL + SPATIAL COVARIATES --------------------------------
# Adds coordinates in addition to climate/land-use variables.
ambientais <- data_raw %>%
  distinct(community_cluster, .keep_all = TRUE) %>%
  dplyr::select(community_cluster,
                HII_por_ano, NDVI_por_ano,
                bio01, bio04, bio12,
                Latitude_media, Longitude_media)

data_agregated <- agregated %>%
  left_join(ambientais, by = "community_cluster")

write_xlsx(data_agregated, "agregated_PL_comunidade_group.xlsx")

# E. STANDARDIZE -------------------------------------------------------------
data <- data_agregated %>%
  mutate(
    NDVI.scaled    = as.numeric(scale(NDVI_por_ano)),
    HII.scaled     = as.numeric(scale(HII_por_ano)),
    bio12.scaled   = as.numeric(scale(bio12)),
    bio04.scaled   = as.numeric(scale(bio04)),
    bio01.scaled   = as.numeric(scale(bio01)),
    q1_aves.scaled = as.numeric(scale(q1_aves))
  )

# F. FORMULA (same common_terms as PART 1 for comparability) ---------------
common_terms <- "
  group +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  bio12.scaled + I(bio12.scaled^2) +
  bio04.scaled +
  bio01.scaled + I(bio01.scaled^2) +
  q1_aves.scaled +
  group:NDVI.scaled +
  group:HII.scaled + group:I(HII.scaled^2) +
  group:bio12.scaled + group:I(bio12.scaled^2) +
  group:bio04.scaled +
  group:bio01.scaled + group:I(bio01.scaled^2) +
  group:q1_aves.scaled +
  (1 | community_cluster)
"

formula_pl <- as.formula(paste("cbind(n_inf_PL, total - n_inf_PL) ~", common_terms))

# G. FIT, DIAGNOSE, EXPORT ----------------------------------------------------------
model_pl3.2 <- glmmTMB(formula_pl, data = data, family = binomial)
params_pl <- parameters::parameters(model_pl3.2)

# AIC for model comparison across parasite models.
cat("\n===== AIC =====\n")
cat("Plasmodium:", AIC(model_pl3.2), "\n\n")

cat("\n===== SINGULARITY =====\n")
print(check_singularity(model_pl3.2))

cat("\n===== OVERDISPERSION =====\n")
print(check_overdispersion(model_pl3.2))

cat("\n===== CONVERGENCY =====\n")
print(check_convergence(model_pl3.2))

cat("\n===== R² =====\n")
r2_pl <- performance::r2(model_pl3.2)
cat("R² marginal   :", round(r2_pl$R2_marginal, 3), "\n")
cat("R² condicional:", round(r2_pl$R2_conditional, 3), "\n")

write_xlsx(
  list(
    data_model  = data,
    params_PL = as.data.frame(params_pl)
  ),
  "3.2-prevalence_PL_params.xlsx"
)

# -------- Coefficient table (same recipe as PART 1) ----------
COEFS <- as.data.frame(parameters::parameters(model_pl3.2))

table_COEFS <- COEFS %>%
  select(
    Termo    = Parameter,
    Estimate = Coefficient,
    SE       = SE,
    IC_low   = CI_low,
    IC_high  = CI_high,
    z        = z,
    p        = p
  ) %>%
  mutate(
    Estimate = round(Estimate, 3),
    SE       = round(SE, 3),
    IC_low   = round(IC_low, 3),
    IC_high  = round(IC_high, 3),
    z        = round(z, 2),
    Signif   = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE      ~ ""
    ),
    p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

write_xlsx(table_COEFS, "COEFS_model_PL.xlsx")

#-------------------------------------------------------------------------------
# PART 3 — PARAHAEMOPROTEUS PREVALENCE MODEL
#-------------------------------------------------------------------------------
# GOAL: Same structure as PART 2, but the response is infection by
#       Parahaemoproteus (identified by "PA" in Genus.for.all).
#-------------------------------------------------------------------------------

# A. LIBRARIES------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(glmmTMB)
library(parameters)
library(performance)

# B. READ AND CLEAN RAW DATA -------------------------------------------------
# PA_bin = 1 if Genus.for.all contains "PA" (Parahaemoproteus), else 0.
data_raw <- read_excel("dataset.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    PA_bin = ifelse(str_detect(Genus.for.all, "PA"), 1, 0),
    group  = ifelse(BirdTree.Order == "Passeriformes",
                    "Passeriformes", "non_Passeriformes"),
    group  = factor(group, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      group == "non_Passeriformes" ~ q1_hosts_nonpass,
      group == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    )
  ) %>%
  mutate(PA_bin = replace_na(PA_bin, 0))

# C. AGGREGATE ---------------------------------------------------------------
agregated <- data_raw %>%
  group_by(community_cluster, group) %>%
  summarise(
    total    = n(),
    n_inf_PA = sum(PA_bin, na.rm = TRUE),
    q1_aves  = first(na.omit(q1_aves)),
    .groups = "drop"
  )

# --- D. MERGE ENVIRONMENTAL + SPATIAL COVARIATES ---------------------------
ambientais <- data_raw %>%
  distinct(community_cluster, .keep_all = TRUE) %>%
  dplyr::select(community_cluster,
                HII_por_ano, NDVI_por_ano,
                bio01, bio04, bio12,
                Latitude_media, Longitude_media)

data_agregated <- agregated %>%
  left_join(ambientais, by = "community_cluster")

write_xlsx(data_agregated, "agregated_PA_comunidade_group.xlsx")

# E. STANDARDIZE -------------------------------------------------------------
data <- data_agregated %>%
  mutate(
    NDVI.scaled    = as.numeric(scale(NDVI_por_ano)),
    HII.scaled     = as.numeric(scale(HII_por_ano)),
    bio12.scaled   = as.numeric(scale(bio12)),
    bio04.scaled   = as.numeric(scale(bio04)),
    bio01.scaled   = as.numeric(scale(bio01)),
    q1_aves.scaled = as.numeric(scale(q1_aves))
  )

# F. FORMULA ----------------------------------------------------------------
common_terms <- "
  group +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  bio12.scaled + I(bio12.scaled^2) +
  bio04.scaled +
  bio01.scaled + I(bio01.scaled^2) +
  q1_aves.scaled +
  group:NDVI.scaled +
  group:HII.scaled + group:I(HII.scaled^2) +
  group:bio12.scaled + group:I(bio12.scaled^2) +
  group:bio04.scaled +
  group:bio01.scaled + group:I(bio01.scaled^2) +
  group:q1_aves.scaled +
  (1 | community_cluster)
"

formula_pa <- as.formula(paste("cbind(n_inf_PA, total - n_inf_PA) ~", common_terms))

# G. FIT, DIAGNOSE, EXPORT ----------------------------------------------------------
model_pa3.2 <- glmmTMB(formula_pa, data = data, family = binomial)
params_pa <- parameters::parameters(model_pa3.2)

cat("\n===== SINGULARITY =====\n")
print(check_singularity(model_pa3.2))

cat("\n===== OVERDISPERSION =====\n")
print(check_overdispersion(model_pa3.2))

cat("\n===== CONVERGENCY =====\n")
print(check_convergence(model_pa3.2))

cat("\n===== R² =====\n")
r2_pa <- performance::r2(model_pa3.2)
cat("R² marginal   :", round(r2_pa$R2_marginal, 3), "\n")
cat("R² condicional:", round(r2_pa$R2_conditional, 3), "\n")

write_xlsx(
  list(
    data_model  = data,
    params_PA = as.data.frame(params_pa)
  ),
  "3.2-prevalence_PA_params.xlsx"
)

# Coefficient table ----------------------------------------------------
COEFS <- as.data.frame(parameters::parameters(model_pa3.2))

table_COEFS <- COEFS %>%
  select(
    Termo    = Parameter,
    Estimate = Coefficient,
    SE       = SE,
    IC_low   = CI_low,
    IC_high  = CI_high,
    z        = z,
    p        = p
  ) %>%
  mutate(
    Estimate = round(Estimate, 3),
    SE       = round(SE, 3),
    IC_low   = round(IC_low, 3),
    IC_high  = round(IC_high, 3),
    z        = round(z, 2),
    Signif   = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE      ~ ""
    ),
    p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

write_xlsx(table_COEFS, "COEFS_model_PA.xlsx")

#----------------------------------------------------------------------------------
# DIVERSITY MODELS
#----------------------------------------------------------------------------------

# GOAL: Model parasite DIVERSITY (q1_parasite, a continuous diversity index)
# as a function of bird order, anthropogenic, environmental, and host diversity,
# with group interactions.
#
# MODEL TYPE: Gaussian GLM with log link (ensures positive predictions),
# fitted at community x group level. NO random effect here —
# unlike prevalence models, this is a plain glm().

# NAVIGATION:
#   A. Libraries
#   B. Clean raw data -> build q1_aves and q1_parasite (order-specific)
#   C. Aggregate to community x group (first non-NA per cell)
#   D. Drop missing / zero-diversity cells; standardize predictors
#   E. Fit Gaussian log-link GLM
#   F. Compute deviance R² and adjusted R²; DHARMa residual diagnostics
#   G. Export parameters + model data

#---------------------------------------------------------------------------------
# PART 4 — ALL GENERA DIVERSITY MODEL (Haemoproteus + Plasmodium + Leucocytozoon)
#---------------------------------------------------------------------------------

# A. LIBRARIES -------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(parameters)
library(DHARMa)      # simulation-based residual diagnostics
library(performance)

# B. READ AND CLEAN RAW DATA ---------------------------------------------------
# q1_parasite is built from order-specific columns:
#   non_Passeriformes -> q1_all_genera_nopass
#   Passeriformes     -> q1_all_genera_pass
# This mirrors the q1_aves construction pattern.
data_raw <- read_excel("dataset.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    group = ifelse(BirdTree.Order == "Passeriformes",
                   "Passeriformes", "non_Passeriformes"),
    group = factor(group, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      group == "non_Passeriformes" ~ q1_hosts_nonpass,
      group == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    ),
    q1_parasite = case_when(
      group == "non_Passeriformes" ~ q1_all_genera_nopass,
      group == "Passeriformes"     ~ q1_all_genera_pass,
      TRUE                         ~ NA_real_
    )
  )

# C. AGGREGATE TO COMMUNITY x GROUP ---------------------------------------------
# Diversity and environment are community-level, so we take the first
# non-missing value per cell (equivalent to distinct() per group).
data_agregated <- data_raw %>%
  group_by(community_cluster, group) %>%
  summarise(
    q1_aves      = first(na.omit(q1_aves)),
    q1_parasite  = first(na.omit(q1_parasite)),
    NDVI_por_ano = first(NDVI_por_ano),
    HII_por_ano  = first(HII_por_ano),
    bio01        = first(bio01),
    bio04        = first(bio04),
    bio12        = first(bio12),
    .groups = "drop"
  )

# D. FILTER AND STANDARDIZE -----------------------------------------------------
# Keep only cells with both diversity measures and positive parasite
# diversity (log link requires q1_parasite > 0 to avoid log(0)).
data <- data_agregated %>%
  filter(!is.na(q1_parasite), !is.na(q1_aves)) %>%
  filter(q1_parasite > 0) %>%
  mutate(
    NDVI.scaled    = as.numeric(scale(NDVI_por_ano)),
    HII.scaled     = as.numeric(scale(HII_por_ano)),
    bio01.scaled   = as.numeric(scale(bio01)),
    bio04.scaled   = as.numeric(scale(bio04)),
    bio12.scaled   = as.numeric(scale(bio12)),
    q1_aves.scaled = as.numeric(scale(q1_aves))
  )

# E. FORMULA AND MODEL ---------------------------------------------------------
# Simpler formula than prevalence models: no bio01/bio04/bio12 here —
# only NDVI, HII (with quadratic), and q1_aves, all interacting with group.
formula_all_div <- q1_parasite ~
  group +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  q1_aves.scaled +
  group:NDVI.scaled +
  group:HII.scaled + group:I(HII.scaled^2) +
  group:q1_aves.scaled

# Gaussian GLM with log link (positive continuous response).
model_all_div <- glm(
  formula_all_div,
  data = data,
  family = gaussian(link = "log")
)

params_all_div <- parameters::parameters(model_all_div)

# F. GOODNESS-OF-FIT AND DIAGNOSTICS -------------------------------------------
# Deviance-based R² (McFadden-style for GLM) and an adjusted version
# correcting for number of parameters k and sample size n.
r2_dev <- with(summary(model_all_div), 1 - deviance/null.deviance)
n <- nobs(model_all_div)
k <- length(coef(model_all_div))
r2_adj <- 1 - (1 - r2_dev) * (n - 1) / (n - k - 1)

# DHARMa: simulate residuals (n = 1000) and check uniformity + dispersion.
sim_all <- simulateResiduals(model_all_div, n = 1000)
plot(sim_all)
print(testDispersion(sim_all))

# G. EXPORT --------------------------------------------------------------------
write_xlsx(
  list(
    params = as.data.frame(params_all_div),
    data      = data
  ),
  "model_diversity_all_genera.xlsx"
)

#------------------------------------------------------------------------------
# PART 5 — PLASMODIUM DIVERSITY MODEL
#------------------------------------------------------------------------------
# GOAL: Same as PART 4 but for Plasmodium diversity (q1_plasmodium_*).
# STRUCTURE: Identical workflow — only the q1_parasite source columns and
# output filenames change.
#------------------------------------------------------------------------------

# A. LIBRARIES ----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(parameters)
library(DHARMa)
library(performance)

# B. READ AND CLEAN ----------------------------------------------------------
data_raw <- read_excel("dataset.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    group = ifelse(BirdTree.Order == "Passeriformes",
                   "Passeriformes", "non_Passeriformes"),
    group = factor(group, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      group == "non_Passeriformes" ~ q1_hosts_nonpass,
      group == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    ),
    q1_parasite = case_when(
      group == "non_Passeriformes" ~ q1_plasmodium_nopass,
      group == "Passeriformes"     ~ q1_plasmodium_pass,
      TRUE                         ~ NA_real_
    )
  )

# C. AGGREGATE ---------------------------------------------------------------
data_agregated <- data_raw %>%
  group_by(community_cluster, group) %>%
  summarise(
    q1_aves      = first(na.omit(q1_aves)),
    q1_parasite  = first(na.omit(q1_parasite)),
    NDVI_por_ano = first(NDVI_por_ano),
    HII_por_ano  = first(HII_por_ano),
    bio01        = first(bio01),
    bio04        = first(bio04),
    bio12        = first(bio12),
    .groups = "drop"
  )

# D. FILTER AND STANDARDIZE -------------------------------------------------
data <- data_agregated %>%
  filter(!is.na(q1_parasite), !is.na(q1_aves)) %>%
  filter(q1_parasite > 0) %>%
  mutate(
    NDVI.scaled    = as.numeric(scale(NDVI_por_ano)),
    HII.scaled     = as.numeric(scale(HII_por_ano)),
    bio01.scaled   = as.numeric(scale(bio01)),
    bio04.scaled   = as.numeric(scale(bio04)),
    bio12.scaled   = as.numeric(scale(bio12)),
    q1_aves.scaled = as.numeric(scale(q1_aves))
  )

# E. FORMULA AND MODEL ------------------------------------------------------
formula_pl_div <- q1_parasite ~
  group +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  q1_aves.scaled +
  group:NDVI.scaled +
  group:HII.scaled + group:I(HII.scaled^2) +
  group:q1_aves.scaled

model_pl_div <- glm(
  formula_pl_div,
  data = data,
  family = gaussian(link = "log")
)

params_pl_div <- parameters::parameters(model_pl_div)

# F. GOODNESS-OF-FIT AND DIAGNOSTICS ------------------------------------------
r2_dev <- with(summary(model_pl_div), 1 - deviance/null.deviance)
n <- nobs(model_pl_div)
k <- length(coef(model_pl_div))
r2_adj <- 1 - (1 - r2_dev) * (n - 1) / (n - k - 1)

sim_pl <- simulateResiduals(model_pl_div, n = 1000)
plot(sim_pl)
print(testDispersion(sim_pl))

# G. EXPORT -------------------------------------------------------------------
write_xlsx(
  list(
    params = as.data.frame(params_pl_div),
    data      = data
  ),
  "model_diversity_PL.xlsx"
)

#------------------------------------------------------------------------------
# PART 6 — PARAHAEMOPROTEUS DIVERSITY MODEL
#------------------------------------------------------------------------------
# GOAL: Same as PARTS 4–5 but for Parahaemoproteus diversity
#       (q1_parahaemoproteus_*).
#------------------------------------------------------------------------------

# A. LIBRARIES ----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(parameters)
library(DHARMa)
library(performance)

# B. READ AND CLEAN ------------------------------------------------------------
data_raw <- read_excel("dataset.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    group = ifelse(BirdTree.Order == "Passeriformes",
                   "Passeriformes", "non_Passeriformes"),
    group = factor(group, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      group == "non_Passeriformes" ~ q1_hosts_nonpass,
      group == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    ),
    q1_parasite = case_when(
      group == "non_Passeriformes" ~ q1_parahaemoproteus_nopass,
      group == "Passeriformes"     ~ q1_parahaemoproteus_pass,
      TRUE                         ~ NA_real_
    )
  )

# C. AGGREGATE ---------------------------------------------------------------
data_agregated <- data_raw %>%
  group_by(community_cluster, group) %>%
  summarise(
    q1_aves      = first(na.omit(q1_aves)),
    q1_parasite  = first(na.omit(q1_parasite)),
    NDVI_por_ano = first(NDVI_por_ano),
    HII_por_ano  = first(HII_por_ano),
    bio01        = first(bio01),
    bio04        = first(bio04),
    bio12        = first(bio12),
    .groups = "drop"
  )

# D. FILTER AND STANDARDIZE ----------------------------------------------------
data <- data_agregated %>%
  filter(!is.na(q1_parasite), !is.na(q1_aves)) %>%
  filter(q1_parasite > 0) %>%
  mutate(
    NDVI.scaled    = as.numeric(scale(NDVI_por_ano)),
    HII.scaled     = as.numeric(scale(HII_por_ano)),
    bio01.scaled   = as.numeric(scale(bio01)),
    bio04.scaled   = as.numeric(scale(bio04)),
    bio12.scaled   = as.numeric(scale(bio12)),
    q1_aves.scaled = as.numeric(scale(q1_aves))
  )

# E. FORMULA AND MODEL ---------------------------------------------------------
formula_pa_div <- q1_parasite ~
  group +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  q1_aves.scaled +
  group:NDVI.scaled +
  group:HII.scaled + group:I(HII.scaled^2) +
  group:q1_aves.scaled

model_pa_div <- glm(
  formula_pa_div,
  data = data,
  family = gaussian(link = "log")
)

params_pa_div <- parameters::parameters(model_pa_div)

# F. GOODNESS-OF-FIT AND DIAGNOSTICS ------------------------------------------
r2_dev <- with(summary(model_pa_div), 1 - deviance/null.deviance)
n <- nobs(model_pa_div)
k <- length(coef(model_pa_div))
r2_adj <- 1 - (1 - r2_dev) * (n - 1) / (n - k - 1)

# Explicit reporting block (see header: extra vs PARTS 4–5).
cat("\n===== R² =====\n")
cat("R² deviance :", round(r2_dev, 3), "\n")
cat("R² ajusted :", round(r2_adj, 3), "\n")
cat("N           :", n, "\n")
cat("k           :", k, "\n")
cat("ratio N/k   :", round(n/k, 2), "\n")

# G. EXPORT ------------------------------------------------------------------
write_xlsx(
  list(
    params = as.data.frame(params_pa_div),
    data      = data
  ),
  "model_diversity_PA.xlsx"
)
