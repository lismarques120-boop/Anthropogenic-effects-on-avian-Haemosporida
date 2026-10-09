#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###ALL GENERA PREVALENCE MODEL####################################################

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(glmmTMB)
library(parameters)
library(performance)

data_raw <- read_excel("Supporting information - Appendix S1.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    inf_bin = ifelse(Diagnosis == "P", 1, 0),   # All genera
    grupo   = ifelse(BirdTree.Order == "Passeriformes",
                     "Passeriformes", "non_Passeriformes"),
    grupo   = factor(grupo, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      grupo == "non_Passeriformes" ~ q1_hosts_nonpass,
      grupo == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    )
  ) %>%
  mutate(inf_bin = replace_na(inf_bin, 0))

agregated <- data_raw %>%
  group_by(community_cluster, grupo) %>%
  summarise(
    total     = n(),
    n_inf_all = sum(inf_bin, na.rm = TRUE),
    q1_aves   = first(na.omit(q1_aves)),
    .groups = "drop"
  )


ambientais <- data_raw %>%
  distinct(community_cluster, .keep_all = TRUE) %>%
  dplyr::select(community_cluster,
                HII_por_ano, NDVI_por_ano,
                bio01, bio04, bio12)


data_agregated <- agregated %>%
  left_join(ambientais, by = "community_cluster")

write_xlsx(data_agregated, "agregated_all_genera.xlsx")


data <- data_agregated %>%
  mutate(
    NDVI.scaled    = as.numeric(scale(NDVI_por_ano)),
    HII.scaled     = as.numeric(scale(HII_por_ano)),
    bio12.scaled   = as.numeric(scale(bio12)),
    bio04.scaled   = as.numeric(scale(bio04)),
    bio01.scaled   = as.numeric(scale(bio01)),
    q1_aves.scaled = as.numeric(scale(q1_aves))
  )


common_terms <- "
  grupo +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  bio12.scaled + I(bio12.scaled^2) +
  bio04.scaled +
  bio01.scaled + I(bio01.scaled^2) +
  q1_aves.scaled +
  grupo:NDVI.scaled +
  grupo:HII.scaled + grupo:I(HII.scaled^2) +
  grupo:bio12.scaled + grupo:I(bio12.scaled^2) +
  grupo:bio04.scaled +
  grupo:bio01.scaled + grupo:I(bio01.scaled^2) +
  grupo:q1_aves.scaled +
  (1 | community_cluster)
"

formula_all <- as.formula(paste("cbind(n_inf_all, total - n_inf_all) ~",
                                common_terms))

model_all3.2 <- glmmTMB(formula_all, data = data, family = binomial)

params_all <- parameters::parameters(model_all3.2)

cat("\n===== SINGULARITY =====\n")
print(check_singularity(model_all3.2))

cat("\n===== OVERDISPERSION =====\n")
print(check_overdispersion(model_all3.2))

cat("\n===== CONVERGENCY =====\n")
print(check_convergence(model_all3.2))

cat("\n===== R² =====\n")
r2_all <- performance::r2(model_all3.2)
cat("R² marginal   :", round(r2_all$R2_marginal, 3), "\n")
cat("R² condicional:", round(r2_all$R2_conditional, 3), "\n")

write_xlsx(
  list(
    data_model = data,
    params   = as.data.frame(params_all)
  ),
  "prevalence_all_genera_params.xlsx"
)

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###PLASMODIUM PREVALENCE MODEL####################################################

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(glmmTMB)
library(parameters)

data_raw <- read_excel("Supporting information - Appendix S1.xlsx")

data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    PL_bin = ifelse(str_detect(Genus.for.all, "PL"), 1, 0),
    grupo  = ifelse(BirdTree.Order == "Passeriformes",
                    "Passeriformes", "non_Passeriformes"),
    grupo  = factor(grupo, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      grupo == "non_Passeriformes" ~ q1_hosts_nonpass,
      grupo == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    )
  ) %>%
  mutate(PL_bin = replace_na(PL_bin, 0))


agregar <- function(df) {
  df %>%
    group_by(community_cluster, grupo) %>%
    summarise(
      total    = n(),
      n_inf_PL = sum(PL_bin, na.rm = TRUE),
      q1_aves  = first(na.omit(q1_aves)),
      .groups = "drop"
    )
}

agregated <- agregar(data_raw)

ambientais <- data_raw %>%
  distinct(community_cluster, .keep_all = TRUE) %>%
  dplyr::select(community_cluster,
                HII_por_ano, NDVI_por_ano,
                bio01, bio04, bio12,
                Latitude_media, Longitude_media)


data_agregated <- agregated %>%
  left_join(ambientais, by = "community_cluster")

write_xlsx(data_agregated, "agregated_PL_comunidade_grupo.xlsx")

data <- data_agregated %>%
  mutate(
    NDVI.scaled    = as.numeric(scale(NDVI_por_ano)),
    HII.scaled     = as.numeric(scale(HII_por_ano)),
    bio12.scaled   = as.numeric(scale(bio12)),
    bio04.scaled   = as.numeric(scale(bio04)),
    bio01.scaled   = as.numeric(scale(bio01)),
    q1_aves.scaled = as.numeric(scale(q1_aves))
  )


common_terms <- "
  grupo +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  bio12.scaled + I(bio12.scaled^2) +
  bio04.scaled +
  bio01.scaled + I(bio01.scaled^2) +
  q1_aves.scaled +
  grupo:NDVI.scaled +
  grupo:HII.scaled + grupo:I(HII.scaled^2) +
  grupo:bio12.scaled + grupo:I(bio12.scaled^2) +
  grupo:bio04.scaled +
  grupo:bio01.scaled + grupo:I(bio01.scaled^2) +
  grupo:q1_aves.scaled +
  (1 | community_cluster)
"

formula_pl <- as.formula(paste("cbind(n_inf_PL, total - n_inf_PL) ~", common_terms))

model_pl3.2 <- glmmTMB(formula_pl, data = data, family = binomial)

params_pl <- parameters::parameters(model_pl3.2)

cat("\n===== AIC =====\n")
cat("Plasmodium:", AIC(model_pl3.2), "\n\n")

cat("===== PARÂMETROS =====\n")
print(params_pl)

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###PARAHAEMOPROTEUS PREVALENCE MODEL####################################################

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(glmmTMB)
library(parameters)
library(performance)
library(ggplot2)
library(ggeffects)
library(patchwork)
library(ggtext)
library(scales)

data_raw <- read_excel("Supporting information - Appendix S1.xlsx")

data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    PA_bin = ifelse(str_detect(Genus.for.all, "PA"), 1, 0),
    grupo  = ifelse(BirdTree.Order == "Passeriformes",
                    "Passeriformes", "non_Passeriformes"),
    grupo  = factor(grupo, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      grupo == "non_Passeriformes" ~ q1_hosts_nonpass,
      grupo == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    )
  ) %>%
  mutate(PA_bin = replace_na(PA_bin, 0))

agregated <- data_raw %>%
  group_by(community_cluster, grupo) %>%
  summarise(
    total    = n(),
    n_inf_PA = sum(PA_bin, na.rm = TRUE),
    q1_aves  = first(na.omit(q1_aves)),
    .groups = "drop"
  )


ambientais <- data_raw %>%
  distinct(community_cluster, .keep_all = TRUE) %>%
  dplyr::select(community_cluster,
                HII_por_ano, NDVI_por_ano,
                bio01, bio04, bio12,
                Latitude_media, Longitude_media)


data_agregated <- agregated %>%
  left_join(ambientais, by = "community_cluster")

write_xlsx(data_agregated, "agregated_PA_comunidade_grupo.xlsx")


data <- data_agregated %>%
  mutate(
    NDVI.scaled    = as.numeric(scale(NDVI_por_ano)),
    HII.scaled     = as.numeric(scale(HII_por_ano)),
    bio12.scaled   = as.numeric(scale(bio12)),
    bio04.scaled   = as.numeric(scale(bio04)),
    bio01.scaled   = as.numeric(scale(bio01)),
    q1_aves.scaled = as.numeric(scale(q1_aves))
  )


common_terms <- "
  grupo +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  bio12.scaled + I(bio12.scaled^2) +
  bio04.scaled +
  bio01.scaled + I(bio01.scaled^2) +
  q1_aves.scaled +
  grupo:NDVI.scaled +
  grupo:HII.scaled + grupo:I(HII.scaled^2) +
  grupo:bio12.scaled + grupo:I(bio12.scaled^2) +
  grupo:bio04.scaled +
  grupo:bio01.scaled + grupo:I(bio01.scaled^2) +
  grupo:q1_aves.scaled +
  (1 | community_cluster)
"

formula_pa <- as.formula(paste("cbind(n_inf_PA, total - n_inf_PA) ~", common_terms))

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###ALL GENERA DIVERSITY MODEL####################################################

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(parameters)
library(DHARMa)
library(performance)

data_raw <- read_excel("Supporting information - Appendix S1.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    grupo = ifelse(BirdTree.Order == "Passeriformes",
                   "Passeriformes", "non_Passeriformes"),
    grupo = factor(grupo, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      grupo == "non_Passeriformes" ~ q1_hosts_nonpass,
      grupo == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    ),
    q1_parasite = case_when(
      grupo == "non_Passeriformes" ~ q1_all_genera_nopass,
      grupo == "Passeriformes"     ~ q1_all_genera_pass,
      TRUE                         ~ NA_real_
    )
  )

data_agregated <- data_raw %>%
  group_by(community_cluster, grupo) %>%
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

formula_all_div <- q1_parasite ~
  grupo +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  q1_aves.scaled +
  grupo:NDVI.scaled +
  grupo:HII.scaled + grupo:I(HII.scaled^2) +
  grupo:q1_aves.scaled


model_all_div <- glm(
  formula_all_div,
  data = data,
  family = gaussian(link = "log")
)

params_all_div <- parameters::parameters(model_all_div)


r2_dev <- with(summary(model_all_div), 1 - deviance/null.deviance)
n <- nobs(model_all_div)
k <- length(coef(model_all_div))
r2_adj <- 1 - (1 - r2_dev) * (n - 1) / (n - k - 1)

sim_all <- simulateResiduals(model_all_div, n = 1000)
plot(sim_all)
print(testDispersion(sim_all))

write_xlsx(
  list(
    params = as.data.frame(params_all_div),
    data      = data
  ),
  "model_diversity_all_genera.xlsx"
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###PLASMODIUM DIVERSITY MODEL####################################################

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(parameters)
library(DHARMa)
library(performance)

data_raw <- read_excel("Supporting information - Appendix S1.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    grupo = ifelse(BirdTree.Order == "Passeriformes",
                   "Passeriformes", "non_Passeriformes"),
    grupo = factor(grupo, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      grupo == "non_Passeriformes" ~ q1_hosts_nonpass,
      grupo == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    ),
    q1_parasite = case_when(
      grupo == "non_Passeriformes" ~ q1_plasmodium_nopass,
      grupo == "Passeriformes"     ~ q1_plasmodium_pass,
      TRUE                         ~ NA_real_
    )
  )

data_agregated <- data_raw %>%
  group_by(community_cluster, grupo) %>%
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

formula_pl_div <- q1_parasite ~
  grupo +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  q1_aves.scaled +
  grupo:NDVI.scaled +
  grupo:HII.scaled + grupo:I(HII.scaled^2) +
  grupo:q1_aves.scaled

model_pl_div <- glm(
  formula_pl_div,
  data = data,
  family = gaussian(link = "log")
)

params_pl_div <- parameters::parameters(model_pl_div)

r2_dev <- with(summary(model_pl_div), 1 - deviance/null.deviance)
n <- nobs(model_pl_div)
k <- length(coef(model_pl_div))
r2_adj <- 1 - (1 - r2_dev) * (n - 1) / (n - k - 1)

sim_pl <- simulateResiduals(model_pl_div, n = 1000)
plot(sim_pl)
print(testDispersion(sim_pl))

write_xlsx(
  list(
    params = as.data.frame(params_pl_div),
    data      = data
  ),
  "model_diversity_PL.xlsx"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###PARAHAEMOPROTEUS DIVERSITY MODEL####################################################

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(parameters)
library(DHARMa)
library(performance)

data_raw <- read_excel("Supporting information - Appendix S1.xlsx")
data_raw$community_cluster <- as.character(data_raw$community_cluster)

data_raw <- data_raw %>%
  mutate(
    grupo = ifelse(BirdTree.Order == "Passeriformes",
                   "Passeriformes", "non_Passeriformes"),
    grupo = factor(grupo, levels = c("Passeriformes", "non_Passeriformes")),
    q1_aves = case_when(
      grupo == "non_Passeriformes" ~ q1_hosts_nonpass,
      grupo == "Passeriformes"     ~ q1_hosts_pass,
      TRUE                         ~ NA_real_
    ),
    q1_parasite = case_when(
      grupo == "non_Passeriformes" ~ q1_parahaemoproteus_nopass,
      grupo == "Passeriformes"     ~ q1_parahaemoproteus_pass,
      TRUE                         ~ NA_real_
    )
  )

data_agregated <- data_raw %>%
  group_by(community_cluster, grupo) %>%
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

formula_pa_div <- q1_parasite ~
  grupo +
  NDVI.scaled +
  HII.scaled + I(HII.scaled^2) +
  q1_aves.scaled +
  grupo:NDVI.scaled +
  grupo:HII.scaled + grupo:I(HII.scaled^2) +
  grupo:q1_aves.scaled

model_pa_div <- glm(
  formula_pa_div,
  data = data,
  family = gaussian(link = "log")
)

params_pa_div <- parameters::parameters(model_pa_div)

r2_dev <- with(summary(model_pa_div), 1 - deviance/null.deviance)
n <- nobs(model_pa_div)
k <- length(coef(model_pa_div))
r2_adj <- 1 - (1 - r2_dev) * (n - 1) / (n - k - 1)

cat("\n===== R² =====\n")
cat("R² deviance :", round(r2_dev, 3), "\n")
cat("R² ajusted :", round(r2_adj, 3), "\n")
cat("N           :", n, "\n")
cat("k           :", k, "\n")
cat("ratio N/k   :", round(n/k, 2), "\n")

sim_pa <- simulateResiduals(model_pa_div, n = 1000)
plot(sim_pa)
print(testDispersion(sim_pa))

write_xlsx(
  list(
    params = as.data.frame(params_pa_div),
    data      = data
  ),
  "model_diversity_PA.xlsx"
)
