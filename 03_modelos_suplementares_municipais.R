# PUBLICATION REPRODUCIBILITY — BREAST CANCER NEOADJUVANT CHEMOTHERAPY
# Source supplementary S2 municipal models.
# Executes locally on authorized DATASUS data; never commit patient-level data.
# 2026-10-15. Historical crosswalks 2024/2026 are not validated for 2023.
# All R scripts should be run from project root; MAMA_RAIZ configures data root.
# NOTE: full end-to-end re-execution was not possible in the delivery environment.
# -------------------------------------------------------------------------

# ============================================================================
# SCRIPT SUPLEMENTAR S2 — ANÁLISES ESTATÍSTICAS
# Projeto: desigualdades raciais e territoriais no início da quimioterapia
#          neoadjuvante para câncer de mama no SUS
# Coorte principal: Brasil, 2023
#
# Este script NÃO repete extração/mineração do DATASUS.
# Ele parte do arquivo produzido pelo script principal:
#   resultados/coorte_episodios_neoadjuvantes_2023_pre_analise.fst
#
# Reproduz:
#   S1) descrição da coorte e raça x deslocamento;
#   S2) modelos Poisson M0-M3;
#   S3) inferência CR2/Satterthwaite e teste global HTZ;
#   S4) riscos ajustados e contrastes raciais/territoriais;
#   S5) sensibilidade abril-dezembro/2023;
#   S6) sensibilidade sem teto de 365 dias;
#   S7) sensibilidade usando raça/cor da APAC índice;
#   S8) leave-one-CNES-out (opcional).
# ============================================================================

# 0. Pacotes -----------------------------------------------------------------
pacotes <- c(
  "dplyr", "readr", "fst", "splines", "sandwich", "emmeans",
  "clubSandwich", "purrr", "tibble"
)

faltantes <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltantes) > 0) {
  stop("Instale os pacotes: ", paste(faltantes, collapse = ", "))
}

library(dplyr)
library(readr)
library(fst)
library(splines)
library(sandwich)
library(emmeans)
library(clubSandwich)
library(purrr)
library(tibble)

# 1. Diretórios ---------------------------------------------------------------
PASTA_BASE <- file.path(Sys.getenv("MAMA_RAIZ",unset=getwd()), "dados_mama_SIA_AQ_2019_2025")
DIR_RESULTADOS <- file.path(PASTA_BASE, "resultados")
DIR_SUPL <- file.path(DIR_RESULTADOS, "material_suplementar")
dir.create(DIR_SUPL, recursive = TRUE, showWarnings = FALSE)

ARQ_COORTE <- file.path(
  DIR_RESULTADOS,
  "coorte_episodios_neoadjuvantes_2023_pre_analise.fst"
)

if (!file.exists(ARQ_COORTE)) stop("Arquivo não encontrado: ", ARQ_COORTE)

# TRUE = executa exclusão sucessiva de cada CNES; pode ser demorado no Brasil.
RODAR_LOO <- FALSE

# 2. Leitura e padronização ---------------------------------------------------
coorte_2023 <- read_fst(ARQ_COORTE, as.data.table = FALSE) %>%
  as_tibble() %>%
  mutate(
    raca_final = factor(raca_final,
      levels = c("Branca", "Parda", "Preta", "Amarela", "Indígena",
                 "Inconsistente", "Sem informação")
    ),
    fora_municipio = factor(fora_municipio,
      levels = c("Mesmo município", "Outro município")
    ),
    estadio = as.character(estadio),
    procedimento = as.character(procedimento),
    uf_res = as.character(uf_res),
    mes_inicio = as.integer(mes_inicio),
    atraso60 = as.integer(tempo_dias > 60)
  )

# 3. Fluxo analítico resumido -------------------------------------------------
fluxo_suplementar <- tibble(
  etapa = c(
    "Episódios índice de 2023 antes dos filtros inferenciais",
    "Datas elegíveis: 0-365 dias e consistência temporal",
    "Estádios I-III",
    "Raça/cor consistente: branca, parda ou preta",
    "Covariáveis essenciais completas"
  ),
  n = c(
    nrow(coorte_2023),
    coorte_2023 %>% filter(problema_data == "Elegível análise temporal") %>% nrow(),
    coorte_2023 %>% filter(
      problema_data == "Elegível análise temporal",
      estadio %in% c("1", "2", "3")
    ) %>% nrow(),
    coorte_2023 %>% filter(
      problema_data == "Elegível análise temporal",
      estadio %in% c("1", "2", "3"),
      raca_final %in% c("Branca", "Parda", "Preta")
    ) %>% nrow(),
    coorte_2023 %>% filter(
      problema_data == "Elegível análise temporal",
      estadio %in% c("1", "2", "3"),
      raca_final %in% c("Branca", "Parda", "Preta"),
      !is.na(fora_municipio), !is.na(idade), !is.na(cnes),
      !is.na(uf_res), !is.na(mes_inicio)
    ) %>% nrow()
  )
)

write_csv(fluxo_suplementar, file.path(DIR_SUPL, "Tabela_S1_fluxo_analitico.csv"))

# Distribuição de estágios antes da restrição I-III
tabela_estadio <- coorte_2023 %>%
  filter(problema_data == "Elegível análise temporal") %>%
  count(estadio, name = "mulheres") %>%
  mutate(percentual = round(100 * mulheres / sum(mulheres), 2))
write_csv(tabela_estadio, file.path(DIR_SUPL, "Tabela_S2_estadiamento_pre_restricao.csv"))

# 4. População analítica principal --------------------------------------------
neo_principal <- coorte_2023 %>%
  filter(
    problema_data == "Elegível análise temporal",
    estadio %in% c("1", "2", "3"),
    raca_final %in% c("Branca", "Parda", "Preta"),
    !is.na(fora_municipio),
    !is.na(idade),
    !is.na(cnes),
    !is.na(uf_res),
    !is.na(mes_inicio)
  ) %>%
  mutate(
    raca_final = factor(raca_final, levels = c("Branca", "Parda", "Preta")),
    fora_municipio = factor(
      fora_municipio,
      levels = c("Mesmo município", "Outro município")
    ),
    estadio = factor(estadio, levels = c("1", "2", "3")),
    procedimento = factor(procedimento),
    uf_res = factor(uf_res),
    mes_inicio = factor(mes_inicio),
    atraso60 = as.integer(tempo_dias > 60)
  ) %>%
  droplevels()

# 5. Descrição ----------------------------------------------------------------
resumo_coorte <- neo_principal %>%
  summarise(
    mulheres = n(),
    brancas = sum(raca_final == "Branca"),
    pardas = sum(raca_final == "Parda"),
    pretas = sum(raca_final == "Preta"),
    atraso60_n = sum(atraso60),
    atraso60_pct = round(100 * mean(atraso60), 1),
    mediana_tempo = median(tempo_dias),
    p25 = quantile(tempo_dias, .25),
    p75 = quantile(tempo_dias, .75),
    n_cnes = n_distinct(cnes),
    n_uf = n_distinct(uf_res)
  )

write_csv(resumo_coorte, file.path(DIR_SUPL, "Tabela_S3_resumo_coorte.csv"))

tabela_raca_desloc <- neo_principal %>%
  group_by(raca_final, fora_municipio) %>%
  summarise(
    mulheres = n(),
    mediana_dias = median(tempo_dias),
    p25 = quantile(tempo_dias, .25),
    p75 = quantile(tempo_dias, .75),
    atraso60_n = sum(atraso60),
    atraso60_pct = round(100 * mean(atraso60), 1),
    .groups = "drop"
  )

write_csv(
  tabela_raca_desloc,
  file.path(DIR_SUPL, "Tabela_S4_raca_x_deslocamento.csv")
)

# 6. Modelos Poisson sequenciais ---------------------------------------------
m0 <- glm(
  atraso60 ~ raca_final * fora_municipio,
  family = poisson(link = "log"), data = neo_principal
)

m1 <- glm(
  atraso60 ~ raca_final * fora_municipio + ns(idade, df = 3),
  family = poisson(link = "log"), data = neo_principal
)

m2 <- glm(
  atraso60 ~ raca_final * fora_municipio + ns(idade, df = 3) +
    estadio + procedimento,
  family = poisson(link = "log"), data = neo_principal
)

m3 <- glm(
  atraso60 ~ raca_final * fora_municipio + ns(idade, df = 3) +
    estadio + procedimento + uf_res + mes_inicio,
  family = poisson(link = "log"), data = neo_principal
)

# 7. CR2/Satterthwaite + HTZ --------------------------------------------------
analisa_cr2 <- function(modelo, dados, nome) {
  V <- vcovCR(modelo, cluster = dados$cnes, type = "CR2")

  tab <- coef_test(modelo, vcov = V, test = "Satterthwaite") %>%
    as.data.frame()

  write_csv(tab, file.path(DIR_SUPL, paste0("Tabela_S5_CR2_", nome, ".csv")))

  idx <- grep("raca_final.*:fora_municipio", names(coef(modelo)))
  global <- Wald_test(
    modelo,
    constraints = constrain_zero(idx),
    vcov = V,
    test = "HTZ"
  ) %>% as.data.frame()

  write_csv(global, file.path(DIR_SUPL, paste0("Tabela_S6_HTZ_", nome, ".csv")))

  list(V = V, coef = tab, global = global)
}

cr2_m0 <- analisa_cr2(m0, neo_principal, "M0_bruto")
cr2_m1 <- analisa_cr2(m1, neo_principal, "M1_idade")
cr2_m2 <- analisa_cr2(m2, neo_principal, "M2_clinico")
cr2_m3 <- analisa_cr2(m3, neo_principal, "M3_UF_mes")

# Tabela enxuta: termos de interação nos modelos sequenciais
extrai_interacoes <- function(obj, modelo) {
  obj$coef %>%
    as_tibble() %>%
    filter(grepl("raca_final.*fora_municipio", Coef)) %>%
    mutate(modelo = modelo, .before = 1)
}

interacoes_modelos <- bind_rows(
  extrai_interacoes(cr2_m0, "M0 bruto"),
  extrai_interacoes(cr2_m1, "M1 + idade"),
  extrai_interacoes(cr2_m2, "M2 + clínico"),
  extrai_interacoes(cr2_m3, "M3 + UF + mês")
)

write_csv(
  interacoes_modelos,
  file.path(DIR_SUPL, "Tabela_S7_interacoes_modelos_sequenciais.csv")
)

# 8. Riscos ajustados e contrastes do modelo final ---------------------------
# Para apresentação: erros robustos clusterizados por CNES.
V3_CL <- vcovCL(m3, cluster = ~ cnes, type = "HC1")

emm3 <- emmeans(
  m3,
  ~ raca_final * fora_municipio,
  vcov. = V3_CL
)

riscos_ajustados <- summary(
  emm3, type = "response", infer = c(TRUE, TRUE)
) %>% as.data.frame()

contraste_raca <- summary(
  pairs(emm3, by = "fora_municipio", adjust = "holm"),
  type = "response", infer = c(TRUE, TRUE)
) %>% as.data.frame()

contraste_desloc <- summary(
  pairs(emm3, by = "raca_final", adjust = "holm"),
  type = "response", infer = c(TRUE, TRUE)
) %>% as.data.frame()

write_csv(riscos_ajustados, file.path(DIR_SUPL, "Tabela_S8_riscos_ajustados_M3.csv"))
write_csv(contraste_raca, file.path(DIR_SUPL, "Tabela_S9_contrastes_raciais_M3.csv"))
write_csv(contraste_desloc, file.path(DIR_SUPL, "Tabela_S10_contrastes_deslocamento_M3.csv"))

# 9. Sensibilidade: abril-dezembro/2023 ---------------------------------------
neo_abr_dez <- neo_principal %>%
  filter(dt_inicio >= as.Date("2023-04-01")) %>%
  droplevels()

m_abr_dez <- glm(
  atraso60 ~ raca_final * fora_municipio + ns(idade, df = 3) +
    estadio + procedimento + uf_res + mes_inicio,
  family = poisson(link = "log"), data = neo_abr_dez
)

cr2_abr_dez <- analisa_cr2(
  m_abr_dez, neo_abr_dez, "sensibilidade_abril_dezembro_2023"
)

# 10. Sensibilidade: retirar teto de 365 dias ---------------------------------
neo_sem_teto <- coorte_2023 %>%
  filter(
    !is.na(dt_diagnostico),
    !is.na(dt_inicio),
    dt_inicio <= comp_fim,
    dt_diagnostico <= dt_inicio,
    estadio %in% c("1", "2", "3"),
    raca_final %in% c("Branca", "Parda", "Preta"),
    !is.na(fora_municipio), !is.na(idade), !is.na(cnes),
    !is.na(uf_res), !is.na(mes_inicio)
  ) %>%
  mutate(
    raca_final = factor(raca_final, levels = c("Branca", "Parda", "Preta")),
    fora_municipio = factor(
      fora_municipio,
      levels = c("Mesmo município", "Outro município")
    ),
    estadio = factor(estadio, levels = c("1", "2", "3")),
    procedimento = factor(procedimento),
    uf_res = factor(uf_res),
    mes_inicio = factor(mes_inicio),
    atraso60 = as.integer(tempo_dias > 60)
  ) %>%
  droplevels()

m_sem_teto <- glm(
  atraso60 ~ raca_final * fora_municipio + ns(idade, df = 3) +
    estadio + procedimento + uf_res + mes_inicio,
  family = poisson(link = "log"), data = neo_sem_teto
)

cr2_sem_teto <- analisa_cr2(
  m_sem_teto, neo_sem_teto, "sensibilidade_sem_teto_365"
)

# 11. Sensibilidade: raça/cor da APAC índice ---------------------------------
recodifica_raca_indice <- function(x) {
  case_when(
    x == "01" ~ "Branca",
    x == "02" ~ "Preta",
    x == "03" ~ "Parda",
    TRUE ~ NA_character_
  )
}

neo_raca_indice <- coorte_2023 %>%
  mutate(raca_indice = recodifica_raca_indice(raca_cor)) %>%
  filter(
    problema_data == "Elegível análise temporal",
    estadio %in% c("1", "2", "3"),
    raca_indice %in% c("Branca", "Parda", "Preta"),
    !is.na(fora_municipio), !is.na(idade), !is.na(cnes),
    !is.na(uf_res), !is.na(mes_inicio)
  ) %>%
  mutate(
    raca_indice = factor(raca_indice, levels = c("Branca", "Parda", "Preta")),
    fora_municipio = factor(
      fora_municipio,
      levels = c("Mesmo município", "Outro município")
    ),
    estadio = factor(estadio, levels = c("1", "2", "3")),
    procedimento = factor(procedimento),
    uf_res = factor(uf_res),
    mes_inicio = factor(mes_inicio),
    atraso60 = as.integer(tempo_dias > 60)
  ) %>%
  droplevels()

m_raca_indice <- glm(
  atraso60 ~ raca_indice * fora_municipio + ns(idade, df = 3) +
    estadio + procedimento + uf_res + mes_inicio,
  family = poisson(link = "log"), data = neo_raca_indice
)

V_raca_indice <- vcovCR(
  m_raca_indice,
  cluster = neo_raca_indice$cnes,
  type = "CR2"
)

write_csv(
  coef_test(m_raca_indice, vcov = V_raca_indice, test = "Satterthwaite") %>%
    as.data.frame(),
  file.path(DIR_SUPL, "Tabela_S11_sensibilidade_raca_APAC_indice.csv")
)

# 12. Leave-one-CNES-out -------------------------------------------------------
if (RODAR_LOO) {
  termo_alvo <- "raca_finalPreta:fora_municipioOutro município"

  loo <- map_dfr(unique(neo_principal$cnes), function(cnes_excluido) {
    d <- neo_principal %>% filter(cnes != cnes_excluido) %>% droplevels()

    mod <- glm(
      atraso60 ~ raca_final * fora_municipio + ns(idade, df = 3) +
        estadio + procedimento + uf_res + mes_inicio,
      family = poisson(link = "log"), data = d
    )

    beta <- unname(coef(mod)[termo_alvo])

    tibble(
      cnes_excluido = cnes_excluido,
      n = nrow(d),
      beta_interacao = beta,
      RR_interacao = exp(beta)
    )
  })

  write_csv(loo, file.path(DIR_SUPL, "Tabela_S12_leave_one_CNES_out.csv"))

  resumo_loo <- loo %>%
    summarise(
      n_cnes = n(),
      RR_min = min(RR_interacao, na.rm = TRUE),
      RR_p25 = quantile(RR_interacao, .25, na.rm = TRUE),
      RR_mediana = median(RR_interacao, na.rm = TRUE),
      RR_p75 = quantile(RR_interacao, .75, na.rm = TRUE),
      RR_max = max(RR_interacao, na.rm = TRUE)
    )

  write_csv(resumo_loo, file.path(DIR_SUPL, "Tabela_S13_resumo_leave_one_CNES_out.csv"))
}

# 13. Diagnóstico dos clusters ------------------------------------------------
diagnostico_cnes <- neo_principal %>%
  count(cnes, sort = TRUE, name = "mulheres")

resumo_cnes <- diagnostico_cnes %>%
  summarise(
    n_cnes = n(),
    minimo = min(mulheres),
    mediana = median(mulheres),
    maximo = max(mulheres)
  )

write_csv(diagnostico_cnes, file.path(DIR_SUPL, "Tabela_S14_distribuicao_CNES.csv"))
write_csv(resumo_cnes, file.path(DIR_SUPL, "Tabela_S15_resumo_clusters_CNES.csv"))

# 14. Resumo no console --------------------------------------------------------
cat("\n================ MATERIAL SUPLEMENTAR ================\n")
cat("Coorte principal: Brasil, 2023\n\n")
print(resumo_coorte)
cat("\nRaça/cor x deslocamento:\n")
print(tabela_raca_desloc)
cat("\nCR2 — modelo final M3:\n")
print(cr2_m3$coef)
cat("\nTeste global HTZ — modelo final M3:\n")
print(cr2_m3$global)
cat("\nArquivos exportados em: ", normalizePath(DIR_SUPL), "\n")
cat("Concluído.\n")
