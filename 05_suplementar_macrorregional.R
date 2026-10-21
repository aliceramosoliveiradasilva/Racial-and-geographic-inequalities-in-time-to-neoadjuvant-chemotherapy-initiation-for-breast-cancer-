# PUBLICATION REPRODUCIBILITY — BREAST CANCER NEOADJUVANT CHEMOTHERAPY
# Macroregional supplementary sensitivity and model outputs.
# Executes locally on authorized DATASUS data; never commit patient-level data.
# 2026-10-15. Historical crosswalks 2024/2026 are not validated for 2023.
# All R scripts should be run from project root; MAMA_RAIZ configures data root.
# NOTE: full end-to-end re-execution was not possible in the delivery environment.
# -------------------------------------------------------------------------

# ============================================================================
# SCRIPT SUPLEMENTAR S2 — ANÁLISES ESTATÍSTICAS E MACRORREGIOES
# REANALISE TERRITORIAL: 14/10/2026.
# Municipal original preservado; novas tabelas em subpasta separada.
# CSV de 2026 NAO comprovado como recorte historico de 2023.
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
# RODADA RAPIDA: aproveita a coorte de 2023 JA EXPORTADA no script 02.
# Nao relê/minera APAC de 2019 a 2025. Nao sobrescreve o suplemento antigo.
RAIZ_MACRO <- Sys.getenv("MAMA_RAIZ", unset = getwd())
PASTA_BASE <- file.path(RAIZ_MACRO, "dados_mama_SIA_AQ_2019_2025")
DIR_RESULTADOS_ORIG <- file.path(PASTA_BASE, "resultados")
DIR_RESULTADOS <- file.path(PASTA_BASE, "resultados_com_macrorregioes")
DIR_SUPL <- file.path(DIR_RESULTADOS, "material_suplementar")
dir.create(DIR_SUPL, recursive = TRUE, showWarnings = FALSE)

# Primeiro procura a coorte original para comparacao independente.
ARQ_COORTE <- file.path(
  DIR_RESULTADOS_ORIG,
  "coorte_episodios_neoadjuvantes_2023_pre_analise.fst"
)
if (!file.exists(ARQ_COORTE)) {
  ARQ_COORTE <- file.path(
    DIR_RESULTADOS, "coorte_episodios_neoadjuvantes_2023_pre_analise.fst")
}
if (!file.exists(ARQ_COORTE)) {
  stop("Coorte pre-analise nao encontrada. Execute primeiro o script 02 original ",
       "ou o 02 com macrorregioes. Esperado: ", ARQ_COORTE)
}
cat("Coorte de entrada para S2: ", ARQ_COORTE, "\n", sep = "")

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

# Evita expandir a grade completa de UF, mes, estadio e procedimento.
emm3 <- emmeans(
  m3,
  ~ raca_final * fora_municipio,
  vcov. = V3_CL,
  nuisance = c("uf_res", "mes_inicio", "estadio", "procedimento"),
  wt.nuis = "proportional"
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


# ---------------------------------------------------------------------------
# 15. Analises suplementares adicionais com macrorregioes de saude
# ---------------------------------------------------------------------------
# ===========================================================================
# EXTENSAO: MACRORREGIOES DE SAUDE (2026) — COORTE ONCOLOGICA DE 2023
# Fonte: CSV do Ministerio da Saude fornecido pela pesquisadora.
# NOTA: a referencia territorial foi atualizada em 2026 e sua validade
# retrospectiva para 2023 precisa ser verificada antes da submissao.
# Mantem desfecho, elegibilidade e inferencia CR2/HTZ do estudo original.
# ===========================================================================

macro_normalizar_mun <- function(x) {
  z <- trimws(as.character(x))
  z <- sub("\\.0+$", "", z)
  z[z %in% c("", "NA", "NaN", "NULL", "000000", "0000000")] <- NA_character_
  valido <- !is.na(z) & grepl("^[0-9]{6,7}$", z)
  z[!valido] <- NA_character_
  substr(z, 1L, 6L)
}

macro_localizar_csv <- function(raiz) {
  opcoes <- c(
    file.path(raiz, "macroregiao_de_saude(1).csv"),
    file.path(raiz, "macroregiao_de_saude.csv")
  )
  existentes <- opcoes[file.exists(opcoes)]
  if (length(existentes) > 0) return(existentes[1L])
  encontrados <- list.files(raiz, pattern = "^macroregiao_de_saude.*\\.csv$",
                            full.names = TRUE, ignore.case = TRUE)
  if (length(encontrados) != 1L) {
    stop("Nao foi possivel identificar UNIVOCAMENTE o CSV de macrorregioes em ",
         raiz, ". Encontrados: ", paste(basename(encontrados), collapse = ", "))
  }
  encontrados[1L]
}

macro_ler_referencia <- function(csv) {
  mapa <- readr::read_delim(
    file = csv, delim = ";",
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE, progress = FALSE,
    locale = readr::locale(encoding = "UTF-8")
  )
  necessarias <- c(
    "cod_municipio", "cod_macrorregiao_de_saude", "macrorregiao_de_saude",
    "cod_regiao_de_saude", "regiao_de_saude", "sg_uf", "co_uf"
  )
  ausentes <- setdiff(necessarias, names(mapa))
  if (length(ausentes)) {
    stop("Colunas ausentes no CSV de macrorregioes: ", paste(ausentes, collapse = ", "))
  }
  mapa <- dplyr::transmute(
    mapa,
    municipio_6 = macro_normalizar_mun(cod_municipio),
    macro_cod = as.character(cod_macrorregiao_de_saude),
    macro_nome = as.character(macrorregiao_de_saude),
    regiao_cod = as.character(cod_regiao_de_saude),
    regiao_nome = as.character(regiao_de_saude),
    sg_uf = as.character(sg_uf),
    cod_uf = as.character(co_uf)
  )
  if (nrow(mapa) == 0 || anyNA(mapa$municipio_6) ||
      anyDuplicated(mapa$municipio_6) > 0 ||
      any(is.na(mapa$macro_cod) | !nzchar(mapa$macro_cod)) ||
      any(is.na(mapa$regiao_cod) | !nzchar(mapa$regiao_cod))) {
    stop("CSV territorial tem chaves duplicadas, ausentes ou vazias; vinculo interrompido.")
  }
  if (any(substr(mapa$municipio_6, 1, 2) != mapa$cod_uf)) {
    stop("Codigos municipais com prefixo de UF inconsistente no CSV.")
  }
  if (any(substr(mapa$macro_cod, 1, 2) != mapa$cod_uf)) {
    stop("Macrorregioes e UF inconsistentes no CSV.")
  }
  if (nrow(mapa) != 5570L) {
    warning("CSV contem ", nrow(mapa), " municipios (arquivo examinado tinha 5.570).")
  }
  if (dplyr::n_distinct(mapa$macro_cod) != 121L) {
    warning("CSV contem ", dplyr::n_distinct(mapa$macro_cod),
            " macrorregioes (arquivo examinado tinha 121).")
  }
  mapa
}

macro_vincular_coorte <- function(d, mapa) {
  requeridas <- c("mun_res", "mun_estab")
  faltam <- setdiff(requeridas, names(d))
  if (length(faltam)) {
    stop("A base nao tem as colunas usadas no script original: ",
         paste(faltam, collapse = ", "),
         ". E necessario usar a COORTE original, nao o CSV municipal.")
  }
  n_antes <- nrow(d)
  d <- tibble::as_tibble(d)
  d$mun_res6_macro <- macro_normalizar_mun(d$mun_res)
  d$mun_trat6_macro <- macro_normalizar_mun(d$mun_estab)
  i_res <- match(d$mun_res6_macro, mapa$municipio_6)
  i_trat <- match(d$mun_trat6_macro, mapa$municipio_6)
  d$macro_res_cod <- mapa$macro_cod[i_res]
  d$macro_res_nome <- mapa$macro_nome[i_res]
  d$macro_trat_cod <- mapa$macro_cod[i_trat]
  d$macro_trat_nome <- mapa$macro_nome[i_trat]
  d$reg_saude_res_cod <- mapa$regiao_cod[i_res]
  d$reg_saude_res_nome <- mapa$regiao_nome[i_res]
  d$reg_saude_trat_cod <- mapa$regiao_cod[i_trat]
  d$reg_saude_trat_nome <- mapa$regiao_nome[i_trat]

  # NA permanece NA, mesmo se os codigos municipais nao estiverem na referencia.
  mesma_macro <- !is.na(d$macro_res_cod) & !is.na(d$macro_trat_cod) &
    d$macro_res_cod == d$macro_trat_cod
  outra_macro <- !is.na(d$macro_res_cod) & !is.na(d$macro_trat_cod) &
    d$macro_res_cod != d$macro_trat_cod
  mesmo_mun <- !is.na(d$mun_res6_macro) & !is.na(d$mun_trat6_macro) &
    d$mun_res6_macro == d$mun_trat6_macro
  outro_mun <- !is.na(d$mun_res6_macro) & !is.na(d$mun_trat6_macro) &
    d$mun_res6_macro != d$mun_trat6_macro
  mesma_uf <- !is.na(d$mun_res6_macro) & !is.na(d$mun_trat6_macro) &
    substr(d$mun_res6_macro, 1, 2) == substr(d$mun_trat6_macro, 1, 2)
  outra_uf <- !is.na(d$mun_res6_macro) & !is.na(d$mun_trat6_macro) &
    substr(d$mun_res6_macro, 1, 2) != substr(d$mun_trat6_macro, 1, 2)
  mesma_regiao <- !is.na(d$reg_saude_res_cod) & !is.na(d$reg_saude_trat_cod) &
    d$reg_saude_res_cod == d$reg_saude_trat_cod
  outra_regiao <- !is.na(d$reg_saude_res_cod) & !is.na(d$reg_saude_trat_cod) &
    d$reg_saude_res_cod != d$reg_saude_trat_cod

  d$fora_macro <- factor(
    dplyr::case_when(
      mesma_macro ~ "Mesma macrorregiao",
      outra_macro ~ "Outra macrorregiao",
      TRUE ~ NA_character_
    ),
    levels = c("Mesma macrorregiao", "Outra macrorregiao")
  )
  # Modelo exploratorio mais parcimonioso: 3 niveis operacionais.
  # Evita separar deslocamentos inter-UF quando ha celulas raciais pequenas.
  d$trajeto_macro3 <- factor(
    dplyr::case_when(
      mesmo_mun & mesma_macro ~ "1 Mesmo municipio",
      outro_mun & mesma_macro ~ "2 Outro municipio na mesma macro",
      outra_macro ~ "3 Fora da macrorregiao",
      TRUE ~ NA_character_
    ),
    levels = c("1 Mesmo municipio", "2 Outro municipio na mesma macro",
               "3 Fora da macrorregiao")
  )
  d$trajeto_macro4 <- factor(
    dplyr::case_when(
      mesmo_mun & mesma_macro ~ "1 Mesmo municipio",
      outro_mun & mesma_macro ~ "2 Outro municipio na mesma macro",
      outra_macro & mesma_uf ~ "3 Outra macro na mesma UF",
      outra_macro & outra_uf ~ "4 Outra UF",
      TRUE ~ NA_character_
    ),
    levels = c("1 Mesmo municipio", "2 Outro municipio na mesma macro",
               "3 Outra macro na mesma UF", "4 Outra UF")
  )
  # Uma regiao de saude no CSV (RO, 11003) atravessa duas macros.
  # Portanto NAO presumir que as regioes de saude sao perfeitamente aninhadas;
  # o fluxo entre macros tem precedencia na classificacao combinada.
  d$trajeto_rede5 <- factor(
    dplyr::case_when(
      mesmo_mun & mesma_macro ~ "1 Mesmo municipio",
      outra_macro & mesma_uf ~ "4 Outra macro na mesma UF",
      outra_macro & outra_uf ~ "5 Outra UF",
      outro_mun & mesma_regiao & mesma_macro ~ "2 Outro municipio na mesma regiao de saude",
      outro_mun & outra_regiao & mesma_macro ~ "3 Outra regiao de saude na mesma macro",
      TRUE ~ NA_character_
    ),
    levels = c("1 Mesmo municipio", "2 Outro municipio na mesma regiao de saude",
               "3 Outra regiao de saude na mesma macro", "4 Outra macro na mesma UF",
               "5 Outra UF")
  )
  if (any(mesmo_mun & outra_macro, na.rm = TRUE)) {
    stop("Inconsistencia territorial: mesmo municipio classificado em macros distintas.")
  }
  if (nrow(d) != n_antes) stop("A vinculacao alterou o numero de mulheres/linhas!")
  d
}

macro_auditoria <- function(d, descricao) {
  tibble::tibble(
    base = descricao,
    linhas = nrow(d),
    municipio_res_valido = sum(!is.na(d$mun_res6_macro)),
    municipio_trat_valido = sum(!is.na(d$mun_trat6_macro)),
    macro_res_encontrada = sum(!is.na(d$macro_res_cod)),
    macro_trat_encontrada = sum(!is.na(d$macro_trat_cod)),
    fluxo_macro_identificado = sum(!is.na(d$fora_macro)),
    fora_macro = sum(d$fora_macro == "Outra macrorregiao", na.rm = TRUE),
    mesma_macro = sum(d$fora_macro == "Mesma macrorregiao", na.rm = TRUE),
    fluxo_3_identificado = sum(!is.na(d$trajeto_macro3)),
    fluxo_4_identificado = sum(!is.na(d$trajeto_macro4)),
    fluxo_5_identificado = sum(!is.na(d$trajeto_rede5)),
    pendentes = sum(is.na(d$fora_macro))
  )
}

macro_preparar_modelo <- function(d) {
  d <- dplyr::filter(d, !is.na(fora_macro))
  d <- dplyr::mutate(
    d,
    fora_macro = factor(fora_macro,
      levels = c("Mesma macrorregiao", "Outra macrorregiao")),
    raca_final = factor(raca_final, levels = c("Branca", "Parda", "Preta")),
    estadio = factor(estadio, levels = c("1", "2", "3")),
    procedimento = factor(procedimento),
    uf_res = factor(uf_res),
    mes_inicio = factor(mes_inicio),
    atraso60 = as.integer(tempo_dias > 60)
  )
  droplevels(d)
}

macro_checar_celulas <- function(d, expo, minimo = 10L) {
  tab <- table(d$raca_final, d[[expo]])
  if (any(tab < minimo) || ncol(tab) < 2L || nrow(tab) < 3L) {
    return(list(ok = FALSE, tabela = as.data.frame(tab)))
  }
  list(ok = TRUE, tabela = as.data.frame(tab))
}

macro_ajustar_modelos <- function(d, exposicao = "fora_macro") {
  f <- function(adicionais = "") {
    stats::as.formula(paste0("atraso60 ~ raca_final * ", exposicao, adicionais))
  }
  list(
    M0 = stats::glm(f(), family = stats::poisson(link = "log"), data = d),
    M1 = stats::glm(f(" + splines::ns(idade, df = 3)"),
                    family = stats::poisson(link = "log"), data = d),
    M2 = stats::glm(f(" + splines::ns(idade, df = 3) + estadio + procedimento"),
                    family = stats::poisson(link = "log"), data = d),
    M3 = macro_ajustar_final(d, exposicao)
  )
}

macro_ajustar_final <- function(d, exposicao = "fora_macro") {
  anos <- if ("ano_inicio" %in% names(d)) d$ano_inicio else rep(2023, nrow(d))
  ajuste_ano <- if (dplyr::n_distinct(anos) > 1L) " + factor(ano_inicio)" else ""
  formula <- stats::as.formula(
    paste0("atraso60 ~ raca_final * ", exposicao,
           " + splines::ns(idade, df = 3) + estadio + procedimento + uf_res + mes_inicio",
           ajuste_ano)
  )
  stats::glm(formula, family = stats::poisson(link = "log"), data = d)
}

macro_cr2_htz <- function(modelo, d, exposicao, nome, destino) {
  if (nrow(stats::model.frame(modelo)) != nrow(d)) {
    stop("Modelo ", nome, ": glm excluiu linhas por NA; cluster CNES pode ficar desalinhado.")
  }
  vc <- clubSandwich::vcovCR(modelo, cluster = d$cnes, type = "CR2")
  coef <- as.data.frame(clubSandwich::coef_test(
    modelo, vcov = vc, test = "Satterthwaite"))
  readr::write_csv(coef, file.path(destino, paste0(nome, "_coef_CR2.csv")))
  indices <- grep(paste0("^raca_final[^:]*:", exposicao), names(stats::coef(modelo)))
  if (length(indices) == 0L) {
    stop("Sem termos de interacao estimaveis no modelo ", nome)
  }
  global <- as.data.frame(clubSandwich::Wald_test(
    modelo, constraints = clubSandwich::constrain_zero(indices),
    vcov = vc, test = "HTZ"))
  readr::write_csv(global, file.path(destino, paste0(nome, "_interacao_global_HTZ.csv")))
  list(V = vc, coef = coef, global = global)
}

macro_estimar_riscos <- function(modelo, V, d, exposicao, destino, prefixo = "M3") {
  # Mesma parametrizacao da analise original (emmeans); aqui CR2 tambem
  # informa os intervalos dos riscos e dos contrastes, evitando HC1 x CR2.
  emm <- emmeans::emmeans(
    modelo,
    specs = stats::as.formula(paste0("~ raca_final * ", exposicao)),
    vcov. = V,
    nuisance = c("uf_res", "mes_inicio", "estadio", "procedimento"),
    wt.nuis = "proportional"
  )
  riscos <- as.data.frame(summary(emm, type = "response", infer = c(TRUE, TRUE)))
  entre_racas <- as.data.frame(summary(
    pairs(emm, by = exposicao, adjust = "holm"),
    type = "response", infer = c(TRUE, TRUE)))
  entre_trajetos <- as.data.frame(summary(
    pairs(emm, by = "raca_final", adjust = "holm"),
    type = "response", infer = c(TRUE, TRUE)))
  preta_branca <- as.data.frame(summary(
    emmeans::contrast(
      emm, method = list("Preta_vs_Branca" = c(-1, 0, 1)),
      by = exposicao, adjust = "none"),
    type = "response", infer = c(TRUE, TRUE)))
  readr::write_csv(riscos, file.path(destino, paste0(prefixo, "_riscos_ajustados.csv")))
  readr::write_csv(entre_racas, file.path(destino, paste0(prefixo, "_contrastes_raciais.csv")))
  readr::write_csv(entre_trajetos, file.path(destino, paste0(prefixo, "_contrastes_territoriais.csv")))
  readr::write_csv(preta_branca, file.path(destino, paste0(prefixo, "_preta_vs_branca_por_fluxo.csv")))
  invisible(list(riscos = riscos, raca = entre_racas,
                 trajeto = entre_trajetos, preta_branca = preta_branca))
}

macro_guardar_proveniencia <- function(mapa, arq_mapa, pasta) {
  readr::write_csv(mapa, file.path(pasta, "dicionario_municipio_macrorregiao_utilizado.csv"))
  # Auditoria de regioes de saude que aparecem em mais de uma macro
  # (na referencia fornecida, regiao 11003, Rondônia).
  regioes_multimacro <- mapa %>%
    dplyr::group_by(regiao_cod, regiao_nome, sg_uf) %>%
    dplyr::summarise(n_macros = dplyr::n_distinct(macro_cod),
                     macros = paste(sort(unique(macro_cod)), collapse = ";"),
                     municipios = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n_macros > 1L)
  readr::write_csv(regioes_multimacro,
                   file.path(pasta, "AUDITORIA_REGIOES_SAUDE_EM_VARIAS_MACROS.csv"))
  writeLines(c(
    "Projeto: tempo ate quimioterapia neoadjuvante para cancer de mama no SUS",
    paste("Gerado em:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("CSV:", normalizePath(arq_mapa)),
    paste("MD5:", unname(tools::md5sum(arq_mapa))),
    paste("Municipios:", nrow(mapa)),
    paste("Macrorregioes distintas:", dplyr::n_distinct(mapa$macro_cod)),
    paste("Regioes de saude distintas:", dplyr::n_distinct(mapa$regiao_cod)),
    "Base do Ministerio da Saude consultada em 13/10/2026.",
    "IMPORTANTE: divisao territorial de 2026 NAO validada para ano da coorte (2023).",
    "Unidade de analise: MULHERES da coorte indice, e nao registros APAC.",
    "Exposicao principal: municipio de residencia e municipio de tratamento",
    "classificados pelas macrorregioes de saude do CSV fornecido.",
    "Outro municipio dentro da mesma macro = fluxo intramacrorregional.",
    "Outra macrorregiao nao implica necessariamente barreira assistencial.",
    "A classificacao NAO mede distancia, tempo de viagem ou TFD.",
    "Macroregiao de saude e diferente das 5 grandes regioes brasileiras.",
    "O CSV pode conter regiao de saude em mais de uma macro; veja auditoria correspondente.",
    "Verificar pactuacoes e possiveis mudancas territoriais antes de submeter."
  ), file.path(pasta, "PROVENIENCIA_E_LIMITACOES.txt"))
}

macro_executar_analises <- function(
  d_principal, coorte_pre, sensibilidades, raiz, destino,
  rodar_loo = FALSE, esperado_n = 15415L
) {
  dir.create(destino, showWarnings = FALSE, recursive = TRUE)
  arq_mapa <- macro_localizar_csv(raiz)
  mapa <- macro_ler_referencia(arq_mapa)
  macro_guardar_proveniencia(mapa, arq_mapa, destino)
  cat("\n>>> Vinculando tabela de macrorregioes: ", basename(arq_mapa), "\n", sep = "")

  # A coorte antes das exclusoes inferenciais permanece a mesma; a tabela
  # municipal so acrescenta atributos, sem alterar registros ou identidades.
  pre <- macro_vincular_coorte(coorte_pre, mapa)
  fst::write_fst(as.data.frame(pre),
    file.path(destino, "coorte_episodios_2023_COM_MACROS_pre_analise.fst"), compress = 50)

  d_p <- macro_vincular_coorte(d_principal, mapa)
  if (nrow(d_p) != nrow(d_principal)) stop("Alteracao de N no join da coorte.")
  aud <- dplyr::bind_rows(
    macro_auditoria(pre, "coorte_2023_pre_analise"),
    macro_auditoria(d_p, "coorte_2023_analitica_original")
  )
  if (!is.null(esperado_n) && nrow(d_p) != esperado_n) {
    warning("Coorte original tem ", nrow(d_p),
            " mulheres; artigo reporta ", esperado_n, ". Verifique versao do banco.")
  }
  if (anyNA(d_p$fora_macro)) {
    warning("Ha ", sum(is.na(d_p$fora_macro)),
            " mulheres sem fluxo macro identificavel; somente essas ficam fora dos modelos macro.")
  }
  readr::write_csv(aud, file.path(destino, "AUDITORIA_COBERTURA_VINCULO.csv"))

  # Auditabilidade: codigos municipais nao encontrados, SEM dados de paciente.
  cod_sem_res <- unique(
    d_p$mun_res6_macro[!is.na(d_p$mun_res6_macro) & is.na(d_p$macro_res_cod)])
  cod_sem_trat <- unique(
    d_p$mun_trat6_macro[!is.na(d_p$mun_trat6_macro) & is.na(d_p$macro_trat_cod)])
  cod_nao_map <- dplyr::bind_rows(
    tibble::tibble(lado = rep("residencia", length(cod_sem_res)), municipio = cod_sem_res),
    tibble::tibble(lado = rep("tratamento", length(cod_sem_trat)), municipio = cod_sem_trat)
  )
  readr::write_csv(cod_nao_map, file.path(destino, "MUNICIPIOS_NAO_VINCULADOS.csv"))

  # Preserva comparabilidade: este arquivo inclui TODAS as mulheres da
  # populacao municipal original, incluindo os eventuais NA de macrorregiao.
  fst::write_fst(as.data.frame(d_p),
    file.path(destino, "coorte_analitica_2023_COM_MACROS.fst"), compress = 50)
  comparacao <- d_p %>%
    dplyr::count(fora_municipio, fora_macro, name = "mulheres") %>%
    dplyr::arrange(fora_municipio, fora_macro)
  readr::write_csv(comparacao, file.path(destino, "MATRIZ_MUNICIPIO_VS_MACRO.csv"))

  d_macro <- macro_preparar_modelo(d_p)
  if (nrow(d_macro) < 100L) {
    stop("Menos de 100 mulheres com macrorregiao identificada: verifique a correspondencia.")
  }
  tab_celulas <- macro_checar_celulas(d_macro, "fora_macro", minimo = 10L)
  readr::write_csv(tab_celulas$tabela,
                   file.path(destino, "CELULAS_RACA_X_FLUXO_MACRO.csv"))
  if (!tab_celulas$ok) {
    stop("Algum estrato raca/cor x macro tem menos de 10 mulheres. ",
         "Modelos inferenciais macro nao foram ajustados; veja as tabelas descritivas.")
  }

  # Descritivos: os totais aqui sao MULHERES, nao autorizacoes repetidas.
  desc <- d_macro %>%
    dplyr::group_by(raca_final, fora_macro) %>%
    dplyr::summarise(
      mulheres = dplyr::n(), atrasos_60 = sum(atraso60),
      percentual_atraso = round(100 * mean(atraso60), 2),
      mediana_dias = median(tempo_dias),
      p25 = as.numeric(stats::quantile(tempo_dias, 0.25)),
      p75 = as.numeric(stats::quantile(tempo_dias, 0.75)),
      n_cnes = dplyr::n_distinct(cnes),
      .groups = "drop"
    )
  readr::write_csv(desc, file.path(destino, "TABELA_RACA_X_MACRO.csv"))

  tab_macro3 <- d_macro %>% dplyr::group_by(raca_final, trajeto_macro3) %>%
    dplyr::summarise(mulheres = dplyr::n(), atrasos_60 = sum(atraso60),
      percentual_atraso = round(100 * mean(atraso60), 2),
      mediana_dias = median(tempo_dias), .groups = "drop")
  readr::write_csv(tab_macro3, file.path(destino, "TABELA_FLUXO_3_NIVEIS.csv"))

  tab_macro4 <- d_macro %>% dplyr::group_by(raca_final, trajeto_macro4) %>%
    dplyr::summarise(mulheres = dplyr::n(), atrasos_60 = sum(atraso60),
      percentual_atraso = round(100 * mean(atraso60), 2),
      mediana_dias = median(tempo_dias), .groups = "drop")
  readr::write_csv(tab_macro4, file.path(destino, "TABELA_FLUXO_4_NIVEIS.csv"))

  tab_rede5 <- d_macro %>% dplyr::group_by(raca_final, trajeto_rede5) %>%
    dplyr::summarise(mulheres = dplyr::n(), atrasos_60 = sum(atraso60),
      percentual_atraso = round(100 * mean(atraso60), 2),
      mediana_dias = median(tempo_dias), .groups = "drop")
  readr::write_csv(tab_rede5, file.path(destino, "TABELA_FLUXO_5_NIVEIS_REGIAO_SAUDE.csv"))

  orig_dest <- d_macro %>%
    dplyr::count(macro_res_cod, macro_res_nome,
                 macro_trat_cod, macro_trat_nome, name = "mulheres") %>%
    dplyr::arrange(dplyr::desc(mulheres))
  readr::write_csv(orig_dest, file.path(destino, "MATRIZ_ORIGEM_DESTINO_MACROS_2023.csv"))

  distrib_macro <- d_macro %>%
    dplyr::group_by(macro_res_cod, macro_res_nome, raca_final) %>%
    dplyr::summarise(mulheres = dplyr::n(),
                     pct_atraso = round(100 * mean(atraso60), 2),
                     .groups = "drop")
  readr::write_csv(distrib_macro,
                   file.path(destino, "DESCRITIVO_POR_MACRO_RESIDENCIA_E_RACA.csv"))

  # Modelos sequenciais M0-M3: preservam ajustes originais, mudando
  # SOMENTE a exposicao municipio -> macrorregiao de saude.
  mod <- macro_ajustar_modelos(d_macro, "fora_macro")
  resultados <- list()
  for (nome in names(mod)) {
    resultados[[nome]] <- macro_cr2_htz(
      mod[[nome]], d_macro, "fora_macro", paste0("MACRO_", nome, "_2023"), destino)
  }
  macro_estimar_riscos(mod$M3, resultados$M3$V, d_macro,
                       "fora_macro", destino, prefixo = "MACRO_M3_2023")

  # O modelo municipal na MESMA subamostra avalia se eventuais falhas de
  # vinculacao mudaram N. Mantem inferencia CR2/HTZ.
  mod_mun_pareado <- macro_ajustar_final(d_macro, "fora_municipio")
  res_mun_pareado <- macro_cr2_htz(mod_mun_pareado, d_macro,
    "fora_municipio", "MUNICIPIO_M3_PAREADO_2023", destino)
  macro_estimar_riscos(mod_mun_pareado, res_mun_pareado$V, d_macro,
    "fora_municipio", destino, prefixo = "MUNICIPIO_M3_PAREADO_2023")

  # Analise secundaria em 3 niveis, com criterio de estabilidade por celula.
  tab3 <- macro_checar_celulas(d_macro, "trajeto_macro3", minimo = 10L)
  readr::write_csv(tab3$tabela,
                   file.path(destino, "CELULAS_RACA_X_TRAJETO_3_NIVEIS.csv"))
  status3 <- "Nao realizado: alguma celula com <10 mulheres ou nivel ausente."
  if (tab3$ok && ncol(table(d_macro$raca_final, d_macro$trajeto_macro3)) == 3L) {
    ajuste3 <- macro_ajustar_final(d_macro, "trajeto_macro3")
    cr3 <- macro_cr2_htz(ajuste3, d_macro, "trajeto_macro3",
                         "TRAJETO3_M3_2023", destino)
    macro_estimar_riscos(ajuste3, cr3$V, d_macro,
      "trajeto_macro3", destino, prefixo = "TRAJETO3_M3_2023")
    status3 <- "Realizado: todos os cruzamentos raca/cor x trajeto >=10 mulheres."
  }
  writeLines(status3, file.path(destino, "STATUS_MODELO_TRAJETO_3_NIVEIS.txt"))

  # Analise secundaria em 4 niveis apenas quando TODOS os cruzamentos
  # raca x trajeto tiverem pelo menos 10 mulheres. Caso contrario, nao
  # estima RR instaveis: publica a tabela descritiva e o motivo.
  tab4 <- macro_checar_celulas(d_macro, "trajeto_macro4", minimo = 10L)
  readr::write_csv(tab4$tabela,
                   file.path(destino, "CELULAS_RACA_X_TRAJETO_4_NIVEIS.csv"))
  status4 <- "Nao realizado: alguma celula tem <10 mulheres ou nivel ausente."
  if (tab4$ok && ncol(table(d_macro$raca_final, d_macro$trajeto_macro4)) == 4L) {
    ajuste4 <- macro_ajustar_final(d_macro, "trajeto_macro4")
    cr4 <- macro_cr2_htz(ajuste4, d_macro, "trajeto_macro4",
                         "TRAJETO4_M3_2023", destino)
    macro_estimar_riscos(ajuste4, cr4$V, d_macro,
      "trajeto_macro4", destino, prefixo = "TRAJETO4_M3_2023")
    status4 <- "Realizado: todos os cruzamentos raca/cor x trajeto >=10 mulheres."
  }
  writeLines(status4, file.path(destino, "STATUS_MODELO_TRAJETO_4_NIVEIS.txt"))

  # Sensibilidades que ja estavam na estrutura original, se disponibilizadas.
  # Principal: 2023; ampliadas 2022-23 e 2023-24 quando o script 02 inteiro roda.
  audit_sens <- list()
  for (nm in names(sensibilidades)) {
    d_s <- macro_vincular_coorte(sensibilidades[[nm]], mapa)
    audit_sens[[nm]] <- macro_auditoria(d_s, nm)
    d_s <- macro_preparar_modelo(d_s)
    if (nrow(d_s) < 100L || !macro_checar_celulas(d_s, "fora_macro")$ok) {
      warning("Sensibilidade ", nm, ": amostra/celulas insuficientes; modelo nao executado.")
      next
    }
    m_s <- macro_ajustar_final(d_s, "fora_macro")
    macro_cr2_htz(m_s, d_s, "fora_macro", paste0("SENS_", nm), destino)
    readr::write_csv(
      dplyr::summarise(d_s, mulheres = dplyr::n(),
          atrasos_60 = sum(atraso60), pct_atraso = 100 * mean(atraso60),
          n_cnes = dplyr::n_distinct(cnes)),
      file.path(destino, paste0("SENS_", nm, "_resumo.csv")))
  }
  if (length(audit_sens)) {
    readr::write_csv(dplyr::bind_rows(audit_sens),
                     file.path(destino, "AUDITORIA_SENSIBILIDADES.csv"))
  }

  # LOO: avalia influencia da exclusao de cada estabelecimento CNES
  # no termo racial preto x deslocamento macrorregional (opcional).
  if (isTRUE(rodar_loo)) {
    termo <- "raca_finalPreta:fora_macroOutra macrorregiao"
    loo <- purrr::map_dfr(unique(d_macro$cnes), function(cnes_fora) {
      ds <- droplevels(d_macro[d_macro$cnes != cnes_fora, , drop = FALSE])
      if (!macro_checar_celulas(ds, "fora_macro")$ok) {
        return(tibble::tibble(cnes_excluido = as.character(cnes_fora),
          mulheres = nrow(ds), rr_interacao = NA_real_))
      }
      ms <- macro_ajustar_final(ds, "fora_macro")
      val <- unname(stats::coef(ms)[termo])
      tibble::tibble(cnes_excluido = as.character(cnes_fora),
                     mulheres = nrow(ds), rr_interacao = exp(val))
    })
    readr::write_csv(loo, file.path(destino, "SENS_LEAVE_ONE_CNES_OUT.csv"))
  }

  cat("\n========= MACRORREGIOES DE SAUDE: ANALISE 2023 =========\n")
  print(aud)
  print(desc)
  cat("\nCR2 modelo M3 (macro):\n")
  print(resultados$M3$coef)
  cat("\nInteracao global HTZ: \n")
  print(resultados$M3$global)
  cat("\nSaidas: ", normalizePath(destino), "\n", sep = "")
  invisible(list(principal = d_macro, modelos = mod, auditoria = aud,
                 resultados = resultados))
}

# A sensibilidade por raca/cor da APAC indice altera explicitamente a
# exposicao racial para raca_indice; caso contrario repetiria a harmonizada.
neo_raca_indice_macro <- neo_raca_indice %>%
  mutate(raca_final = factor(raca_indice,
         levels = c("Branca", "Parda", "Preta")))

RODAR_LOO_MACRO <- FALSE
analises_macro_s2 <- macro_executar_analises(
  d_principal = neo_principal,
  coorte_pre = coorte_2023,
  sensibilidades = list(
    abril_dezembro_2023 = neo_abr_dez,
    sem_teto_365_2023 = neo_sem_teto,
    raca_APAC_indice_2023 = neo_raca_indice_macro
  ),
  raiz = RAIZ_MACRO,
  destino = file.path(DIR_SUPL, "ANALISE_MACRORREGIOES_S2"),
  rodar_loo = RODAR_LOO_MACRO,
  esperado_n = 15415L
)
