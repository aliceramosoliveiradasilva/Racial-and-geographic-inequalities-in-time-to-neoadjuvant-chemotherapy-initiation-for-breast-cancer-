# PUBLICATION REPRODUCIBILITY — BREAST CANCER NEOADJUVANT CHEMOTHERAPY
# Macroregion extension. Executes separately; outputs isolated.
# Executes locally on authorized DATASUS data; never commit patient-level data.
# 2026-10-15. Historical crosswalks 2024/2026 are not validated for 2023.
# All R scripts should be run from project root; MAMA_RAIZ configures data root.
# NOTE: full end-to-end re-execution was not possible in the delivery environment.
# -------------------------------------------------------------------------

# ============================================================================
# 02_analisar_mama_neoadjuvancia_2019_2025_ATUALIZADO.R
# Projeto: desigualdades raciais e territoriais no tempo até quimioterapia
#          neoadjuvante para câncer de mama no SUS
#
# REANALISE MACRORREGIONAL - 14/10/2026
# Preserva integralmente o script municipal e acrescenta modelos territoriais.
# Os resultados vao para "resultados_com_macrorregioes"; originais intactos.
# A classificacao territorial do CSV e de 2026: verificar validade em 2023.
# METODOLOGIA MUNICIPAL DE REFERENCIA v1 - 25/09/2026
#
# Análise principal:
#   - coorte nacional de episódios de quimioterapia neoadjuvante iniciados em 2023;
#   - 2019-2022 como lookback para identificar neoadjuvância prévia;
#   - estádios I-III;
#   - raça/cor harmonizada longitudinalmente usando registros até 2023;
#   - branca, parda e preta mantidas separadamente;
#   - exposição principal: mesmo município vs outro município;
#   - desfecho principal: início >60 dias após identificação patológica;
#   - interação principal pré-especificada: Preta x tratamento fora do município;
#   - Poisson com link log e inferência clusterizada por CNES;
#   - CR2/Satterthwaite e teste conjunto HTZ.
#
# Sensibilidades pré-especificadas:
#   1) abril-dezembro de 2023;
#   2) 2022-2023;
#   3) 2023-2024, reconhecendo perda de ID em 2024;
#   4) raça/cor da APAC índice;
#   5) retirada do teto de 365 dias, mantendo datas formalmente consistentes;
#   6) leave-one-CNES-out opcional.
#
# 2025 é mantido para auditoria de qualidade, mas NÃO integra a análise principal.
# ============================================================================

# -------------------------------
# 0. Pacotes
# -------------------------------
pacotes <- c(
  "dplyr", "stringr", "lubridate", "readr", "purrr", "tibble", "fst",
  "splines", "sandwich", "lmtest", "emmeans", "clubSandwich"
)

faltantes <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltantes) > 0) {
  stop(
    "Instale os pacotes antes de continuar: ",
    paste(faltantes, collapse = ", ")
  )
}

library(dplyr)
library(stringr)
library(lubridate)
library(readr)
library(purrr)
library(tibble)
library(fst)
library(splines)
library(sandwich)
library(lmtest)
library(emmeans)
library(clubSandwich)

# -------------------------------
# 1. Parâmetros e diretórios
# -------------------------------
# O script funciona se for executado:
#   a) a partir da pasta que contém dados_mama_SIA_AQ_2019_2025; OU
#   b) de dentro da própria pasta dados_mama_SIA_AQ_2019_2025.
# Caminhos absolutos confirmados na maquina Fedora da pesquisadora.
# A pasta de saida e NOVA para nao sobrescrever resultados do artigo original.
RAIZ_MACRO <- Sys.getenv("MAMA_RAIZ",unset=getwd())
PASTA_BASE <- file.path(RAIZ_MACRO, "dados_mama_SIA_AQ_2019_2025")
DIR_DADOS <- file.path(PASTA_BASE, "dados")
DIR_RESULTADOS <- file.path(PASTA_BASE, "resultados_com_macrorregioes")
dir.create(DIR_RESULTADOS, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(DIR_DADOS)) {
  stop("Pasta de dados não encontrada: ", DIR_DADOS)
}

ANO_PRINCIPAL <- 2023L
ANOS_LOOKBACK <- 2019:2022
ANOS_SENS_AMP <- 2022:2023
ANOS_SENS_CONTEMP <- 2023:2024
DATA_INICIO_SENS_ABRIL_2023 <- as.Date("2023-04-01")
DATA_FIM_2023 <- as.Date("2023-12-31")

# Leave-one-CNES-out nacional pode ser demorado.
RODAR_LOO <- FALSE

# Códigos considerados candidatos à quimioterapia neoadjuvante de mama.
# O script exporta a frequência observada para auditoria; códigos ausentes não
# geram registros artificialmente.
COD_NEO_QUIMIO <- c(
  "0304040029",
  "0304040185",
  "0304040215",
  "0304040223"
)

RACAS_VALIDAS <- c("01", "02", "03", "04", "05")
RACAS_MODELO <- c("Branca", "Parda", "Preta")

# -------------------------------
# 2. Funções auxiliares
# -------------------------------
limpa_cid <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[^A-Z0-9]", "", x)
  x[x == ""] <- NA_character_
  x
}

limpa_proc <- function(x) {
  x <- as.character(x)
  x <- gsub("[^0-9]", "", x)
  x[x == ""] <- NA_character_
  x
}

limpa_mun <- function(x) {
  x <- as.character(x)
  x <- gsub("[^0-9]", "", x)
  x[nchar(x) != 6] <- NA_character_
  x
}

parse_data_sia <- function(x) {
  x <- as.character(x)
  x_ok <- ifelse(
    str_detect(x, "^(19|20)[0-9]{6}$"),
    x,
    NA_character_
  )
  suppressWarnings(ymd(x_ok, quiet = TRUE))
}

recodifica_raca <- function(x) {
  case_when(
    x == "01" ~ "Branca",
    x == "02" ~ "Preta",
    x == "03" ~ "Parda",
    x == "04" ~ "Amarela",
    x == "05" ~ "Indígena",
    x == "INCONSISTENTE" ~ "Inconsistente",
    TRUE ~ "Sem informação"
  )
}

regiao_por_uf_ibge <- function(mun) {
  uf <- substr(mun, 1, 2)
  case_when(
    uf %in% c("11", "12", "13", "14", "15", "16", "17") ~ "Norte",
    uf %in% c("21", "22", "23", "24", "25", "26", "27", "28", "29") ~ "Nordeste",
    uf %in% c("31", "32", "33", "35") ~ "Sudeste",
    uf %in% c("41", "42", "43") ~ "Sul",
    uf %in% c("50", "51", "52", "53") ~ "Centro-Oeste",
    TRUE ~ NA_character_
  )
}

# Harmoniza raça/cor usando apenas registros até o ano definido.
# Isso impede que a classificação da coorte de 2023 dependa de registros futuros.
harmoniza_raca_ate <- function(dados_mama, ano_fim) {
  dados_mama %>%
    filter(
      ano_comp <= ano_fim,
      !is.na(id_paciente)
    ) %>%
    group_by(id_paciente) %>%
    summarise(
      n_racas_validas = n_distinct(raca_cor[raca_cor %in% RACAS_VALIDAS]),
      raca_final_cod = case_when(
        n_racas_validas == 0 ~ NA_character_,
        n_racas_validas == 1 ~ first(raca_cor[raca_cor %in% RACAS_VALIDAS]),
        TRUE ~ "INCONSISTENTE"
      ),
      .groups = "drop"
    ) %>%
    mutate(raca_final = recodifica_raca(raca_final_cod))
}

# Prepara a primeira APAC representativa de cada episódio de início.
# ordem_episodio permite identificar neoadjuvância prévia no lookback.
cria_episodios <- function(neo_registros) {
  neo_registros %>%
    filter(!is.na(id_paciente), !is.na(dt_inicio)) %>%
    arrange(id_paciente, dt_inicio, competencia) %>%
    distinct(id_paciente, dt_inicio, .keep_all = TRUE) %>%
    group_by(id_paciente) %>%
    arrange(dt_inicio, .by_group = TRUE) %>%
    mutate(
      ordem_episodio = row_number(),
      tem_neo_previa_observada = ordem_episodio > 1
    ) %>%
    ungroup()
}

# Seleciona um episódio índice por mulher dentro da janela especificada.
# Por padrão, exige que seja o primeiro episódio neoadjuvante observado desde 2019.
seleciona_coorte <- function(
    episodios,
    dados_mama,
    data_inicio,
    data_fim,
    ano_fim_raca,
    excluir_neo_previa = TRUE) {

  d <- episodios %>%
    filter(
      dt_inicio >= as.Date(data_inicio),
      dt_inicio <= as.Date(data_fim)
    ) %>%
    arrange(id_paciente, dt_inicio, competencia) %>%
    group_by(id_paciente) %>%
    slice_min(dt_inicio, n = 1, with_ties = FALSE) %>%
    ungroup()

  if (excluir_neo_previa) {
    d <- d %>% filter(!tem_neo_previa_observada)
  }

  raca_long <- harmoniza_raca_ate(dados_mama, ano_fim = ano_fim_raca)

  d %>%
    left_join(
      raca_long %>% select(id_paciente, n_racas_validas, raca_final),
      by = "id_paciente"
    ) %>%
    mutate(
      ano_inicio = year(dt_inicio),
      mes_inicio = month(dt_inicio),
      atraso60 = as.integer(tempo_dias > 60),
      fora_municipio = case_when(
        is.na(mun_res) | is.na(mun_estab) ~ NA_character_,
        mun_res == mun_estab ~ "Mesmo município",
        TRUE ~ "Outro município"
      ),
      fora_municipio = factor(
        fora_municipio,
        levels = c("Mesmo município", "Outro município")
      ),
      raca_final = factor(
        raca_final,
        levels = c(
          "Branca", "Parda", "Preta", "Amarela", "Indígena",
          "Inconsistente", "Sem informação"
        )
      )
    )
}

# População inferencial principal: datas consistentes, estádios I-III e
# branca/parda/preta com raça longitudinalmente consistente.
prepara_populacao <- function(d, retirar_teto_365 = FALSE) {

  if (retirar_teto_365) {
    d <- d %>%
      filter(
        !is.na(dt_diagnostico),
        !is.na(dt_inicio),
        dt_inicio <= comp_fim,
        dt_diagnostico <= dt_inicio
      )
  } else {
    d <- d %>%
      filter(problema_data == "Elegível análise temporal")
  }

  d %>%
    filter(
      estadio %in% c("1", "2", "3"),
      raca_final %in% RACAS_MODELO,
      !is.na(fora_municipio),
      !is.na(idade),
      !is.na(cnes),
      !is.na(uf_res),
      !is.na(mes_inicio)
    ) %>%
    mutate(
      raca_final = factor(raca_final, levels = RACAS_MODELO),
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
}

# Ajusta os quatro modelos sequenciais previstos na metodologia.
ajusta_modelos_principais <- function(dados) {
  m0 <- glm(
    atraso60 ~ raca_final * fora_municipio,
    family = poisson(link = "log"),
    data = dados
  )

  m1 <- glm(
    atraso60 ~
      raca_final * fora_municipio +
      ns(idade, df = 3),
    family = poisson(link = "log"),
    data = dados
  )

  m2 <- glm(
    atraso60 ~
      raca_final * fora_municipio +
      ns(idade, df = 3) +
      estadio +
      procedimento,
    family = poisson(link = "log"),
    data = dados
  )

  m3 <- glm(
    atraso60 ~
      raca_final * fora_municipio +
      ns(idade, df = 3) +
      estadio +
      procedimento +
      uf_res +
      mes_inicio,
    family = poisson(link = "log"),
    data = dados
  )

  list(M0 = m0, M1 = m1, M2 = m2, M3 = m3)
}

# Modelo final para sensibilidades multiano. Ajusta também por ano quando
# houver mais de um ano de início na população.
ajusta_modelo_final <- function(dados) {
  if (n_distinct(dados$ano_inicio) > 1) {
    formula_final <- as.formula(
      "atraso60 ~ raca_final * fora_municipio + ns(idade, df = 3) + estadio + procedimento + factor(ano_inicio) + uf_res + mes_inicio"
    )
  } else {
    formula_final <- as.formula(
      "atraso60 ~ raca_final * fora_municipio + ns(idade, df = 3) + estadio + procedimento + uf_res + mes_inicio"
    )
  }

  glm(
    formula_final,
    family = poisson(link = "log"),
    data = dados
  )
}

# CR2/Satterthwaite + teste global HTZ.
analisa_cr2 <- function(modelo, dados, nome_modelo) {
  V_CR2 <- vcovCR(
    modelo,
    cluster = dados$cnes,
    type = "CR2"
  )

  teste <- coef_test(
    modelo,
    vcov = V_CR2,
    test = "Satterthwaite"
  )

  teste_df <- as.data.frame(teste)
  write_csv(
    teste_df,
    file.path(DIR_RESULTADOS, paste0("coeficientes_CR2_", nome_modelo, ".csv"))
  )

  idx_interacao <- grep(
    "raca_final.*:fora_municipio",
    names(coef(modelo))
  )

  global <- Wald_test(
    modelo,
    constraints = constrain_zero(idx_interacao),
    vcov = V_CR2,
    test = "HTZ"
  )

  write_csv(
    as.data.frame(global),
    file.path(DIR_RESULTADOS, paste0("teste_global_interacao_HTZ_", nome_modelo, ".csv"))
  )

  list(V = V_CR2, coef = teste, global = global)
}

# -------------------------------
# 3. Leitura dos arquivos anuais 2019-2025
# -------------------------------
arquivos <- list.files(
  DIR_DADOS,
  pattern = "^SIA_AQ_MAMA_BRASIL_20(19|20|21|22|23|24|25)\\.fst$",
  full.names = TRUE
)

if (length(arquivos) == 0) {
  stop("Nenhum arquivo FST anual encontrado em: ", DIR_DADOS)
}

if (length(arquivos) < 7) {
  warning(
    "Foram encontrados ", length(arquivos),
    " arquivos anuais; esperavam-se 7 arquivos (2019-2025)."
  )
}

cat("Arquivos encontrados: ", length(arquivos), "\n")
print(basename(arquivos))

mama_raw <- map_dfr(
  arquivos,
  ~ read_fst(.x, as.data.table = FALSE) %>% as_tibble()
)

if (!"source" %in% names(mama_raw)) mama_raw$source <- NA_character_

# -------------------------------
# 4. Auditoria do identificador longitudinal antes da deduplicação
# -------------------------------
qualidade_id_base <- mama_raw %>%
  mutate(
    ano_comp = suppressWarnings(as.integer(substr(AP_MVM, 1, 4))),
    procedimento_audit = limpa_proc(AP_PRIPAL),
    id_ausente = is.na(ID_PACIENTE_HEX) |
      coalesce(trimws(as.character(ID_PACIENTE_HEX)) == "", FALSE),
    neo_audit = procedimento_audit %in% COD_NEO_QUIMIO,
    uf_estab_audit = substr(limpa_mun(AP_UFMUN), 1, 2)
  )

qualidade_id_anual <- qualidade_id_base %>%
  group_by(ano_comp) %>%
  summarise(
    registros_mama = n(),
    sem_id_mama = sum(id_ausente),
    pct_sem_id_mama = round(100 * mean(id_ausente), 2),
    registros_neo = sum(neo_audit),
    sem_id_neo = sum(id_ausente & neo_audit),
    pct_sem_id_neo = round(
      100 * sum(id_ausente & neo_audit) / sum(neo_audit),
      2
    ),
    .groups = "drop"
  )

write_csv(
  qualidade_id_anual,
  file.path(DIR_RESULTADOS, "qualidade_ID_anual_mama_e_neoadjuvancia_2019_2025.csv")
)

qualidade_id_mensal_neo <- qualidade_id_base %>%
  filter(neo_audit) %>%
  group_by(AP_MVM) %>%
  summarise(
    registros_neo = n(),
    sem_id = sum(id_ausente),
    pct_sem_id = round(100 * mean(id_ausente), 2),
    .groups = "drop"
  )

write_csv(
  qualidade_id_mensal_neo,
  file.path(DIR_RESULTADOS, "qualidade_ID_mensal_neoadjuvancia_2019_2025.csv")
)

qualidade_id_uf_2024_2025 <- qualidade_id_base %>%
  filter(
    ano_comp %in% c(2024, 2025),
    neo_audit
  ) %>%
  group_by(ano_comp, uf_estab_audit) %>%
  summarise(
    registros_neo = n(),
    sem_id = sum(id_ausente),
    pct_sem_id = round(100 * mean(id_ausente), 2),
    .groups = "drop"
  )

write_csv(
  qualidade_id_uf_2024_2025,
  file.path(DIR_RESULTADOS, "qualidade_ID_neoadjuvancia_por_UF_2024_2025.csv")
)

rm(qualidade_id_base)
invisible(gc())

# -------------------------------
# 5. Padronização da base de câncer de mama
# -------------------------------
mama <- mama_raw %>%
  transmute(
    id_paciente = case_when(
      is.na(ID_PACIENTE_HEX) ~ NA_character_,
      trimws(as.character(ID_PACIENTE_HEX)) == "" ~ NA_character_,
      TRUE ~ as.character(ID_PACIENTE_HEX)
    ),
    competencia = AP_MVM,
    ano_comp = suppressWarnings(as.integer(substr(AP_MVM, 1, 4))),
    autorizacao = AP_AUTORIZ,
    cnes = AP_CODUNI,
    procedimento = limpa_proc(AP_PRIPAL),

    codigo_idade = AP_COIDADE,
    idade = suppressWarnings(as.integer(AP_NUIDADE)),
    sexo = AP_SEXO,
    raca_cor = AP_RACACOR,

    mun_res = limpa_mun(AP_MUNPCN),
    mun_estab = limpa_mun(AP_UFMUN),
    uf_diferente = AP_UFDIF,
    mun_diferente_sistema = AP_MNDIF,

    tipo_apac = AP_TPAPAC,
    dt_inicio_apac = parse_data_sia(AP_DTINIC),
    dt_fim_apac = parse_data_sia(AP_DTFIM),
    motivo_saida = AP_MOTSAI,
    obito = AP_OBITO,

    cid_principal = limpa_cid(AP_CIDPRI),
    cid_oncologico = limpa_cid(AQ_CID10),
    estadio = as.character(AQ_ESTADI),
    linfonodos = AQ_LINFIN,
    grau_histologico = AQ_GRAHIS,

    tratamento_anterior = AQ_TRANTE,
    continuidade = AQ_CONTTR,
    cid_ini1 = limpa_cid(AQ_CIDINI1),
    dt_ini1 = parse_data_sia(AQ_DTINI1),
    cid_ini2 = limpa_cid(AQ_CIDINI2),
    dt_ini2 = parse_data_sia(AQ_DTINI2),
    cid_ini3 = limpa_cid(AQ_CIDINI3),
    dt_ini3 = parse_data_sia(AQ_DTINI3),

    dt_diag_raw = AQ_DTIDEN,
    dt_inicio_raw = AQ_DTINTR,
    dt_diagnostico = parse_data_sia(AQ_DTIDEN),
    dt_inicio = parse_data_sia(AQ_DTINTR),

    source = source
  ) %>%
  filter(
    sexo == "F",
    str_detect(cid_principal, "^C50") |
      str_detect(cid_oncologico, "^C50")
  ) %>%
  mutate(
    comp_inicio = ymd(paste0(competencia, "01")),
    comp_fim = ceiling_date(comp_inicio, "month") - days(1),
    uf_res = substr(mun_res, 1, 2),
    regiao = regiao_por_uf_ibge(mun_res)
  )

rm(mama_raw)
invisible(gc())

# -------------------------------
# 6. Auditoria da raça/cor até 2023 e até 2024
# -------------------------------
raca_ate_2023 <- harmoniza_raca_ate(mama, 2023)
raca_ate_2024 <- harmoniza_raca_ate(mama, 2024)

qualidade_raca_ate_2023 <- raca_ate_2023 %>%
  count(n_racas_validas, raca_final, name = "pacientes") %>%
  mutate(percentual = round(100 * pacientes / sum(pacientes), 2))

qualidade_raca_ate_2024 <- raca_ate_2024 %>%
  count(n_racas_validas, raca_final, name = "pacientes") %>%
  mutate(percentual = round(100 * pacientes / sum(pacientes), 2))

write_csv(
  qualidade_raca_ate_2023,
  file.path(DIR_RESULTADOS, "qualidade_raca_longitudinal_ate_2023.csv")
)
write_csv(
  qualidade_raca_ate_2024,
  file.path(DIR_RESULTADOS, "qualidade_raca_longitudinal_ate_2024.csv")
)

# -------------------------------
# 7. Identificação dos registros e episódios neoadjuvantes
# -------------------------------
neo_registros <- mama %>%
  filter(procedimento %in% COD_NEO_QUIMIO) %>%
  mutate(
    tempo_dias = as.integer(dt_inicio - dt_diagnostico),
    problema_data = case_when(
      is.na(dt_diagnostico) ~ "Data diagnóstico inválida",
      is.na(dt_inicio) ~ "Data início inválida",
      dt_inicio > comp_fim ~ "Início posterior à competência",
      dt_diagnostico > dt_inicio ~ "Diagnóstico após tratamento",
      tempo_dias > 365 ~ "Intervalo >365 dias",
      TRUE ~ "Elegível análise temporal"
    )
  )

# Frequência real dos procedimentos por ano para documentação.
procedimentos_neo_observados <- neo_registros %>%
  count(ano_comp, procedimento, name = "registros") %>%
  arrange(ano_comp, desc(registros))

write_csv(
  procedimentos_neo_observados,
  file.path(DIR_RESULTADOS, "procedimentos_neoadjuvantes_observados_por_ano.csv")
)

episodios_neo <- cria_episodios(neo_registros)

# Auditoria de episódios com tratamento prévio no lookback.
auditoria_episodios <- episodios_neo %>%
  mutate(ano_inicio = year(dt_inicio)) %>%
  group_by(ano_inicio) %>%
  summarise(
    episodios = n(),
    mulheres = n_distinct(id_paciente),
    com_neo_previa_n = sum(tem_neo_previa_observada),
    com_neo_previa_pct = round(100 * mean(tem_neo_previa_observada), 2),
    .groups = "drop"
  )

write_csv(
  auditoria_episodios,
  file.path(DIR_RESULTADOS, "auditoria_episodios_neoadjuvantes_e_lookback.csv")
)

# -------------------------------
# 8. Construção das coortes pré-especificadas
# -------------------------------
coorte_2023 <- seleciona_coorte(
  episodios = episodios_neo,
  dados_mama = mama,
  data_inicio = "2023-01-01",
  data_fim = "2023-12-31",
  ano_fim_raca = 2023,
  excluir_neo_previa = TRUE
)

coorte_abr_dez_2023 <- seleciona_coorte(
  episodios = episodios_neo,
  dados_mama = mama,
  data_inicio = DATA_INICIO_SENS_ABRIL_2023,
  data_fim = DATA_FIM_2023,
  ano_fim_raca = 2023,
  excluir_neo_previa = TRUE
)

coorte_2022_2023 <- seleciona_coorte(
  episodios = episodios_neo,
  dados_mama = mama,
  data_inicio = "2022-01-01",
  data_fim = "2023-12-31",
  ano_fim_raca = 2023,
  excluir_neo_previa = TRUE
)

coorte_2023_2024 <- seleciona_coorte(
  episodios = episodios_neo,
  dados_mama = mama,
  data_inicio = "2023-01-01",
  data_fim = "2024-12-31",
  ano_fim_raca = 2024,
  excluir_neo_previa = TRUE
)

# Salva a coorte de episódios 2023 antes das exclusões inferenciais.
write_fst(
  as.data.frame(coorte_2023),
  file.path(DIR_RESULTADOS, "coorte_episodios_neoadjuvantes_2023_pre_analise.fst"),
  compress = 50
)

# -------------------------------
# 9. Fluxo de seleção da coorte principal 2023
# -------------------------------
fluxo_2023 <- tibble(
  etapa = c(
    "Episódios neoadjuvantes com início em 2023 e ID disponível",
    "Primeiro episódio neoadjuvante observado desde 2019",
    "Datas elegíveis (0-365 dias e consistência temporal)",
    "Estádios I-III",
    "Raça/cor consistente: branca, parda ou preta",
    "Município, idade, CNES, UF e mês disponíveis"
  ),
  n = c(
    episodios_neo %>%
      filter(dt_inicio >= as.Date("2023-01-01"), dt_inicio <= as.Date("2023-12-31")) %>%
      summarise(n = n_distinct(id_paciente)) %>% pull(n),
    nrow(coorte_2023),
    coorte_2023 %>%
      filter(problema_data == "Elegível análise temporal") %>% nrow(),
    coorte_2023 %>%
      filter(
        problema_data == "Elegível análise temporal",
        estadio %in% c("1", "2", "3")
      ) %>% nrow(),
    coorte_2023 %>%
      filter(
        problema_data == "Elegível análise temporal",
        estadio %in% c("1", "2", "3"),
        raca_final %in% RACAS_MODELO
      ) %>% nrow(),
    coorte_2023 %>%
      filter(
        problema_data == "Elegível análise temporal",
        estadio %in% c("1", "2", "3"),
        raca_final %in% RACAS_MODELO,
        !is.na(mun_res), !is.na(mun_estab), !is.na(idade),
        !is.na(cnes), !is.na(uf_res), !is.na(mes_inicio)
      ) %>% nrow()
  )
)

write_csv(
  fluxo_2023,
  file.path(DIR_RESULTADOS, "fluxo_selecao_coorte_principal_2023.csv")
)

# -------------------------------
# 10. Populações analíticas
# -------------------------------
neo_principal <- prepara_populacao(coorte_2023)
neo_sens_abr_dez_2023 <- prepara_populacao(coorte_abr_dez_2023)
neo_sens_2022_2023 <- prepara_populacao(coorte_2022_2023)
neo_sens_2023_2024 <- prepara_populacao(coorte_2023_2024)
neo_sens_sem_teto_365 <- prepara_populacao(coorte_2023, retirar_teto_365 = TRUE)

if (nrow(neo_principal) == 0) {
  stop("A população analítica principal de 2023 ficou vazia. Revise os filtros.")
}

# -------------------------------
# 11. Tabelas descritivas da análise principal 2023
# -------------------------------
resumo_principal <- neo_principal %>%
  summarise(
    mulheres = n(),
    brancas = sum(raca_final == "Branca"),
    pardas = sum(raca_final == "Parda"),
    pretas = sum(raca_final == "Preta"),
    atraso60_n = sum(atraso60),
    atraso60_pct = round(100 * mean(atraso60), 1),
    mediana_tempo = median(tempo_dias),
    p25_tempo = quantile(tempo_dias, .25),
    p75_tempo = quantile(tempo_dias, .75),
    n_cnes = n_distinct(cnes),
    n_uf_res = n_distinct(uf_res)
  )

write_csv(
  resumo_principal,
  file.path(DIR_RESULTADOS, "resumo_coorte_principal_2023.csv")
)

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
  file.path(DIR_RESULTADOS, "tabela_raca_x_deslocamento_2023.csv")
)

# Tabela de composição clínica/demográfica por raça e deslocamento.
tabela_composicao <- neo_principal %>%
  group_by(raca_final, fora_municipio) %>%
  summarise(
    mulheres = n(),
    idade_mediana = median(idade, na.rm = TRUE),
    idade_p25 = quantile(idade, .25, na.rm = TRUE),
    idade_p75 = quantile(idade, .75, na.rm = TRUE),
    estadio_I_pct = round(100 * mean(estadio == "1"), 1),
    estadio_II_pct = round(100 * mean(estadio == "2"), 1),
    estadio_III_pct = round(100 * mean(estadio == "3"), 1),
    .groups = "drop"
  )

write_csv(
  tabela_composicao,
  file.path(DIR_RESULTADOS, "composicao_clinica_por_raca_deslocamento_2023.csv")
)

# -------------------------------
# 12. Modelos sequenciais da análise principal
# -------------------------------
modelos <- ajusta_modelos_principais(neo_principal)

cr2_m0 <- analisa_cr2(modelos$M0, neo_principal, "M0_2023_bruto")
cr2_m1 <- analisa_cr2(modelos$M1, neo_principal, "M1_2023_idade")
cr2_m2 <- analisa_cr2(modelos$M2, neo_principal, "M2_2023_clinico")
cr2_m3 <- analisa_cr2(modelos$M3, neo_principal, "M3_2023_UF_mes")

# -------------------------------
# 13. Riscos ajustados e contrastes do modelo final M3
# -------------------------------
# A inferência principal dos termos de interação permanece CR2/HTZ.
# Emmeans é usado aqui para apresentação interpretável de riscos/contrastes.
V3_CL <- vcovCL(
  modelos$M3,
  cluster = ~ cnes,
  type = "HC1"
)

# -------------------------------
# 13. Riscos ajustados e contrastes do modelo final M3
# -------------------------------

V3_CL <- vcovCL(
  modelos$M3,
  cluster = ~ cnes,
  type = "HC1"
)

# Fatores de ajuste que não participam de interações
emm3 <- emmeans(
  modelos$M3,
  ~ raca_final * fora_municipio,
  vcov. = V3_CL,
  nuisance = c(
    "uf_res",
    "mes_inicio",
    "estadio",
    "procedimento"
  ),
  wt.nuis = "proportional"
)

# Riscos ajustados
riscos_ajustados <- as.data.frame(
  summary(
    emm3,
    type = "response",
    infer = c(TRUE, TRUE)
  )
)

riscos_ajustados

contraste_raca_todos <- as.data.frame(
  summary(
    pairs(emm3, by = "fora_municipio", adjust = "holm"),
    type = "response",
    infer = c(TRUE, TRUE)
  )
)

contraste_desloc_todos <- as.data.frame(
  summary(
    pairs(emm3, by = "raca_final", adjust = "holm"),
    type = "response",
    infer = c(TRUE, TRUE)
  )
)

# ============================================================
# CONTRASTE PRIMÁRIO: Preta vs Branca dentro de cada situação
# territorial, usando o emm3 já construído
# ============================================================

# ============================================================
# CONTRASTE PRIMÁRIO: PRETA vs BRANCA
# dentro de cada situação territorial
# ============================================================

contraste_preta_branca <- as.data.frame(
  summary(
    contrast(
      emm3,
      method = list(
        "Preta_vs_Branca" = c(-1, 0, 1)
      ),
      by = "fora_municipio",
      adjust = "none"
    ),
    type = "response",
    infer = c(TRUE, TRUE)
  )
)

print(contraste_preta_branca)


write_csv(
  riscos_ajustados,
  file.path(DIR_RESULTADOS, "riscos_ajustados_M3_2023.csv")
)
write_csv(
  contraste_raca_todos,
  file.path(DIR_RESULTADOS, "contrastes_raca_dentro_deslocamento_M3_2023.csv")
)
write_csv(
  contraste_desloc_todos,
  file.path(DIR_RESULTADOS, "contrastes_deslocamento_dentro_raca_M3_2023.csv")
)
write_csv(
  contraste_preta_branca,
  file.path(
    DIR_RESULTADOS,
    "contraste_PRIMARIO_preta_vs_branca_por_deslocamento_2023.csv"
  )
)

# -------------------------------
# 14. Sensibilidades temporais e de definição
# -------------------------------
roda_sensibilidade <- function(dados, nome) {
  if (nrow(dados) == 0) {
    warning("Sensibilidade vazia: ", nome)
    return(NULL)
  }

  modelo <- ajusta_modelo_final(dados)
  res <- analisa_cr2(modelo, dados, nome)

  resumo <- dados %>%
    summarise(
      mulheres = n(),
      atraso60_n = sum(atraso60),
      atraso60_pct = round(100 * mean(atraso60), 1),
      mediana_tempo = median(tempo_dias),
      n_cnes = n_distinct(cnes),
      n_uf = n_distinct(uf_res)
    )

  write_csv(
    resumo,
    file.path(DIR_RESULTADOS, paste0("resumo_", nome, ".csv"))
  )

  list(modelo = modelo, cr2 = res, resumo = resumo)
}

sens_abr_dez_2023 <- roda_sensibilidade(
  neo_sens_abr_dez_2023,
  "sensibilidade_abr_dez_2023"
)

sens_2022_2023 <- roda_sensibilidade(
  neo_sens_2022_2023,
  "sensibilidade_2022_2023"
)

sens_2023_2024 <- roda_sensibilidade(
  neo_sens_2023_2024,
  "sensibilidade_2023_2024"
)

sens_sem_teto_365 <- roda_sensibilidade(
  neo_sens_sem_teto_365,
  "sensibilidade_sem_teto_365_dias_2023"
)

# -------------------------------
# 15. Sensibilidade: raça/cor da APAC índice
# -------------------------------
neo_raca_indice <- coorte_2023 %>%
  mutate(
    raca_indice = recodifica_raca(raca_cor),
    raca_final = factor(raca_indice, levels = RACAS_MODELO)
  ) %>%
  prepara_populacao()

if (nrow(neo_raca_indice) > 0) {
  m_raca_indice <- ajusta_modelo_final(neo_raca_indice)
  cr2_raca_indice <- analisa_cr2(
    m_raca_indice,
    neo_raca_indice,
    "sensibilidade_raca_APAC_indice_2023"
  )
}

# -------------------------------
# 16. Diagnóstico dos clusters CNES
# -------------------------------
diagnostico_cnes <- neo_principal %>%
  count(cnes, sort = TRUE, name = "mulheres")

write_csv(
  diagnostico_cnes,
  file.path(DIR_RESULTADOS, "distribuicao_mulheres_por_CNES_2023.csv")
)

resumo_cnes <- tibble(
  n_cnes = n_distinct(neo_principal$cnes),
  min_por_cnes = min(diagnostico_cnes$mulheres),
  mediana_por_cnes = median(diagnostico_cnes$mulheres),
  max_por_cnes = max(diagnostico_cnes$mulheres)
)

write_csv(
  resumo_cnes,
  file.path(DIR_RESULTADOS, "resumo_clusters_CNES_2023.csv")
)

# -------------------------------
# 17. Leave-one-CNES-out opcional
# -------------------------------
if (RODAR_LOO) {
  cnes_lista <- unique(neo_principal$cnes)
  termo_alvo <- "raca_finalPreta:fora_municipioOutro município"

  loo <- map_dfr(
    cnes_lista,
    function(cnes_excluido) {
      d <- neo_principal %>% filter(cnes != cnes_excluido)

      mod <- ajusta_modelo_final(d)
      beta <- unname(coef(mod)[termo_alvo])

      tibble(
        cnes_excluido = cnes_excluido,
        n = nrow(d),
        beta = beta,
        RR_interacao = exp(beta)
      )
    }
  )

  write_csv(
    loo,
    file.path(DIR_RESULTADOS, "leave_one_CNES_out_2023.csv")
  )
}

# -------------------------------
# 18. Auditoria específica da perda de ID em 2024-2025
# -------------------------------
# Mantém documentado por que 2025 não compõe a análise principal.
audit_id_2024_2025 <- qualidade_id_anual %>%
  filter(ano_comp %in% c(2024, 2025))

write_csv(
  audit_id_2024_2025,
  file.path(DIR_RESULTADOS, "justificativa_exclusao_2025_qualidade_ID.csv")
)

# -------------------------------
# 19. Resumo no console
# -------------------------------
cat("\n================ ANÁLISE PRINCIPAL 2023 ================\n")
print(resumo_principal)

cat("\nRaça/cor x deslocamento - 2023:\n")
print(tabela_raca_desloc)

cat("\nCR2 - M3 principal 2023:\n")
print(cr2_m3$coef)

cat("\nTeste global HTZ - M3 principal 2023:\n")
print(cr2_m3$global)

cat("\nQualidade do ID por ano:\n")
print(qualidade_id_anual)

cat("\nResultados exportados em: ", normalizePath(DIR_RESULTADOS), "\n")
cat("Concluído.\n")

print(cr2_m3$coef)
print(cr2_m3$global)

write_csv(
  as.data.frame(cr2_m3$coef),
  file.path(DIR_RESULTADOS, "coeficientes_CR2_M3_PRINCIPAL_2023.csv")
)

write_csv(
  as.data.frame(cr2_m3$global),
  file.path(DIR_RESULTADOS, "teste_global_HTZ_M3_PRINCIPAL_2023.csv")
)
neo_principal_outro <- neo_principal %>%
  mutate(
    fora_municipio = relevel(
      fora_municipio,
      ref = "Outro município"
    )
  )

m3_outro <- glm(
  atraso60 ~
    raca_final * fora_municipio +
    ns(idade, df = 3) +
    estadio +
    procedimento +
    uf_res +
    mes_inicio,
  family = poisson(link = "log"),
  data = neo_principal_outro
)

V3_outro_CR2 <- vcovCR(
  m3_outro,
  cluster = neo_principal_outro$cnes,
  type = "CR2"
)

teste_outro_CR2 <- coef_test(
  m3_outro,
  vcov = V3_outro_CR2,
  test = "Satterthwaite"
)

teste_outro_CR2[
  grepl("raca_finalPreta", teste_outro_CR2$Coef),
]

# ---------------------------------------------------------------------------
# 20. Analises adicionais: macrorregioes de saude
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

# Usa a populacao analitica JA CONSTRUIDA e os mesmos criterios de elegibilidade.
# Nao reexecuta a mineracao na etapa territorial nem agrega APAC por municipio.
RODAR_LOO_MACRO <- FALSE

analises_macro <- macro_executar_analises(
  d_principal = neo_principal,
  coorte_pre = coorte_2023,
  sensibilidades = list(
    abril_dezembro_2023 = neo_sens_abr_dez_2023,
    anos_2022_2023 = neo_sens_2022_2023,
    anos_2023_2024 = neo_sens_2023_2024,
    sem_teto_365_2023 = neo_sens_sem_teto_365,
    raca_APAC_indice_2023 = neo_raca_indice
  ),
  raiz = RAIZ_MACRO,
  destino = file.path(DIR_RESULTADOS, "ANALISE_MACRORREGIOES"),
  rodar_loo = RODAR_LOO_MACRO,
  esperado_n = 15415L
)
