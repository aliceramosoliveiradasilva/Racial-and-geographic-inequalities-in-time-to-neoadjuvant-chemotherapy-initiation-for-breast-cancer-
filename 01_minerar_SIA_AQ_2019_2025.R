# PUBLICATION REPRODUCIBILITY — BREAST CANCER NEOADJUVANT CHEMOTHERAPY
# Original SIA-AQ mining (microdatasus); paths parameterized.
# Executes locally on authorized DATASUS data; never commit patient-level data.
# 2026-10-15. Historical crosswalks 2024/2026 are not validated for 2023.
# All R scripts should be run from project root; MAMA_RAIZ configures data root.
# NOTE: full end-to-end re-execution was not possible in the delivery environment.
# -------------------------------------------------------------------------

# ============================================================================
# 01_extrair_mama_SIA_AQ_2019_2025.R
# Projeto: desigualdades raciais e territoriais no acesso à quimioterapia
#          neoadjuvante para câncer de mama no SUS
# Objetivo deste script:
#   1) baixar SIA-AQ nacional, 2019-2025, ano a ano;
#   2) manter somente mulheres com câncer de mama (C50.x);
#   3) converter o identificador criptografado do usuário para hexadecimal ASCII;
#   4) exportar um arquivo FST por ano;
#   5) auditar completude/consistência de raça/cor e variáveis centrais.
#
# Observação: NÃO aplicar enc2utf8() a AP_CNSPCN. O campo é criptografado e pode
# chegar marcado como bytes. Ele é convertido para uma representação hexadecimal
# estável antes da exportação.
# ============================================================================

# -------------------------------
# 0. Pacotes
# -------------------------------
pacotes <- c(
  "microdatasus", "dplyr", "stringr", "lubridate",
  "readr", "purrr", "tibble", "fst"
)

faltantes <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltantes) > 0) {
  stop(
    "Instale os pacotes antes de continuar: ",
    paste(faltantes, collapse = ", ")
  )
}

library(microdatasus)
library(dplyr)
library(stringr)
library(lubridate)
library(readr)
library(purrr)
library(tibble)
library(fst)

cat("microdatasus: ", as.character(packageVersion("microdatasus")), "\n")

# -------------------------------
# 1. Parâmetros
# -------------------------------
ANOS <- 2019:2025
PASTA_SAIDA <- file.path(Sys.getenv("MAMA_RAIZ", unset = getwd()), "dados_mama_SIA_AQ_2019_2025")
SOBRESCREVER <- FALSE
TIMEOUT <- 900

DIR_DADOS <- file.path(PASTA_SAIDA, "dados")
DIR_QUALIDADE <- file.path(PASTA_SAIDA, "qualidade")
DIR_LOG <- file.path(PASTA_SAIDA, "logs")

dir.create(PASTA_SAIDA, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_DADOS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_QUALIDADE, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_LOG, showWarnings = FALSE, recursive = TRUE)

# Campos necessários para a análise atual + auditoria e futuras sensibilidades.
# Mantemos dados de mama em geral, e NÃO somente neoadjuvância.
VARS_AQ <- c(
  "AP_MVM", "AP_CODUNI", "AP_AUTORIZ", "AP_CNSPCN",
  "AP_COIDADE", "AP_NUIDADE", "AP_SEXO", "AP_RACACOR",
  "AP_MUNPCN", "AP_UFDIF", "AP_MNDIF", "AP_UFMUN",
  "AP_DTINIC", "AP_DTFIM", "AP_TPAPAC", "AP_MOTSAI", "AP_OBITO",
  "AP_PRIPAL", "AP_CIDPRI", "AP_CIDSEC",
  "AQ_CID10", "AQ_LINFIN", "AQ_ESTADI", "AQ_GRAHIS", "AQ_DTIDEN",
  "AQ_TRANTE", "AQ_CIDINI1", "AQ_DTINI1", "AQ_CIDINI2", "AQ_DTINI2",
  "AQ_CIDINI3", "AQ_DTINI3", "AQ_CONTTR", "AQ_DTINTR"
)

# -------------------------------
# 2. Funções auxiliares
# -------------------------------

# Converte bytes do identificador criptografado em hexadecimal ASCII.
id_para_hex <- function(x) {
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x) & x != ""

  out[ok] <- vapply(
    x[ok],
    function(z) paste(sprintf("%02x", as.integer(charToRaw(z))), collapse = ""),
    FUN.VALUE = character(1)
  )

  out
}

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

# Parser estrito. Não tenta corrigir datas digitadas incorretamente.
parse_data_sia <- function(x) {
  x <- as.character(x)
  x_ok <- ifelse(
    str_detect(x, "^(19|20)[0-9]{6}$"),
    x,
    NA_character_
  )
  suppressWarnings(ymd(x_ok, quiet = TRUE))
}

pct <- function(x) {
  if (length(x) == 0) return(NA_real_)
  round(100 * mean(x, na.rm = TRUE), 2)
}

# Auditoria anual por REGISTRO.
calcula_qualidade_anual <- function(d, ano) {
  limite_sup <- as.Date(sprintf("%d-12-31", ano))

  d2 <- d %>%
    mutate(
      idade_num = suppressWarnings(as.integer(AP_NUIDADE)),
      dt_diag = parse_data_sia(AQ_DTIDEN),
      dt_inicio = parse_data_sia(AQ_DTINTR),
      raca_valida = AP_RACACOR %in% c("01", "02", "03", "04", "05"),
      estadio_preenchido = !is.na(AQ_ESTADI) & AQ_ESTADI != "",
      estadio_0a4 = AQ_ESTADI %in% as.character(0:4),
      idade_anos_valida = AP_COIDADE == "4" & !is.na(idade_num) &
        idade_num >= 18 & idade_num <= 120,
      mun_res_valido = str_detect(as.character(AP_MUNPCN), "^[0-9]{6}$"),
      mun_estab_valido = str_detect(as.character(AP_UFMUN), "^[0-9]{6}$"),
      cnes_valido = str_detect(as.character(AP_CODUNI), "^[0-9]{7}$"),
      diag_formato_valido = !is.na(dt_diag),
      inicio_formato_valido = !is.na(dt_inicio),
      diag_plausivel = !is.na(dt_diag) & dt_diag >= as.Date("1990-01-01") &
        dt_diag <= limite_sup,
      inicio_plausivel = !is.na(dt_inicio) & dt_inicio >= as.Date("1990-01-01") &
        dt_inicio <= limite_sup
    )

  tibble(
    ano = ano,
    registros_mama = nrow(d2),
    pacientes = n_distinct(d2$ID_PACIENTE_HEX, na.rm = TRUE),
    raca_valida_pct = pct(d2$raca_valida),
    estadio_preenchido_pct = pct(d2$estadio_preenchido),
    estadio_0a4_pct = pct(d2$estadio_0a4),
    idade_anos_valida_pct = pct(d2$idade_anos_valida),
    mun_res_valido_pct = pct(d2$mun_res_valido),
    mun_estab_valido_pct = pct(d2$mun_estab_valido),
    cnes_valido_pct = pct(d2$cnes_valido),
    dt_diagnostico_formato_valido_pct = pct(d2$diag_formato_valido),
    dt_inicio_formato_valido_pct = pct(d2$inicio_formato_valido),
    dt_diagnostico_plausivel_pct = pct(d2$diag_plausivel),
    dt_inicio_plausivel_pct = pct(d2$inicio_plausivel)
  )
}

# Auditoria mensal, sobretudo útil para documentar a mudança de completude de raça/cor.
calcula_qualidade_mensal <- function(d) {
  d %>%
    mutate(
      competencia = as.character(AP_MVM),
      raca_valida = AP_RACACOR %in% c("01", "02", "03", "04", "05"),
      estadio_preenchido = !is.na(AQ_ESTADI) & AQ_ESTADI != "",
      dt_diag_valida = !is.na(parse_data_sia(AQ_DTIDEN)),
      dt_inicio_valida = !is.na(parse_data_sia(AQ_DTINTR))
    ) %>%
    group_by(competencia) %>%
    summarise(
      registros = n(),
      pacientes = n_distinct(ID_PACIENTE_HEX),
      raca_valida_pct = pct(raca_valida),
      estadio_preenchido_pct = pct(estadio_preenchido),
      dt_diagnostico_valida_pct = pct(dt_diag_valida),
      dt_inicio_valida_pct = pct(dt_inicio_valida),
      .groups = "drop"
    )
}

# Qualidade da raça/cor por paciente dentro de cada ano.
calcula_raca_paciente_ano <- function(d, ano) {
  d %>%
    group_by(ID_PACIENTE_HEX) %>%
    summarise(
      n_racas_validas = n_distinct(
        AP_RACACOR[AP_RACACOR %in% c("01", "02", "03", "04", "05")]
      ),
      .groups = "drop"
    ) %>%
    count(n_racas_validas, name = "pacientes") %>%
    mutate(
      ano = ano,
      classificacao = case_when(
        n_racas_validas == 0 ~ "Sem raça/cor válida",
        n_racas_validas == 1 ~ "Uma raça/cor válida",
        n_racas_validas > 1 ~ "Raça/cor inconsistente",
        TRUE ~ NA_character_
      ),
      percentual = round(100 * pacientes / sum(pacientes), 2)
    ) %>%
    select(ano, classificacao, n_racas_validas, pacientes, percentual)
}

# -------------------------------
# 3. Download, filtro e exportação
# -------------------------------
qualidade_anual_lista <- list()
qualidade_mensal_lista <- list()
raca_paciente_lista <- list()
fontes_lista <- list()
procedimentos_lista <- list()

for (ano in ANOS) {

  cat("\n============================================================\n")
  cat("ANO: ", ano, "\n")
  cat("============================================================\n")

  arq_saida <- file.path(
    DIR_DADOS,
    sprintf("SIA_AQ_MAMA_BRASIL_%d.fst", ano)
  )

  if (file.exists(arq_saida) && !SOBRESCREVER) {
    cat("Arquivo já existe. Reutilizando: ", arq_saida, "\n")

    mama_ano <- read_fst(
      arq_saida,
      as.data.table = FALSE
    ) %>%
      as_tibble()

  } else {

    cat("Baixando SIA-AQ nacional...\n")

    dados_raw <- fetch_datasus(
      year_start = ano,
      month_start = 1,
      year_end = ano,
      month_end = 12,
      uf = "all",
      information_system = "SIA-AQ",
      vars = VARS_AQ,
      track_source = TRUE,
      stop_on_error = FALSE,
      timeout = TIMEOUT,
      quiet = FALSE
    )

    if (is.null(dados_raw) || nrow(dados_raw) == 0) {
      warning("Nenhum registro retornado para ", ano, ". Ano ignorado.")
      next
    }

    # Registra quais DBCs efetivamente contribuíram com dados antes do filtro.
    if ("source" %in% names(dados_raw)) {
      fontes_lista[[as.character(ano)]] <- dados_raw %>%
        distinct(source) %>%
        mutate(ano_solicitado = ano, .before = 1)
    }

    # Filtra antes de converter o identificador; isso reduz bastante o trabalho.
    mama_ano <- dados_raw %>%
      mutate(
        CID_PRI_LIMPO = limpa_cid(AP_CIDPRI),
        CID_AQ_LIMPO = limpa_cid(AQ_CID10)
      ) %>%
      filter(
        AP_SEXO == "F",
        str_detect(CID_PRI_LIMPO, "^C50") |
          str_detect(CID_AQ_LIMPO, "^C50")
      ) %>%
      mutate(
        ID_PACIENTE_HEX = id_para_hex(AP_CNSPCN),
        ANO_ARQUIVO = as.character(ano),
        AP_PRIPAL = limpa_proc(AP_PRIPAL)
      ) %>%
      select(
        ANO_ARQUIVO,
        ID_PACIENTE_HEX,
        everything(),
        -AP_CNSPCN,
        -CID_PRI_LIMPO,
        -CID_AQ_LIMPO
      )

    cat("Registros de mama selecionados: ", nrow(mama_ano), "\n")
    cat("Pacientes distintos: ", n_distinct(mama_ano$ID_PACIENTE_HEX), "\n")

    write_fst(as.data.frame(mama_ano), arq_saida, compress = 50)
    cat("Exportado: ", arq_saida, "\n")

    rm(dados_raw)
    invisible(gc())
  }

  # Se o arquivo reutilizado já traz source, aproveita também no log.
  if (is.null(fontes_lista[[as.character(ano)]]) && "source" %in% names(mama_ano)) {
    fontes_lista[[as.character(ano)]] <- mama_ano %>%
      distinct(source) %>%
      mutate(ano_solicitado = ano, .before = 1)
  }

  qualidade_anual_lista[[as.character(ano)]] <- calcula_qualidade_anual(mama_ano, ano)
  qualidade_mensal_lista[[as.character(ano)]] <- calcula_qualidade_mensal(mama_ano)
  raca_paciente_lista[[as.character(ano)]] <- calcula_raca_paciente_ano(mama_ano, ano)

  procedimentos_lista[[as.character(ano)]] <- mama_ano %>%
    count(AP_PRIPAL, sort = TRUE, name = "registros") %>%
    mutate(ano = ano, .before = 1)

  rm(mama_ano)
  invisible(gc())
}

# -------------------------------
# 4. Consolidação dos relatórios de qualidade
# -------------------------------
qualidade_anual <- bind_rows(qualidade_anual_lista)
qualidade_mensal <- bind_rows(qualidade_mensal_lista)
qualidade_raca_paciente <- bind_rows(raca_paciente_lista)
fontes_carregadas <- bind_rows(fontes_lista)
procedimentos_mama <- bind_rows(procedimentos_lista)

write_csv(
  qualidade_anual,
  file.path(DIR_QUALIDADE, "qualidade_completude_anual_2019_2025.csv")
)

write_csv(
  qualidade_mensal,
  file.path(DIR_QUALIDADE, "qualidade_completude_mensal_2019_2025.csv")
)

write_csv(
  qualidade_raca_paciente,
  file.path(DIR_QUALIDADE, "qualidade_raca_por_paciente_ano_2019_2025.csv")
)

write_csv(
  fontes_carregadas,
  file.path(DIR_LOG, "fontes_DBC_carregadas_2019_2025.csv")
)

write_csv(
  procedimentos_mama,
  file.path(DIR_QUALIDADE, "procedimentos_mama_por_ano_2019_2025.csv")
)

# -------------------------------
# 5. Resumo no console
# -------------------------------
cat("\n================ RESUMO FINAL ================\n")
print(qualidade_anual)

cat("\nCompletude mensal de raça/cor em 2023:\n")
print(
  qualidade_mensal %>%
    filter(substr(competencia, 1, 4) == "2023") %>%
    select(competencia, registros, pacientes, raca_valida_pct)
)

cat("\nArquivos salvos em: ", normalizePath(PASTA_SAIDA), "\n")
cat("Concluído.\n")
