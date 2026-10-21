# PUBLICATION REPRODUCIBILITY — BREAST CANCER NEOADJUVANT CHEMOTHERAPY
# Finite-degree CR2/HTZ and influence tests.
# Executes locally on authorized DATASUS data; never commit patient-level data.
# 2026-10-15. Historical crosswalks 2024/2026 are not validated for 2023.
# All R scripts should be run from project root; MAMA_RAIZ configures data root.
# NOTE: full end-to-end re-execution was not possible in the delivery environment.
# -------------------------------------------------------------------------

# ==============================================================================
# 03 - VALIDACAO FINAL: CANCER DE MAMA, APAC, MACRORREGIOES DE SAUDE
# Versao: 14/10/2026 | Autocontida, sem reler os FST anuais de producao
#
# Objetivos:
# 1) Reestimar M3 macrorregional e municipal na MESMA amostra, com CR2/HTZ;
# 2) calcular IC95% com graus de liberdade finitos para TODOS os contrastes;
# 3) registrar familia de testes e ajustes de Holm SEM reinterpretar interacao;
# 4) avaliar concentracao por CNES e leave-one-CNES-out (LOO);
# 5) auditar selecao da coorte, perdas e extremos do intervalo temporal;
# 6) testar (se fornecida) classificacao territorial DOCUMENTADA de 2023;
# 7) produzir tabelas de auditoria, sem sobrescrever nenhum resultado anterior.
#
# LIMITACOES QUE O SCRIPT NAO PODE CORRIGIR SOZINHO:
# - a tabela de outubro/2026 nao comprova o mapa vigente em 2023;
# - deslocamento territorial nao mede distancia, tempo de viagem ou TFD;
# - a coorte inclui somente mulheres que iniciaram a quimioterapia;
# - variaveis nao observadas nao sao corrigidas por ajuste estatistico.
#
# COMO USAR: abra este arquivo no RStudio, altere somente o bloco CONFIGURACAO
# se necessario e clique Source. Requer os FST de coorte gerados pelo script 02.
# ==============================================================================

# 0. CONFIGURACAO --------------------------------------------------------------
RAIZ <- Sys.getenv("MAMA_RAIZ", unset = getwd())
BASE <- file.path(RAIZ, "dados_mama_SIA_AQ_2019_2025")
DIR_MACROS <- file.path(BASE, "resultados_com_macrorregioes", "ANALISE_MACRORREGIOES")
ARQ_ANALITICA <- file.path(DIR_MACROS, "coorte_analitica_2023_COM_MACROS.fst")
ARQ_PRE <- file.path(DIR_MACROS, "coorte_episodios_2023_COM_MACROS_pre_analise.fst")

# Deixe em branco se AINDA NAO tiver correspondencia validada para 2023.
# Aceita CSV com as colunas: cod_municipio e cod_macrorregiao_de_saude.
# NAO utilize a tabela de 2026 renomeada como se fosse de 2023.
ARQ_MAPA_2023 <- file.path(RAIZ, "municipio_macrorregiao_REFERENCIA_2023.csv")
FONTE_DOCUMENTAL_2023 <- ""  # ex.: resolucao CIB + URL oficial e data de referencia
MAPA_2023_VALIDADO <- FALSE  # TRUE somente apos verificacao documental
ARQ_MAPA_2026 <- file.path(RAIZ, "macroregiao_de_saude(1).csv")

SAIDA <- file.path(DIR_MACROS, "VALIDACAO_FINAL_20261014")
RODAR_LOO_CNES <- TRUE        # reestima coeficientes apos remover cada CNES
RODAR_HTZ_LOO_TODOS <- FALSE  # TRUE: calcula CR2/HTZ para TODOS os LOO;
                              # FALSE: recalcula CR2/HTZ em cenarios extremos
QTD_EXTREMOS_LOO <- 10L      # extremos + maiores CNES + maiores clusters pretos

# 1. PACOTES E UTILITARIOS -----------------------------------------------------
pacotes <- c("dplyr", "readr", "fst", "emmeans", "clubSandwich", "splines", "tibble")
faltam <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltam)) stop("Instale os pacotes ausentes: install.packages(c('",
                        paste(faltam, collapse = "','"), "'))")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(fst); library(emmeans)
  library(clubSandwich); library(splines); library(tibble)
})
if (!file.exists(ARQ_ANALITICA)) stop("Coorte analitica nao encontrada: ", ARQ_ANALITICA)
if (!file.exists(ARQ_PRE)) stop("Coorte pre-analise nao encontrada: ", ARQ_PRE)
dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)

escrever <- function(obj, nome) {
  readr::write_csv(as.data.frame(obj), file.path(SAIDA, paste0(nome, ".csv")), na = "")
}
registrar <- function(...) cat(..., "\n", sep = "")
checagem_colunas <- function(d, nomes, origem) {
  ausentes <- setdiff(nomes, names(d))
  if (length(ausentes)) stop(origem, ": colunas ausentes: ", paste(ausentes, collapse = ", "))
}
normalizar_mun6 <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  x[x %in% c("", "NA", "NaN", "NULL", "000000", "0000000")] <- NA_character_
  x[!is.na(x) & !grepl("^[0-9]{6,7}$", x)] <- NA_character_
  substr(x, 1L, 6L)
}
seguro_div <- function(a, b) ifelse(b > 0, a / b, NA_real_)

registrar("Iniciando validacao. Diretorio de saida: ", SAIDA)
registrar("R version: ", R.version.string)

# 2. COORTES -------------------------------------------------------------------
d_todas <- tibble::as_tibble(fst::read_fst(ARQ_ANALITICA))
pre <- tibble::as_tibble(fst::read_fst(ARQ_PRE))
checagem_colunas(d_todas, c("raca_final", "fora_macro", "fora_municipio",
 "idade", "estadio", "procedimento", "uf_res", "mes_inicio",
 "cnes", "atraso60", "tempo_dias", "mun_res", "mun_estab"), "Coorte analitica")
checagem_colunas(pre, c("problema_data", "estadio", "raca_final", "fora_municipio",
 "idade", "cnes", "uf_res", "mes_inicio", "tempo_dias", "fora_macro"), "Coorte pre-analise")

n_coorte <- nrow(d_todas)
if (n_coorte != 15415L) warning("N da coorte original: ", n_coorte, " (esperado 15.415).")
registrar("Coorte analitica original: ", n_coorte, " registros")

# O unico NA na macrorregiao (se persistir) e excluido APENAS da analise macro.
d <- d_todas %>%
  filter(!is.na(fora_macro)) %>%
  mutate(
    raca_final = factor(as.character(raca_final), levels = c("Branca", "Parda", "Preta")),
    fora_macro = factor(as.character(fora_macro),
                        levels = c("Mesma macrorregiao", "Outra macrorregiao")),
    fora_municipio = factor(as.character(fora_municipio),
                            levels = c("Mesmo município", "Outro município")),
    estadio = factor(as.character(estadio)),
    procedimento = factor(as.character(procedimento)),
    uf_res = factor(as.character(uf_res)),
    mes_inicio = factor(as.character(mes_inicio)),
    cnes = as.character(cnes),
    atraso60 = as.integer(atraso60)
  ) %>% droplevels()

model_vars <- c("raca_final", "fora_macro", "fora_municipio", "idade", "estadio",
                "procedimento", "uf_res", "mes_inicio", "cnes", "atraso60", "tempo_dias")
if (any(!stats::complete.cases(d[, model_vars]))) {
  faltante <- colSums(is.na(d[, model_vars]))
  print(faltante[faltante > 0]); stop("Ha variaveis ausentes na amostra modelo. Corrija antes de ajustar.")
}
if (any(!d$atraso60 %in% c(0L, 1L))) stop("atraso60 nao e binaria 0/1.")
if (nrow(d) < 10000L) stop("Coorte muito menor do que a esperada: verificar FST.")
if (length(unique(d$cnes)) < 20L) stop("Numero de CNES insuficiente/inesperado.")
if (nlevels(d$raca_final) != 3L || nlevels(d$fora_macro) != 2L ||
    nlevels(d$fora_municipio) != 2L) stop("Categorias nao correspondem ao projeto.")

# Spline definida UMA vez e reutilizada no LOO: exclusoes nao alteram os nos.
mat_spline <- splines::ns(d$idade, df = 3)
col_ns <- paste0("idade_ns", seq_len(ncol(mat_spline)))
for (j in seq_along(col_ns)) d[[col_ns[j]]] <- as.numeric(mat_spline[, j])
escrever(tibble(parametro = c("limite_inferior", "limite_superior",
                            paste0("no_interno_", seq_along(attr(mat_spline, "knots")))),
                 valor = c(attr(mat_spline, "Boundary.knots"), attr(mat_spline, "knots"))),
         "SPLINE_NOS_FIXOS_LOO")

# 3. SELECAO E EXTREMOS -------------------------------------------------------
pre <- pre %>% mutate(
  raca_chr = ifelse(is.na(raca_final), "NA", as.character(raca_final)),
  macro_chr = ifelse(is.na(fora_macro), "Nao classificada", as.character(fora_macro)),
  municipio_chr = ifelse(is.na(fora_municipio), "Nao identificado", as.character(fora_municipio)),
  etapa_data = !is.na(problema_data) & problema_data == "Elegível análise temporal",
  etapa_estadio = as.character(estadio) %in% c("1", "2", "3"),
  etapa_raca = as.character(raca_final) %in% c("Branca", "Parda", "Preta"),
  etapa_completa = !is.na(fora_municipio) & !is.na(idade) & !is.na(cnes) &
                   !is.na(uf_res) & !is.na(mes_inicio),
  causa_saida = case_when(
    !etapa_data ~ "01 Data nao elegivel",
    !etapa_estadio ~ "02 Fora dos estadios I-III",
    !etapa_raca ~ "03 Raca/cor fora das categorias elegiveis ou inconsistente",
    !etapa_completa ~ "04 Covariavel principal ausente",
    TRUE ~ "05 Elegivel para analise municipal"
  ),
  elegivel_municipal = causa_saida == "05 Elegivel para analise municipal"
)

fluxo <- tibble(
  etapa = c("Episodios pre-analise", "Datas elegiveis", "Datas + estadios I-III",
            "Datas + estadios + raca/cor", "Datas + estadios + raca/cor + covariaveis"),
  mulheres = c(nrow(pre), sum(pre$etapa_data),
              sum(pre$etapa_data & pre$etapa_estadio),
              sum(pre$etapa_data & pre$etapa_estadio & pre$etapa_raca),
              sum(pre$elegivel_municipal))
)
fluxo$pct_etapa_inicial <- round(100 * fluxo$mulheres / nrow(pre), 2)
escrever(fluxo, "AUDITORIA_FLUXO_SEQUENCIAL")
escrever(pre %>% count(causa_saida, name = "mulheres"), "AUDITORIA_MOTIVO_PRINCIPAL_EXCLUSAO")
escrever(pre %>% group_by(raca_chr, macro_chr) %>%
  summarise(total = n(), incluidos = sum(elegivel_municipal),
            excluidos = total - incluidos, pct_exclusao = 100 * excluidos / total,
            .groups = "drop"), "AUDITORIA_EXCLUSAO_POR_RACA_E_MACRO")
escrever(pre %>% group_by(raca_chr, municipio_chr, macro_chr, causa_saida) %>%
  summarise(mulheres = n(), .groups = "drop"), "AUDITORIA_EXCLUSAO_RACA_FLUXO_MOTIVO")
if (tail(fluxo$mulheres, 1L) != n_coorte) {
  warning("Selecao sequencial N=", tail(fluxo$mulheres, 1L),
          " difere de coorte analitica N=", n_coorte,
          ". Investigue possivel ausencia em procedimento/desfecho.")
}

escrever(d_todas %>%
  mutate(fora_macro = ifelse(is.na(fora_macro), "Nao identificada", as.character(fora_macro))) %>%
  group_by(raca_final, fora_municipio, fora_macro) %>%
  summarise(mulheres = n(), min_dias = min(tempo_dias, na.rm = TRUE),
            p25 = as.numeric(quantile(tempo_dias, .25, na.rm = TRUE)),
            mediana = median(tempo_dias, na.rm = TRUE),
            p75 = as.numeric(quantile(tempo_dias, .75, na.rm = TRUE)),
            max_dias = max(tempo_dias, na.rm = TRUE),
            n_tempo_zero = sum(tempo_dias == 0, na.rm = TRUE),
            n_atraso_60 = sum(atraso60 == 1L, na.rm = TRUE),
            .groups = "drop"), "AUDITORIA_INTERVALO_MIN_MAX_IQR")

# Sem microdados de pacientes nas saídas de auditoria.
pendentes <- bind_rows(
  d_todas %>% filter(is.na(fora_macro)) %>% count(mun_res, name = "mulheres") %>%
    transmute(campo = "residencia", cod_municipio = as.character(mun_res), mulheres),
  d_todas %>% filter(is.na(fora_macro)) %>% count(mun_estab, name = "mulheres") %>%
    transmute(campo = "estabelecimento", cod_municipio = as.character(mun_estab), mulheres)
)
escrever(pendentes, "AUDITORIA_CODIGOS_REGISTROS_SEM_MACRO")

# 4. MODELOS CR2/HTZ NA MESMA AMOSTRA ----------------------------------------
f_m3 <- function(expo) {
  as.formula(paste0("atraso60 ~ raca_final * ", expo,
    " + ", paste(col_ns, collapse = " + "),
    " + estadio + procedimento + uf_res + mes_inicio"))
}

ajustar <- function(dat, expo) {
  fit <- stats::glm(f_m3(expo), data = dat, family = stats::poisson(link = "log"),
                    control = stats::glm.control(maxit = 50))
  if (nrow(stats::model.frame(fit)) != nrow(dat)) stop("GLM excluiu observacoes; verifique NA.")
  if (!isTRUE(fit$converged) || anyNA(stats::coef(fit)))
    stop("Modelo nao convergiu ou ha coeficiente nao estimavel.")
  fit
}

# CR2 do clubSandwich e HTZ para termos conjuntos; Satterthwaite por coeficiente.
calc_global <- function(fit, dat, expo, prefixo, exportar_coef = TRUE) {
  V <- clubSandwich::vcovCR(fit, cluster = dat$cnes, type = "CR2")
  if (isTRUE(exportar_coef)) {
    escrever(as.data.frame(clubSandwich::coef_test(fit, vcov = V, test = "Satterthwaite")),
             paste0(prefixo, "_COEFICIENTES_CR2"))
  }
  idx <- grep(paste0("^raca_final[^:]*:", expo), names(stats::coef(fit)))
  if (length(idx) != 2L) stop("Esperados 2 parametros de interacao; encontrados: ", length(idx))
  w <- as.data.frame(clubSandwich::Wald_test(
     fit, constraints = clubSandwich::constrain_zero(idx), vcov = V, test = "HTZ"))
  w$amostra <- nrow(dat)
  w$cnes_distintos <- dplyr::n_distinct(dat$cnes)
  escrever(w, paste0(prefixo, "_INTERACAO_GLOBAL_HTZ"))
  list(V = V, global = w)
}

# Funcao generica: usa a matriz CR2 e df finitos HTZ para cada contraste.
contraste_cr2 <- function(fit, V, emm_obj, method, by, family_label) {
  obj <- emmeans::contrast(emm_obj, method = method, by = by, adjust = "none")
  lab <- as.data.frame(summary(obj))
  L <- obj@linfct
  b <- stats::coef(fit)
  if (!all(names(b) %in% colnames(L))) stop("Linhas de contraste sem nomes completos de coeficientes.")
  L <- L[, names(b), drop = FALSE]
  ans <- lapply(seq_len(nrow(L)), function(i) {
    li <- matrix(as.numeric(L[i, ]), nrow = 1)
    colnames(li) <- names(b)
    est <- as.numeric(li %*% b)
    se <- sqrt(as.numeric(li %*% V %*% t(li)))
    w <- as.data.frame(clubSandwich::Wald_test(fit, constraints = li,
                                                vcov = V, test = "HTZ"))
    gl <- as.numeric(w$df_denom[1])
    if (!is.finite(gl) || gl <= 0) stop("gl HTZ invalidos no contraste ", i)
    crit <- stats::qt(.975, df = gl)
    tibble(comparacao = as.character(lab$contrast[i]),
           estrato = as.character(lab[[by]][i]),
           grupo_estratificacao = by, familia_teste = family_label,
           RR = exp(est), IC95_inf = exp(est - crit * se),
           IC95_sup = exp(est + crit * se),
           beta_log = est, EP_CR2 = se, gl_HTZ = gl,
           F_HTZ = as.numeric(w$Fstat[1]), p_HTZ = as.numeric(w$p_val[1]))
  })
  dplyr::bind_rows(ans)
}

# Riscos marginais sao exportados com EP CR2; intervalos assintoticos do emmeans
# sao rotulados separadamente e NAO chamados de IC95% HTZ.
riscos <- function(fit, V, expo, prefixo) {
  emm <- emmeans::emmeans(fit,
    specs = stats::as.formula(paste0("~ raca_final * ", expo)),
    vcov. = V, nuisance = c("uf_res", "mes_inicio", "estadio", "procedimento"),
    wt.nuis = "proportional")
  tab_r <- as.data.frame(summary(emm, type = "response", infer = c(TRUE, FALSE)))
  tab_r$tipo_IC <- "CR2 + aproximacao de grande amostra do emmeans; nao IC HTZ"
  escrever(tab_r, paste0(prefixo, "_RISCOS_AJUSTADOS"))
  cr <- contraste_cr2(fit, V, emm,
    method = list(Parda_vs_Branca = c(-1, 1, 0),
                  Preta_vs_Branca = c(-1, 0, 1),
                  Preta_vs_Parda = c(0, -1, 1)),
    by = expo, family_label = "Raca/cor por territorio: familia de 6 testes exploratorios")
  cr$p_Holm_6 <- stats::p.adjust(cr$p_HTZ, "holm")
  # Ajuste opcional entre os DOIS contrastes preta/branca; relatar
  # o ajuste de seis comparacoes se a familia for exploratoria completa.
  cr$p_Holm_2_PretaBranca <- NA_real_
  j <- which(cr$comparacao == "Preta_vs_Branca")
  if (length(j) == 2L) cr$p_Holm_2_PretaBranca[j] <- stats::p.adjust(cr$p_HTZ[j], "holm")
  escrever(cr, paste0(prefixo, "_CONTRASTES_RACIAIS_CR2_HTZ"))
  # Um contraste territorial dentro de cada uma das tres categorias raciais.
  levels_expo <- levels(emm@grid[[expo]])
  if (length(levels_expo) != 2L) stop("Exposicao ", expo, " nao tem dois niveis.")
  ct <- contraste_cr2(fit, V, emm,
    method = list(Segundo_vs_primeiro_territorio = c(-1, 1)),
    by = "raca_final", family_label = "Territorio por raca: familia de 3 testes exploratorios")
  ct$p_Holm_3 <- stats::p.adjust(ct$p_HTZ, "holm")
  ct$exposicao <- expo
  ct$primeiro_nivel <- levels_expo[1]
  ct$segundo_nivel <- levels_expo[2]
  escrever(ct, paste0(prefixo, "_CONTRASTES_TERRITORIAIS_CR2_HTZ"))
  invisible(list(riscos = tab_r, raciais = cr, territoriais = ct))
}

registrar("Ajustando M3 macrorregional e M3 municipal na mesma amostra de ", nrow(d), " mulheres.")
fit_macro <- ajustar(d, "fora_macro")
res_macro <- calc_global(fit_macro, d, "fora_macro", "M3_MACRO")
res_cr_macro <- riscos(fit_macro, res_macro$V, "fora_macro", "M3_MACRO")
fit_mun <- ajustar(d, "fora_municipio")
res_mun <- calc_global(fit_mun, d, "fora_municipio", "M3_MUNICIPIO_PAREADO")
riscos(fit_mun, res_mun$V, "fora_municipio", "M3_MUNICIPIO_PAREADO")

# Auditoria de equivalencia do spline fixo e formula original com ns(idade,df=3)
fit_convencional <- stats::glm(
  atraso60 ~ raca_final * fora_macro + splines::ns(idade, df = 3) +
    estadio + procedimento + uf_res + mes_inicio,
  data = d, family = stats::poisson(link = "log"))
escrever(tibble(n = nrow(d),
 max_diff_risco_ajustado = max(abs(stats::fitted(fit_macro) - stats::fitted(fit_convencional))),
 risco_medio_modelo_fixo = mean(stats::fitted(fit_macro)),
 risco_medio_modelo_convencional = mean(stats::fitted(fit_convencional))),
 "AUDITORIA_EQUIVALENCIA_SPLINE_FIXO")

# 5. CONCENTRACAO E INFLUENCIA DE CNES ---------------------------------------
clusters <- d %>% group_by(cnes) %>% summarise(
  mulheres = n(), atrasos = sum(atraso60),
  pretas_fora = sum(raca_final == "Preta" & fora_macro == "Outra macrorregiao"),
  pardas_fora = sum(raca_final == "Parda" & fora_macro == "Outra macrorregiao"),
  brancas_fora = sum(raca_final == "Branca" & fora_macro == "Outra macrorregiao"),
  .groups = "drop") %>% arrange(desc(mulheres))
clusters$pct_da_coorte <- 100 * clusters$mulheres / nrow(d)
clusters$pct_das_pretas_fora <- 100 * clusters$pretas_fora /
  sum(d$raca_final == "Preta" & d$fora_macro == "Outra macrorregiao")
escrever(clusters, "CNES_DISTRIBUICAO_E_CONCENTRACAO")

celulas <- d %>% group_by(raca_final, fora_macro) %>% summarise(
  mulheres = n(), cnes_distintos = n_distinct(cnes),
  max_pacientes_por_CNES = max(as.numeric(table(cnes))),
  n_cnes_efetivo = {
    ct <- as.numeric(table(cnes))
    sum(ct)^2 / sum(ct^2)
  }, .groups = "drop")
escrever(celulas, "AUDITORIA_CN_ES_EFETIVOS_POR_ESTRATO")

cn_full <- names(stats::coef(fit_macro))
i_preta <- grep("^raca_finalPreta:fora_macro", cn_full)
i_parda <- grep("^raca_finalParda:fora_macro", cn_full)
if (length(i_preta) != 1L || length(i_parda) != 1L) stop("Nao localizados termos de interacao esperados.")
beta_full_preta <- unname(stats::coef(fit_macro)[i_preta])
beta_full_parda <- unname(stats::coef(fit_macro)[i_parda])

loo <- tibble()
if (isTRUE(RODAR_LOO_CNES)) {
  registrar("Iniciando LOO CNES: ", nrow(clusters), " estabelecimentos (coeficientes pontuais).")
  loo_list <- vector("list", nrow(clusters))
  for (k in seq_len(nrow(clusters))) {
    exc <- clusters$cnes[k]
    dd <- droplevels(d[d$cnes != exc, , drop = FALSE])
    row <- tryCatch({
      mod <- ajustar(dd, "fora_macro")
      b <- stats::coef(mod)
      tibble(cnes_omitido = exc, status = "ok", n_restante = nrow(dd),
        n_omitido = nrow(d) - nrow(dd),
        beta_preta = as.numeric(b[cn_full[i_preta]]),
        beta_parda = as.numeric(b[cn_full[i_parda]]),
        RRR_preta = exp(as.numeric(b[cn_full[i_preta]])),
        RRR_parda = exp(as.numeric(b[cn_full[i_parda]])),
        delta_beta_preta = as.numeric(b[cn_full[i_preta]]) - beta_full_preta,
        delta_beta_parda = as.numeric(b[cn_full[i_parda]]) - beta_full_parda)
    }, error = function(e) tibble(cnes_omitido = exc, status = paste0("falha: ", conditionMessage(e)),
      n_restante = nrow(dd), n_omitido = nrow(d) - nrow(dd),
      beta_preta = NA_real_, beta_parda = NA_real_,
      RRR_preta = NA_real_, RRR_parda = NA_real_,
      delta_beta_preta = NA_real_, delta_beta_parda = NA_real_))
    loo_list[[k]] <- row
    if (k %% 25L == 0L || k == nrow(clusters)) registrar("  CNES ", k, "/", nrow(clusters))
  }
  loo <- bind_rows(loo_list) %>% left_join(clusters, by = c("cnes_omitido" = "cnes"))
  escrever(loo, "LOO_CNES_TODOS_COEFICIENTES")
  # Casos selecionados SEM usar p-valor: maior influencia absoluta, maior tamanho
  # e maior concentracao das pacientes pretas deslocadas para fora da macro.
  extremos <- loo %>% filter(status == "ok", !is.na(delta_beta_preta))
  if (nrow(extremos)) {
    if (isTRUE(RODAR_HTZ_LOO_TODOS)) {
      escolhidos <- extremos$cnes_omitido
    } else {
      esc_delta <- extremos %>% arrange(desc(abs(delta_beta_preta))) %>%
        slice_head(n = QTD_EXTREMOS_LOO) %>% pull(cnes_omitido)
      esc_tamanho <- extremos %>% arrange(desc(n_omitido)) %>%
        slice_head(n = QTD_EXTREMOS_LOO) %>% pull(cnes_omitido)
      esc_pretas <- extremos %>% arrange(desc(pretas_fora)) %>%
        slice_head(n = QTD_EXTREMOS_LOO) %>% pull(cnes_omitido)
      escolhidos <- unique(c(esc_delta, esc_tamanho, esc_pretas))
    }
    registrar("HTZ reestimado apos omissao de ", length(escolhidos), " CNES selecionados.")
    tab_htz <- lapply(escolhidos, function(exc) {
      dd <- droplevels(d[d$cnes != exc, , drop = FALSE])
      tryCatch({
        mm <- ajustar(dd, "fora_macro")
        VV <- clubSandwich::vcovCR(mm, cluster = dd$cnes, type = "CR2")
        ix <- grep("^raca_final[^:]*:fora_macro", names(stats::coef(mm)))
        ww <- as.data.frame(clubSandwich::Wald_test(mm,
          constraints = clubSandwich::constrain_zero(ix), vcov = VV, test = "HTZ"))
        tibble(cnes_omitido = exc, status = "ok", p_interacao_HTZ = as.numeric(ww$p_val[1]),
          F_HTZ = as.numeric(ww$Fstat[1]), gl_denominador_HTZ = as.numeric(ww$df_denom[1]))
      }, error = function(e) tibble(cnes_omitido = exc, status = paste0("falha: ", conditionMessage(e)),
         p_interacao_HTZ = NA_real_, F_HTZ = NA_real_, gl_denominador_HTZ = NA_real_))
    }) %>% bind_rows()
    escrever(tab_htz, "LOO_CNES_HTZ_SELECIONADOS")
    stats <- extremos %>% summarise(
      modelos_reestimados = n(),
      min_RRR_preta = min(RRR_preta), max_RRR_preta = max(RRR_preta),
      min_RRR_parda = min(RRR_parda), max_RRR_parda = max(RRR_parda),
      min_delta_preta = min(delta_beta_preta), max_delta_preta = max(delta_beta_preta))
    escrever(stats, "LOO_CNES_RESUMO_ESTABILIDADE")
  }
}

# 6. VALIDACAO TEMPORAL (NAO INVENTAR UMA MALHA DE 2023) -----------------------
ler_mapa <- function(arq) {
  cabecalho <- readLines(arq, n = 1L, warn = FALSE, encoding = "UTF-8")
  separador <- if (length(cabecalho) && lengths(regmatches(cabecalho, gregexpr(";", cabecalho, fixed = TRUE))) >
                   lengths(regmatches(cabecalho, gregexpr(",", cabecalho, fixed = TRUE)))) ";" else ","
  raw <- readr::read_delim(arq, delim = separador,
    col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE,
    locale = readr::locale(encoding = "UTF-8"), progress = FALSE)
  checagem_colunas(raw, c("cod_municipio", "cod_macrorregiao_de_saude"), basename(arq))
  raw <- raw %>% transmute(
     mun6 = normalizar_mun6(cod_municipio),
     macro = trimws(as.character(cod_macrorregiao_de_saude)))
  raw$macro[!is.na(raw$macro) & raw$macro == ""] <- NA_character_
  if (anyNA(raw$mun6) || anyDuplicated(raw$mun6))
    stop("CSV territorial tem codigos municipais nulos ou duplicados: ", arq)
  raw
}

status_map <- character()
if (file.exists(ARQ_MAPA_2026)) {
  map26 <- ler_mapa(ARQ_MAPA_2026)
  escrever(tibble(ano_referencia = "2026", linhas = nrow(map26),
                 macros_distintas = n_distinct(map26$macro),
                 arquivo = basename(ARQ_MAPA_2026),
                 md5 = as.character(tools::md5sum(ARQ_MAPA_2026))),
           "AUDITORIA_MAPA_2026")
  if (!file.exists(ARQ_MAPA_2023)) {
    modelo <- tibble(cod_municipio = map26$mun6,
       cod_macrorregiao_de_saude = rep(NA_character_, nrow(map26)),
       fonte_documental_2023 = rep(NA_character_, nrow(map26)),
       data_vigencia = rep(NA_character_, nrow(map26)))
    escrever(modelo, "MODELO_A_PREENCHER_MAPA_MUNICIPAL_2023_SEM_IMPUTACAO")
  }
} else {
  warning("CSV 2026 nao encontrado. Os modelos utilizam o mapa JA VINCULADO no FST.")
  status_map <- c(status_map, "CSV 2026 nao localizado; arquivos FST ja enriquecidos foram utilizados.")
}

if (!file.exists(ARQ_MAPA_2023)) {
  status_map <- c(status_map,
    "SEM MAPA 2023: arquivo municipio_macrorregiao_REFERENCIA_2023.csv nao encontrado.",
    "As analises com macro do FST permanecem EXPLORATORIAS, codificacao disponivel em outubro/2026.",
    "Documentar correspondencia de 2023 com atos estaduais/CIB e referencia territorial do MS.",
    "NAO afirmar que os territorios de 2026 eram os mesmos de 2023.")
} else {
  mapa23 <- ler_mapa(ARQ_MAPA_2023)
  if (nrow(mapa23) < 5400L || sum(!is.na(mapa23$macro)) < 5400L)
    stop("Mapa 2023 incompleto: nao e seguro gerar os modelos historicos.")
  status_map <- c(status_map,
      paste("Arquivo candidato 2023 encontrado:", ARQ_MAPA_2023),
      paste("Fonte informada:", if (nzchar(FONTE_DOCUMENTAL_2023)) FONTE_DOCUMENTAL_2023 else "NAO INFORMADA"),
      paste("Validacao documental:", if (MAPA_2023_VALIDADO && nzchar(FONTE_DOCUMENTAL_2023))
            "DECLARADA (conferir atos)" else "PENDENTE - APENAS EXPLORATORIO"))
  escrever(tibble(referencia = "2023_candidata", n_municipios = nrow(mapa23),
                 n_macros = n_distinct(mapa23$macro),
                 fonte_documental = FONTE_DOCUMENTAL_2023,
                 validado = MAPA_2023_VALIDADO && nzchar(FONTE_DOCUMENTAL_2023)),
           "AUDITORIA_MAPA_2023_CANDIDATO")
  if (exists("map26")) {
    mm <- map26 %>% rename(macro_2026 = macro) %>%
      full_join(mapa23 %>% rename(macro_2023 = macro), by = "mun6") %>%
      mutate(codigo_igual = !is.na(macro_2026) & !is.na(macro_2023) &
               macro_2026 == macro_2023)
    escrever(mm, "COMPARACAO_CODIGOS_MACRO_2023_2026_POR_MUNICIPIO")
    # Mudanca de CODIGO nao equivale necessariamente a mudanca territorial.
  }
  mres23 <- mapa23$macro[match(normalizar_mun6(d_todas$mun_res), mapa23$mun6)]
  mtrat23 <- mapa23$macro[match(normalizar_mun6(d_todas$mun_estab), mapa23$mun6)]
  novo_fluxo <- ifelse(is.na(mres23) | is.na(mtrat23), NA_character_,
                   ifelse(mres23 == mtrat23, "Mesma macrorregiao", "Outra macrorregiao"))
  pares <- tibble(fora_macro_2026 = ifelse(is.na(d_todas$fora_macro), "NA", as.character(d_todas$fora_macro)),
                 fora_macro_2023 = ifelse(is.na(novo_fluxo), "NA", novo_fluxo))
  escrever(pares %>% count(fora_macro_2026, fora_macro_2023, name = "mulheres"),
           "MATRIZ_RECLASSIFICACAO_COORTE_2023_VS_2026")
  # Comparacao justa: mesmAS pacientes com ambos os mapas completos.
  histor <- d_todas
  histor$fora_macro_hist <- novo_fluxo
  histor <- histor %>% filter(!is.na(fora_macro), !is.na(fora_macro_hist))
  # Reaproveitar as mesmas bases spline: nao alterar demais covariaveis.
  hist_idx <- which(!is.na(d_todas$fora_macro) & !is.na(novo_fluxo))
  # Gerar splines fixos na subamostra comum para AMBAS especificacoes.
  compar <- d_todas[hist_idx, , drop = FALSE] %>% mutate(
    raca_final = factor(as.character(raca_final), levels = c("Branca", "Parda", "Preta")),
    fora_municipio = factor(as.character(fora_municipio)),
    estadio = factor(as.character(estadio)),
    procedimento = factor(as.character(procedimento)),
    uf_res = factor(as.character(uf_res)), mes_inicio = factor(as.character(mes_inicio)),
    cnes = as.character(cnes), atraso60 = as.integer(atraso60))
  # Para comparabilidade, splines recomputadas UMA vez e reutilizadas.
  mat_hist <- splines::ns(compar$idade, df = 3)
  for (j in seq_along(col_ns)) {
    compar[[col_ns[j]]] <- as.numeric(mat_hist[, j])
  }
  for (qual in c("ATUAL_2026", "CANDIDATO_2023")) {
    compar$fora_macro <- factor(if (qual == "ATUAL_2026")
      as.character(compar$fora_macro) else novo_fluxo[hist_idx],
      levels = c("Mesma macrorregiao", "Outra macrorregiao"))
    if (any(table(compar$raca_final, compar$fora_macro) < 10L)) {
       warning("Celulas raciais pequenas para modelo historico: ", qual)
       next
    }
    ff <- ajustar(compar, "fora_macro")
    rr <- calc_global(ff, compar, "fora_macro", paste0("HIST_COMPARAVEL_", qual))
    riscos(ff, rr$V, "fora_macro", paste0("HIST_COMPARAVEL_", qual))
  }
}
writeLines(status_map, file.path(SAIDA, "STATUS_E_PROVENIENCIA_MACRORREGIOES_2023.txt"), useBytes = TRUE)

# 7. RESUMO E PROVENIENCIA -----------------------------------------------------
p_global_macro <- as.numeric(res_macro$global$p_val[1])
p_global_mun <- as.numeric(res_mun$global$p_val[1])
ca <- res_cr_macro$raciais
peb <- ca %>% filter(comparacao == "Preta_vs_Branca")
resumo <- c(
  paste("Coorte municipal original:", n_coorte),
  paste("Coorte pareada com macro identificada:", nrow(d)),
  paste("Estabelecimentos CNES:", dplyr::n_distinct(d$cnes)),
  sprintf("Interacao global raca x macro CR2/HTZ: p = %.6f", p_global_macro),
  sprintf("Interacao global raca x municipio (mesma amostra) CR2/HTZ: p = %.6f", p_global_mun),
  "Contrastes Preta/Branca com CR2 e IC HTZ (apenas resultados condicionais):"
)
for (i in seq_len(nrow(peb))) {
  rr <- peb[i, ]
  resumo <- c(resumo, sprintf("  %s: RR %.4f; IC95%% %.4f - %.4f; p_HTZ %.6f; p_Holm_2 %.6f",
    rr$estrato, rr$RR, rr$IC95_inf, rr$IC95_sup, rr$p_HTZ, rr$p_Holm_2_PretaBranca))
}
resumo <- c(resumo,
  "Importante: p distintos por estrato NAO provam interacao. Dar precedencia ao teste global.",
  "RR, IC e contrastes estimados entre mulheres que INICIARAM quimioterapia, nao entre todas as elegiveis.",
  "Mapa 2026 nao demonstra organizacao assistencial vigente em 2023.")
writeLines(resumo, file.path(SAIDA, "RESUMO_PARA_MANUSCRITO_E_VERIFICACAO.txt"))
writeLines(capture.output(utils::sessionInfo()), file.path(SAIDA, "SESSION_INFO_R.txt"))
md5 <- tools::md5sum(c(ARQ_ANALITICA, ARQ_PRE))
escrever(tibble(arquivo = basename(names(md5)), md5 = unname(as.character(md5))),
         "PROVENIENCIA_FST_MD5")
registrar("\nVALIDACAO CONCLUIDA. Confira: ", SAIDA)
registrar(paste(resumo, collapse = "\n"))
