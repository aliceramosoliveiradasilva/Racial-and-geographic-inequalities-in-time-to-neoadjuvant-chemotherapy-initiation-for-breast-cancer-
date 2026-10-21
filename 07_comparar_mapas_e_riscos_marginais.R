# PUBLICATION REPRODUCIBILITY — BREAST CANCER NEOADJUVANT CHEMOTHERAPY
# 2024 vs 2026 crosswalks and Poisson margins.
# Executes locally on authorized DATASUS data; never commit patient-level data.
# 2026-10-15. Historical crosswalks 2024/2026 are not validated for 2023.
# All R scripts should be run from project root; MAMA_RAIZ configures data root.
# NOTE: full end-to-end re-execution was not possible in the delivery environment.
# -------------------------------------------------------------------------

# ============================================================================
# 04 - RESOLVER PENDENCIAS DO ARTIGO: MACRORREGIOES 2023 E RISCOS MARGINAIS
# Coorte: quimioterapia neoadjuvante, SUS, 2023
# Versao de trabalho: 14/10/2026
#
# ETAPA A: obtem uma REFERENCIA DOCUMENTADA DE MARCO/2024 (NAO DE 2023)
#          na Portaria GM/MS 3.233/2024 e compara com o mapa usado de 2026.
# ETAPA B: recalcula TODOS os riscos absolutos por padronizacao marginal
#          sobre as MESMAS 15.414 mulheres, mantendo spline e covariaveis.
# ETAPA C: testa interacoes e contrastes CR2/HTZ para ambos os mapas em
#          AMOSTRA PAREADA; gera tabela de reclassificacao por paciente.
# ETAPA D: permite incluir um CSV de 2023 quando VALIDADO DOCUMENTALMENTE.
#
# OBS: a Portaria de 2024 NAO comprova, sozinha, o mapa vigente em 2023.
# O programa nao refaz a selecao individual e nao altera arquivos originais.
# ============================================================================

# 0. CONFIGURACAO -------------------------------------------------------------
RAIZ <- Sys.getenv("MAMA_RAIZ", unset = getwd())
DIR_MACROS <- file.path(RAIZ, 'dados_mama_SIA_AQ_2019_2025',
                       'resultados_com_macrorregioes', 'ANALISE_MACRORREGIOES')
ARQ_COORTE <- file.path(DIR_MACROS, 'coorte_analitica_2023_COM_MACROS.fst')
MAPA_2026 <- file.path(RAIZ, 'macroregiao_de_saude(1).csv')

# O programa obtem esta tabela diretamente da Portaria federal, se necessario.
MAPA_2024 <- file.path(RAIZ, 'municipio_macrorregiao_PORTARIA_MAR2024.csv')
HTML_2024_LOCAL <- ''    # opcional: caminho de HTML salvo manualmente no navegador
URL_2024 <- paste0('https://bvsms.saude.gov.br/bvs/saudelegis/gm/2024/',
                   'prt3233_06_03_2024.html')

# Se houver um mapa real de 2023, coloque o arquivo e INFORME a fonte precisa.
# Cabecalhos aceitos: cod_municipio;cod_macrorregiao_de_saude
MAPA_2023 <- file.path(RAIZ, 'municipio_macrorregiao_REFERENCIA_2023.csv')
FONTE_OFICIAL_2023 <- ''  # Ex.: URL + identificador do ato, data de vigencia

SAIDA <- file.path(DIR_MACROS, 'RESOLUCAO_PENDENCIAS_1_2_20261014')
dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)

# 1. PACOTES ------------------------------------------------------------------
necessarios <- c('fst', 'data.table', 'clubSandwich', 'splines', 'rvest', 'xml2')
faltantes <- necessarios[!vapply(necessarios, requireNamespace, logical(1), quietly=TRUE)]
if (length(faltantes)) stop('Instale os pacotes: install.packages(c(',
                           paste(sprintf('"%s"', faltantes), collapse=', '), '))')
suppressPackageStartupMessages(library(data.table))

escrever <- function(d, nome) data.table::fwrite(as.data.table(d),
  file.path(SAIDA, paste0(nome, '.csv')), sep=';', na='')
registrar <- function(...) cat(..., '\n', sep='')
normalizar_mun6 <- function(x) {
  s <- trimws(as.character(x))
  s <- sub('\\.0+$', '', s)
  s[s %in% c('', 'NA', 'NaN', 'NULL', '000000', '0000000')] <- NA_character_
  s[!is.na(s) & !grepl('^[0-9]{6,7}$', s)] <- NA_character_
  substr(s, 1, 6)
}
checar <- function(d, cols, nome) {
  aus <- setdiff(cols, names(d))
  if (length(aus)) stop(nome, ': faltam colunas: ', paste(aus, collapse=', '))
}

# 2. MAPA DE 2026 -------------------------------------------------------------
ler_mapa <- function(arq, origem) {
  d <- data.table::fread(arq, sep='auto', colClasses='character',
                         encoding='UTF-8', showProgress=FALSE)
  possiveis <- list(
    codigo_mun = c('cod_municipio', 'municipio6', 'codigo_ibge', 'mun6'),
    codigo_macro = c('cod_macrorregiao_de_saude', 'macro', 'codigo_macro', 'macro_cod'))
  pegar <- function(alternativas) {
    achou <- intersect(alternativas, names(d))
    if (!length(achou)) stop('Falta campo entre: ', paste(alternativas, collapse='/'),
                             ' no arquivo ', arq)
    achou[1]
  }
  cm <- pegar(possiveis$codigo_mun)
  cg <- pegar(possiveis$codigo_macro)
  out <- unique(data.table(mun6=normalizar_mun6(d[[cm]]),
                macro=trimws(as.character(d[[cg]]))))
  out[macro == '', macro := NA_character_]
  if (anyNA(out$mun6) || anyNA(out$macro))
    stop('O mapa ', origem, ' tem codigos nulos. Corrija antes de executar.')
  diverg <- out[, .(macros=data.table::uniqueN(macro)), by=mun6][macros > 1]
  if (nrow(diverg)) stop('Municipios com mais de uma macro em ', origem,
                         ': ', paste(head(diverg$mun6), collapse=', '))
  setkey(out, mun6)
  out <- unique(out, by='mun6')
  out[, fonte := origem]
  out
}
if (!file.exists(ARQ_COORTE)) stop('FST nao encontrado: ', ARQ_COORTE)
if (!file.exists(MAPA_2026)) stop('Tabela oficial de 2026 nao encontrada: ', MAPA_2026)
map26 <- ler_mapa(MAPA_2026, 'Portal Dados Abertos SUS; extraido outubro/2026')
registrar('Mapa 2026: ', nrow(map26), ' municipios, ', uniqueN(map26$macro), ' macros.')

# 3. CRIAR MAPA FEDERAL DE MARCO/2024 -----------------------------------------
# A portaria traz linhas IBGE/municipio/macrorregiao, incluindo subtotais.
# Extraimos somente linhas contendo codigo IBGE com 6 e codigo macro com 4 digitos.
extrair_portaria_2024 <- function() {
  registrar('Buscando Portaria GM/MS 3.233/2024 (NAO mapa 2023)...')
  doc <- if (nzchar(HTML_2024_LOCAL) && file.exists(HTML_2024_LOCAL)) {
    xml2::read_html(HTML_2024_LOCAL)
  } else {
    xml2::read_html(URL_2024)
  }
  trs <- rvest::html_elements(doc, 'tr')
  if (!length(trs)) stop('Nao achei linhas <tr> na portaria. Salve HTML manualmente.')
  itens <- lapply(trs, function(tr) {
    cel <- trimws(rvest::html_text2(rvest::html_elements(tr, 'th,td')))
    ind_mun <- which(grepl('^[0-9]{6}$', cel))
    if (!length(ind_mun)) return(NULL)
    ind_mun <- ind_mun[1]
    mun <- cel[ind_mun]
    ind_macro <- which(grepl('^[0-9]{4}$', cel) &
                       substr(cel, 1, 2) == substr(mun, 1, 2))
    if (!length(ind_macro)) return(NULL)
    macro <- cel[ind_macro[1]]
    nome <- if (ind_mun < length(cel)) cel[ind_mun+1] else ''
    data.table(mun6=mun, macro=macro, nome_municipio=nome)
  })
  saida <- unique(rbindlist(itens, fill=TRUE))
  if (!nrow(saida)) stop('Nao foi possivel extrair codigo do municipio e macro.')
  conflitos <- saida[, .(nmacros=data.table::uniqueN(macro)), by=mun6][nmacros > 1]
  if (nrow(conflitos)) stop('Ha municipios com macro conflitante na portaria: ',
                          paste(head(conflitos$mun6), collapse=', '))
  setorder(saida, mun6)
  saida <- unique(saida, by='mun6')
  if (nrow(saida) < 5000L) stop('A portaria forneceu so ', nrow(saida),
    ' municipios; conferir estrutura/versao da pagina antes de usar.')
  if (nrow(saida) != 5570L) warning('Portaria de 2024 contem ', nrow(saida),
    ' municipios (e nao 5570). Analise de sensibilidade usara apenas pares completos.')
  escrever(saida, 'MAPA_FEDERAL_MAR2024_EXTRAIDO_BRUTO')
  data.table::fwrite(saida[, .(cod_municipio=mun6,
                                cod_macrorregiao_de_saude=macro)],
                      MAPA_2024, sep=';')
  registrar('Portaria extraida: ', nrow(saida), ' municipios, ',
            uniqueN(saida$macro), ' macrorregioes em marco/2024.')
  invisible(saida)
}
if (!file.exists(MAPA_2024)) {
  tryCatch(extrair_portaria_2024(), error=function(e) {
    registrar('Aviso: nao consegui baixar/ler Portaria 2024: ', conditionMessage(e))
    registrar('Se necessario, salve a pagina HTML no navegador e use HTML_2024_LOCAL.')
  })
}
map24 <- if (file.exists(MAPA_2024)) {
  tryCatch(ler_mapa(MAPA_2024, 'Portaria GM/MS 3.233/2024 (marco/2024)'),
           error=function(e) {registrar('Mapa 2024 invalido: ', conditionMessage(e)); NULL})
} else NULL

# Mapa 2023 so sera rotulado validado se houver REFERENCIA DOCUMENTAL.
map23 <- NULL
if (file.exists(MAPA_2023) && nzchar(trimws(FONTE_OFICIAL_2023))) {
  map23 <- ler_mapa(MAPA_2023, paste('Referencia documentada de 2023:', FONTE_OFICIAL_2023))
  registrar('Mapa 2023 documentado recebido: ', nrow(map23), ' municipios, ',
            uniqueN(map23$macro), ' macros. Conferir consistencia do ato e vigencia.')
}

# 4. VINCULAR MAPAS A MESMA COORTE -------------------------------------------
d0 <- as.data.table(fst::read_fst(ARQ_COORTE))
necessario_coorte <- c('mun_res', 'mun_estab', 'raca_final', 'idade', 'estadio',
 'procedimento', 'uf_res', 'mes_inicio', 'cnes', 'atraso60', 'fora_municipio')
checar(d0, necessario_coorte, 'coorte analitica')
if (nrow(d0) != 15415L) warning('Coorte com ', nrow(d0), ' linhas (esperado 15415).')
d0[, `:=`(mun_res6=normalizar_mun6(mun_res),
          mun_estab6=normalizar_mun6(mun_estab))]

anexar_mapa <- function(d, m, suf) {
  ir <- match(d$mun_res6, m$mun6)
  it <- match(d$mun_estab6, m$mun6)
  cr <- paste0('macro_res_', suf); ct <- paste0('macro_trat_', suf)
  fl <- paste0('fora_macro_', suf)
  d[, (cr) := m$macro[ir]]
  d[, (ct) := m$macro[it]]
  ok <- !is.na(d[[cr]]) & !is.na(d[[ct]])
  d[, (fl) := fifelse(ok, fifelse(get(cr)==get(ct),
        'Mesma macrorregiao', 'Outra macrorregiao'), NA_character_)]
  invisible(d)
}
anexar_mapa(d0, map26, '2026')
if (!is.null(map24)) anexar_mapa(d0, map24, '2024')
if (!is.null(map23)) anexar_mapa(d0, map23, '2023')

# A auditoria so publica informacoes AGREGADAS e codigos municipais, nao CNS.
aud_mun <- data.table(mun6=map26$mun6, macro_2026=map26$macro)
if (!is.null(map24)) aud_mun[, macro_2024 := map24$macro[match(mun6, map24$mun6)]]
if (!is.null(map23)) aud_mun[, macro_2023 := map23$macro[match(mun6, map23$mun6)]]
escrever(aud_mun, 'AUDITORIA_CORRESPONDENCIA_MUNICIPAL_2024_2026')

# A mudanca de CODIGO nao necessariamente equivale a alteracao de FRONTEIRA.
# Por isso, a analise central compara o status de deslocamento da MESMA mulher.
if (!is.null(map24)) {
  tab24 <- d0[!is.na(fora_macro_2024) & !is.na(fora_macro_2026),
                .N, by=.(raca_final, fora_macro_2024, fora_macro_2026)]
  tab24[, mudou := fora_macro_2024 != fora_macro_2026]
  escrever(tab24, 'RECLASSIFICACAO_FLUXO_MAR2024_VS_2026_POR_RACA')
  escrever(d0[, .(n_coorte=.N,  n_com_mapa_2024=sum(!is.na(fora_macro_2024)),
                 n_com_mapa_2026=sum(!is.na(fora_macro_2026)),
                 n_em_ambos=sum(!is.na(fora_macro_2024)&!is.na(fora_macro_2026)),
                 n_fluxo_alterado=sum(fora_macro_2024 != fora_macro_2026, na.rm=TRUE))],
           'AUDITORIA_COBERTURA_E_RECLASSIFICACAO_2024_2026')
}
if (!is.null(map23)) {
  escrever(d0[!is.na(fora_macro_2023)&!is.na(fora_macro_2026),
     .N, by=.(raca_final, fora_macro_2023, fora_macro_2026)],
     'RECLASSIFICACAO_FLUXO_2023_VS_2026_POR_RACA')
}

# 5. AJUSTAR MODELO M3 E PADRONIZAR SOBRE TODA A AMOSTRA -----------------------
# Cada paciente contribui com sua idade, estadio, procedimento, UF e mes reais.
# Em cada celula, variamos somente categoria racial e exposicao territorial.
# E uma padronizacao DESCRITIVA; jamais interpretar mudanca de raca como
# intervencao causal e jamais usar linguagem de efeito causal identificado.
preparar <- function(d, coluna) {
  d <- as.data.frame(copy(d))
  d$raca_final <- factor(as.character(d$raca_final),
                         levels=c('Branca','Parda','Preta'))
  d$exposicao <- factor(as.character(d[[coluna]]),
                         levels=c('Mesma macrorregiao','Outra macrorregiao'))
  d$estadio <- factor(as.character(d$estadio))
  d$procedimento <- factor(as.character(d$procedimento))
  d$uf_res <- factor(as.character(d$uf_res))
  d$mes_inicio <- factor(as.character(d$mes_inicio))
  d$cnes <- as.character(d$cnes)
  d$idade <- as.numeric(d$idade)
  d$atraso60 <- as.integer(d$atraso60)
  cols <- c('raca_final','exposicao','estadio','procedimento','uf_res',
            'mes_inicio','cnes','idade','atraso60')
  d <- d[stats::complete.cases(d[, cols]), , drop=FALSE]
  if (length(unique(d$cnes)) < 20) stop('Poucos clusters CNES; conferir coorte.')
  if (any(!d$atraso60 %in% c(0L, 1L))) stop('atraso60 nao binaria.')
  if (any(table(d$raca_final, d$exposicao) < 10))
    stop('Pelo menos um estrato raca x territorio tem menos de 10 observacoes.')
  # Reconstroi o spline UMA VEZ, sobre a amostra especificada.
  # Ao construir previsoes contrafactuais, estas colunas permanecem inalteradas.
  s <- splines::ns(d$idade, df=3)
  for (j in seq_len(ncol(s))) d[[paste0('idade_ns', j)]] <- as.numeric(s[,j])
  attr(d, 'spline_knots') <- attr(s, 'knots')
  attr(d, 'spline_boundary') <- attr(s, 'Boundary.knots')
  d
}

ajustar_m3 <- function(d) {
  f <- atraso60 ~ raca_final * exposicao + idade_ns1 + idade_ns2 + idade_ns3 +
    estadio + procedimento + uf_res + mes_inicio
  fit <- stats::glm(f, family=stats::poisson(link='log'), data=d,
                    control=stats::glm.control(maxit=50))
  if (!isTRUE(fit$converged) || anyNA(stats::coef(fit)))
    stop('M3 nao convergiu ou contem coeficientes nao estimaveis.')
  if (nrow(stats::model.frame(fit)) != nrow(d)) stop('M3 excluiu observacoes.')
  V <- clubSandwich::vcovCR(fit, cluster=d$cnes, type='CR2')
  ind <- grep('^raca_final[^:]*:exposicao', names(coef(fit)))
  if (length(ind) != 2L) stop('Interacao raca x territorio nao tem 2 termos.')
  HTZ <- as.data.frame(clubSandwich::Wald_test(fit,
    constraints=clubSandwich::constrain_zero(ind), vcov=V, test='HTZ'))
  list(fit=fit, V=V, global=HTZ)
}

# Graus de liberdade para delta de uma funcao nao linear via gradiente
# de primeira ordem + matriz CR2. Aproximacao t de Satterthwaite.
gl_delta_CR2 <- function(m, gradiente) {
  L <- matrix(unname(gradiente), nrow=1,
              dimnames=list(NULL, names(coef(m$fit))))
  tryCatch({
    w <- as.data.frame(clubSandwich::Wald_test(m$fit,
         constraints=L, vcov=m$V, test='HTZ'))
    df <- as.numeric(w$df_denom[1])
    if (!is.finite(df) || df<=0) NA_real_ else df
  }, error=function(e) NA_real_)
}

padronizar <- function(m, dados, nome_analise) {
  beta <- stats::coef(m$fit)
  V <- as.matrix(m$V)
  saidas <- list(); gr <- list(); k <- 0L
  for (territorio in levels(dados$exposicao)) {
    for (raca in levels(dados$raca_final)) {
      novo <- dados
      novo$raca_final <- factor(rep(raca, nrow(dados)), levels=levels(dados$raca_final))
      novo$exposicao <- factor(rep(territorio, nrow(dados)),
                               levels=levels(dados$exposicao))
      X <- stats::model.matrix(stats::delete.response(stats::terms(m$fit)),
                                data=novo, contrasts.arg=m$fit$contrasts)
      if (!setequal(colnames(X), names(beta))) stop('Matriz de previsao incompleta.')
      X <- X[, names(beta), drop=FALSE]
      pred <- drop(exp(X %*% beta))
      if (any(!is.finite(pred))) stop('Previsoes nao finitas em ', nome_analise)
      risco <- mean(pred)
      grad <- drop(crossprod(X, pred)/length(pred))
      names(grad) <- names(beta)
      ep <- sqrt(max(0, drop(t(grad) %*% V %*% grad)))
      # IC95% delta/CR2 na escala logito, com gl de Satterthwaite
      # aproximados a partir do gradiente de primeira ordem.
      # Nao equivale a um teste HTZ exato para funcao nao linear.
      if (risco <= 0 || risco >= 1) stop('Risco padronizado fora (0,1).')
      glogit <- grad / (risco * (1-risco))
      ep_logit <- sqrt(max(0, drop(t(glogit) %*% V %*% glogit)))
      df_aprox <- gl_delta_CR2(m, glogit)
      crit <- if (is.finite(df_aprox)) stats::qt(.975, df=df_aprox)
              else stats::qnorm(.975)
      intervalo <- stats::plogis(stats::qlogis(risco) + c(-1,1)*crit*ep_logit)
      k <- k+1L
      chave <- paste(raca, territorio, sep='::')
      saidas[[k]] <- data.table(analise=nome_analise, n_padronizacao=nrow(dados),
        raca=raca, territorio=territorio, risco=risco,
        risco_percentual=100*risco, IC95_inf=intervalo[1], IC95_sup=intervalo[2],
        EP_CR2_delta=ep, df_Satterthwaite_aprox=df_aprox,
        n_previsoes_maiores_1=sum(pred>1),
        metodo_IC='delta logito + CR2 + t-Satterthwaite aproximado; nao HTZ exato')
      gr[[chave]] <- list(risco=risco, grad=grad)
    }
  }
  tab <- rbindlist(saidas)
  escrever(tab, paste0(nome_analise, '_RISCOS_MARGINAIS_PADRONIZADOS'))
  # Diferencas absolutas (pontos percentuais) usando o mesmo alvo de padronizacao
  comp <- list(); k <- 0L
  for (territorio in levels(dados$exposicao)) {
    for (rac_comp in c('Parda','Preta')) {
      a <- gr[[paste(rac_comp, territorio, sep='::')]]
      b <- gr[[paste('Branca', territorio, sep='::')]]
      gd <- a$grad - b$grad
      rd <- a$risco - b$risco
      ep <- sqrt(max(0, drop(t(gd) %*% V %*% gd)))
      df_rd <- gl_delta_CR2(m, gd)
      crit_rd <- if (is.finite(df_rd)) stats::qt(.975, df=df_rd)
                 else stats::qnorm(.975)
      k <- k+1L
      comp[[k]] <- data.table(analise=nome_analise,
         contraste=paste(rac_comp, 'menos Branca'), territorio=territorio,
         diferenca_pp=100*rd, IC95_inf_pp=100*(rd-crit_rd*ep),
         IC95_sup_pp=100*(rd+crit_rd*ep),
         EP_CR2_delta_pp=100*ep, df_Satterthwaite_aprox=df_rd)
    }
  }
  escrever(rbindlist(comp), paste0(nome_analise, '_DIFERENCAS_RISCO_PP'))
  tab
}

contrastes_rr_cr2 <- function(m, dados, nome_analise) {
  # Tres comparacoes raciais em cada estrato (familia exploratoria de 6).
  # RR e IC CR2/HTZ com gl finitos; nao inferir interacao de p estratificados.
  b <- coef(m$fit); V <- as.matrix(m$V)
  nms <- names(b); out <- list(); k <- 0L
  levels_expo <- levels(dados$exposicao)
  montar_L <- function(raca, territorio) {
    L <- setNames(rep(0, length(b)), nms)
    if (raca != 'Branca') {
      co <- paste0('raca_final',raca)
      if (!(co %in% nms)) stop('Coeficiente ausente: ', co)
      L[co] <- 1
      if (territorio == levels_expo[2]) {
        nome_int <- paste0(co, ':exposicao', levels_expo[2])
        if (!(nome_int %in% nms)) stop('Interacao ausente: ', nome_int)
        L[nome_int] <- 1
      }
    }
    L
  }
  pares <- list(c('Parda','Branca'), c('Preta','Branca'), c('Preta','Parda'))
  for (territorio in levels_expo) {
    for (par in pares) {
      L <- montar_L(par[1], territorio) - montar_L(par[2], territorio)
      est <- sum(L*b)
      EP <- sqrt(max(0, drop(t(L) %*% V %*% L)))
      mat <- matrix(unname(L), nrow=1L, dimnames=list(NULL, names(L)))
      w <- as.data.frame(clubSandwich::Wald_test(m$fit,
         constraints=mat, vcov=m$V, test='HTZ'))
      gl <- as.numeric(w$df_denom[1]); p <- as.numeric(w$p_val[1])
      se_mult <- stats::qt(.975, df=gl)*EP
      k <- k+1L
      out[[k]] <- data.table(analise=nome_analise, territorio=territorio,
       contraste=paste0(par[1],'_vs_',par[2]), RR=exp(est),
       IC95_inf=exp(est-se_mult), IC95_sup=exp(est+se_mult),
       EP_logRR_CR2=EP, gl_HTZ=gl, p_HTZ=p)
    }
  }
  z <- rbindlist(out)
  z[, p_Holm_6 := p.adjust(p_HTZ, method='holm')]
  escrever(z, paste0(nome_analise,'_CONTRASTES_RR_CR2_HTZ'))
  z
}

rodar <- function(base, exposicao, nome) {
  base <- preparar(base, exposicao)
  fit <- ajustar_m3(base)
  marg <- padronizar(fit, base, nome)
  rr <- contrastes_rr_cr2(fit, base, nome)
  glob <- as.data.table(fit$global)
  glob[, `:=`(analise=nome, n=nrow(base), cnes_distintos=uniqueN(base$cnes))]
  escrever(glob, paste0(nome,'_INTERACAO_CR2_HTZ'))
  auditoria <- data.table(analise=nome, n=nrow(base),
    cnes=uniqueN(base$cnes), atraso_observado=mean(base$atraso60),
    media_previsao_original=mean(fitted(fit$fit)),
    max_diferenca_modelo_original= max(abs(fitted(fit$fit)-
                                      predict(fit$fit, type='response'))),
    risco_max_individual_modelo=max(fitted(fit$fit)))
  escrever(auditoria, paste0(nome,'_AUDITORIA_PREVISOES'))
  list(global=glob, riscos=marg, RR=rr, auditoria=auditoria)
}

# Amostra da analise principal anterior: dados completos com macro de 2026.
d26 <- copy(d0[!is.na(fora_macro_2026)])
if (nrow(d26) != 15414L) warning('N 2026=', nrow(d26), ' (esperado 15414).')
main <- rodar(d26, 'fora_macro_2026', 'M3_MAPA_2026_N15414')

# Comparacao territorial somente onde AMBOS os mapas cobrem resid+tratamento.
res_compare <- list()
comparar_cenarios <- function(col_hist, etiqueta) {
  pareada <- copy(d0[!is.na(get(col_hist)) & !is.na(fora_macro_2026)])
  if (nrow(pareada)<10000L) {
     registrar('Mapa ', etiqueta, ': amostra pareada muito pequena ou incompleta; sem ajuste.')
     return(invisible(NULL))
  }
  a <- rodar(pareada, 'fora_macro_2026', paste0('M3_PAREADO_2026_VS_',etiqueta))
  b <- rodar(pareada, col_hist, paste0('M3_PAREADO_',etiqueta,'_VS_2026'))
  tot <- data.table(mapa=c('2026 pareado',etiqueta),
     n=c(nrow(pareada), nrow(pareada)),
     n_fora_macro=c(sum(pareada$fora_macro_2026 == 'Outra macrorregiao'),
       sum(pareada[[col_hist]] == 'Outra macrorregiao')),
     p_interacao_HTZ=c(a$global$p_val[1], b$global$p_val[1]))
  escrever(tot, paste0('COMPARACAO_PAREADA_',etiqueta,'_VS_2026'))
  registrar('Comparacao ',etiqueta,' vs 2026: n=',nrow(pareada),
           ', p_HTZ=',signif(tot$p_interacao_HTZ[2],4),' vs ',
           signif(tot$p_interacao_HTZ[1],4))
  invisible(tot)
}
if (!is.null(map24)) comparar_cenarios('fora_macro_2024', 'MAR2024')
if (!is.null(map23)) comparar_cenarios('fora_macro_2023', 'ANO2023_DOCUMENTADO')

# 6. PROVENIENCIA E NOTAS PARA MANUSCRITO -------------------------------------
if (file.exists(MAPA_2024)) {
  arq_fonte <- c(ARQ_COORTE, MAPA_2026, MAPA_2024)
} else arq_fonte <- c(ARQ_COORTE, MAPA_2026)
if (!is.null(map23)) arq_fonte <- c(arq_fonte, MAPA_2023)
md5 <- tools::md5sum(arq_fonte)
escrever(data.table(arquivo=names(md5), md5=as.character(md5)), 'ARQUIVOS_MD5')
proveniencia <- c(
  'Coorte: mulheres que iniciaram quimioterapia, nao todas as elegiveis.',
  'Mapas: 2026=portal de dados abertos SUS; marco/2024=Portaria GM/MS 3.233/2024.',
  paste('Portaria 2024:', URL_2024),
  'Portaria de 2024 NAO deve ser rotulada como delimitacao de 2023.',
  if (!is.null(map23)) paste('Fonte indicada para 2023:', FONTE_OFICIAL_2023)
    else 'Mapa de 2023 ainda sem correspondencia nacional documentada; obter de MS/CIB.',
  'Riscos absolutos: media das previsoes contrafactuais nas covariaveis observadas',
  'das mulheres da AMOSTRA ANALITICA; ajuste por idade spline fixa (3 df),',
  'estadio, procedimento, UF de residencia e mes de inicio.',
  'IC de riscos absolutos: delta logito, CR2, gl de Satterthwaite aproximados via gradiente.',
  'RR estratificados e testes de interacao: CR2/HTZ com graus de liberdade finitos.',
  'A classificacao territorial e os contrastes raciais sao ASSOCIACOES, nao efeitos causais.',
  'Nao usar p estratificados diferentes como prova de interacao.',
  'Familia de seis RR raciais exploratorios: ajuste de Holm nas seis comparacoes.',
  paste('Sessao R:', R.version.string)
)
writeLines(proveniencia, file.path(SAIDA, 'LEIA_ANTES_DE_INTERPRETAR.txt'), useBytes=TRUE)
writeLines(capture.output(utils::sessionInfo()),file.path(SAIDA,'SESSION_INFO.txt'))
registrar('\nCONCLUIDO. Resultados em: ', SAIDA)
registrar('Confira os arquivos M3_MAPA_2026_N15414_RISCOS_MARGINAIS_PADRONIZADOS.csv')
registrar('e COMPARACAO_PAREADA_MAR2024_VS_2026.csv (se a Portaria foi lida).')
