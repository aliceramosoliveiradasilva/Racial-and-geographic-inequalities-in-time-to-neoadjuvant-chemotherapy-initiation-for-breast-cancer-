# PUBLICATION REPRODUCIBILITY — BREAST CANCER NEOADJUVANT CHEMOTHERAPY
# Complete historical map 2024 and logistic predictive margins.
# Executes locally on authorized DATASUS data; never commit patient-level data.
# 2026-10-15. Historical crosswalks 2024/2026 are not validated for 2023.
# All R scripts should be run from project root; MAMA_RAIZ configures data root.
# NOTE: full end-to-end re-execution was not possible in the delivery environment.
# -------------------------------------------------------------------------

# ============================================================================
# 05 - AUDITORIA DO MAPA MAR/2024 + PADRONIZACAO LOGISTICA DOS RISCOS ABSOLUTOS
# Cancer de mama / quimioterapia neoadjuvante, coorte de 2023
# Elaborado em 14/10/2026
#
# PROPÓSITOS
# (1) Corrigir omissao na extracao da Portaria GM/MS 3.233/2024:
#     44 municipios da macro 2310 (Fortaleza, CE), documentados na versao
#     https://bvs.saude.gov.br/bvs/saudelegis/gm/2024/prt3233_04_03_2024.html
#     A versao _06_03_2024 usada anteriormente nao os trouxe ao CSV.
# (2) Comparar 2024 e 2026 na MESMA coorte de pacientes.
# (3) Manter Poisson modificado + CR2/HTZ para RR e interacao principal.
# (4) Calcular riscos absolutos ajustados com regressao LOGISTICA e
#     padronizacao marginal, garantindo previsoes entre 0 e 1.
#
# IMPORTANTE: portaria de 2024 NAO equivale a delimitacao validada de 2023.
# Os modelos sao associativos; nao sustentam afirmacoes causais.
# As saidas novas NAO sobrescrevem resultados de scripts anteriores.
# ============================================================================

RAIZ <- Sys.getenv("MAMA_RAIZ", unset = getwd())
DIR_MACROS <- file.path(RAIZ, 'dados_mama_SIA_AQ_2019_2025',
                       'resultados_com_macrorregioes', 'ANALISE_MACRORREGIOES')
ARQ_COORTE <- file.path(DIR_MACROS, 'coorte_analitica_2023_COM_MACROS.fst')
MAPA_2026 <- file.path(RAIZ, 'macroregiao_de_saude.csv')
# Opcao de compatibilidade caso o CSV esteja com sufixo (1).
if (!file.exists(MAPA_2026)) {
  outro_csv <- file.path(RAIZ, 'macroregiao_de_saude(1).csv')
  if (file.exists(outro_csv)) MAPA_2026 <- outro_csv
}
MAPA_2024_INCOMPLETO <- file.path(RAIZ, 'municipio_macrorregiao_PORTARIA_MAR2024.csv')
SAIDA <- file.path(DIR_MACROS, 'RESOLUCAO_PENDENCIAS_1_2_20261014',
                   'VALIDACAO_COMPLEMENTAR_05')
dir.create(SAIDA, showWarnings=FALSE, recursive=TRUE)

pacotes <- c('data.table','fst','clubSandwich','splines')
ausentes <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly=TRUE)]
if (length(ausentes)) stop('Instale os pacotes: ', paste(ausentes, collapse=', '))
suppressPackageStartupMessages(library(data.table))
registrar <- function(...) cat(..., '\n', sep='')
salvar <- function(x, nome) {
  data.table::fwrite(as.data.table(x), file.path(SAIDA, paste0(nome,'.csv')),
                     sep=';', na='')
}
mun6 <- function(x) {
  x <- trimws(as.character(x))
  x <- sub('\\.0+$', '', x)
  x[!grepl('^[0-9]{6,7}$', x)] <- NA_character_
  substr(x, 1, 6)
}
ler_mapa <- function(caminho) {
  if (!file.exists(caminho)) stop('Arquivo nao encontrado: ', caminho)
  x <- data.table::fread(caminho, colClasses='character', encoding='UTF-8')
  campos_mun <- intersect(c('cod_municipio','mun6','codigo_ibge'), names(x))
  campos_macro <- intersect(c('cod_macrorregiao_de_saude','macro','codigo_macro'), names(x))
  if (!length(campos_mun) || !length(campos_macro))
    stop('Cabecalhos municipais/macrorregionais ausentes em ', caminho)
  out <- unique(data.table(mun6=mun6(x[[campos_mun[1]]]),
                           macro=trimws(x[[campos_macro[1]]])))
  if (anyNA(out$mun6) || anyNA(out$macro) || any(!nzchar(out$macro)))
    stop('Chaves ou codigos nulos em ', caminho)
  if (anyDuplicated(out$mun6)) stop('Municipio repetido com macro conflitante em ', caminho)
  out
}

# 1. COMPLETAR A TABELA DE MARCO DE 2024, COM RASTREABILIDADE ------------------
map24 <- ler_mapa(MAPA_2024_INCOMPLETO)
map26 <- ler_mapa(MAPA_2026)
if (nrow(map26) != 5570L) stop('Mapa 2026 nao possui 5570 municipios.')

# Codigos IBGE conferidos na lista oficial da 1a Macro de Fortaleza da Portaria
# GM/MS 3.233 de 01/03/2024, Anexo II (versao _04_03_2024).
ce_2310_oficial <- c(
  '230015','230075','230090','230100','230120','230140','230195','230210',
  '230220','230290','230350','230370','230395','230428','230440','230460',
  '230495','230510','230523','230625','230630','230640','230650','230765',
  '230770','230837','230910','230945','230960','230970','230980','231010',
  '231020','231025','231070','231085','231160','231240','231260','231335',
  '231350','231355','231375','231380')

faltando <- setdiff(map26$mun6, map24$mun6)
nao_oficiais <- setdiff(faltando, ce_2310_oficial)
if (length(nao_oficiais)) stop('Ha municipios ausentes NAO validados no Anexo II: ',
                               paste(head(nao_oficiais, 20), collapse=', '))
if (nrow(map24) != 5570L && length(faltando) != 44L)
  stop('O mapa parcial de 2024 tem uma omissao diferente das 44 esperadas; auditar manualmente.')
if (length(faltando)) {
  macro26_faltantes <- map26$macro[match(faltando, map26$mun6)]
  if (!all(macro26_faltantes == '2310'))
    stop('Os municipios omitidos nao correspondem a macro 2310 no mapa 2026.')
  map24 <- rbindlist(list(map24, data.table(mun6=faltando, macro='2310')))
}
setorder(map24, mun6)
stopifnot(nrow(map24)==5570L, !anyDuplicated(map24$mun6),
          setequal(map24$mun6,map26$mun6))

salvar(map24[, .(cod_municipio=mun6,cod_macrorregiao_de_saude=macro)],
       'MAPA_MAR2024_COMPLETO_DOCUMENTADO')
salvar(data.table(cod_municipio=faltando, macro_marc2024='2310',
                  fonte='Portaria GM/MS 3233/2024, Anexo II, versao _04_03_2024'),
       'CORRECAO_44_MUNICIPIOS_CEARA')
registrar('Mapa MAR/2024 completado: ', nrow(map24), ' municipios; ',
          uniqueN(map24$macro), ' macrorregioes.')
registrar('Foram acrescentados ', length(faltando), ' municipios oficiais da macro 2310.')

# 2. VINCULACAO E AUDITORIA DE RECLASSIFICACAO -------------------------------
if (!file.exists(ARQ_COORTE)) stop('Coorte FST nao encontrada: ', ARQ_COORTE)
d <- as.data.table(fst::read_fst(ARQ_COORTE))
obrigatorios <- c('mun_res','mun_estab','raca_final','idade','estadio',
                  'procedimento','uf_res','mes_inicio','cnes','atraso60')
if (length(setdiff(obrigatorios,names(d))))
  stop('Faltam colunas: ', paste(setdiff(obrigatorios,names(d)),collapse=', '))
if (nrow(d) != 15415L) warning('Numero de mulheres difere da coorte original: ',nrow(d))
d[, `:=`(mun_res6=mun6(mun_res), mun_estab6=mun6(mun_estab))]
classificar <- function(mapa, etiqueta) {
  res <- mapa$macro[match(d$mun_res6,mapa$mun6)]
  trt <- mapa$macro[match(d$mun_estab6,mapa$mun6)]
  d[, (paste0('macro_res_',etiqueta)) := res]
  d[, (paste0('macro_trat_',etiqueta)) := trt]
  d[, (paste0('fora_',etiqueta)) := fifelse(!is.na(res)&!is.na(trt),
                   fifelse(res==trt,'Mesma macrorregiao','Outra macrorregiao'),
                   NA_character_)]
}
classificar(map24,'2024'); classificar(map26,'2026')
pareada <- copy(d[!is.na(fora_2024) & !is.na(fora_2026)])
registrar('Coorte original: ',nrow(d), '; cobertas em ambos os mapas: ',nrow(pareada))
if (nrow(pareada) != 15414L)
  warning('Amostra pareada difere da esperada (15414); verificar exclusoes.')
pareada[, mudou := fora_2024 != fora_2026]
salvar(pareada[, .N, by=.(raca_final, fora_2024, fora_2026, mudou)],
       'RECLASSIFICACAO_2024_VS_2026_POR_RACA')
salvar(data.table(coorte_original=nrow(d), completa_em_ambas=nrow(pareada),
                  mudaram_de_categoria=sum(pareada$mudou),
                  pct_reclassificacao=100*mean(pareada$mudou),
                  municipios_mapa_2024=nrow(map24),
                  macros_mapa_2024=uniqueN(map24$macro),
                  municipios_mapa_2026=nrow(map26),
                  macros_mapa_2026=uniqueN(map26$macro)), 'RESUMO_MAPAS_E_COBERTURA')

# 3. PREPARAR MESMA AMOSTRA E SPLINE PARA OS DOIS CENARIOS -------------------
b <- as.data.frame(pareada)
b$raca_final <- factor(as.character(b$raca_final),levels=c('Branca','Parda','Preta'))
for (coluna in c('estadio','procedimento','uf_res','mes_inicio'))
  b[[coluna]] <- factor(as.character(b[[coluna]]))
b$cnes <- as.character(b$cnes)
b$idade <- as.numeric(b$idade)
b$atraso60 <- as.integer(b$atraso60)
completas <- complete.cases(b[,c('raca_final','estadio','procedimento','uf_res',
                                  'mes_inicio','cnes','idade','atraso60','fora_2024','fora_2026')])
if (!all(completas)) {
  warning('Foram excluidas ',sum(!completas),' mulheres por covariaveis ausentes.')
  b <- b[completas,,drop=FALSE]
}
if (any(!b$atraso60 %in% 0:1)) stop('Desfecho nao binario.')
s <- splines::ns(b$idade,df=3)
for (j in seq_len(ncol(s))) b[[paste0('idade_ns',j)]] <- as.numeric(s[,j])
FORMULA <- atraso60 ~ raca_final * exposicao + idade_ns1 + idade_ns2 + idade_ns3 +
  estadio + procedimento + uf_res + mes_inicio

# Estimar gl para gradiente da padronizacao (aproximacao de primeira ordem).
gl_delta <- function(fit, V, grad) {
  # Wald_test exige o objeto CR2 original, com classe 'clubSandwich'.
  # Nunca passar aqui as.matrix(V): isso remove os atributos necessarios.
  if (!inherits(V, 'clubSandwich'))
    stop('gl_delta: V precisa ser o objeto original de clubSandwich::vcovCR().')
  mat <- matrix(unname(grad), nrow=1,
                dimnames=list(NULL,names(stats::coef(fit))))
  w <- tryCatch(
    as.data.frame(clubSandwich::Wald_test(fit, constraints=mat,
                                         vcov=V, test='HTZ')),
    error=function(e) stop('Falha ao estimar gl CR2 para o gradiente: ',
                           conditionMessage(e), call.=FALSE)
  )
  df <- as.numeric(w$df_denom[1])
  if (!is.finite(df) || df <= 0)
    stop('gl_delta: graus de liberdade nao finitos/nao positivos; conferir clusters.')
  df
}

padronizar_logistica <- function(fit,V,base,etiqueta) {
  beta <- stats::coef(fit)
  Vmat <- as.matrix(V)  # Somente algebra; preserva V (classe clubSandwich) para HTZ
  riscos <- list(); gradientes <- list()
  k <- 0L
  for (territorio in levels(base$exposicao)) for (raca in levels(base$raca_final)) {
    novo <- base
    novo$raca_final <- factor(rep(raca,nrow(base)),levels=levels(base$raca_final))
    novo$exposicao <- factor(rep(territorio,nrow(base)),levels=levels(base$exposicao))
    X <- stats::model.matrix(stats::delete.response(stats::terms(fit)),
                             data=novo,contrasts.arg=fit$contrasts)
    if (!setequal(colnames(X),names(beta))) stop('Matriz de previsao nao corresponde ao modelo.')
    X <- X[,names(beta),drop=FALSE]
    p <- stats::plogis(drop(X%*%beta))
    stopifnot(all(is.finite(p)),all(p >= 0),all(p <= 1))
    risco <- mean(p)
    grad <- drop(crossprod(X,p*(1-p))/length(p))
    names(grad) <- names(beta)
    glogit <- grad/(risco*(1-risco))
    ep_l <- sqrt(max(0,drop(t(glogit)%*%Vmat%*%glogit)))
    df <- gl_delta(fit,V,glogit)
    crit <- if (is.finite(df) && df>0) stats::qt(.975,df) else stats::qnorm(.975)
    IC <- stats::plogis(stats::qlogis(risco)+c(-1,1)*crit*ep_l)
    k <- k+1L
    riscos[[k]] <- data.table(mapa=etiqueta,raca=raca,territorio=territorio,
       n_padronizacao=nrow(base), risco_percentual=100*risco,
       IC95_inf_percentual=100*IC[1],IC95_sup_percentual=100*IC[2],
       df_aprox=df,max_predicao=max(p),n_pred_maior_1=sum(p>1),
       modelo='logistico M3; padronizacao marginal; delta CR2')
    gradientes[[paste(raca,territorio,sep='|')]] <- list(risco=risco,grad=grad)
  }
  salvar(rbindlist(riscos),paste0('MAPA_',etiqueta,'_RISCOS_MARGINAIS_LOGISTICOS'))
  out_dif <- list(); k <- 0L
  for (territorio in levels(base$exposicao)) for (raca in c('Parda','Preta')) {
    a <- gradientes[[paste(raca,territorio,sep='|')]]
    z <- gradientes[[paste('Branca',territorio,sep='|')]]
    grad <- a$grad-z$grad
    rd <- a$risco-z$risco
    ep <- sqrt(max(0,drop(t(grad)%*%Vmat%*%grad)))
    df <- gl_delta(fit,V,grad)
    crit <- if (is.finite(df)&&df>0) stats::qt(.975,df) else stats::qnorm(.975)
    k <- k+1L
    out_dif[[k]] <- data.table(mapa=etiqueta,territorio=territorio,
      contraste=paste0(raca,' menos Branca'),diferenca_pp=100*rd,
      IC95_inf_pp=100*(rd-crit*ep),IC95_sup_pp=100*(rd+crit*ep),df_aprox=df)
  }
  salvar(rbindlist(out_dif),paste0('MAPA_',etiqueta,'_DIFERENCAS_RISCO_LOGISTICAS_PP'))
  rbindlist(riscos)
}

testes_poisson <- function(fit,V,base,etiqueta) {
  nm <- names(stats::coef(fit))
  ind <- grep('^raca_final[^:]*:exposicao',nm)
  if (length(ind)!=2L) stop('A interacao nao possui dois coeficientes.')
  w <- as.data.table(as.data.frame(clubSandwich::Wald_test(fit,
          constraints=clubSandwich::constrain_zero(ind),vcov=V,test='HTZ')))
  w[, `:=`(mapa=etiqueta,n=nrow(base),cnes=uniqueN(base$cnes))]
  salvar(w,paste0('MAPA_',etiqueta,'_INTERACAO_POISSON_CR2_HTZ'))
  bcoef <- stats::coef(fit)
  Vmat <- as.matrix(V)  # So para calculo de EP; HTZ recebe V original
  out <- list(); k <- 0L
  L_raca <- function(raca,territorio) {
    L <- setNames(rep(0,length(bcoef)),nm)
    if (raca != 'Branca') {
      co <- paste0('raca_final',raca)
      L[co] <- 1
      if (territorio == levels(base$exposicao)[2]) {
        int <- paste0(co,':exposicao',levels(base$exposicao)[2])
        if (!(int %in% nm)) stop('Termo de interacao ausente: ',int)
        L[int] <- 1
      }
    }
    L
  }
  for (territorio in levels(base$exposicao)) {
    for (par in list(c('Parda','Branca'),c('Preta','Branca'),c('Preta','Parda'))) {
      L <- L_raca(par[1],territorio)-L_raca(par[2],territorio)
      est <- sum(L*bcoef)
      ep <- sqrt(max(0,drop(t(L)%*%Vmat%*%L)))
      mat <- matrix(unname(L),nrow=1,dimnames=list(NULL,nm))
      zz <- as.data.frame(clubSandwich::Wald_test(fit,constraints=mat,vcov=V,test='HTZ'))
      df <- as.numeric(zz$df_denom[1]); p <- as.numeric(zz$p_val[1])
      limite <- stats::qt(.975,df)*ep
      k <- k+1L
      out[[k]] <- data.table(mapa=etiqueta,territorio=territorio,
        contraste=paste0(par[1],'_vs_',par[2]),
        RR=exp(est),IC95_inf=exp(est-limite),IC95_sup=exp(est+limite),
        df_HTZ=df,p_HTZ=p)
    }
  }
  ct <- rbindlist(out)
  ct[, p_Holm_6 := stats::p.adjust(p_HTZ,method='holm')]
  salvar(ct,paste0('MAPA_',etiqueta,'_RR_POISSON_CR2_HTZ'))
  w
}

rodar_mapa <- function(etiqueta,exposicao_col) {
  base <- b
  base$exposicao <- factor(as.character(base[[exposicao_col]]),
                   levels=c('Mesma macrorregiao','Outra macrorregiao'))
  if (anyNA(base$exposicao) || any(table(base$raca_final,base$exposicao)<10L))
    stop('Falta exposicao ou existem estratos com menos de 10 mulheres.')
  registrar('Ajustando modelos ',etiqueta,'; N=',nrow(base))
  fit_p <- stats::glm(FORMULA,family=stats::poisson(link='log'),data=base)
  fit_l <- stats::glm(FORMULA,family=stats::binomial(link='logit'),data=base)
  if (!fit_p$converged||!fit_l$converged||anyNA(coef(fit_p))||anyNA(coef(fit_l)))
    stop('Modelo nao convergiu ou contem coeficientes nao estimaveis em ',etiqueta)
  Vp <- clubSandwich::vcovCR(fit_p,cluster=base$cnes,type='CR2')
  Vl <- clubSandwich::vcovCR(fit_l,cluster=base$cnes,type='CR2')
  w <- testes_poisson(fit_p,Vp,base,etiqueta)
  marg <- padronizar_logistica(fit_l,Vl,base,etiqueta)
  salvar(data.table(mapa=etiqueta,n=nrow(base),n_cnes=uniqueN(base$cnes),
            maior_predicao_poisson_observada=max(stats::fitted(fit_p)),
            n_poisson_observado_pred_maior_1=sum(stats::fitted(fit_p)>1),
            maior_predicao_logistica_observada=max(stats::fitted(fit_l)),
            maior_risco_logistico_padronizado=max(marg$risco_percentual),
            p_interacao_poisson_HTZ=w$p_val[1]),
         paste0('MAPA_',etiqueta,'_AUDITORIA'))
  invisible(list(teste=w,riscos=marg))
}
res24 <- rodar_mapa('MAR2024','fora_2024')
res26 <- rodar_mapa('OUT2026','fora_2026')

salvar(data.table(mapa=c('MAR2024','OUT2026'),n=nrow(b),
                  p_interacao_HTZ=c(res24$teste$p_val[1],res26$teste$p_val[1]),
                  maior_risco_logistico_ajustado_pct=c(
                    max(res24$riscos$risco_percentual),
                    max(res26$riscos$risco_percentual))),
       'COMPARACAO_M3_PAREADA_2024_2026_FINAL')
writeLines(c(
 'LEIA ANTES DE UTILIZAR NO MANUSCRITO:',
 '1. Base 2024 completada pela Portaria GM/MS 3.233/2024, Anexo II,',
 '   https://bvs.saude.gov.br/bvs/saudelegis/gm/2024/prt3233_04_03_2024.html',
 '   para os 44 municipios da macrorregiao 2310 omitidos na extracao anterior.',
 '2. A divisao de 2024 NAO comprova a vigente em 2023.',
 '3. RRs e testes globais: Poisson modificado com CR2/HTZ (escala multiplicativa).',
 '4. Riscos absolutos e diferencas de risco: logistica M3 padronizada marginalmente,',
 '   mesma distribuicao de covariaveis para as seis celulas, IC delta-CR2.',
 '5. Riscos logisticos e RRs Poisson sao estimandos de modelos diferentes;',
 '   nao apresentar como se viessem do mesmo modelo.',
 '6. Comparacoes raciais exploratorias Poisson: correcao de Holm para seis testes.',
 '7. Cuidado com extrapolacao nos estratos pouco frequentes (preta/outra macro).',
 '8. Nao interpretar variacao no p estratificado como prova de interacao.',
 paste('R version:',R.version.string)),
 file.path(SAIDA,'LEIA_ANTES_DE_PUBLICAR.txt'),useBytes=TRUE)
writeLines(capture.output(utils::sessionInfo()),file.path(SAIDA,'SESSION_INFO.txt'))
registrar('CONCLUIDO. Confira os CSVs em ',SAIDA)
