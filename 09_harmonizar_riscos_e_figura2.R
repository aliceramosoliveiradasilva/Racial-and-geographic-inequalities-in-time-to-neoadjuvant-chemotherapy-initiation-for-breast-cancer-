# PUBLICATION REPRODUCIBILITY — BREAST CANCER NEOADJUVANT CHEMOTHERAPY
# Harmonize municipal/macroregional margins. Figures PT and numeric CSV.
# Executes locally on authorized DATASUS data; never commit patient-level data.
# 2026-10-15. Historical crosswalks 2024/2026 are not validated for 2023.
# All R scripts should be run from project root; MAMA_RAIZ configures data root.
# NOTE: full end-to-end re-execution was not possible in the delivery environment.
# -------------------------------------------------------------------------

# =============================================================================
# 06 - HARMONIZAR FIGURA 2: MUNICIPIO VERSUS MACRORREGIAO DE SAUDE
# Cancer de mama e quimioterapia neoadjuvante no SUS, coorte de 2023
# Objetivo: riscos absolutos por logistica + padronizacao marginal (mesma regra
# usada no script 05), preservando Poisson/CR2/HTZ para a analise de RR.
#
# AMOSTRA PRINCIPAL MUNICIPAL: 15.415 mulheres, se todas as covariaveis completas.
# AMOSTRA PARA PAINEIS COMPARAVEIS: 15.414 mulheres com ambos os municipios
# vinculados ao mapa 2026, como no script 05. Sao resultados distintos.
# Os modelos sao ASSOCIATIVOS; nao estimam o efeito causal da raca/cor.
# Nao modifica nenhum arquivo de entrada e NAO substitui o mapa de 2023.
# =============================================================================

RAIZ <- Sys.getenv("MAMA_RAIZ", unset = getwd())
DIR_MACROS <- file.path(RAIZ, 'dados_mama_SIA_AQ_2019_2025',
                        'resultados_com_macrorregioes', 'ANALISE_MACRORREGIOES')
ARQ_COORTE <- file.path(DIR_MACROS, 'coorte_analitica_2023_COM_MACROS.fst')
MAPA_2026 <- file.path(RAIZ, 'macroregiao_de_saude.csv')
if (!file.exists(MAPA_2026)) {
  MAPA_2026 <- file.path(RAIZ, 'macroregiao_de_saude(1).csv')
}
ARQ_MACRO_LOGISTICA <- file.path(
  DIR_MACROS, 'RESOLUCAO_PENDENCIAS_1_2_20261014', 'VALIDACAO_COMPLEMENTAR_05',
  'MAPA_OUT2026_RISCOS_MARGINAIS_LOGISTICOS.csv')
SAIDA <- file.path(DIR_MACROS, 'FIGURA2_PADRONIZACAO_HARMONIZADA')
dir.create(SAIDA, recursive=TRUE, showWarnings=FALSE)

pacotes <- c('fst','data.table','clubSandwich','splines','ggplot2')
faltantes <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly=TRUE)]
if (length(faltantes)) {
  stop('Instale os pacotes ausentes: install.packages(c(',
       paste(sprintf('"%s"', faltantes), collapse=', '), '))')
}
suppressPackageStartupMessages(library(data.table))

salvar <- function(x,nome) data.table::fwrite(
  as.data.table(x), file.path(SAIDA,paste0(nome,'.csv')), sep=';',na='')
mun6 <- function(x) {
  s <- trimws(as.character(x))
  s <- sub('\\.0+$','',s)
  s[!grepl('^[0-9]{6,7}$',s)] <- NA_character_
  substr(s,1,6)
}
obrigar <- function(x,colunas,origem) {
  faltam <- setdiff(colunas,names(x))
  if(length(faltam)) stop(origem, ': colunas ausentes: ',paste(faltam,collapse=', '))
}
if (!file.exists(ARQ_COORTE)) stop('Coorte nao encontrada: ',ARQ_COORTE)
if (!file.exists(MAPA_2026)) stop('Mapa 2026 nao encontrado: ',MAPA_2026)

d <- as.data.table(fst::read_fst(ARQ_COORTE))
campos <- c('mun_res','mun_estab','fora_municipio','raca_final','idade',
            'estadio','procedimento','uf_res','mes_inicio','cnes','atraso60')
obrigar(d,campos,'Coorte de 2023')
if (nrow(d)!=15415L) warning('Coorte difere das 15.415 mulheres originais: ',nrow(d))

# Reproduzir exatamente o criterio territorial do script 05 para obter N pareado.
map <- data.table::fread(MAPA_2026,colClasses='character',showProgress=FALSE)
obrigar(map,c('cod_municipio','cod_macrorregiao_de_saude'),'Mapa 2026')
map[, mun6 := mun6(cod_municipio)]
if (anyNA(map$mun6) || anyDuplicated(map$mun6)) {
  stop('Mapa 2026 contem municipios sem codigo ou duplicados.')
}
i_res <- match(mun6(d$mun_res),map$mun6)
i_trat <- match(mun6(d$mun_estab),map$mun6)
com_macro <- !is.na(i_res) & !is.na(i_trat) &
  !is.na(map$cod_macrorregiao_de_saude[i_res]) &
  !is.na(map$cod_macrorregiao_de_saude[i_trat])
if (sum(com_macro)!=15414L) {
  warning('Amostra pareada N=',sum(com_macro),', esperado N=15414. '
          ,'Auditar diferencas em relacao ao script 05.')
}
salvar(data.table(N_original=nrow(d), N_pareado=sum(com_macro),
                  sem_vinculo=sum(!com_macro),
                  N_municipio_mesmo=sum(as.character(d$fora_municipio)=='Mesmo município',na.rm=TRUE),
                  N_municipio_outro=sum(as.character(d$fora_municipio)=='Outro município',na.rm=TRUE)),
       'AUDITORIA_AMOSTRAS')

# A normalizacao acima eh usada somente para formar a amostra pareada.
# O indicador municipal ja consta da coorte; sua semantica nao eh modificada.
preparar <- function(df) {
  b <- as.data.frame(data.table::copy(df))
  b$raca_final <- factor(as.character(b$raca_final), levels=c('Branca','Parda','Preta'))
  b$exposicao <- factor(as.character(b$fora_municipio),
                       levels=c('Mesmo município','Outro município'))
  for (c in c('estadio','procedimento','uf_res','mes_inicio')) {
    b[[c]] <- factor(as.character(b[[c]]))
  }
  b$cnes <- as.character(b$cnes)
  b$idade <- suppressWarnings(as.numeric(b$idade))
  b$atraso60 <- suppressWarnings(as.integer(b$atraso60))
  modelos <- c('raca_final','exposicao','estadio','procedimento',
              'uf_res','mes_inicio','cnes','idade','atraso60')
  n_pre <- nrow(b)
  b <- b[stats::complete.cases(b[,modelos,drop=FALSE]),,drop=FALSE]
  if (nrow(b)!=n_pre) stop('Existem ',n_pre-nrow(b),
                          ' mulheres com dados ausentes. Auditar antes do ajuste.')
  if (any(!b$atraso60 %in% c(0L,1L))) stop('atraso60 nao eh binaria.')
  if (any(table(b$raca_final,b$exposicao)<10)) stop('Celula com menos de 10 mulheres.')
  ns <- splines::ns(b$idade,df=3)
  for (j in seq_len(ncol(ns))) b[[paste0('idade_ns',j)]] <- as.numeric(ns[,j])
  b
}

f_m3 <- atraso60 ~ raca_final * exposicao + idade_ns1 + idade_ns2 + idade_ns3 +
  estadio + procedimento + uf_res + mes_inicio

# Graus de liberdade aproximados para um gradiente de funcao nao linear.
# Crucial: Wald_test() recebe a matriz V CR2 ORIGINAL, nao as.matrix(V).
df_cr2 <- function(fit,V,g) {
  C <- matrix(unname(g),nrow=1,dimnames=list(NULL,names(coef(fit))))
  w <- as.data.frame(clubSandwich::Wald_test(fit,
       constraints=C,vcov=V,test='HTZ'))
  ans <- as.numeric(w$df_denom[1])
  if (!is.finite(ans) || ans <= 0) stop('Graus de liberdade CR2 invalidos.')
  ans
}

padronizar <- function(fit,V,b,rotulo) {
  beta <- stats::coef(fit)
  Vm <- as.matrix(V)
  rows <- list(); gradients <- list(); k <- 0L
  for (tr in levels(b$exposicao)) for (r in levels(b$raca_final)) {
    novo <- b
    novo$raca_final <- factor(rep(r,nrow(b)),levels=levels(b$raca_final))
    novo$exposicao <- factor(rep(tr,nrow(b)),levels=levels(b$exposicao))
    X <- stats::model.matrix(stats::delete.response(stats::terms(fit)),
       data=novo,contrasts.arg=fit$contrasts)
    if (!setequal(colnames(X),names(beta))) {
      stop('Matriz de previsao nao coincide com os coeficientes.')
    }
    X <- X[,names(beta),drop=FALSE]
    pi <- stats::plogis(drop(X %*% beta))
    if (any(!is.finite(pi)) || any(pi < 0) || any(pi > 1)) {
      stop('Previsoes logisticamente invalidas.')
    }
    p <- mean(pi)
    grad <- drop(crossprod(X,pi*(1-pi))/length(pi))
    names(grad) <- names(beta)
    gl <- grad/(p*(1-p))
    ep <- sqrt(max(0,drop(t(gl) %*% Vm %*% gl)))
    df <- df_cr2(fit,V,gl)
    crit <- stats::qt(.975,df)
    ic <- stats::plogis(stats::qlogis(p)+c(-1,1)*crit*ep)
    k <- k+1L
    rows[[k]] <- data.table(analise=rotulo,raca=r,territorio=tr,
      n_padronizacao=nrow(b),risco_percentual=100*p,
      IC95_inf_percentual=100*ic[1],IC95_sup_percentual=100*ic[2],
      df_CR2_aproximado=df,max_probabilidade_individual=max(pi),
      metodo='Logistica M3; marginalizacao empirica; IC delta/CR2')
    gradients[[paste(r,tr,sep='|')]] <- list(p=p,g=grad)
  }
  res <- rbindlist(rows)
  salvar(res,paste0(rotulo,'_RISCOS_LOGISTICOS_MARGINAIS'))
  # Diferencas absolutas raciais em pontos percentuais.
  saidas <- list();k <- 0L
  for (tr in levels(b$exposicao)) for (r in c('Parda','Preta')) {
    a <- gradients[[paste(r,tr,sep='|')]]
    ref <- gradients[[paste('Branca',tr,sep='|')]]
    gd <- a$g-ref$g
    rd <- a$p-ref$p
    ep <- sqrt(max(0,drop(t(gd)%*%Vm%*%gd)))
    df <- df_cr2(fit,V,gd)
    c95 <- stats::qt(.975,df)
    k <- k+1L
    saidas[[k]] <- data.table(analise=rotulo,territorio=tr,
      contraste=paste(r,'menos Branca'),
      diferenca_pp=100*rd,IC95_inf_pp=100*(rd-c95*ep),
      IC95_sup_pp=100*(rd+c95*ep),df_CR2_aproximado=df)
  }
  salvar(rbindlist(saidas),paste0(rotulo,'_DIFERENCAS_RISCO_PP'))
  res
}

rodar <- function(df,rotulo) {
  b <- preparar(df)
  cat('Ajustando ',rotulo,', N=',nrow(b),', CNES=',length(unique(b$cnes)),'\n',sep='')
  fit_l <- stats::glm(f_m3, family=stats::binomial(link='logit'),data=b,
                      control=stats::glm.control(maxit=50))
  fit_p <- stats::glm(f_m3, family=stats::poisson(link='log'),data=b,
                      control=stats::glm.control(maxit=50))
  if (!fit_l$converged || !fit_p$converged || anyNA(coef(fit_l)) || anyNA(coef(fit_p))) {
    stop('Modelo nao convergiu ou teve coeficiente nao estimavel.')
  }
  V_l <- clubSandwich::vcovCR(fit_l,cluster=b$cnes,type='CR2')
  V_p <- clubSandwich::vcovCR(fit_p,cluster=b$cnes,type='CR2')
  marg <- padronizar(fit_l,V_l,b,rotulo)
  idx <- grep('^raca_final[^:]*:exposicao',names(coef(fit_p)))
  if (length(idx)!=2L) stop('Esperava dois termos de interacao racial.')
  ht <- as.data.table(as.data.frame(clubSandwich::Wald_test(
    fit_p,constraints=clubSandwich::constrain_zero(idx),vcov=V_p,test='HTZ')))
  ht[, `:=`(analise=rotulo,n=nrow(b),cnes=uniqueN(b$cnes))]
  salvar(ht,paste0(rotulo,'_INTERACAO_POISSON_CR2_HTZ'))
  list(riscos=marg,teste=ht)
}

# Primeiro, o modelo comparavel (mesma coorte do script 05).
pareado <- rodar(d[com_macro], 'MUNICIPIO_PAREADO_N15414')
# Depois, a estimativa municipal na coorte total original.
completo <- rodar(d, 'MUNICIPIO_ORIGINAL_N15415')

# Opcional: reutilizar a estimativa MACRO, ja calculada no script 05,
# para gerar Figura 2 com paineis comparaveis (N=15.414).
if (file.exists(ARQ_MACRO_LOGISTICA)) {
  macro <- data.table::fread(ARQ_MACRO_LOGISTICA)
  obrigar(macro,c('raca','territorio','risco_percentual',
      'IC95_inf_percentual','IC95_sup_percentual'),'Resultados logísticos 05')
  if (any(macro$n_padronizacao!=nrow(d[com_macro]))) {
    stop('N de padronizacao da macro nao corresponde ao N pareado.')
  }
  mun <- data.table::copy(pareado$riscos)
  mun[, `:=`(painel='A. Município',exposicao_fig=as.character(territorio))]
  macro[, `:=`(painel='B. Macrorregião',exposicao_fig=as.character(territorio))]
  mun[exposicao_fig=='Mesmo município',exposicao_fig:='Mesmo território']
  mun[exposicao_fig=='Outro município',exposicao_fig:='Outro território']
  macro[exposicao_fig=='Mesma macrorregiao',exposicao_fig:='Mesmo território']
  macro[exposicao_fig=='Outra macrorregiao',exposicao_fig:='Outro território']
  cols_fig <- c('painel','raca','exposicao_fig','risco_percentual',
                'IC95_inf_percentual','IC95_sup_percentual')
  fig <- rbindlist(list(mun[,..cols_fig],macro[,..cols_fig]))
  fig[, raca := factor(raca,levels=c('Branca','Parda','Preta'))]
  fig[, exposicao_fig := factor(exposicao_fig,
                               levels=c('Mesmo território','Outro território'))]
  salvar(fig,'FIGURA2_DADOS_PADRONIZADOS_PAREADOS')
  pos <- ggplot2::position_dodge(width=0.45)
  grafico <- ggplot2::ggplot(fig,ggplot2::aes(x=raca,y=risco_percentual,
                color=exposicao_fig,group=exposicao_fig))+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=IC95_inf_percentual,
                    ymax=IC95_sup_percentual),width=0.10,position=pos,linewidth=0.75)+
    ggplot2::geom_point(position=pos,size=2.7)+
    ggplot2::facet_wrap(~painel,nrow=1)+
    ggplot2::scale_y_continuous(breaks=seq(0,100,20),limits=c(0,100),
                                labels=function(x) paste0(x,'%'))+
    ggplot2::labs(x='Raça/cor registrada',y='Risco ajustado de início após 60 dias',
       color='Local do tratamento',
       caption='N = 15.414 em ambos os painéis; logística M3 com padronização marginal; IC95% delta/CR2.')+
    ggplot2::theme_bw(base_size=11)+
    ggplot2::theme(legend.position='bottom',panel.grid.minor=ggplot2::element_blank(),
       plot.caption=ggplot2::element_text(size=8,hjust=0))
  for (ext in c('png','pdf','jpeg')) {
    nome <- file.path(SAIDA,paste0('FIGURA2_MUNICIPIO_VS_MACRO_N15414.',ext))
    ggplot2::ggsave(nome,plot=grafico,width=11,height=5.4,dpi=500,bg='white')
  }
} else {
  warning('CSV do script 05 nao encontrado. Foram criados os riscos municipais; ',
          'rode o script 05 para obter o grafico com dois paineis.')
}

writeLines(c(
 'AUDITORIA FIGURA 2:',
 paste('Coorte FST:',ARQ_COORTE),
 paste('Mapa apenas para formar amostra pareada:',MAPA_2026),
 paste('Analise municipal original N=',nrow(d)),
 paste('Analise municipal pareada N=',sum(com_macro)),
 'Modelo de riscos absolutos: logistica M3, spline natural idade (df=3),',
 'estadio, procedimento, UF de residencia, mes e raca/cor x deslocamento.',
 'Padronizacao empirica dentro de cada amostra, para todas as seis celulas.',
 'IC95%: metodo delta com CR2 e graus de liberdade de primeira ordem.',
 'RR/teste global: Poisson modificado com CR2/HTZ (estimando distinto).',
 'Modelo nao ajusta por distancia real, custo ou disponibilidade do servico.',
 'Figura com dois paineis usa exclusivamente N pareado=15414.',
 'Classificacao 2026 usada somente como referencia territorial; NAO prova 2023.',
 paste('R:',R.version.string)),
 file.path(SAIDA,'LEIA_ANTES_DA_FIGURA2.txt'),useBytes=TRUE)
writeLines(capture.output(utils::sessionInfo()),
           file.path(SAIDA,'SESSION_INFO.txt'))
cat('CONCLUIDO. Resultados em: ',SAIDA,'\n',sep='')
