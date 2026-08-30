###############################################################################
# credito.R -- VaR de credito por migracao de rating multiperiodo
#
# O QUE O MODELO FAZ
# -----------------
# Cada ativo de credito e projetado ano a ano ate a sua duration, limitada a 10
# anos (horizonte maximo das curvas de PD acumulada publicadas). A cada periodo
# o rating do ativo e RESSORTEADO por risk_discrete() a partir da linha da
# matriz de transicao correspondente ao rating vigente. A perda esperada
# mantem a forma classica
#
#     EL_i = Exposicao_i * PD_i * LGD_i
#
# so que a PD deixa de ser um numero de tabela e passa a ser a PD acumulada ao
# longo do caminho de rating efetivamente sorteado:
#
#     PD_i(caminho) = 1 - prod_{t=1..T_i} ( 1 - pd(R_i(t-1)) )
#
# Repetindo isso em N cenarios sai a distribuicao empirica da perda agregada.
# VaR e TVaR sao percentis dessa distribuicao -- sem ajuste parametrico
# intermediario, basta ordenar e cortar no percentil desejado (99,5%).
#
# OS DOIS MODOS
# -------------
# modo = "perda_esperada"  (padrao)
#     A migracao e sorteada da matriz CONDICIONAL A SOBREVIVENCIA (a linha
#     renormalizada sem a coluna D) e a inadimplencia entra pelo valor
#     esperado, via pd() do rating vigente. O ativo nunca "pula" para D: ele
#     acumula EL periodo a periodo. Formato de perda esperada preservado.
#
# modo = "default"
#     A migracao e sorteada da matriz COMPLETA, com D absorvente. Quando o
#     caminho entra em D o ativo perde Exposicao * LGD naquele periodo. E o
#     modelo default-mode classico (CreditMetrics).
#
# Os dois nao sao alternativas soltas: E[perda] e IDENTICA nos dois modos, por
# construcao. Escrevendo Mc para a matriz condicional a sobrevivencia e q_r
# para a PD de um periodo do rating r,
#
#     prod_t M[R_{t-1}, R_t]  =  prod_t Mc[R_{t-1}, R_t] * (1 - q_{R_{t-1}})
#
# logo P(sobreviver ate T) = E_Mc[ prod_t (1 - q) ], ou seja
#
#     E_Mc[ PD(caminho) ] = 1 - P(sobreviver) = P(default ate T) = E[1{default}]
#
# O modo "perda_esperada" e portanto a versao Rao-Blackwellizada do modo
# "default": mesma media, variancia menor, porque integra analiticamente o
# ruido idiossincratico do evento de default e deixa na distribuicao so o
# risco de MIGRACAO. tests/regressao.R verifica essa igualdade numericamente.
#
# O FATOR SISTEMICO
# -----------------
# Sorteio independente por ativo faz a lei dos grandes numeros achatar a cauda:
# o percentil 99,5% de uma carteira grande fica colado na media e o numero
# perde sentido economico. O sorteio aqui e feito na representacao de limiares
# (CreditMetrics): o rating de destino e determinado por onde cai
#
#     X_{i,t} = sqrt(rho_i) * Z_t + sqrt(1 - rho_i) * eps_{i,t}
#
# com Z_t ~ N(0,1) COMUM a todos os ativos naquele periodo/cenario e eps
# idiossincratico. rho = 0 devolve exatamente o sorteio multinomial
# independente da linha da matriz; rho > 0 correlaciona as migracoes. O padrao
# usa a correlacao de ativos do IRB de Basileia, calculada por ativo a partir
# da sua PD de um ano.
#
# DEPENDENCIAS: nenhuma. R base, >= 3.5.
###############################################################################

# Ordem canonica dos estados: do PIOR para o MELHOR. E a ordem do eixo normal
# na representacao de limiares -- D ocupa a cauda esquerda. Todos os indices de
# rating no codigo se referem a esta ordem.
ESTADOS <- c("D", "Caa_C", "B", "Ba", "Baa", "A", "Aa", "Aaa")

# ---------------------------------------------------------------------------
# PREMISSAS
# ---------------------------------------------------------------------------

# Le a matriz de transicao publicada (em %, com coluna WR), remove a coluna WR
# redistribuindo-a proporcionalmente entre os estados observados, acrescenta a
# linha de D como estado absorvente e reordena tudo para ESTADOS.
normalizar_matriz <- function(bruta) {
  rn <- bruta$rating
  m  <- as.matrix(bruta[, setdiff(names(bruta), c("rating", "WR")), drop = FALSE])
  rownames(m) <- rn
  m <- m / rowSums(m)                       # tira WR e renormaliza a linha

  M <- matrix(0, length(ESTADOS), length(ESTADOS),
              dimnames = list(ESTADOS, ESTADOS))
  M["D", "D"] <- 1                          # default e absorvente
  M[rownames(m), colnames(m)] <- m
  M <- M[ESTADOS, ESTADOS, drop = FALSE]

  stopifnot(all(abs(rowSums(M) - 1) < 1e-9), all(M >= 0))
  M
}

# Limiares normais no estilo CreditMetrics. lim[r, k] e o ponto do eixo normal
# padrao abaixo do qual o ativo de rating r termina em um estado de indice <= k.
montar_limiares <- function(M) {
  lim <- t(apply(M, 1, function(p) stats::qnorm(cumsum(p))))
  lim[, ncol(lim)] <- Inf                   # protege contra cumsum < 1 por FP
  dimnames(lim) <- dimnames(M)
  lim
}

ler_premissas <- function(dir = "dados") {
  ler <- function(f) utils::read.csv(file.path(dir, f), stringsAsFactors = FALSE,
                                     check.names = FALSE)

  bruta <- ler("matriz_transicao_moodys.csv")
  M     <- normalizar_matriz(bruta)

  pd <- ler("pd_acumulada_moodys.csv")
  pd_acum <- as.matrix(pd[, paste0("a", 1:10)]) / 100
  rownames(pd_acum) <- pd$rating
  # D ja inadimplente: PD acumulada 1 em qualquer horizonte.
  pd_acum <- rbind(pd_acum, D = rep(1, 10))[ESTADOS, , drop = FALSE]

  list(
    matriz_bruta = bruta,
    matriz       = M,
    limiares     = montar_limiares(M),
    pd_acum      = pd_acum,
    pd_1a        = M[, "D"],                # PD de 1 ano implicita na matriz
    lgd          = ler("lgd_basileia.csv"),
    equiv        = ler("equivalencia_ratings.csv")
  )
}

# Correlacao de ativos do IRB de Basileia para exposicoes corporativas.
rho_basileia <- function(pd) {
  pd <- pmin(pmax(pd, 1e-6), 1)
  k  <- (1 - exp(-50 * pd)) / (1 - exp(-50))
  0.12 * k + 0.24 * (1 - k)
}

# Traduz o rating da carteira para a escala do provedor das premissas.
mapear_rating <- function(rating, equiv) {
  r <- trimws(as.character(rating))
  i <- match(r, equiv$rating_origem)
  if (anyNA(i)) {
    stop("rating sem equivalencia em dados/equivalencia_ratings.csv: ",
         paste(unique(r[is.na(i)]), collapse = ", "))
  }
  equiv$rating_moodys[i]
}

# ---------------------------------------------------------------------------
# SORTEIO DO NOVO RATING
# ---------------------------------------------------------------------------
# Um passo da cadeia, vetorizado. `rating` e um vetor de indices em ESTADOS.
#
#   condicional = FALSE -> sorteia da linha completa da matriz; D pode sair.
#   condicional = TRUE  -> sorteia da linha condicional a sobrevivencia; D nunca
#                          sai, e a PD do periodo volta em $pd para entrar na EL.
#
# Devolve sempre $pd = probabilidade de default do periodo, condicional ao
# fator sistemico z. Com rho = 0 ela e simplesmente M[rating, "D"].
risk_discrete <- function(rating, lim, rho = 0, z = 0, condicional = FALSE,
                          u = NULL) {
  n <- length(rating)
  if (length(rho) == 1L) rho <- rep.int(rho, n)
  if (length(z)   == 1L) z   <- rep.int(z, n)

  sr <- sqrt(rho)
  s1 <- sqrt(1 - rho)
  L  <- lim[rating, , drop = FALSE]             # n x K
  q  <- stats::pnorm((L[, 1L] - sr * z) / s1)   # PD do periodo dado z

  if (is.null(u)) u <- stats::runif(n)
  w <- if (condicional) q + u * (1 - q) else u
  w <- pmin(pmax(w, 1e-15), 1 - 1e-15)
  x <- sr * z + s1 * stats::qnorm(w)            # variavel latente do ativo

  novo <- rowSums(x > L[, -ncol(L), drop = FALSE]) + 1L
  list(rating = novo, pd = q)
}

# ---------------------------------------------------------------------------
# PREPARO DA CARTEIRA
# ---------------------------------------------------------------------------
# Entrada minima: data.frame com ativo, rating, duration, posicao.
# Colunas opcionais: id, lgd (sobrepoe a tabela), perfil ("bullet"/"linear").
preparar_carteira <- function(carteira, prem, horizonte_max = 10,
                              fonte_lgd = c("basileia", "moodys")) {
  fonte_lgd <- match.arg(fonte_lgd)
  obrig <- c("ativo", "rating", "duration", "posicao")
  falta <- setdiff(obrig, names(carteira))
  if (length(falta)) stop("carteira sem as colunas: ", paste(falta, collapse = ", "))

  d <- carteira
  if (is.null(d$id))     d$id     <- seq_len(nrow(d))
  if (is.null(d$perfil)) d$perfil <- "bullet"

  d$rating_moodys <- mapear_rating(d$rating, prem$equiv)
  d$rating_idx    <- match(d$rating_moodys, ESTADOS)

  # Duration vira numero inteiro de periodos anuais, no minimo 1 e no maximo o
  # horizonte das curvas de PD. Ativo sem duration contratual (0 ou NA) entra
  # pelo teto -- premissa conservadora, herdada do modelo original.
  dur <- ifelse(is.na(d$duration) | d$duration <= 0, horizonte_max, round(d$duration))
  d$periodos <- pmin(pmax(1, dur), horizonte_max)

  if (is.null(d$lgd)) {
    col <- if (fonte_lgd == "basileia") "lgd_basileia" else "lgd_moodys"
    i   <- match(d$ativo, prem$lgd$ativo)
    if (anyNA(i)) stop("ativo sem LGD em dados/lgd_basileia.csv: ",
                       paste(unique(d$ativo[is.na(i)]), collapse = ", "))
    d$lgd <- prem$lgd[[col]][i]
  }
  d
}

# Matriz T x n com a exposicao viva de cada ativo em cada periodo.
# bullet: principal no vencimento, exposicao cheia ate T_i.
# linear: amortizacao constante, exposicao media do periodo.
matriz_exposicao <- function(d, Tmax) {
  E <- matrix(0, Tmax, nrow(d))
  for (j in seq_len(nrow(d))) {
    Tj <- d$periodos[j]
    E[seq_len(Tj), j] <- if (identical(d$perfil[j], "linear")) {
      d$posicao[j] * (Tj - seq_len(Tj) + 1) / Tj
    } else {
      d$posicao[j]
    }
  }
  E
}

# ---------------------------------------------------------------------------
# SIMULACAO
# ---------------------------------------------------------------------------
simular_carteira <- function(carteira, prem,
                             n_cen         = 10000,
                             modo          = c("perda_esperada", "default"),
                             rho           = "basileia",
                             horizonte_max = 10,
                             taxa_desconto = 0,
                             fonte_lgd     = c("basileia", "moodys"),
                             bloco         = 500,
                             seed          = NULL,
                             guardar_fluxo = FALSE) {
  modo <- match.arg(modo)
  if (!is.null(seed)) set.seed(seed)

  d0   <- preparar_carteira(carteira, prem, horizonte_max, match.arg(fonte_lgd))
  nat  <- nrow(d0)
  Tmax <- max(d0$periodos)
  lim  <- prem$limiares

  # Ordena por prazo decrescente. Assim, no periodo t, os ativos ainda vivos
  # sao exatamente as PRIMEIRAS n_vivo[t] colunas e o laco encolhe conforme a
  # carteira vence, em vez de simular ativos ja liquidados.
  ord  <- order(d0$periodos, decreasing = TRUE)
  d    <- d0[ord, , drop = FALSE]
  n_vivo <- vapply(seq_len(Tmax), function(t) sum(d$periodos >= t), integer(1))

  rho_ativo <- if (identical(rho, "basileia")) {
                 unname(rho_basileia(prem$pd_1a[d$rating_idx]))
               } else {
                 rep_len(as.numeric(rho), nat)
               }
  rho_ativo <- pmin(pmax(rho_ativo, 0), 0.999)

  Emat <- matriz_exposicao(d, Tmax)
  desc <- 1 / (1 + taxa_desconto)^seq_len(Tmax)

  perdas   <- numeric(n_cen)
  fluxo_ac <- numeric(Tmax)                            # soma da EL por periodo
  fluxo_q  <- if (guardar_fluxo) matrix(0, n_cen, Tmax) else NULL
  el_ativo <- numeric(nat)                             # soma da EL por ativo
  pd_ativo <- numeric(nat)                             # soma da PD do caminho

  cortes <- split(seq_len(n_cen), ceiling(seq_len(n_cen) / bloco))
  for (ce in cortes) {
    nb <- length(ce)
    R  <- matrix(d$rating_idx, nb, nat, byrow = TRUE)
    S  <- matrix(1, nb, nat)                           # sobrevivencia acumulada
    ok <- matrix(TRUE, nb, nat)                        # ainda nao deu default
    el_bloco <- matrix(0, nb, Tmax)

    for (t in seq_len(Tmax)) {
      cols <- seq_len(n_vivo[t])
      Rt   <- R[, cols, drop = FALSE]
      RHO  <- matrix(rho_ativo[cols], nb, length(cols), byrow = TRUE)
      Z    <- matrix(stats::rnorm(nb), nb, length(cols))  # fator sistemico do ano
      EXP  <- matrix(Emat[t, cols],    nb, length(cols), byrow = TRUE)
      LGD  <- matrix(d$lgd[cols],      nb, length(cols), byrow = TRUE)

      st <- risk_discrete(as.vector(Rt), lim, as.vector(RHO), as.vector(Z),
                          condicional = (modo == "perda_esperada"))
      Rn <- matrix(st$rating, nb, length(cols))

      if (modo == "perda_esperada") {
        q  <- matrix(st$pd, nb, length(cols))
        el <- S[, cols, drop = FALSE] * q * EXP * LGD * desc[t]
        S[, cols] <- S[, cols, drop = FALSE] * (1 - q)
      } else {
        novo_d <- ok[, cols, drop = FALSE] & (Rn == 1L)   # entrou em D neste ano
        el     <- novo_d * EXP * LGD * desc[t]
        ok[, cols] <- ok[, cols, drop = FALSE] & !(Rn == 1L)
      }

      el_bloco[, t]  <- rowSums(el)
      el_ativo[cols] <- el_ativo[cols] + colSums(el)
      R[, cols] <- Rn
    }

    perdas[ce] <- rowSums(el_bloco)
    fluxo_ac   <- fluxo_ac + colSums(el_bloco)
    if (guardar_fluxo) fluxo_q[ce, ] <- el_bloco
    pd_ativo <- pd_ativo + if (modo == "perda_esperada") colSums(1 - S)
                           else colSums(!ok)
  }

  inv <- order(ord)                                    # devolve na ordem original
  list(
    perdas      = perdas,
    fluxo_medio = fluxo_ac / n_cen,
    fluxo       = fluxo_q,
    el_ativo    = (el_ativo / n_cen)[inv],
    pd_caminho  = (pd_ativo / n_cen)[inv],
    carteira    = d0,
    rho         = rho_ativo[inv],
    modo        = modo,
    n_cen       = n_cen,
    horizonte   = Tmax,
    exposicao   = sum(d0$posicao)
  )
}

# ---------------------------------------------------------------------------
# METRICAS
# ---------------------------------------------------------------------------
# VaR: percentil da distribuicao empirica das perdas.
# TVaR: media das perdas que excedem o VaR.
# Capital: perda inesperada, VaR - EL, na convencao de Basileia/Solvencia II.
metricas_risco <- function(perdas, niveis = c(0.95, 0.99, 0.995),
                           exposicao = NA_real_) {
  el <- mean(perdas)
  out <- do.call(rbind, lapply(niveis, function(a) {
    v <- unname(stats::quantile(perdas, a, type = 7))
    cauda <- perdas[perdas > v]
    if (!length(cauda)) cauda <- perdas[perdas >= v]
    data.frame(nivel = a, VaR = v, TVaR = mean(cauda), capital = v - el)
  }))
  data.frame(out,
             EL           = el,
             desvio       = stats::sd(perdas),
             VaR_pct_expo = out$VaR / exposicao,
             row.names    = NULL)
}

# Benchmark deterministico: EL com a PD acumulada de tabela no rating inicial.
# E o calculo do modelo antigo, sem simulacao. Serve de ancora de sanidade.
el_analitico <- function(carteira, prem, horizonte_max = 10,
                         fonte_lgd = c("basileia", "moodys")) {
  d <- preparar_carteira(carteira, prem, horizonte_max, match.arg(fonte_lgd))
  d$pd_tabela <- prem$pd_acum[cbind(d$rating_idx, d$periodos)]
  d$EL        <- d$posicao * d$pd_tabela * d$lgd
  d
}

# ---------------------------------------------------------------------------
# CARTEIRA SINTETICA -- so para validar o modelo sem dado proprietario
# ---------------------------------------------------------------------------
carteira_aleatoria <- function(n = 200, seed = 42) {
  set.seed(seed)
  ativos <- c("Debenture", "Debenture_Subordinada", "Nota_Comercial", "Bond",
              "Letra_Financeira", "Letra_Financeira_Sub", "CRI", "CRA",
              "FIDC_Senior", "FIDC_Mezanino", "FIDC_Subordinada")
  peso_ativo <- c(30, 4, 6, 5, 18, 3, 9, 8, 9, 4, 4)

  ratings <- c("AAA", "AA+", "AA", "AA-", "A+", "A", "A-", "BBB+", "BBB",
               "BBB-", "BB+", "BB", "BB-", "B+", "B", "B-", "CCC", "ND")
  peso_rating <- c(3, 4, 8, 9, 11, 13, 12, 10, 8, 6, 4, 3, 3, 2, 1.5, 1, 1.5, 1)

  ativo <- sample(ativos, n, TRUE, prob = peso_ativo)
  data.frame(
    id       = sprintf("ATV%04d", seq_len(n)),
    ativo    = ativo,
    rating   = sample(ratings, n, TRUE, prob = peso_rating),
    # duration ate 12 anos de proposito: exercita o teto de 10.
    duration = ifelse(grepl("^FIDC", ativo) & runif(n) < 0.25, 0,
                      pmin(12, round(rgamma(n, shape = 2.2, scale = 1.7), 1))),
    posicao  = round(rlnorm(n, meanlog = 16.1, sdlog = 1.0), 2),
    perfil   = ifelse(grepl("^(CRI|CRA|FIDC)", ativo), "linear", "bullet"),
    stringsAsFactors = FALSE
  )
}

# Formatacao em reais, para os relatorios.
brl <- function(x, dig = 0) {
  paste0("R$ ", formatC(x, format = "f", digits = dig, big.mark = ".",
                        decimal.mark = ","))
}
pct <- function(x, dig = 2) {
  paste0(formatC(100 * x, format = "f", digits = dig, decimal.mark = ","), "%")
}
num <- function(x, dig = 2) {
  formatC(x, format = "f", digits = dig, big.mark = ".", decimal.mark = ",")
}
