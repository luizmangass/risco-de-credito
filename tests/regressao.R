###############################################################################
# tests/regressao.R -- validacao do modelo de risco de credito
#
# Rodar da raiz do repositorio:
#     Rscript tests/regressao.R
#
# Sai com status 1 se qualquer verificacao falhar. Nao depende de pacote algum.
#
# Os testes 3 a 8 sao os que importam: verificam propriedades que o modelo TEM
# de ter por construcao, comparando a simulacao com o valor teorico. Se o
# sorteio de rating ou a contabilizacao da perda estiverem errados, eles
# quebram.
###############################################################################

raiz <- if (file.exists("credito.R")) "." else ".."
source(file.path(raiz, "credito.R"))
prem <- ler_premissas(file.path(raiz, "dados"))

falhas <- 0L
ok <- function(nome, cond, detalhe = "") {
  cond <- isTRUE(cond)
  cat(sprintf("%-4s %-58s %s\n", if (cond) "[ok]" else "[XX]", nome, detalhe))
  if (!cond) falhas <<- falhas + 1L
  invisible(cond)
}
cabec <- function(x) cat("\n== ", x, " ", strrep("=", max(0, 60 - nchar(x))), "\n", sep = "")

M   <- prem$matriz
lim <- prem$limiares
K   <- length(ESTADOS)

# ---------------------------------------------------------------------------
cabec("1. Matriz de transicao")
ok("linhas somam 1",            max(abs(rowSums(M) - 1)) < 1e-12)
ok("probabilidades em [0,1]",   all(M >= 0 & M <= 1))
ok("D e absorvente",            M["D", "D"] == 1 && sum(M["D", -1]) == 0)
ok("PD de 1 ano e monotona no rating",
   !is.unsorted(prem$pd_1a[c("Aaa", "Aa", "A", "Baa", "Ba", "B", "Caa_C")]),
   paste(sprintf("%.3f%%", 100 * prem$pd_1a[c("Aaa", "Baa", "Caa_C")]), collapse = " "))
ok("massa na diagonal > 60% em todo rating",
   all(diag(M)[-1] > 0.60), sprintf("min = %.3f", min(diag(M)[-1])))

cabec("2. Limiares de CreditMetrics")
ok("limiares crescentes dentro da linha",
   all(apply(lim, 1, function(x) !is.unsorted(x))))
ok("ultimo limiar = Inf",       all(is.infinite(lim[, K]) & lim[, K] > 0))
ok("limiar de default = qnorm(PD)",
   max(abs(lim[-1, 1] - qnorm(M[-1, "D"])), na.rm = TRUE) < 1e-12)

# ---------------------------------------------------------------------------
# O teste central do sorteio: simular muitas vezes a partir de cada rating e
# conferir se a frequencia observada dos ratings de destino reproduz a linha da
# matriz. Tolerancia = 4 desvios-padrao binomiais (falso alarme ~ 6e-5 por
# celula), com piso para nao acusar celulas de probabilidade minuscula.
# ---------------------------------------------------------------------------
cabec("3. risk_discrete reproduz a matriz (rho = 0)")
set.seed(20240101)
n <- 200000L
for (r in ESTADOS[-1]) {
  i    <- match(r, ESTADOS)
  novo <- risk_discrete(rep.int(i, n), lim, rho = 0)$rating
  freq <- tabulate(novo, K) / n
  alvo <- M[i, ]
  tol  <- pmax(4 * sqrt(alvo * (1 - alvo) / n), 2e-4)
  ok(sprintf("linha %-6s reproduzida", r), all(abs(freq - alvo) <= tol),
     sprintf("max desvio = %.5f", max(abs(freq - alvo))))
}

# ---------------------------------------------------------------------------
# Com rho > 0 e o fator sistemico Z integrado, a distribuicao MARGINAL do
# rating de destino tem de continuar sendo a mesma linha da matriz. E o que
# garante que introduzir correlacao muda a cauda sem mexer na media.
# ---------------------------------------------------------------------------
cabec("4. Correlacao nao distorce a distribuicao marginal")
set.seed(20240102)
for (r in c("Baa", "Ba", "B")) {
  i    <- match(r, ESTADOS)
  z    <- rnorm(n)                       # um Z por sorteio = Z integrado
  novo <- risk_discrete(rep.int(i, n), lim, rho = 0.24, z = z)$rating
  freq <- tabulate(novo, K) / n
  alvo <- M[i, ]
  ok(sprintf("linha %-6s com rho = 0.24", r),
     all(abs(freq - alvo) <= pmax(4 * sqrt(alvo * (1 - alvo) / n), 2e-4)),
     sprintf("max desvio = %.5f", max(abs(freq - alvo))))
}

cabec("5. Sorteio condicional a sobrevivencia")
set.seed(20240103)
i    <- match("B", ESTADOS)
alvo <- M[i, ]; alvo[1] <- 0; alvo <- alvo / sum(alvo)

# Sem fator sistemico, o sorteio condicional tem de reproduzir exatamente a
# linha da matriz renormalizada sem a coluna D.
s0   <- risk_discrete(rep.int(i, n), lim, rho = 0, condicional = TRUE)
freq <- tabulate(s0$rating, K) / n
ok("nunca sorteia D",           !any(s0$rating == 1L))
ok("reproduz a linha renormalizada (rho = 0)",
   all(abs(freq - alvo) <= pmax(4 * sqrt(alvo * (1 - alvo) / n), 2e-4)),
   sprintf("max desvio = %.5f", max(abs(freq - alvo))))

# Com rho > 0 a marginal do sorteio condicional NAO e a linha renormalizada, e
# nem deveria ser: a renormalizacao e por 1 - q(Z), que depende do cenario, e
# essa covariancia com Z desloca a marginal. O que tem de valer -- e e o que a
# igualdade das medias entre os dois modos usa -- e que a PD do periodo tenha
# media igual a PD incondicional da matriz.
st <- risk_discrete(rep.int(i, n), lim, rho = 0.2, z = rnorm(n), condicional = TRUE)
ok("nunca sorteia D (rho = 0,2)", !any(st$rating == 1L))
ok("E[PD do periodo] = PD incondicional da matriz",
   abs(mean(st$pd) - M[i, "D"]) < 4 * sd(st$pd) / sqrt(n),
   sprintf("%.5f vs %.5f", mean(st$pd), M[i, "D"]))

# ---------------------------------------------------------------------------
cabec("6. Contabilidade da perda")
ct <- carteira_aleatoria(120, seed = 7)
s  <- simular_carteira(ct, prem, n_cen = 3000, seed = 1, bloco = 1000,
                       guardar_fluxo = TRUE)
ok("fluxo por periodo soma a perda do cenario",
   max(abs(rowSums(s$fluxo) - s$perdas)) < 1e-6)
ok("media do fluxo bate com a EL",
   abs(sum(s$fluxo_medio) - mean(s$perdas)) < 1e-6)
ok("EL por ativo soma a EL da carteira",
   abs(sum(s$el_ativo) - mean(s$perdas)) < 1e-6)
ok("perda nunca excede exposicao x LGD",
   max(s$perdas) <= sum(s$carteira$posicao * s$carteira$lgd) + 1e-6)
ok("perda sempre nao negativa", all(s$perdas >= 0))

# ---------------------------------------------------------------------------
# O invariante mais forte do modelo: os dois modos tem a MESMA perda esperada.
# Ver a demonstracao no cabecalho de credito.R.
# ---------------------------------------------------------------------------
cabec("7. Equivalencia das medias entre os dois modos")
nc <- 20000
a  <- simular_carteira(ct, prem, n_cen = nc, modo = "perda_esperada",
                       seed = 99, bloco = 2000)
b  <- simular_carteira(ct, prem, n_cen = nc, modo = "default",
                       seed = 99, bloco = 2000)
ep <- sqrt(var(a$perdas) / nc + var(b$perdas) / nc)
ok("E[perda] igual nos dois modos (3 erros-padrao)",
   abs(mean(a$perdas) - mean(b$perdas)) < 3 * ep,
   sprintf("dif = %.2f%% da EL", 100 * (mean(b$perdas) / mean(a$perdas) - 1)))
ok("modo perda_esperada tem variancia menor",
   sd(a$perdas) < sd(b$perdas),
   sprintf("sd %.3g vs %.3g", sd(a$perdas), sd(b$perdas)))
ok("PD media do caminho tambem coincide",
   abs(mean(a$pd_caminho) - mean(b$pd_caminho)) < 0.01,
   sprintf("%.4f vs %.4f", mean(a$pd_caminho), mean(b$pd_caminho)))

# ---------------------------------------------------------------------------
cabec("8. Efeito da correlacao")
r0 <- simular_carteira(ct, prem, n_cen = nc, rho = 0.00, seed = 5, bloco = 2000)
r2 <- simular_carteira(ct, prem, n_cen = nc, rho = 0.20, seed = 5, bloco = 2000)
r4 <- simular_carteira(ct, prem, n_cen = nc, rho = 0.40, seed = 5, bloco = 2000)
epc <- sqrt(var(r0$perdas) / nc + var(r4$perdas) / nc)
ok("EL nao depende de rho (3 erros-padrao)",
   abs(mean(r0$perdas) - mean(r4$perdas)) < 3 * epc,
   sprintf("dif = %.2f%%", 100 * (mean(r4$perdas) / mean(r0$perdas) - 1)))
q <- function(x) unname(quantile(x$perdas, 0.995))
ok("VaR 99,5% cresce com rho",
   q(r0) < q(r2) && q(r2) < q(r4),
   sprintf("%.3g < %.3g < %.3g", q(r0), q(r2), q(r4)))
ok("capital economico cresce com rho",
   (q(r0) - mean(r0$perdas)) < (q(r4) - mean(r4$perdas)))

# ---------------------------------------------------------------------------
cabec("9. Horizonte, duration e amortizacao")
base <- data.frame(ativo = "Debenture", rating = "BBB", duration = 5,
                   posicao = 1e6, stringsAsFactors = FALSE)
p15 <- preparar_carteira(transform(base, duration = 15), prem)
p10 <- preparar_carteira(transform(base, duration = 10), prem)
p0  <- preparar_carteira(transform(base, duration = 0),  prem)
ok("duration 15 e truncada em 10 anos", p15$periodos == 10 && p10$periodos == 10)
ok("duration 0 assume o teto de 10 anos", p0$periodos == 10)
ok("duration 0,4 vira 1 periodo",
   preparar_carteira(transform(base, duration = 0.4), prem)$periodos == 1)
e5  <- mean(simular_carteira(base, prem, n_cen = 8000, seed = 3)$perdas)
e10 <- mean(simular_carteira(transform(base, duration = 10), prem,
                             n_cen = 8000, seed = 3)$perdas)
ok("prazo maior gera perda esperada maior", e10 > e5,
   sprintf("%.0f vs %.0f", e5, e10))
ok("amortizacao linear reduz a EX vs bullet",
   mean(simular_carteira(transform(base, duration = 10, perfil = "linear"),
                         prem, n_cen = 8000, seed = 3)$perdas) < e10)
ok("desconto reduz a EL",
   mean(simular_carteira(transform(base, duration = 10), prem, n_cen = 8000,
                         seed = 3, taxa_desconto = 0.10)$perdas) < e10)

cabec("10. Ativo ja inadimplente")
dd <- data.frame(ativo = "Debenture", rating = "D", duration = 5,
                 posicao = 1e6, stringsAsFactors = FALSE)
pd_lgd <- prem$lgd$lgd_basileia[prem$lgd$ativo == "Debenture"] * 1e6
ok("rating D perde exposicao x LGD (modo perda_esperada)",
   all(abs(simular_carteira(dd, prem, n_cen = 200, seed = 2)$perdas - pd_lgd) < 1e-6))
ok("rating D perde exposicao x LGD (modo default)",
   all(abs(simular_carteira(dd, prem, n_cen = 200, seed = 2,
                            modo = "default")$perdas - pd_lgd) < 1e-6))

# ---------------------------------------------------------------------------
cabec("11. PD simulada x curva publicada da Moody's")
rr <- c("Aa", "A", "Baa", "Ba", "B", "Caa_C")
cv <- data.frame(ativo = "Debenture", rating = rr, duration = 1, posicao = 1,
                 stringsAsFactors = FALSE)
sv <- simular_carteira(cv, prem, n_cen = 40000, rho = 0, seed = 13, bloco = 4000)
tab1 <- prem$pd_acum[match(rr, ESTADOS), 1]
raz  <- sv$pd_caminho / tab1
ok("PD de 1 ano dentro de 15% da tabela publicada",
   all(abs(raz - 1) < 0.15),
   sprintf("razao %.2f a %.2f", min(raz), max(raz)))
# Em 10 anos a cadeia de Markov e a tabela de coortes divergem por construcao
# (a matriz de 1 ano nao reproduz a dinamica de longo prazo dos ratings).
# Aqui so exigimos que a ordenacao por rating continue de pe.
cv10 <- transform(cv, duration = 10)
s10  <- simular_carteira(cv10, prem, n_cen = 20000, rho = 0, seed = 14, bloco = 4000)
ok("PD de 10 anos monotona no rating", !is.unsorted(s10$pd_caminho),
   paste(sprintf("%.1f%%", 100 * s10$pd_caminho), collapse = " "))
ok("PD de 10 anos > PD de 1 ano em todo rating", all(s10$pd_caminho > sv$pd_caminho))

cabec("12. Metricas e reprodutibilidade")
m <- metricas_risco(a$perdas, exposicao = a$exposicao)
ok("VaR crescente no nivel de confianca", !is.unsorted(m$VaR))
ok("TVaR > VaR em todo nivel",  all(m$TVaR > m$VaR))
ok("VaR 95% > EL",              m$VaR[1] > m$EL[1])
ok("capital = VaR - EL",        max(abs(m$capital - (m$VaR - m$EL))) < 1e-9)
x1 <- simular_carteira(ct, prem, n_cen = 1000, seed = 4242, bloco = 250)$perdas
x2 <- simular_carteira(ct, prem, n_cen = 1000, seed = 4242, bloco = 250)$perdas
ok("mesma semente, mesmo resultado", identical(x1, x2))
ok("blocos nao alteram o resultado esperado",
   abs(mean(x1) - mean(simular_carteira(ct, prem, n_cen = 20000, seed = 9,
                                        bloco = 5000)$perdas)) <
   4 * sd(x1) / sqrt(1000))

cabec("13. Correlacao de Basileia")
pds <- c(0.0003, 0.003, 0.03, 0.15)
rb  <- rho_basileia(pds)
ok("rho decresce com a PD",     !is.unsorted(rev(rb)))
ok("rho dentro de [0,12; 0,24]", all(rb >= 0.12 - 1e-9 & rb <= 0.24 + 1e-9),
   paste(sprintf("%.3f", rb), collapse = " "))

cabec("14. Erros esperados")
ok("rating desconhecido gera erro",
   inherits(try(preparar_carteira(transform(base, rating = "ZZZ"), prem),
                silent = TRUE), "try-error"))
ok("ativo sem LGD gera erro",
   inherits(try(preparar_carteira(transform(base, ativo = "Nao_Existe"), prem),
                silent = TRUE), "try-error"))
ok("coluna faltando gera erro",
   inherits(try(preparar_carteira(base[, c("ativo", "rating")], prem),
                silent = TRUE), "try-error"))

cat("\n", strrep("-", 68), "\n", sep = "")
if (falhas == 0L) {
  cat("TODOS OS TESTES PASSARAM\n")
} else {
  cat(sprintf("%d TESTE(S) FALHARAM\n", falhas))
  quit(status = 1L)
}
