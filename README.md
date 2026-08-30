# Cálculo de Risco de Crédito

Modelo de **VaR e TVaR de crédito** em R: cada ativo da carteira é projetado ano
a ano até a sua *duration* (teto de 10 anos), o rating é ressorteado a cada
período pela matriz de transição, a perda esperada é acumulada ao longo do
caminho percorrido e o percentil 99,5% sai direto da distribuição empírica das
perdas agregadas.

📄 **[Metodologia completa, com resultados](Metodologia.md)** — o documento
executável, com as premissas, a formulação, os gráficos e a validação.

## O que o modelo faz

1. **Projeta o fluxo de cada ativo** em períodos anuais até
   `min(duration, 10)`. Dez anos é o teto porque é até onde vai a curva de PD
   acumulada publicada, sem extrapolação. Ativo sem *duration* contratual entra
   pelo teto.
2. **Ressorteia o rating a cada período** com `risk_discrete()`, a partir da
   linha da matriz de transição do rating **vigente** — não do rating original.
3. **Mantém a forma da perda esperada**, `EL = Exposição × PD × LGD`, com a PD
   passando a ser a acumulada ao longo do caminho de rating sorteado:
   `PD(caminho) = 1 − Π(1 − q(R_t))`.
4. **Simula N cenários e ordena.** VaR e TVaR são percentis da distribuição
   empírica. Não há ajuste de distribuição teórica no meio do caminho.

## Instalação e uso

Não há dependências: o motor é R base (>= 3.5). Só o documento precisa de
`knitr`.

```r
source("credito.R")
prem     <- ler_premissas("dados")
carteira <- carteira_aleatoria(n = 250, seed = 42)   # ou a sua carteira

sim <- simular_carteira(carteira, prem, n_cen = 20000, seed = 2024)
metricas_risco(sim$perdas, niveis = c(0.95, 0.99, 0.995),
               exposicao = sim$exposicao)
```

A carteira é um `data.frame` com as colunas obrigatórias `ativo`, `rating`,
`duration` e `posicao`, e as opcionais `id`, `lgd` e `perfil`
(`"bullet"` ou `"linear"`).

### Principais argumentos de `simular_carteira()`

| Argumento | Padrão | O que faz |
| --- | --- | --- |
| `n_cen` | 10000 | número de cenários |
| `modo` | `"perda_esperada"` | ver abaixo |
| `rho` | `"basileia"` | correlação de ativos: `"basileia"` (por ativo, via PD) ou um número |
| `horizonte_max` | 10 | teto de períodos de projeção |
| `taxa_desconto` | 0 | desconta a perda de cada ano |
| `fonte_lgd` | `"basileia"` | `"basileia"` ou `"moodys"` |
| `seed` | `NULL` | reprodutibilidade |
| `guardar_fluxo` | `FALSE` | devolve a matriz cenário × ano |

### Os dois modos

- **`"perda_esperada"`** (padrão): a migração é sorteada da matriz condicional
  à sobrevivência e a inadimplência entra pelo valor esperado. O ativo acumula
  EL período a período. É o formato de perda esperada preservado.
- **`"default"`**: a migração é sorteada da matriz completa, com D absorvente;
  quando o caminho entra em D o ativo perde `Exposição × LGD`. É o
  *default-mode* clássico do CreditMetrics.

A perda esperada é **idêntica** nos dois modos — não por acaso, por construção
(a demonstração está no cabeçalho de [`credito.R`](credito.R) e a verificação
numérica é o teste 7). O modo padrão é a versão Rao-Blackwellizada do outro:
mesma média, variância menor, porque integra analiticamente o ruído do evento
de default e deixa na distribuição apenas o risco de migração.

## Premissas

Todas em [`dados/`](dados), como CSV. Trocar de fonte é editar arquivo —
nenhuma linha de código muda.

| Arquivo | Conteúdo | Fonte |
| --- | --- | --- |
| `matriz_transicao_moodys.csv` | migração média anual por letra cheia | Moody's, *Default and Recovery Rates of Corporate Bond Issuers, 1920–2004*, Exhibit 31 |
| `pd_acumulada_moodys.csv` | PD acumulada, 1 a 10 anos | idem, Exhibit 17 |
| `lgd_basileia.csv` | LGD supervisória por tipo de ativo, com a alternativa empírica da Moody's | BCBS d424 (Basel III), CRE32; recuperação da Moody's, Exhibit 27 |
| `equivalencia_ratings.csv` | escala S&P/Fitch → escala Moody's | — |

## Validação

```bash
Rscript tests/regressao.R
```

45 verificações, sem dependência de pacote. As que de fato prendem o modelo:

- 200 mil sorteios por rating reproduzem a linha da matriz **célula a célula**,
  dentro de 4 desvios binomiais;
- com correlação positiva e o fator sistêmico integrado, a distribuição
  marginal continua sendo a mesma linha;
- a soma do fluxo anual é igual à perda do cenário, ao centavo (identidade de
  telescopagem);
- a perda esperada dos dois modos coincide dentro de 3 erros-padrão;
- a EL independe de ρ, mas o VaR 99,5% cresce com ρ;
- a PD simulada de 1 ano fica a menos de 15% da curva publicada.

## Gerar o documento

```bash
Rscript render.R
```

Produz [`Metodologia.md`](Metodologia.md) a partir de `Metodologia.Rmd` via
`knitr::knit()`, que não precisa de pandoc — o GitHub renderiza o Markdown
direto. Quem tiver pandoc pode gerar HTML com
`rmarkdown::render("Metodologia.Rmd")`.

## Limitações conhecidas

A matriz de um ano não reproduz a curva de default de 10 anos observada: a
migração real tem memória e a hipótese markoviana subestima o default do grau
de investimento alto e superestima o do grau especulativo. A correlação é o
parâmetro dominante e é uma premissa regulatória, não uma medição. A LGD é
determinística, sem correlação com a PD. Não há risco de *spread*. A seção 6
da [metodologia](Metodologia.md) detalha cada ponto.

## Arquivos

| Arquivo | O que é |
| --- | --- |
| [`credito.R`](credito.R) | o modelo inteiro, em R base |
| [`Metodologia.Rmd`](Metodologia.Rmd) | documento executável |
| [`Metodologia.md`](Metodologia.md) | saída renderizada |
| [`render.R`](render.R) | gera o `.md` |
| [`tests/regressao.R`](tests/regressao.R) | bateria de validação |
| [`dados/`](dados) | premissas em CSV |
