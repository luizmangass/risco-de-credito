# Gera Metodologia.md a partir de Metodologia.Rmd.
#
#     Rscript render.R
#
# Usa knitr::knit, que nao precisa de pandoc: a saida ja e Markdown, que o
# GitHub renderiza direto. Quem tiver pandoc pode gerar HTML com
# rmarkdown::render("Metodologia.Rmd").

if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("instale o knitr: install.packages(\"knitr\")")
}

knitr::knit("Metodologia.Rmd", output = "Metodologia.md", quiet = FALSE)

# Troca o cabecalho YAML por titulo e subtitulo em Markdown puro. O GitHub
# renderiza YAML de topo como uma tabela, o que ficaria feio.
linhas <- readLines("Metodologia.md", encoding = "UTF-8", warn = FALSE)
if (linhas[1] == "---") {
  fim <- which(linhas == "---")[2]
  linhas <- c("# Cálculo de Risco de Crédito",
              "",
              "**VaR e TVaR de crédito por migração de rating multiperíodo**",
              linhas[-seq_len(fim)])
}
writeLines(linhas, "Metodologia.md", useBytes = TRUE)
cat("Metodologia.md gerado com", length(linhas), "linhas\n")
