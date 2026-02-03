library(quarto)
library(knitr)

# render whole book
quarto::quarto_render(output_format = "html")

# render single chapters
quarto::quarto_render("index.qmd", output_format = "html")

