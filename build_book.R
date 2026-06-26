library(quarto)
library(knitr)

# render whole book
quarto::quarto_render(output_format = "html", cache_refresh = TRUE)

# render single chapters
quarto::quarto_render("index.qmd", output_format = "html")
quarto::quarto_render("Habitat.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("Functions.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("GrowthRegime.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("IBM_simulation.qmd", output_format = "html", cache_refresh = TRUE)
