library(quarto)
library(knitr)

# render whole book
quarto::quarto_render(output_format = "html")

# render single chapters
quarto::quarto_render("index.qmd", output_format = "html")
quarto::quarto_render("DataViz.qmd", output_format = "html")
quarto::quarto_render("Breakpoints.qmd", output_format = "html")
quarto::quarto_render("FormatData.qmd", output_format = "html")
quarto::quarto_render("ModelTemp-Letcher.qmd", output_format = "html")
quarto::quarto_render("ModelTemp-EcoDrought.qmd", output_format = "html")
# quarto::quarto_render("ModelTemp-EcoDrought-StdBasin.qmd", output_format = "html")
# quarto::quarto_render("ModelTemp-EcoDrought-withCovs.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("ModelTemp-EcoDrought-ByMonth.qmd", output_format = "html", cache_refresh = TRUE)
