library(quarto)
library(knitr)

# render whole book
quarto::quarto_render(output_format = "html", cache_refresh = FALSE)

# render single chapters
quarto::quarto_render("index.qmd", output_format = "html")

quarto::quarto_render("Habitat.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("Functions.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("SimulationLoop.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("PlotResults.qmd", output_format = "html", cache_refresh = TRUE)

quarto::quarto_render("PreCalGrowth.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("GrowthRegime.qmd", output_format = "html", cache_refresh = TRUE)

quarto::quarto_render("ScenarioDefs.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("SensitivityAnalysis.qmd", output_format = "html", cache_refresh = FALSE)

quarto::quarto_render("ViewSimResults.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("CompareOutcomes.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("MechanisticDrivers.qmd", output_format = "html", cache_refresh = TRUE)
quarto::quarto_render("SurplusProduction.qmd", output_format = "html", cache_refresh = TRUE)



