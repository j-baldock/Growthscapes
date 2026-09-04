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
quarto::quarto_render("CompareFixedHab.qmd", output_format = "html", cache_refresh = TRUE)


# quarto::quarto_render("Scenario_1_null_cold.qmd", output_format = "html", cache_refresh = TRUE)
# quarto::quarto_render("Scenario_2_temp_mult.qmd", output_format = "html", cache_refresh = TRUE)
# quarto::quarto_render("Scenario_3_temp_offset.qmd", output_format = "html", cache_refresh = TRUE)
# quarto::quarto_render("Scenario_4_temp_offset_diffP.qmd", output_format = "html", cache_refresh = TRUE)



