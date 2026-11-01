# Renders the server's own outputs for the tabs a headless screenshot cannot click (Sites, simulation), via
# shiny::testServer on the real app. Rscript tests/render_tabs.R  -> docs/screenshots/*.png, sim result text.
app <- shiny::shinyAppDir("app")
OBS <- jsonlite::fromJSON("app/defaults.json")$observed
shiny::testServer(app, args = list(), {
  OBS <- get("OBS", globalenv())
  session$setInputs(alpha = "0.05", power = 0.80, rrr = 0.20, prev = round(OBS$prevalence, 4),
                    p_hyper = round(OBS$p_hyper, 4), p_hypo = round(OBS$p_hypo, 4), sens = 0.85, spec = 0.85,
                    site = "restricted_setB", scoreable = 0.8, n_sims = 1000)
  for (id in c("curve_plot", "site_plot")) {
    src <- output[[id]]$src
    writeBin(jsonlite::base64_dec(sub("^data:image/png;base64,", "", src)),
             file.path("..", "docs", "screenshots", paste0(id, ".png")))
  }
  session$setInputs(run_sim = 1)
  cat(gsub("<[^>]+>", " ", output$sim_result$html), "\n")
})
