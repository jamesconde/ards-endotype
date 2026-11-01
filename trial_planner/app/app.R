# Enrichment trial planner: how many patients a prognostically enriched trial needs, and how many a site must
# screen to enrol them, given the classifier's accuracy and how many patients the site can actually score.
# Calculations: R/enrichment.R (closed form, parity-tested against the ARDS project's Python) and R/simulate.R
# (patient-level Monte Carlo). Shiny sources R/ automatically.

library(shiny)
library(bslib)
library(ggplot2)
library(markdown)   # used inside includeMarkdown(); named here so rsconnect's dependency scan deploys it

DEF <- jsonlite::fromJSON("defaults.json", simplifyVector = FALSE)
OBS <- DEF$observed
PRESETS <- c(setNames(names(DEF$labels), unlist(DEF$labels)), "Custom" = "custom")

INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8984"; GRID <- "#e6e5e0"; SURFACE <- "#fcfcfb"
SERIES <- "#2a78d6"; RECESSIVE <- "#c9c8c2"
SIM_BUDGET <- 3e7          # screened patients per simulation run (~3 s); sims are capped to stay under it

theme_planner <- function() {
  theme_minimal(base_size = 13) +
    theme(text = element_text(colour = INK2), axis.text = element_text(colour = INK2),
          plot.title = element_text(colour = INK, face = "bold", size = 14),
          plot.subtitle = element_text(colour = INK2), panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.4),
          plot.background = element_rect(fill = SURFACE, colour = NA), plot.title.position = "plot")
}

pct <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")
num <- function(x) ifelse(is.finite(x), scales::comma(round(x)), "not runnable")

ui <- page_sidebar(
  title = "Enrichment trial planner",
  theme = bs_theme(version = 5, bg = SURFACE, fg = INK, primary = SERIES),
  sidebar = sidebar(
    width = 330,
    accordion(
      open = c("Classifier", "Site"),
      accordion_panel("Trial design",
        selectInput("alpha", "Two-sided alpha", c(0.01, 0.025, 0.05, 0.10), selected = 0.05),
        sliderInput("power", "Power", 0.70, 0.95, 0.80, step = 0.01),
        sliderInput("rrr", "Relative risk reduction (assumed)", 0.05, 0.50, 0.20, step = 0.01)),
      accordion_panel("Cohort",
        numericInput("prev", "Share of patients in the high-risk class", round(OBS$prevalence, 4), 0.01, 0.99, 0.005),
        numericInput("p_hyper", "Event rate, high-risk class", round(OBS$p_hyper, 4), 0.001, 0.999, 0.005),
        numericInput("p_hypo", "Event rate, low-risk class", round(OBS$p_hypo, 4), 0.001, 0.999, 0.005)),
      accordion_panel("Classifier",
        sliderInput("sens", "Sensitivity", 0.50, 1.00, DEF$constants$acc_assumed, step = 0.01),
        sliderInput("spec", "Specificity", 0.50, 1.00, DEF$constants$acc_assumed, step = 0.01)),
      accordion_panel("Site",
        selectInput("site", "Where you screen", PRESETS, selected = "restricted_setB"),
        conditionalPanel("input.site == 'custom'",
          sliderInput("scoreable", "Share of screened patients the classifier can score", 0.01, 1, 0.8,
                      step = 0.01)),
        uiOutput("site_note"))
    ),
    actionButton("reset", "Reset to the study cohort", class = "btn-outline-secondary btn-sm")
  ),
  layout_column_wrap(
    width = 1 / 4, fill = FALSE,
    value_box("Patients per arm", textOutput("vb_n"), p(textOutput("vb_n_sub", inline = TRUE))),
    value_box("Patients to screen", textOutput("vb_screen"), p(textOutput("vb_screen_sub", inline = TRUE))),
    value_box("Screened per enrolled", textOutput("vb_ratio"), p(textOutput("vb_ratio_sub", inline = TRUE))),
    value_box("Enrolled event rate", textOutput("vb_rate"), p(textOutput("vb_rate_sub", inline = TRUE)))
  ),
  navset_card_tab(
    nav_panel("Accuracy", plotOutput("curve_plot", height = "420px")),
    nav_panel("Sites", plotOutput("site_plot", height = "420px")),
    nav_panel("Table",
      p("Required size as classifier accuracy falls (sensitivity = specificity), at full scoring coverage."),
      tableOutput("curve_table"), downloadButton("dl_curve", "Download CSV", class = "btn-sm")),
    nav_panel("Check by simulation",
      p("Simulates whole trials patient by patient: screening, classification, enrolment, randomisation and the",
        "final test. It never uses the formulas above, so agreement is evidence that the sample size is right."),
      layout_columns(col_widths = c(4, 8),
        div(numericInput("n_sims", "Simulated trials", 1000, 100, 5000, 100),
            actionButton("run_sim", "Run simulation", class = "btn-primary")),
        div(uiOutput("sim_result")))),
    nav_panel("About", includeMarkdown("about.md"),
      h5("Where the default numbers come from"), tableOutput("provenance"))
  )
)

# Inputs display the cohort values to 4 decimals; while an input still shows its default, use the exact value so
# the defaults reproduce the ARDS project's numbers (3,596 enrolled, 12,494 screened at system B).
exact <- function(x, full) if (isTRUE(abs(x - round(full, 4)) < 1e-12)) full else x

server <- function(input, output, session) {
  prev <- reactive(exact(input$prev, OBS$prevalence))
  p_hyper <- reactive(exact(input$p_hyper, OBS$p_hyper))
  p_hypo <- reactive(exact(input$p_hypo, OBS$p_hypo))
  scoreable <- reactive({
    if (input$site == "custom") input$scoreable else DEF$scoreable[[input$site]]
  })

  inputs_ok <- reactive({
    validate(
      need(is.finite(prev()) && prev() > 0 && prev() < 1, "Share in the high-risk class must be in (0, 1)."),
      need(is.finite(p_hyper()) && p_hyper() > 0 && p_hyper() < 1, "High-risk event rate must be in (0, 1)."),
      need(is.finite(p_hypo()) && p_hypo() > 0 && p_hypo() < 1, "Low-risk event rate must be in (0, 1)."),
      need(p_hyper() > p_hypo(),
           "The high-risk class must have the higher event rate, or enriching on it would be backwards."))
    TRUE
  })

  d <- reactive({
    inputs_ok()
    design(prev(), p_hyper(), p_hypo(), input$sens, input$spec, scoreable(), rrr = input$rrr,
           alpha = as.numeric(input$alpha), power = input$power)
  })

  output$vb_n <- renderText(num(d()$n_per_arm))
  output$vb_n_sub <- renderText(sprintf("%s fewer than unstratified (%s)", pct(1 - d()$vs_unstratified),
                                        num(d()$n_unstrat)))
  output$vb_screen <- renderText(num(d()$patients_to_screen))
  output$vb_screen_sub <- renderText(sprintf("to enrol %s, scoring %s of patients", num(d()$n_total),
                                             pct(scoreable())))
  output$vb_ratio <- renderText(formatC(d()$screen_per_enrolled, format = "f", digits = 2))
  output$vb_ratio_sub <- renderText(sprintf("%s screen positive", pct(d()$frac_positive)))
  output$vb_rate <- renderText(pct(d()$p_enrolled))
  output$vb_rate_sub <- renderText(sprintf("PPV %s; unstratified %s", pct(d()$ppv),
    pct(prev() * p_hyper() + (1 - prev()) * p_hypo())))

  output$site_note <- renderUI({
    if (input$site == "custom") return(NULL)
    p(class = "small text-muted", sprintf("%s of patients at this site have the variables the model needs.",
                                          pct(scoreable())))
  })

  curve <- reactive({
    inputs_ok()
    enrichment_curve(prev(), p_hyper(), p_hypo(), acc = seq(0.5, 1, by = 0.01), rrr = input$rrr,
                     alpha = as.numeric(input$alpha), power = input$power)
  })

  output$curve_plot <- renderPlot({
    cv <- curve(); dd <- d()
    g <- ggplot(cv, aes(classifier_accuracy, n_per_arm)) +
      geom_hline(yintercept = dd$n_unstrat, colour = MUTED, linewidth = 0.6, linetype = "22") +
      annotate("text", x = 1, y = dd$n_unstrat, label = paste("Unstratified", num(dd$n_unstrat)),
               hjust = 1, vjust = -0.6, colour = INK2, size = 3.8) +
      geom_line(colour = SERIES, linewidth = 1) +
      scale_x_continuous("Classifier accuracy (sensitivity = specificity)", labels = scales::percent) +
      scale_y_continuous("Patients per arm", labels = scales::comma, expand = expansion(mult = c(0.05, 0.12))) +
      labs(title = "Required trial size as the classifier gets worse",
           subtitle = "Misclassified patients dilute the enrolled event rate toward the low-risk class") +
      theme_planner()
    if (abs(input$sens - input$spec) < 1e-9) {
      g <- g + geom_point(data = data.frame(x = input$sens, y = dd$n_per_arm), aes(x, y), colour = SERIES,
                          fill = SURFACE, size = 3.5, shape = 21, stroke = 1.5) +
        annotate("text", x = input$sens, y = dd$n_per_arm, label = paste("Your design", num(dd$n_per_arm)),
                 hjust = 1.15, vjust = -1.2, colour = INK, size = 3.8)
    } else {
      g <- g + labs(caption = sprintf("Your classifier is asymmetric (sens %s, spec %s): %s per arm.",
                                      pct(input$sens, 0), pct(input$spec, 0), num(dd$n_per_arm)))
    }
    g
  }, bg = SURFACE)

  sites <- reactive({
    inputs_ok()
    keys <- names(DEF$labels)
    fs <- vapply(keys, function(k) DEF$scoreable[[k]], numeric(1))
    lab <- unlist(DEF$labels)
    if (input$site == "custom") { keys <- c(keys, "custom"); fs <- c(fs, input$scoreable); lab <- c(lab, "Custom site") }
    scr <- vapply(fs, function(f) design(prev(), p_hyper(), p_hypo(), input$sens, input$spec, f,
                                         input$rrr, as.numeric(input$alpha), input$power)$patients_to_screen,
                  numeric(1))
    data.frame(key = keys, label = factor(lab, levels = lab[order(scr)]), fs = fs, screen = scr,
               selected = keys == input$site)
  })

  output$site_plot <- renderPlot({
    s <- sites()
    ggplot(s, aes(screen, label)) +
      geom_col(aes(fill = selected), width = 0.62) +
      geom_text(aes(label = sprintf("%s  (scores %s)", num(screen), pct(fs))), hjust = -0.06, colour = INK,
                size = 3.8) +
      scale_fill_manual(values = c(`TRUE` = SERIES, `FALSE` = RECESSIVE), guide = "none") +
      scale_x_log10("Patients to screen (log scale)", labels = scales::comma) +
      coord_cartesian(xlim = c(min(s$screen) / 4, max(s$screen) * 40)) +
      labs(y = NULL, title = "Same trial, different hospitals",
           subtitle = "A patient the classifier cannot score cannot be enrolled") +
      theme_planner() + theme(panel.grid.major.y = element_blank())
  }, bg = SURFACE)

  curve_display <- reactive({
    cv <- round_curve(enrichment_curve(prev(), p_hyper(), p_hypo(), rrr = input$rrr,
                                       alpha = as.numeric(input$alpha), power = input$power))
    names(cv) <- c("Accuracy", "PPV", "Enrolled event rate", "N per arm", "N total", "Screened per enrolled",
                   "Screened total", "N vs unstratified")
    cv
  })
  output$curve_table <- renderTable(curve_display(), digits = 3)
  output$dl_curve <- downloadHandler("enrichment_curve.csv", function(f) write.csv(curve_display(), f, row.names = FALSE))

  sim <- eventReactive(input$run_sim, {
    dd <- d()
    per <- expected_screened(dd$n_per_arm, prev(), input$sens, input$spec, scoreable())
    n_sims <- max(50, min(input$n_sims, floor(SIM_BUDGET / per)))
    withProgress(message = "Simulating trials", value = 0.3, {
      r <- simulate_trials(dd$n_per_arm, prev(), p_hyper(), p_hypo(), input$sens, input$spec,
                           scoreable(), input$rrr, as.numeric(input$alpha), n_sims = n_sims)
    })
    c(r, list(capped = n_sims < input$n_sims, target = input$power, formula_ratio = dd$screen_per_enrolled,
              formula_ppv = dd$ppv))
  })

  output$sim_result <- renderUI({
    r <- sim()
    ok <- r$target >= r$power_ci[1] && r$target <= r$power_ci[2]
    tagList(
      tags$table(class = "table table-sm",
        tags$thead(tags$tr(tags$th(""), tags$th("Formula"), tags$th("Simulation"))),
        tags$tbody(
          tags$tr(tags$td("Power"), tags$td(pct(r$target)),
                  tags$td(sprintf("%s (95%% CI %s to %s)", pct(r$power), pct(r$power_ci[1]), pct(r$power_ci[2])))),
          tags$tr(tags$td("Screened per enrolled"), tags$td(formatC(r$formula_ratio, format = "f", digits = 3)),
                  tags$td(formatC(r$screen_per_enrolled, format = "f", digits = 3))),
          tags$tr(tags$td("Enrolled who are high-risk"), tags$td(pct(r$formula_ppv)), tags$td(pct(r$ppv))))),
      p(class = if (ok) "text-success" else "text-danger",
        if (ok) "Nominal power is inside the simulation's 95% interval."
        else "Nominal power is outside the simulation's 95% interval: check the inputs or run more trials."),
      p(class = "small text-muted", sprintf("%s trials of %s patients per arm.%s", scales::comma(r$n_sims),
        scales::comma(r$n_per_arm_simulated),
        if (r$capped) " Capped to keep the run short: this site needs many screens per enrolment." else ""))
    )
  })

  observeEvent(input$reset, {
    updateSelectInput(session, "alpha", selected = 0.05)
    updateSliderInput(session, "power", value = 0.80); updateSliderInput(session, "rrr", value = 0.20)
    updateNumericInput(session, "prev", value = round(OBS$prevalence, 4))
    updateNumericInput(session, "p_hyper", value = round(OBS$p_hyper, 4))
    updateNumericInput(session, "p_hypo", value = round(OBS$p_hypo, 4))
    updateSliderInput(session, "sens", value = DEF$constants$acc_assumed)
    updateSliderInput(session, "spec", value = DEF$constants$acc_assumed)
    updateSelectInput(session, "site", selected = "restricted_setB")
  })

  output$provenance <- renderTable({
    pr <- DEF$`_provenance`
    data.frame(Field = c("Source", "Source commit", "Exported", "Evidence tier", "Basis"),
               Value = c(paste(pr$source_repo, "09_enrichment_simulation.py"), pr$source_commit,
                         pr$created_utc, pr$evidence_tier, pr$basis))
  })
}

shinyApp(ui, server)
