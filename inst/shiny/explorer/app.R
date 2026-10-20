library(shiny)
library(LongOmicsTraj)

data("lot_glucold_ots", package = "LongOmicsTraj")
lot <- lot_glucold_ots
ots <- as.data.frame(lot_ots(lot))
stopifnot(all(c("object_id", "group", "assay", "estimator", "topology_label") %in% names(ots)))

ui <- fluidPage(
  titlePanel("LongOmicsTraj Explorer"),
  p("Explore Omics Trajectory Signatures. Topology matches are descriptive, not significance tests."),
  sidebarLayout(
    sidebarPanel(
      selectInput("assay", "Assay", sort(unique(as.character(ots$assay)))),
      selectInput("estimator", "Estimator", sort(unique(as.character(ots$estimator)))),
      uiOutput("group_controls"),
      selectInput("template", "Biological query template",
                  c("Custom", "Shared increase", "Opposing trajectories", "Stable reference group")),
      actionButton("run", "Run query", class = "btn-primary"),
      hr(),
      downloadButton("csv", "Export CSV"),
      downloadButton("code", "Export R query")
    ),
    mainPanel(
      h4(textOutput("count")),
      DT::DTOutput("results"),
      h4("Feature trajectory"),
      selectInput("feature", "Matching feature", choices = character(0)),
      plotOutput("trajectory", height = "340px")
    )
  )
)

server <- function(input, output, session) {
  scope <- reactive({
    req(input$assay, input$estimator)
    ots[ots$assay == input$assay & ots$estimator == input$estimator, , drop = FALSE]
  })
  arms <- reactive(sort(unique(as.character(scope()$group))))
  output$group_controls <- renderUI({
    labels <- sort(unique(as.character(scope()$topology_label)))
    tagList(lapply(seq_along(arms()), function(i)
      selectInput(paste0("arm_", i), arms()[i],
                  choices = c("Any" = "__any__", stats::setNames(labels, labels)))))
  })
  observeEvent(input$template, {
    g <- arms()
    if (length(g) < 2L || input$template == "Custom") return()
    chosen <- switch(input$template,
      "Shared increase" = rep("up_up", length(g)),
      "Opposing trajectories" = c("up_up", "down_down", rep("__any__", max(0, length(g)-2))),
      "Stable reference group" = c(rep("up_up", length(g)-1), "flat_flat"))
    for (i in seq_along(g)) updateSelectInput(session, paste0("arm_", i), selected = chosen[i])
  })
  conditions <- reactive({
    selected <- lapply(seq_along(arms()), function(i) input[[paste0("arm_", i)]])
    names(selected) <- arms()
    Filter(function(x) !is.null(x) && x != "__any__", selected)
  })
  code_text <- reactive({
    cs <- conditions()
    if (!length(cs)) return("# Select at least one group-specific OTS label")
    items <- vapply(names(cs), function(n)
      paste0("    ", encodeString(n, quote = '"'), " = ",
             encodeString(cs[[n]], quote = '"')), character(1))
    paste0("matches <- lot_query(\n  object = lot,\n  assay = ",
           encodeString(input$assay, quote = '"'), ",\n  estimator = ",
           encodeString(input$estimator, quote = '"'), ",\n  conditions = list(\n",
           paste(items, collapse = ",\n"), "\n  )\n)")
  })
  matches <- eventReactive(input$run, {
    cs <- conditions()
    validate(need(length(cs) > 0L, "Choose at least one OTS label."))
    x <- scope()
    validate(need(!anyDuplicated(x[c("object_id", "group")]),
                  "Duplicate feature/group combinations; select a different assay or estimator."))
    selected <- lot
    selected@ots <- S4Vectors::DataFrame(x)
    as.data.frame(lot_query(selected, conditions = cs))
  })
  output$count <- renderText({
    x <- matches()
    paste(length(unique(x$object_id)), "matching features (", nrow(x), "OTS rows)")
  })
  output$results <- DT::renderDT({
    DT::datatable(matches(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  observeEvent(matches(), {
    ids <- sort(unique(as.character(matches()$object_id)))
    updateSelectInput(session, "feature", choices = ids,
                      selected = if (length(ids)) ids[1] else character(0))
  })
  output$trajectory <- renderPlot({
    req(input$feature)
    tryCatch(print(lot_plot_trajectory(object = lot, features = input$feature,
                                       assay = input$assay)),
             error = function(e) {
               plot.new()
               text(.5, .5, paste("Plot unavailable:", conditionMessage(e)), cex = .8)
             })
  })
  output$csv <- downloadHandler(
    filename = function() "longomicstraj_matches.csv",
    content = function(file) utils::write.csv(matches(), file, row.names = FALSE)
  )
  output$code <- downloadHandler(
    filename = function() "longomicstraj_query.R",
    content = function(file) writeLines(code_text(), file)
  )
}
shinyApp(ui, server)
