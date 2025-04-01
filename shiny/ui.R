
css <- htmltools::HTML(
    "#summarytable > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody {
        transform:rotateX(180deg);
    }
    #summarytable > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody table{
        transform:rotateX(180deg);
    }"
)

ui = dashboardPage(
  header = dashboardHeader(title="2025-swpo-mls"),
  sidebar = dashboardSidebar(
    br(),
    br(),
    sidebarMenu(id="sidebarmenu",
      menuItem("Introduction", tabName="introduction"),
      menuItem("Summary table", tabName="table"),
      menuItem("Time series plots", tabName="ts_plots"),
      menuItem("Biology plots", tabName="bio_plots"),
      menuItem("Catch plots", tabName="catch_plots"),
      menuItem("Selectivity plots", tabName="selex_plots"),
      menuItem("Catch fit plots", tabName="catch_fit_plots"),
      menuItem("CPUE plots", tabName="cpue_plots")
    ),

    # Only show these on the plotting tabs - not Introduction and Summary table tabs
    conditionalPanel(condition="input.sidebarmenu == 'ts_plots'",
      # category
      awesomeCheckboxGroup(
      inputId = "category",
      label = "Metric", 
        choices = c("SSB", "F", "Recruitment", "Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"),
      selected = "SSB"),
      # scaled
      switchInput(
      inputId = "se",  
      label = "Show uncertainty",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger")
    ),
    conditionalPanel(condition="input.sidebarmenu == 'bio_plots'",
      # Growth variability
      switchInput(
      inputId = "growth_var",  
      label = "Show growth variability",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger"),
      # growth sex
      awesomeCheckboxGroup(
      inputId = "bio_sex",
      label = "Sex", 
      choices = c("Male", "Female"),
      selected = c("Male", "Female"))
    ),
    conditionalPanel(condition="input.sidebarmenu == 'catch_plots'",
      # catch chart
      awesomeRadio(
      inputId = "catch_chart",  
      label = "Plot type",
      choices=c("Bar","Area"),
      selected = "Bar"),
      # catch units
      awesomeRadio(
      inputId = "catch_units",  
      label = "Catch units",
      choices=c("Numbers","Biomass"),
      selected = "Numbers"),
      # catch by fleet
      switchInput(
      inputId = "catch_by_fleet",  
      label = "Aggregate catch",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger"),
      switchInput(
      inputId = "catch_rescale",  
      label = "Re-scale catch",
      value=FALSE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger"),
      # catch fill
      awesomeRadio(
        inputId = "catch_fill",
        label = "Fill type", 
          choices = c("None", "Fleet", "Units"),
        selected = "None"
      )
    ),
    conditionalPanel(condition="input.sidebarmenu == 'cpue_plots'",
      # show uncertainty 
      switchInput(
      inputId = "cpue_se",  
      label = "Show uncertainty",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger"),
      # show expected
      switchInput(
      inputId = "cpue_fit",  
      label = "Show fit",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger"),
      # log-scale
      switchInput(
      inputId = "cpue_log",  
      label = "Log-scale?",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger"),
      # var adj
      switchInput(
      inputId = "cpue_varadj",  
      label = "Var. adjust?",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger"),
      # lambdas
      switchInput(
      inputId = "cpue_lambda",  
      label = "Plot used?",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger")
    ),
    conditionalPanel(condition="input.sidebarmenu == 'catch_fit_plots'",
      # Show fitted lines
      switchInput(
        inputId = "catch_fit_show_lines",  
        label = "Show fitted lines",
        value=TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      # Free Y scale
      switchInput(
        inputId = "catch_fit_free_y",  
        label = "Free Y-axis scales",
        value=TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success", 
        offStatus = "danger"
      ),
      # Number of columns
      sliderInput(
        inputId = "catch_fit_n_col",
        label = "Number of columns",
        min = 1,
        max = 4,
        value = 1,
        step = 1
      ),
      # Fleet selector (will be populated in server)
      uiOutput("catch_fit_fleet_selector")
    ),
    conditionalPanel(condition="input.sidebarmenu == 'selex_plots'",
      # Number of columns
      sliderInput(
        inputId = "selex_n_col",
        label = "Number of columns",
        min = 1,
        max = 6,
        value = 3,
        step = 1
      )
    ),
    br(),
    br(),
    tags$footer(
      div(style="text-align:center",
        tags$p("version 0.0.1"),
        tags$p(paste("Copyright", format(Sys.time(),"%Y"), "NOAA Fisheries, PIFSC SAP"))
      )
    )
  ), # End of sidebar

  body = dashboardBody(
    tags$head(tags$style(HTML('.wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}') )),
    tags$head(tags$style(css)),
    # Start of main tab stuff
    tabItems(
      # **** Introduction ****
      tabItem(tabName="introduction", h2("Introduction"),
        fluidRow(column(12, includeMarkdown(file.path("introduction_index.md"))))
      ), # End of introduction tab

      # **** Summary table ****
      tabItem(tabName="table", h2("Summary table"),
        fluidRow(box(title="Model metrics", collapsed=FALSE, solidHeader=TRUE, collapsible=TRUE, status="primary", width=12,
         DT::dataTableOutput("summarytable")))
      ), # End of table tab

      # **** Time series plots ****
      tabItem(tabName="ts_plots", h2("Time series plots"),
        fluidRow(
          box(title="Time series plots", solidHeader=TRUE, collapsible=TRUE, collapsed=FALSE, status="primary", width=12,
            p("Select at least one model."),
            plotOutput("ts_plots", height="auto"))
        )
      ), # End of ts_plots tab

      # **** Biology plots ****
      tabItem(tabName="bio_plots", h2("Biology plots: Growth, M & L-W"),
        fluidRow(
          box(title="Growth", solidHeader=TRUE, collapsible=TRUE, collapsed=FALSE, status="primary", width=12,
            p("Select at least one model."),
            plotOutput("growth_plot", height="auto")),
          box(title="Natural mortality (M)", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, status="primary", width=12,
            p("Select at least one model."),
            plotOutput("m_plot", height="auto")),
          box(title="Length-weight", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, status="primary", width=12,
            p("Select at least one model."),
            plotOutput("lw_plot", height="auto"))
        )
      ), # End of bio_plots tab
      
      # **** Catch plots ****
      tabItem(tabName="catch_plots", h2("Catch plots"),
        fluidRow(
          box(title="Total catch", solidHeader=TRUE, collapsible=TRUE, collapsed=FALSE, status="primary", width=12,
            p("Select at least one model."),
            plotOutput("catch_total_plot", height="auto"))
          # box(title="Catch residuals", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, status="primary", width=12,
          #   p("Select at least one model."),
          #   plotOutput("catch_residuals_plot", height="auto")),
          # box(title="Catch: Observed & expected", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, status="primary", width=12,
          #   p("Select at least one model."),
          #   plotOutput("catch_obs_exp_plot", height="auto"))
        )
      ), # End of catch_plots tab

      tabItem(tabName="selex_plots", h2("Selectivity Plots"),
        fluidRow(
          box(title="Length-based Selectivity", solidHeader=TRUE, collapsible=TRUE, collapsed=FALSE, status="primary", width=12,
            p("Select at least one model to compare selectivity patterns."),
            plotOutput("selex_plot", height="auto")
          )
        )
      ),

      tabItem(tabName="catch_fit_plots", h2("Catch Fit Comparison"),
        fluidRow(
          box(title="Catch Observation vs. Fitted Values", solidHeader=TRUE, collapsible=TRUE, collapsed=FALSE, status="primary", width=12,
            p("Select at least one model to compare catch fits."),
            plotOutput("catch_fit_plot", height="auto")
          )
        )
      ),

      # **** CPUE plots ****
      tabItem(tabName="cpue_plots", h2("CPUE plots"),
        fluidRow(
          box(title="Observed CPUE", solidHeader=TRUE, collapsible=TRUE, collapsed=FALSE, status="primary", width=12,
            p("Select at least one model. Note that if multiple models are selected, any variance adjustment to the standard errors will be based off of the first model."),
            plotOutput("cpue_obs", height="auto"))
          # box(title="Catch residuals", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, status="primary", width=12,
          #   p("Select at least one model."),
          #   plotOutput("catch_residuals_plot", height="auto")),
          # box(title="Catch: Observed & expected", solidHeader=TRUE, collapsible=TRUE, collapsed=TRUE, status="primary", width=12,
          #   p("Select at least one model."),
          #   plotOutput("catch_obs_exp_plot", height="auto"))
        )
      ) # End of cpue_plots tab
    ) # End of tabItems
  ) # End of dashboardBody
)
