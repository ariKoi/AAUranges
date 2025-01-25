
# install.packages("pacman")
pacman::p_load(shiny, tidyverse, rms, gt, bs4Dash, shinyjs, bslib, shinybusy, miceRanger)

mw_dat <- data.frame(compound = c("Glycine","Alanine","Valine","Cysteine","Taurine","Isoleucine","Leucine","Glutamine","Glutamic acid","Methionine","Arginine"), mw = c(75.06,89.09,117.15,121.16,125.14,131.17,131.17,146.14,147.13,149.21,174.2))

final_models <- readRDS("./models/models.rds") 
ids_in <- readRDS("./models/ids_in.rds")

meta_data_train <- readxl::read_xlsx("./data/meta_data1.xlsx", col_types = c("text", "text", "numeric", "numeric", "text", "text", "numeric", "numeric"), na = "NA")

CRIQ <- 0.8 
CQ2 <- 0.6 
CQ3 <- 0.4 
extremeThres <- 0.01

predict_f <- function(mod, new_dat, conf_int = CRIQ, conf_type) {
  new_pred <- predict(mod, newdata = new_dat, conf.int = conf_int, conf.type = conf_type)
  low_uppper <- round(exp(c(new_pred$lower[[1]], new_pred$upper[[1]])), 2)
  names(low_uppper) <- dplyr::case_when(
    conf_int == CRIQ ~ paste0(c("lower","upper"), "CRIQ"),
    conf_int == CQ2 ~ paste0(c("lower","upper"), "CQ2"),
    conf_int == CQ3 ~ paste0(c("lower","upper"), "CQ3"),
    .default = paste0(c("lower","upper"), round(100*conf_int))
  )
  low_uppper
}

f_range <- function(.dat, model_input, confInt = 0.99999, confType = "mean") {

  amino_acid <- sub(pattern = "log", x = as.list(model_input$call$formula)[[2]], replacement = "")[2]
  
  tibble(`Amino Acid` = amino_acid) |> bind_cols(
    .dat |> tidyr::nest(.by = id) |> 
      mutate(preds = map(data, function(df) predict_f(model_input, new_dat = df, conf_int = confInt, conf_type = confType))) |> 
      unnest_wider(col = preds) |> 
      mutate(across(3:4, ~ round(., digits = 0))) |> 
      mutate(across(3:4, ~ format(.x, big.mark=",", trim=TRUE), .names = "{.col}_chr")) |> 
      tidyr::unite(col = "Range", 5, 6, sep = " - ", remove = TRUE) |> 
      relocate(Range, .after = data)
  ) |> 
    select(-c(2:3))
}

rangeGGplot_quantile_scale <- function(idx, res = NULL, rangesCRIQ, plot_vertical_mm = -10.5) {
  
  idx <- as.numeric(idx)
  cqs <- c(1, CRIQ, CQ2, CQ3) 
  
  idx_dat <- rangesCRIQ[idx,]
  
  p1 <- idx_dat |> 
    ggplot() + 
    geom_segment(aes(x = ifelse(lowerCRIQ != 0, 1-cqs[1], (1-cqs[2])/2), y = 0, xend = cqs[1], yend = 0), linewidth = 10, lineend = "round", col = "pink2", alpha = 0.7) + 
    guides(y = "none") +
    labs(x = NULL, y = NULL) +
    geom_segment(aes(x = (1-cqs[2])/2, y = 0, xend = 1 - (1-cqs[2])/2, yend = 0), linewidth = 10, lineend = "round", col = "yellow2", alpha = 0.7) +
    geom_segment(aes(x = (1-cqs[3])/2, y = 0, xend = 1 - (1-cqs[3])/2, yend = 0), linewidth = 10, lineend = "round", col = "green3", alpha = 0.4) +
    geom_segment(aes(x = (1-cqs[4])/2, y = 0, xend = 1 - (1-cqs[4])/2, yend = 0), linewidth = 10, lineend = "round", col = "green3", alpha = 0.7) + 
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks=element_blank()) + 
    theme(plot.margin = grid::unit(c(0,0,plot_vertical_mm,0), "mm"))

    if (res != "NA") {
      res <- as.numeric(res)
      if (idx_dat$lowerCRIQ == 0) { 
        idx_dat$lowerCRIQ <- idx_dat$lowerCRIQ + 1
      }
      mean_lnorm <- apply(log(idx_dat[,3:4]), 1, mean)
      sd_lnorm <- ( log(idx_dat$upperCRIQ) - mean_lnorm ) / qnorm(1 - (1-cqs[2])/2)
      res_quantile <- plnorm(res, mean = mean_lnorm, sd = sd_lnorm)
      p1 + geom_point(aes(x,y), data = data.frame(x=res_quantile, y=0), size = 9, shape = 16, col = scales::alpha("black", alpha = 1))
    } else {
      p1 
    }
}


ggplot_image2 <- function (plot_object, height = 100, aspect_ratio = 1) { 
  rlang::check_installed("ggplot2", "to use the `ggplot_image()` function.")
  if (is.numeric(height)) {
    height <- paste0(height, "px")
  }
  if (inherits(plot_object, "gg")) {
    plot_object <- list(plot_object)
  }
  vapply(seq_along(plot_object), FUN.VALUE = character(1), 
         USE.NAMES = FALSE, FUN = function(x) {
           filename <- paste0("temp_ggplot_", formatC(x, width = 4, flag = "0"), ".png")
           ggplot2::ggsave(filename = filename, plot = plot_object[[x]], 
                           device = "png", dpi = 100, width = 5 * aspect_ratio, 
                           height = 3, scale = 0.25) 
           on.exit(file.remove(filename))
           local_image(filename = filename, height = height)
         })
}

lnorm_valToQuantile <- function(val, lowerCRIQ, upperCRIQ) {
 
  if (lowerCRIQ == 0) {
    lowerCRIQ <- lowerCRIQ + 1
  }
  mean_lnorm <- mean(log(c(lowerCRIQ, upperCRIQ)))
  sd_lnorm <- ( log(upperCRIQ) - mean_lnorm ) / qnorm(1 - (1-CRIQ)/2)
  plnorm(val, mean = mean_lnorm, sd = sd_lnorm)
}

#---- ui ----

ui <- dashboardPage(title = "AminoAcidRanges", dark = NULL, fullscreen = TRUE, scrollToTop = TRUE, help = FALSE,  
    dashboardHeader(HTML('<img src = "theta_biomarkers_logo.png", width = "20%", height = "auto">'), title = "Amino Acid Ranges", border = FALSE), 
    
    dashboardSidebar(elevation = 5, collapsed = FALSE, width = "400px", skin = "light", minified = FALSE,
                     sidebarMenu(compact = TRUE,
                       menuItem(
                         strong("PPM Concentrations"),
                         uiOutput("in_concentrations"),
                         tabName = "in_concentrations",
                         icon = icon("sliders"),
                         startExpanded = TRUE
                       ),
                       menuItem(
                         strong("Meta Data"),
                         uiOutput("in_meta_data"),
                         tabName = "in_meta_data",
                         icon = icon("id-card"),
                         startExpanded = TRUE
                       ),
                       menuItem(
                         strong("ID Selection"),
                         uiOutput("id_selection"),
                         tabName = "id_selection",
                         icon = icon("watchman-monitoring"),
                         startExpanded = TRUE
                       ),
                       menuItem(
                         strong("Creatinine in Urine (mg/dL)"),
                         uiOutput("in_creatinine"),
                         tabName = "in_creatinine",
                         icon = icon("signal"),
                         startExpanded = TRUE
                       ),
                       menuItem(
                         strong("Report Language"),
                         uiOutput("in_language"),
                         tabName = "in_language",
                         icon = icon("edit", lib = "glyphicon"),
                         startExpanded = TRUE
                       ),
                       menuItem(
                         strong("Logos Used"),
                         uiOutput("in_logos"),
                         tabName = "in_logos",
                         icon = icon("record", lib = "glyphicon"),
                         startExpanded = TRUE
                       ),
                       actionButton(inputId = "run_report", label = strong("Generate Report"))
                     )
  ),
  dashboardBody(
    useShinyjs(), 
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom2.css")
    ),
    uiOutput("report"),
    downloadButton('export'),
    add_busy_spinner(spin = "hollow-dots", position = 'top-left', margins = c(325, 760), color = "#007bfffa")
  ),
  dashboardControlbar(overlay = FALSE, skin = "light", pinned = FALSE, width = 300, controlbarMenu(type = "pills", vertical = FALSE, controlbarItem(title = strong("Insert Notes"), icon = icon("comments"), textAreaInput(inputId = "writeup_notes", label = NULL, placeholder = "Specific Notes", height = "600px", value = "None"))
  )) 
)

#---- server ----

server <- function(input, output, session) {
  
  concentrations_in <- reactiveValues()
  meta_data_in <- reactiveValues()
  
  messageImpactNAs <- reactive({
    if (length(subject_info$colsNA) > 0) {
      " Report creation takes a bit longer (i.e. 30 to 60 seconds) when there are missing values in the subject's meta data input. The exact missing values are mentioned in the report. Please wait a bit for report generation to complete."
    } else {
      NULL
    }
  })
  
  output$in_concentrations <- renderUI({
      box(fileInput("concentrations", label = NULL, accept = c('text/csv','text/comma-separated-values','text/plain','.csv'), multiple = FALSE, buttonLabel = icon("database"), placeholder = "No file uploaded yet"), title = actionLink(inputId = "reset1App", label = strong("Upload Concentration Data"), icon = icon("redo")), collapsible = FALSE, id = "box_concentrations", label = boxLabel(text = "?", status = "primary", tooltip = "To upload a concentrations dataset please click on the `table` icon. When the app is already in use and a new dataset is required, please click on the `redo` icon to the left of the `Upload Concentration Data` text. This will clear any current calculation and reset app usage from the start."), width = 12, height = 100)
  })
  
  observeEvent(input$reset1App, {
    runjs("history.go(0)")
  })
  
  output$in_meta_data <- renderUI({
    box(fileInput("meta_data", label = NULL, accept = 'xlsx', multiple = FALSE, buttonLabel = icon("database"), placeholder = "No file uploaded yet"), title = actionLink(inputId = "reset2App", label = strong("Upload Subject ID Data"), icon = icon("redo")), collapsible = FALSE, id = "box_meta_data", label = boxLabel(text = "?", status = "primary", tooltip = "To upload a meta-data dataset please click on the `table` icon. When the app is already in use and a new dataset is required, please click on the `redo` icon to the left of the `Upload Subject Data` text. This will clear any current calculation and reset app usage from the start."), width = 12, height = 100)
  })
  
  observeEvent(input$reset2App, {
    runjs("history.go(0)")
  })
  
  observeEvent(input$concentrations, {
    
    if(!is.null(input$concentrations)) {
      ext <- tools::file_ext(input$concentrations$datapath)
      if(ext != "txt") {
        shinyjs::alert("Please upload a txt file for the concentrations.")
        shinyjs::reset("concentrations")
        shiny::validate(need(ext == "txt", message = FALSE))
      }
    }

    concentrations_in$df <- read.table(input$concentrations$datapath, header = FALSE, sep = "\t", fill = TRUE, skip = 3, comment.char = "") 
    
    if (length(colnames(concentrations_in$df)) != 13) {
      shinyjs::alert("The number of columns of the Concentrations Data input table is not 13. Please check the data table of Concentrations Data input file and reload.")
      shinyjs::reset("concentrations")
      shiny::validate(need(length(colnames(concentrations_in$df)) == 13, message = FALSE))
    }
    
    concentrations_in$df <- concentrations_in$df |> dplyr::select(V1, V3, V11)
    
    compounds <- trimws(substring(concentrations_in$df$V1[grep("Compound", concentrations_in$df$V1)], 14))
    
    concentrations_in$df <- as_tibble(concentrations_in$df) |> 
      filter(!(row_number() %in% grep("Compound", V1))) |> 
      filter(V1 != "") |> 
      mutate(compound = rep(compounds, each = n_distinct(V1))) |> 
      dplyr::select(-1) |> 
      rename(id = V3, ppm = V11) |> 
      mutate(ppm = as.numeric(ppm)) 
    
    compounds <- compounds[-grep(c("IS|C13|N15"), compounds)] 
    
    shiny::validate(need(length((compounds) == 11L), "The number of Amino Acids found is not 11. Please amend and reload the app."))
    
    concentrations_in$df <- concentrations_in$df |> filter(compound %in% compounds)
    concentrations_in$df[concentrations_in$df$compound == "leucine",]$compound <- "Leucine"
    concentrations_in$df <- concentrations_in$df |> 
      filter(!(row_number() %in% grep(paste(c("STD","QC"), collapse="|"), id))) |> 
      mutate(id = sub(".*_","",id))

    concentrations_in$df <- concentrations_in$df |> left_join(mw_dat, by = "compound") |> 
      mutate(micro_mol_per_liter = 1e3*(ppm*4/mw)) 
  }, once = TRUE) 

  output$id_selection <- renderUI({
    box(selectInput("sel_id", label = NULL, choices = NULL, multiple = FALSE, selected = NULL),
        width = 12, title = strong("Subject ID selection"), height = 100, id = "box_sel_id", collapsible = FALSE)
  })

  meta_data_in_req_colNames <- c("id","age","weight_kg","sex","height_m","chronic_disease", "chronic_disease_v3","meds","supplements","aerobics_times_per_week","aerobics_duration_min","resistance_times_per_week","resistance_duration_min")
  meta_data_in_req_genderNames <- c("Female","Male")
  meta_data_in_req_chronic_disease_v3Names <- c("Blood","Depression","Heart","None","Other","Pain","Stomach","Thyroidism")
 
  meta_data_in_age_range <- c(5,95) 
  meta_data_in_weight_kg_range <- c(30,200)
  meta_data_in_height_m_range <- c(1,2.5)
  meta_data_in_aerobics_times_per_week_range <- c(0,14) 
  meta_data_in_aerobics_duration_min_range <- c(0,240) 
  meta_data_in_resistance_times_per_week_range <- c(0,14) 
  meta_data_in_resistance_duration_min_range <- c(0,240)
  
  observeEvent(input$meta_data, {
    
    if(!is.null(input$meta_data)) {
      ext <- tools::file_ext(input$meta_data$datapath)
      if(ext != "xlsx") {
        shinyjs::alert("Please upload an xlsx file for the meta data.")
        shinyjs::reset("meta_data")
        shiny::validate(need(ext == "xlsx", message = FALSE))
      }
    }

    meta_data_in$df <- readxl::read_excel(input$meta_data$datapath)
    
    if (length(colnames(meta_data_in$df)) != 13) {
      shinyjs::alert("The number of columns of the Meta Data input table is not 13. Please check the data table of Meta Data input file and reload.")
      shinyjs::reset("meta_data")
      shiny::validate(need(length(colnames(meta_data_in$df)) == 13, message = FALSE))
    }
    
    if(!all(colnames(meta_data_in$df) %in% meta_data_in_req_colNames)) {
      shinyjs::alert(paste0("The Meta Data column names need to be as follows: ", paste0(meta_data_in_req_colNames, collapse = ", "), ". Please correct and reload."))
      shinyjs::reset("meta_data")
      shiny::validate(need(all(colnames(meta_data_in$df) %in% meta_data_in_req_colNames), message = FALSE))
    }
    
    if(!is.null(meta_data_in$df$sex) & !all(meta_data_in$df$sex %in% meta_data_in_req_genderNames, na.rm = TRUE)) {
      shinyjs::alert(paste0("Valid Meta Data `sex` column values are the following: ", paste0(meta_data_in_req_genderNames, collapse = ", "), ". Please correct and reload."))
    shinyjs::reset("meta_data")
    shiny::validate(need(all(meta_data_in$df$sex %in% meta_data_in_req_genderNames, na.rm = TRUE), message = FALSE))
    }
    
    if(!is.null(meta_data_in$df$chronic_disease_v3) & !all(meta_data_in$df$chronic_disease_v3 %in% meta_data_in_req_chronic_disease_v3Names, na.rm = TRUE)) {
      shinyjs::alert(paste0("Valid Meta Data `chronic_disease_v3` column values are the following: ", paste0(meta_data_in_req_chronic_disease_v3Names, collapse = ", "), ". Please correct and reload."))
    shinyjs::reset("meta_data")
    shiny::validate(need(all(meta_data_in$df$chronic_disease_v3 %in% meta_data_in_req_chronic_disease_v3Names, na.rm = TRUE), message = FALSE))
    }
    
    if(!is.null(meta_data_in$df$age) & !all(all(meta_data_in$df$age >= meta_data_in_age_range[1], na.rm = TRUE) & all(meta_data_in$df$age <= meta_data_in_age_range[2], na.rm = TRUE))) {
      shinyjs::alert(paste0("Valid Meta Data `age` column values are in the range: ", paste0(meta_data_in_age_range, collapse = "-"), ". Please correct and reload."))
    shinyjs::reset("meta_data")
    shiny::validate(need(all(all(meta_data_in$df$age >= meta_data_in_age_range[1], na.rm = TRUE) & all(meta_data_in$df$age <= meta_data_in_age_range[2], na.rm = TRUE)), message = FALSE))
    }
    
    if(!is.null(meta_data_in$df$weight_kg) & !all(all(meta_data_in$df$weight_kg >= meta_data_in_weight_kg_range[1], na.rm = TRUE) & all(meta_data_in$df$weight_kg <= meta_data_in_weight_kg_range[2], na.rm = TRUE))) {
      shinyjs::alert(paste0("Valid Meta Data `weight_kg` column values are in the range: ", paste0(meta_data_in_weight_kg_range, collapse = "-"), ". Please correct and reload."))
    shinyjs::reset("meta_data")
    shiny::validate(need(all(all(meta_data_in$df$weight_kg >= meta_data_in_weight_kg_range[1], na.rm = TRUE) & all(meta_data_in$df$weight_kg <= meta_data_in_weight_kg_range[2], na.rm = TRUE)), message = FALSE))
    }
    
    if(!is.null(meta_data_in$df$height_m) & !all(all(meta_data_in$df$height_m >= meta_data_in_height_m_range[1], na.rm = TRUE) & all(meta_data_in$df$height_m <= meta_data_in_height_m_range[2], na.rm = TRUE))) {
      shinyjs::alert(paste0("Valid Meta Data `height_m` column values are in the range: ", paste0(meta_data_in_height_m_range, collapse = "-"), ". Please correct and reload."))
    shinyjs::reset("meta_data")
    shiny::validate(need(all(all(meta_data_in$df$height_m >= meta_data_in_height_m_range[1], na.rm = TRUE) & all(meta_data_in$df$height_m <= meta_data_in_height_m_range[2], na.rm = TRUE)), message = FALSE))
    }
    
    if(!is.null(meta_data_in$df$aerobics_times_per_week) & !all(all(meta_data_in$df$aerobics_times_per_week >= meta_data_in_aerobics_times_per_week_range[1], na.rm = TRUE) & all(meta_data_in$df$aerobics_times_per_week <= meta_data_in_aerobics_times_per_week_range[2], na.rm = TRUE))) {
      shinyjs::alert(paste0("Valid Meta Data `aerobics_times_per_week` column values are in the range: ", paste0(meta_data_in_aerobics_times_per_week_range, collapse = "-"), ". Please correct and reload."))
    shinyjs::reset("meta_data")
    shiny::validate(need(all(all(meta_data_in$df$aerobics_times_per_week >= meta_data_in_aerobics_times_per_week_range[1], na.rm = TRUE) & all(meta_data_in$df$aerobics_times_per_week <= meta_data_in_aerobics_times_per_week_range[2], na.rm = TRUE)), message = FALSE))
    }
    
    if(!is.null(meta_data_in$df$aerobics_duration_min) & !all(all(meta_data_in$df$aerobics_duration_min >= meta_data_in_aerobics_duration_min_range[1], na.rm = TRUE) & all(meta_data_in$df$aerobics_duration_min <= meta_data_in_aerobics_duration_min_range[2], na.rm = TRUE))) {
      shinyjs::alert(paste0("Valid Meta Data `aerobics_duration_min` column values are in the range: ", paste0(meta_data_in_aerobics_duration_min_range, collapse = "-"), ". Please correct and reload."))
    shinyjs::reset("meta_data")
    shiny::validate(need(all(all(meta_data_in$df$aerobics_duration_min >= meta_data_in_aerobics_duration_min_range[1], na.rm = TRUE) & all(meta_data_in$df$aerobics_duration_min <= meta_data_in_aerobics_duration_min_range[2]), na.rm = TRUE), message = FALSE))
    }
    
    if(!is.null(meta_data_in$df$resistance_times_per_week) & !all(all(meta_data_in$df$resistance_times_per_week >= meta_data_in_resistance_times_per_week_range[1], na.rm = TRUE) & all(meta_data_in$df$resistance_times_per_week <= meta_data_in_resistance_times_per_week_range[2], na.rm = TRUE))) {
      shinyjs::alert(paste0("Valid Meta Data `resistance_times_per_week` column values are in the range: ", paste0(meta_data_in_resistance_times_per_week_range, collapse = "-"), ". Please correct and reload."))
    shinyjs::reset("meta_data")
    shiny::validate(need(all(all(meta_data_in$df$resistance_times_per_week >= meta_data_in_resistance_times_per_week_range[1], na.rm = TRUE) & all(meta_data_in$df$resistance_times_per_week <= meta_data_in_resistance_times_per_week_range[2]), na.rm = TRUE), message = FALSE))
    }
    
    if(!is.null(meta_data_in$df$resistance_duration_min) & !all(all(meta_data_in$df$resistance_duration_min >= meta_data_in_resistance_duration_min_range[1], na.rm = TRUE) & all(meta_data_in$df$resistance_duration_min <= meta_data_in_resistance_duration_min_range[2], na.rm = TRUE))) {
      shinyjs::alert(paste0("Valid Meta Data `resistance_duration_min` column values are in the range: ", paste0(meta_data_in_resistance_duration_min_range, collapse = "-"), ". Please correct and reload."))
    shinyjs::reset("meta_data")
    shiny::validate(need(all(all(meta_data_in$df$resistance_duration_min >= meta_data_in_resistance_duration_min_range[1], na.rm = TRUE) & all(meta_data_in$df$resistance_duration_min <= meta_data_in_resistance_duration_min_range[2]), na.rm = TRUE), message = FALSE))
    }
    
    meta_data_in$df <- meta_data_in$df |> 
      mutate(bmi = weight_kg/height_m^2, 
             supplements_flag = ifelse(is.na(supplements), "No", "Yes")) |> 
      tidyr::replace_na(list(aerobics_times_per_week = 0, aerobics_duration_min = 0, resistance_times_per_week = 0, resistance_duration_min = 0)) |> 
      mutate(aerobics_activity = aerobics_times_per_week * aerobics_duration_min, 
             resistance_activity = resistance_times_per_week * resistance_duration_min) |> 
      select(id, sex, age, bmi, chronic_disease_v3, supplements_flag, aerobics_activity, resistance_activity)
    
    updateSelectInput(session, inputId = "sel_id", choices = meta_data_in$df$id) 
  }, once = FALSE)
  
  output$in_creatinine <- renderUI({
    box(width = 12, height = 100, title = strong("Enter value as per subject ID:"), collapsible = FALSE, id = "box_creatinine",
        numericInput("creatinine_val",
                     label = NULL,
                     value = NULL,
                     min = 0,
                     max = 1000)
    )
  })
  
  output$in_language <- renderUI({
    box(width = 12, height = 100, title = strong("Enter the report language:"), collapsible = FALSE, id = "box_language",
        selectInput("language_choice",
                     label = NULL,
                     choices = c("Greek","English"),
                     selected = "English")
    )
  })
  
  output$in_logos <- renderUI({
    box(width = 12, height = 100, title = strong("Enter the number of logos:"), collapsible = FALSE, id = "box_logos",
        selectInput("logo_choice",
                     label = NULL,
                     choices = 1:2,
                     selected = 1)
    )
  })
  
  subject_info <- reactiveValues(colsNA = NULL, id = NULL)
  reportOut <- reactiveValues()
  profileNotes <- reactiveValues(notes = NULL)
  
  observeEvent(input$sel_id, {
    shinyjs::reset(id = "writeup_notes")
    shinyjs::reset(id = "creatinine_val")
    shinyjs::reset(id = "language_choice")
    shinyjs::reset(id = "run_report")
    })
  
  observeEvent(input$creatinine_val, {
    shinyjs::reset(id = "run_report")
  })
  observeEvent(input$run_report, {
    if (is.na(input$creatinine_val)) shinyjs::alert("Please enter the creatinine value.")
  })
  
  observeEvent(input$language_choice, {
    if (input$language_choice == "Greek") {
      updateTextAreaInput(inputId = "writeup_notes", placeholder = "Ειδικές Σημειώσεις", value = "Χωρίς")
    } else {
      updateTextAreaInput(inputId = "writeup_notes", placeholder = "Specific Notes", value = "None")
    }
  })
  
    output$report <- renderUI({
      
      creatinineVal <- isolate(input$creatinine_val) 
      selID <- isolate(input$sel_id) 
      langChoice <- isolate(input$language_choice)
      logoChoice <- isolate(input$logo_choice)
      
      req(input$run_report, creatinineVal, langChoice, logoChoice)

      if (creatinineVal < 10 | creatinineVal > 400) {
        shinyjs::alert("The creatinine value is expected to be between 10 mg/dL and 400 mg/dL. Please make sure it is a valid value.")
      }
      
      if(length(unique(concentrations_in$df$id)) != length(meta_data_in$df$id)) {
        shinyjs::alert("Subject ids are not the same number in the PPM Concentration Data and the Meta Data. Please check the id consistency between the two input data sets and re-run the app.")
        runjs("history.go(0)")
        shiny::validate(need(length(unique(concentrations_in$df$id)) == length(meta_data_in$df$id), message = FALSE))
      }
     
      id_check <- all(meta_data_in$df$id %in% unique(concentrations_in$df$id))
      in_meta_out_conc <- meta_data_in$df$id[which(meta_data_in$df$id %in% unique(concentrations_in$df$id) == FALSE)]
      in_conc_out_meta <- unique(concentrations_in$df$id)[which(unique(concentrations_in$df$id) %in% meta_data_in$df$id == FALSE)]
      if (!id_check) {
        shinyjs::alert(paste0("Subject ids are not the same in the PPM Concentration Data and the Meta Data. Please check the id consistency between the two input data sets and re-run the app. The non matching ids are in Meta Data: ", in_meta_out_conc, " and in Concentration Data: ", in_conc_out_meta))
        runjs("history.go(0)")
        shiny::validate(need(id_check, message = FALSE)) 
      }
     
    concentrations_in_subj_id <- concentrations_in$df |> 
      filter(id == selID) |> 
      mutate(micro_mol_per_gram = micro_mol_per_liter/(creatinineVal/100)) 
    
    concentrations_in_subj_id$compound[9] <- "Glutamate"
    
    concentrations_in_subj_id <- tibble(id = selID, ppm = NA, compound = "Glutamine/Glutamate", mw = NA, micro_mol_per_liter = NA, micro_mol_per_gram = concentrations_in_subj_id$micro_mol_per_gram[8]/concentrations_in_subj_id$micro_mol_per_gram[9]) |> 
      bind_rows(concentrations_in_subj_id) 
    
    subject_dat <- meta_data_in$df |> filter(id == selID) 
    
    subject_info$colsNA <- subject_dat |> select_if(is.na) |> colnames()
    
    subject_info$id <- subject_dat[1]

    ref_ranges <- map(final_models, ~ f_range(subject_dat, model_input = ., confInt = CRIQ, confType = "individual")) |> list_rbind()
    
    ref_ranges_CQ2 <- map(final_models, ~ f_range(subject_dat, model_input = ., confInt = CQ2, confType = "individual")) |> list_rbind() 
    
    ref_ranges_CQ3 <- map(final_models, ~ f_range(subject_dat, model_input = ., confInt = CQ3, confType = "individual")) |> list_rbind()
    
    compounds_idx_NA <- which(is.na(ref_ranges$lowerCRIQ))
    
    if (length(subject_info$colsNA) > 0) {
      
      showModal(modalDialog(messageImpactNAs()))
      
      result_agg_imp <- map(compounds_idx_NA, function(idx){
     
        vars_in_model <- names(final_models[[idx]]$na.action$nmiss[-1])
      
        data_in_model <- meta_data_train |> 
          filter(id %in% ids_in[[idx]]) |> 
          select(all_of(vars_in_model)) |> 
          filter(!(row_number() %in% final_models[[idx]]$na.action$omit)) # taken out from refit
        
        if (ncol(data_in_model) == 2) {
          data_in_model <- data_in_model |> mutate(rfCol = 1L)
          s_info <- subject_dat |> select(all_of(vars_in_model)) |> mutate(rfCol = 1L)
          } else {
            s_info <- subject_dat |> select(all_of(vars_in_model))
        }
        
        mR <- miceRanger::miceRanger(data_in_model |> bind_rows(s_info), m = 15, maxiter = 5, num.trees = 45, verbose = FALSE) 
        
        imp_dataList <- miceRanger::completeData(mR)
        
        if (ncol(data_in_model) == 2) {
        imp_dataList <- map(imp_dataList, ~ subject_dat[1] |> bind_cols(as_tibble(tail(., n=1)) |> select(-rfCol)))
        } else {
          imp_dataList <- map(imp_dataList, ~ subject_dat[1] |> bind_cols(as_tibble(tail(., n=1))))
        }
        
        ref_ranges_imp_single_compound <- map(imp_dataList, ~ f_range(., model_input = final_models[[idx]], confInt = CRIQ, confType = "individual")) |> list_rbind() 
        
        ref_ranges_imp_single_compound_final <- ref_ranges_imp_single_compound |> 
          group_by(`Amino Acid`) |> 
          summarise(lowerCRIQ = round(mean(lowerCRIQ)), upperCRIQ = round(mean(upperCRIQ))) |> 
          mutate(Range = paste0(format(lowerCRIQ, big.mark=",", trim=TRUE), " - ", format(upperCRIQ, big.mark=",", trim=TRUE))) |> 
          relocate(Range, .after = `Amino Acid`)
        
        ref_ranges_CQ2_imp_single_compound_final <- map(imp_dataList, ~ f_range(., model_input = final_models[[idx]], confInt = CQ2, confType = "individual")) |> 
          list_rbind() 
        
        ref_ranges_CQ2_imp_single_compound_final <- ref_ranges_CQ2_imp_single_compound_final |> 
          group_by(`Amino Acid`) |> 
          summarise(lowerCQ2 = round(mean(lowerCQ2)), upperCQ2 = round(mean(upperCQ2))) |> 
          mutate(Range = paste0(format(lowerCQ2, big.mark=",", trim=TRUE), " - ", format(upperCQ2, big.mark=",", trim=TRUE))) |> 
          relocate(Range, .after = `Amino Acid`)
        
        ref_ranges_CQ3_imp_single_compound_final <- map(imp_dataList, ~ f_range(., model_input = final_models[[idx]], confInt = CQ3, confType = "individual")) |> 
          list_rbind()
        
        ref_ranges_CQ3_imp_single_compound_final <- ref_ranges_CQ3_imp_single_compound_final |> 
          group_by(`Amino Acid`) |> 
          summarise(lowerCQ3 = round(mean(lowerCQ3)), upperCQ3 = round(mean(upperCQ3))) |> 
          mutate(Range = paste0(format(lowerCQ3, big.mark=",", trim=TRUE), " - ", format(upperCQ3, big.mark=",", trim=TRUE))) |> 
          relocate(Range, .after = `Amino Acid`)

        ref_ranges_imp_single_compound_final |>
          bind_rows(ref_ranges_CQ2_imp_single_compound_final) |>
          bind_rows(ref_ranges_CQ3_imp_single_compound_final)
        
      }) |> list_rbind()
    
    ref_ranges <- ref_ranges |> rowid_to_column() |> na.omit() |> bind_rows(
      result_agg_imp |> select(1:4) |> na.omit() |> mutate(rowid = compounds_idx_NA) 
    ) |> arrange(rowid) |> select(-rowid)
    
    ref_ranges_CQ2 <- ref_ranges_CQ2 |> rowid_to_column() |> na.omit() |> bind_rows(
      result_agg_imp |> select(1:2,5:6) |> na.omit() |> mutate(rowid = compounds_idx_NA) 
    ) |> arrange(rowid) |> select(-rowid)
    
    ref_ranges_CQ3 <- ref_ranges_CQ3 |> rowid_to_column() |> na.omit() |> bind_rows(
      result_agg_imp |> select(1:2,7:8) |> na.omit() |> mutate(rowid = compounds_idx_NA) 
    ) |> arrange(rowid) |> select(-rowid)
    
    removeModal()
    } 
    
    if (langChoice == "English") {
      
    age_comment <- ifelse("age" %in% subject_info$colsNA, "(missing Age)", paste0(subject_dat$age, " year old"))
    sex_comment <- ifelse("sex" %in% subject_info$colsNA, "(missing Sex)", tolower(subject_dat$sex))
    bmi_comment <- ifelse("bmi" %in% subject_info$colsNA, "(missing BMI)", paste0("with a Body Mass Index of ", round(subject_dat$bmi,2)))
    chronic_disease_comment <- ifelse(subject_dat$chronic_disease_v3 %in% c("None","Other"), paste0(ifelse(subject_dat$chronic_disease_v3 == "None", "no", "a"), " chronic condition,"), paste0(subject_dat$chronic_disease_v3, " type chronic condition,"))
    
    footnote_subj_meta_title <- gt::md("**Reference Range**")
    
    footnote_subj_meta <- paste("Dynamic - model based on following subject characteristics:" , age_comment, sex_comment, bmi_comment, "who has", chronic_disease_comment, "who takes", ifelse(subject_dat$supplements_flag == "Yes", "", "no"), "supplements, and who on a weekly basis does", ifelse(subject_dat$aerobics_activity == 0, "no", paste0(subject_dat$aerobics_activity, " minutes of")), "aerobics and", ifelse(subject_dat$resistance_activity == 0, "no", paste0(subject_dat$resistance_activity, " minutes of")), "resistance training.")
    
    } else {
      
      age_comment <- ifelse("age" %in% subject_info$colsNA, "(άγνωστη ηλικία)", paste0(subject_dat$age, " χρονών"))
      sex_comment <- ifelse("sex" %in% subject_info$colsNA, "(άγνωστο φύλλο)", ifelse(tolower(subject_dat$sex) == "male", "άνδρας", "γυναίκα"))
      bmi_comment <- ifelse("bmi" %in% subject_info$colsNA, "(άνγωστο ΔΜΣ)", paste0("με δείκτη μάζας σώματος (ΔΜΣ) ", round(subject_dat$bmi,2), ","))
      chronic_disease_comment <- ifelse(subject_dat$chronic_disease_v3 %in% c("None","Other"), paste0(ifelse(subject_dat$chronic_disease_v3 == "None", "χωρίς κάποια", "με κάποια (αόριστη)"), " χρόνια πάθηση,"), paste0(" που πάσχει απο χρόνια πάθηση σχετιζόμενη με ", dplyr::case_when(
        subject_dat$chronic_disease_v3 == "Thyroidism" ~ "τον Θυρεοειδισμό",
        subject_dat$chronic_disease_v3 == "Blood" ~ "το αίμα",
        subject_dat$chronic_disease_v3 == "Heart" ~ "την καρδιά",
        subject_dat$chronic_disease_v3 == "Stomach" ~ "το στομάχι",
        subject_dat$chronic_disease_v3 == "Depression" ~ "την κατάθλιψη",
        subject_dat$chronic_disease_v3 == "Pain" ~ "τον πόνο"
      ), ","))
      
      footnote_subj_meta_title <- gt::md("**Εύρος Αναφοράς**")
        
      footnote_subj_meta <- paste("Δυναμικό μοντέλο βασιζόμενο στα παρακάτω χαρακτηριστικά:" , age_comment, sex_comment, bmi_comment, chronic_disease_comment, "που", ifelse(subject_dat$supplements_flag == "Yes", "", "δεν"), "παίρνει διατροφικά συμπληρώματα, και που σε εβδομαδιαία βάση", ifelse(subject_dat$aerobics_activity == 0, "δεν κάνει αερόβια άσκηση", paste0("κάνει ", subject_dat$aerobics_activity, " λεπτά αερόβιας άσκησης")), "και", ifelse(subject_dat$resistance_activity == 0, "δεν κάνει αναερόβια άσκηση.", paste0("κάνει ", subject_dat$resistance_activity, " λεπτά αναερόβιας άσκησης.")))
    }
 
    concentrations_in_subj_id_res <- concentrations_in_subj_id |> 
      select(compound, micro_mol_per_gram) |> 
      rename(`Amino Acid` = compound, Result = micro_mol_per_gram)
    
    concentrations_in_subj_id_creatinine <- tibble(`Amino Acid` = "Creatinine", Result = creatinineVal, Range = "45 - 225", lowerCRIQ = 45, upperCRIQ = 225) 
    
    ref_ranges$`Amino Acid`[1] <- "Glutamine/Glutamate"
    ref_ranges_CQ2$`Amino Acid`[1] <- "Glutamine/Glutamate"
    ref_ranges_CQ3$`Amino Acid`[1] <- "Glutamine/Glutamate"
    
    concentrations_in_subj_id_res <- concentrations_in_subj_id_creatinine |>
      bind_rows(
        concentrations_in_subj_id_res |> inner_join(ref_ranges, by = c("Amino Acid"))
        )
    
    notes_comments1 <- if (langChoice == "English") {md("**Profile Notes**")} else {md("**Ειδικές Σημειώσεις**")}
    
    notes_comments2 <- isolate(input$writeup_notes)
   
    notes_comments3 <- if (langChoice == "English") {md("**General Information**")} else {md("**Γενικές Πληροφορίες**")}
      
    notes_comments4 <- if (langChoice == "English") {"Amino acids are the building blocks of proteins, which serve to build our tissues and regulate our metabolism (as enzymes and hormones). Since amino acids are involved in all biological processes, deficiencies or imbalances in them can lead to disruptions in the normal functioning of the body. Urine amino acid evaluation can help identify such deficiencies or imbalances and allow for their restoration. Because many vitamins and elements serve as cofactors in amino acid metabolism, problems with amino acids can also indicate problems with these cofactors."} else {
      "Τα αμινοξέα είναι τα δομικά στοιχεία των πρωτεϊνών, που χρησιμεύουν για την κατασκευή των ιστών μας και τη ρύθμιση του μεταβολισμού μας (ως ένζυμα και ορμόνες). Δεδομένου ότι τα αμινοξέα εμπλέκονται σε όλες τις βιολογικές διεργασίες, οι ελλείψεις ή οι ανισορροπίες σε αυτά μπορεί να οδηγήσουν σε διαταραχές στη φυσιολογική λειτουργία του σώματος. Η αξιολόγηση αμινοξέων στα ούρα μπορεί να βοηθήσει στον εντοπισμό τέτοιων ελλείψεων ή ανισορροπιών και να επιτρέψει την αποκατάστασή τους. Επειδή πολλές βιταμίνες και στοιχεία χρησιμεύουν ως συμπαράγοντες στο μεταβολισμό των αμινοξέων, τα προβλήματα με τα αμινοξέα μπορεί επίσης να υποδηλώνουν προβλήματα με αυτούς τους συμπαράγοντες."
    }
    
    notes_comments5 <- if (langChoice == "English") {"The restoration of amino acid problems is primarily achieved through diet. Foods rich in proteins include meat of any origin (red meat, poultry, fish), processed meats, liver, eggs, milk, dairy products, legumes, nuts, and cereals. The proteins in these foods vary in their distribution of amino acids. Therefore, identifying any deficiencies or imbalances with this analysis can guide us in the intake of specific foods or supplements to correct them."} else {
      "Η αποκατάσταση των προβλημάτων αμινοξέων επιτυγχάνεται κατά κύριο λόγο μέσω της διατροφής. Τροφές πλούσιες σε πρωτεΐνες περιλαμβάνουν κρέας οποιασδήποτε προέλευσης (κόκκινο κρέας, πουλερικά, ψάρι), επεξεργασμένα κρέατα, συκώτι, αυγά, γάλα, γαλακτοκομικά προϊόντα, όσπρια, ξηρούς καρπούς και δημητριακά. Οι πρωτεΐνες σε αυτά τα τρόφιμα ποικίλλουν ως προς την κατανομή των αμινοξέων τους. Επομένως, ο εντοπισμός τυχόν ελλείψεων ή ανισορροπιών με αυτήν την ανάλυση μπορεί να μας καθοδηγήσει στην πρόσληψη συγκεκριμένων τροφών ή συμπληρωμάτων για τη διόρθωσή τους."
    }
    
    ref_ranges <- concentrations_in_subj_id_creatinine |> select(-Result) |> bind_rows(ref_ranges)

    mean_lnorm_creatinine <- apply(log(concentrations_in_subj_id_creatinine[,4:5]), 1, mean)
    sd_lnorm_creatinine <- ( log(concentrations_in_subj_id_creatinine$upperCRIQ) - mean_lnorm_creatinine ) / qnorm(1 - (1-CRIQ)/2)
    
    ref_ranges_CQ2 <- concentrations_in_subj_id_creatinine |> 
      select(-Result, -Range,-lowerCRIQ, -upperCRIQ) |> 
      mutate(lowerCQ2 = round(qlnorm((1-CQ2)/2, mean = mean_lnorm_creatinine, sd = sd_lnorm_creatinine)), 
             upperCQ2 = round(qlnorm(1 - (1-CQ2)/2, mean = mean_lnorm_creatinine, sd = sd_lnorm_creatinine))) |> 
      mutate(Range = paste0(format(lowerCQ2, big.mark=",", trim=TRUE), " - ", format(upperCQ2, big.mark=",", trim=TRUE))) |> 
      relocate(Range, .after = `Amino Acid`) |>
      bind_rows(ref_ranges_CQ2)

    ref_ranges_CQ3 <- concentrations_in_subj_id_creatinine |> 
      select(-Result, -Range, -lowerCRIQ, -upperCRIQ) |> 
      mutate(lowerCQ3 = round(qlnorm((1-CQ3)/2, mean = mean_lnorm_creatinine, sd = sd_lnorm_creatinine)), 
             upperCQ3 = round(qlnorm(1 - (1-CQ3)/2, mean = mean_lnorm_creatinine, sd = sd_lnorm_creatinine))) |> 
      mutate(Range = paste0(format(lowerCQ3, big.mark=",", trim=TRUE), " - ", format(upperCQ3, big.mark=",", trim=TRUE))) |> 
      relocate(Range, .after = `Amino Acid`) |>
      bind_rows(ref_ranges_CQ3)

    rep <- concentrations_in_subj_id_res |> 
      mutate(idx = as.character(1:13)) |> 
      group_by(idx) |> 
      mutate(Distribution = list(c(idx, Result))) |> 
      ungroup() |> 
      select(1:3,7) |> 
      rename(Indices = `Amino Acid`) 

    rep_digits <- rep |> 
      mutate(chck1 = round(Result - floor(Result), digits = 2)) |> 
      mutate(chck2 = round(Result - floor(Result), digits = 1)) |> 
      mutate(digits = case_when( 
        Result < 1 & chck1 != 0 ~ 2,
        Result < 1 & chck1 == 0 ~ 0,
        (Result >= 1 & Result <= 1e3) & chck2 != 0 ~ 1,
        (Result >= 1 & Result <= 1e3) & chck2 == 0 ~ 0,
        Result > 1e3 ~ 0,
        is.na(Result) ~ 1 
      )) |> pull(digits)
    
    if (langChoice == "Greek") {
      
      rep[[1]] <- c("Κρεατινίνη","Γλουταμίνη/Γλουταμικό Οξύ","Γλυκίνη","Αλανίνη","Βαλίνη","Κυστεΐνη","Ταυρίνη","Ισολευκίνη","Λευκίνη","Γλουταμίνη","Γλουταμικό Οξύ","Μεθιονίνη","Αργινίνη")
      
      date_title <- "Ημερομηνία Παραλαβής: "
      result_info_footnote <- "LOQ: Όριο κατώτατης ποσοτικής αποτίμησης"
      notes_comments6 <- "Η παρούσα αναφορά δεν προορίζεται για τη διάγνωση, τη θεραπεία ή την πρόληψη παθολογικών καταστάσεων, ούτε προορίζεται να αντικαταστήσει μια ιατρική διάγνωση. Μέθοδος ανάλυσης: UPLC-MS/MS."
      notes_comments7 <- "Για περισσότερες πληροφορίες παρακαλώ επισκεφθείτε τον ιστότοπο"
      notes_comments8 <- "Προδιαγραφή χρώματος: Κοινό (πράσινο), Κάπως Κοινό (πορτοκαλί), Ασυνήθιστο (ροζ)"
      notes_comments9 <- "Αμινοξέα Ούρων" 
      notes_comments10 <- "**Αμινοξέα**"
      notes_comments11 <- "**Ειδικοί**"
      notes_comments12 <- "*μmol/g* *κρεατινίνης*" 
      cols_label_fix1_1 <- "**Δείκτης Δείγματος**"
      cols_label_fix1_2 <- "**Βιοδείκτης**"
      cols_label_fix2 <- "**Αποτέλεσμα**"
      cols_label_fix3 <- "**Εύρος Αναφοράς**"
    } else {
      date_title <- "Delivery Date: "
      result_info_footnote <- "LOQ: Limit of quantification"
      notes_comments6 <- "This report is not intended for the diagnosis, treatment, or prevention of pathological conditions, nor is it meant to replace a medical diagnosis. Method of analysis: UPLC-MS/MS."
      notes_comments7 <- "For more information please visit"
      notes_comments8 <- "Color specification: Common (green), Somewhat Common (orange), Uncommon (pink)"
      notes_comments9 <- "Urine Amino Acids"
      notes_comments10 <- "**Amino Acids**"
      notes_comments11 <- "**Specific**"
      notes_comments12 <- "*μmol/g* *creatinine*"
      cols_label_fix1_1 <- "**Specimen Marker**"
      cols_label_fix1_2 <- "**Biomarker**"
      cols_label_fix2 <- "**Result**"
      cols_label_fix3 <- "**Reference Range**"
    }

#---- Table 1 ---- 
    if (logoChoice == 1) {

    gtbl1 <- head(rep, 2) |> 
      rename(`Specimen Marker` = Indices) |> 
      gt(id = "specimenMarkerTable") |> 
      tab_header( # ?tab_header
        title = htmltools::tagList(
          htmltools::tags$div(
            style = htmltools::css(
              `display` = "flex", 
              `align-items` = "center" 
            ),
            htmltools::HTML("&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;"), 
            htmltools::tags$div(
              htmltools::tags$strong(notes_comments9)
            )
          )
        )
      ) |> 
      tab_spanner(
        label = htmltools::HTML(
          local_image("./www/theta_biomarkers_logo.png", height = 86.5) 
        ), columns = 1:2, level = 1, id = "logo"
      ) |> 
      tab_spanner(
        label = htmltools::HTML("<p style = 'margin-left: 23px'; align = 'left'>", "<br>", paste0("ID: ", selID), "<br><br>", paste0(date_title, format(Sys.Date(),format="%d %B %Y")), "<br><br>", "</p>"), columns = 3:4, level = 1, id = "id_date"
      ) |>
      opt_css(css = "#specimenMarkerTable .gt_column_spanner {border-bottom-style: hidden;}") |>
      cols_align('center', columns = everything()) |> 
      fmt_number(columns = Result, rows = 1, decimals = rep_digits[1], sep_mark = ",", pattern = paste0("{x}", " mg/dL")) |> 
      fmt_number(columns = Result, rows = 2, decimals = rep_digits[2], sep_mark = ",") |> 
      sub_missing(
        columns = Result,
        missing_text = "LOQ"
      ) |> 
      text_transform(
        locations = cells_body(columns = 'Distribution'),
        fn = function(column) {
          map(column, ~ stringr::str_split_1(., ', ')) |>
            map(~ rangeGGplot_quantile_scale(.[1], .[2], rangesCRIQ = ref_ranges, plot_vertical_mm = -4.5)) |> 
            ggplot_image2(height = 30, aspect_ratio = 3) 
        }
      ) |> 
      cols_label(`Specimen Marker` = cols_label_fix1_1, 
                 Result = cols_label_fix2,
                 Range = cols_label_fix3, 
                 Distribution = "", .fn = md) |> 
      cols_width(`Specimen Marker` ~ px(150)) |> 
      cols_width(c(Result,Range) ~ px(120)) |> 
      opt_table_font(font = google_font("Roboto Slab")) |> 
      tab_options(table.font.size = 12, table.width = pct(80))  
    } else {

      gtbl1 <- head(rep, 2) |> 
        rename(`Specimen Marker` = Indices) |> 
        gt(id = "specimenMarkerTable") |> 
        tab_header( 
          title = htmltools::tagList(
            htmltools::tags$div(
              style = htmltools::css(
                `display` = "flex", 
                `align-items` = "center" 
              ),
              htmltools::HTML("&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;"), 
              htmltools::tags$div(
                htmltools::tags$strong(notes_comments9)
              )
            )
          )
        ) |> 
        tab_spanner(
          label = htmltools::HTML(
            local_image("./www/theta_biomarkers_logo.png", height = 57) 
          ), columns = 1, level = 1, id = "logo1"
        ) |> 
        tab_spanner(
          label = htmltools::HTML(
            local_image("./www/analysis_diagnostic_logo.png", height = 57) 
          ), columns = 2, level = 1, id = "logo2"
        ) |> 
        tab_spanner(
          label = htmltools::HTML("<p style = 'margin-left: 23px'; align = 'left'>", "<br>", paste0("ID: ", selID), "<br><br>", paste0(date_title, format(Sys.Date(),format="%d %B %Y")), "<br><br>", "</p>"), columns = 3:4, level = 1, id = "id_date"
        ) |>
        opt_css(css = "#specimenMarkerTable .gt_column_spanner {border-bottom-style: hidden;}") |>
        cols_align('center', columns = everything()) |> 
        fmt_number(columns = Result, rows = 1, decimals = rep_digits[1], sep_mark = ",", pattern = paste0("{x}", " mg/dL")) |> 
        fmt_number(columns = Result, rows = 2, decimals = rep_digits[2], sep_mark = ",") |> 
        sub_missing(
          columns = Result,
          missing_text = "LOQ"
        ) |> 
        text_transform(
          locations = cells_body(columns = 'Distribution'),
          fn = function(column) {
            map(column, ~ stringr::str_split_1(., ', ')) |>
              map(~ rangeGGplot_quantile_scale(.[1], .[2], rangesCRIQ = ref_ranges, plot_vertical_mm = -4.5)) |>
              ggplot_image2(height = 30, aspect_ratio = 3) 
          }
        ) |> 
        cols_label(`Specimen Marker` = cols_label_fix1_1, 
                   Result = cols_label_fix2,
                   Range = cols_label_fix3, 
                   Distribution = "", .fn = md) |> 
        cols_width(`Specimen Marker` ~ px(150)) |> 
        cols_width(c(Result,Range) ~ px(120)) |> 
        opt_table_font(font = google_font("Roboto Slab")) |> 
        tab_options(table.font.size = 12, table.width = pct(80))  
    }
    
    textBold_gtbl1 <- head(concentrations_in_subj_id_res, 2) |> 
      group_by(`Amino Acid`) |> 
      mutate(q = lnorm_valToQuantile(Result, lowerCRIQ, upperCRIQ)) |> 
      mutate(boldCheck = ifelse(!is.na(Result) & (q < extremeThres | q > (1-extremeThres)), "bold", "normal")) |> pull(boldCheck)
    
    for (i in 1:nrow(head(concentrations_in_subj_id_res, 2))) {
      gtbl1 <- gtbl1 |> 
        tab_style(
          style = list(cell_text(weight = textBold_gtbl1[i])),
          locations = cells_body(columns = Result, rows = i)
        ) 
    }

#---- Table 2 ---- 
    gtbl2 <- tail(rep, -2) |> 
      rename(Biomarkers = Indices) |> 
      gt(id = "aminoAcidTable") |> 
      cols_align('center', columns = everything()) |> 
      tab_footnote(
        footnote = result_info_footnote,
        locations = cells_column_labels(columns = Result)
      ) |> 
      fmt_number(columns = Result, rows = 1, decimals = rep_digits[3], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 2, decimals = rep_digits[4], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 3, decimals = rep_digits[5], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 4, decimals = rep_digits[6], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 5, decimals = rep_digits[7], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 6, decimals = rep_digits[8], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 7, decimals = rep_digits[9], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 8, decimals = rep_digits[10], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 9, decimals = rep_digits[11], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 10, decimals = rep_digits[12], sep_mark = ",") |> 
      fmt_number(columns = Result, rows = 11, decimals = rep_digits[13], sep_mark = ",") |> 
      sub_missing(
        columns = Result,
        missing_text = "LOQ"
      ) |> 
      text_transform(
        locations = cells_body(columns = 'Distribution'),
        fn = function(column) {
          map(column, ~ stringr::str_split_1(., ', ')) |>
            map(~ rangeGGplot_quantile_scale(.[1], .[2], rangesCRIQ = ref_ranges, plot_vertical_mm = -5.5)) |> 
            ggplot_image2(height = 30, aspect_ratio = 3) 
        }
      ) |> 
      cols_label(Biomarkers = cols_label_fix1_2, 
                 Result = paste0(cols_label_fix2, "<br>", notes_comments12), 
                 Range = paste0(cols_label_fix3, "<br>", notes_comments12),
                 Distribution = "", .fn = md) |> 
      cols_width(Biomarkers ~ px(150)) |> 
      cols_width(c(Result,Range) ~ px(120)) |> 
      opt_table_font(font = google_font("Roboto Slab")) |> 
      tab_options(table.font.size = 12, table.width = pct(80), column_labels.border.top.color = "white")
      
    textBold_gtbl2 <- tail(concentrations_in_subj_id_res, -2) |> 
      group_by(`Amino Acid`) |> 
      mutate(q = lnorm_valToQuantile(Result, lowerCRIQ, upperCRIQ)) |> 
      mutate(boldCheck = ifelse(!is.na(Result) & (q < extremeThres | q > (1-extremeThres)), "bold", "normal")) |> pull(boldCheck)
    
    for (i in 1:nrow(tail(concentrations_in_subj_id_res, -2))) {
      gtbl2 <- gtbl2 |> 
        tab_style(
          style = list(cell_text(weight = textBold_gtbl2[i])),
          locations = cells_body(columns = Result, rows = i)
        ) 
    }
    
#---- Table 3 ----
    if (logoChoice == 1) {
      
    gtbl3 <- tibble(Indices = "", Result = "", Range = "", Distribution = "") |> 
      gt(id = "allNotes") |> 
      tab_header( # ?tab_header
        title = htmltools::tagList(
          htmltools::tags$div(
            style = htmltools::css(
              `display` = "flex", 
              `align-items` = "center" 
            ),
            htmltools::HTML("&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;"), 
            htmltools::tags$div(
              htmltools::tags$strong(notes_comments9)
            )
          )
        )
      ) |> 
      tab_spanner(
        label = htmltools::HTML(
          local_image("./www/theta_biomarkers_logo.png", height = 86.5)
        ), columns = 1:2, level = 1, id = "logo"
      ) |> 
      tab_spanner(
        label = htmltools::HTML("<p style = 'margin-left: 23px'; align = 'left'>", "<br>", paste0("ID: ", selID), "<br><br>", paste0(date_title, format(Sys.Date(),format="%d %B %Y")), "<br><br>", "</p>"), columns = 3:4, level = 1, id = "id_date"
      ) |>
      opt_css(css = "#allNotes .gt_column_spanner {border-bottom-style: hidden;}") |>
      cols_align('center', columns = everything()) |> 
      cols_width(Indices ~ px(150)) |> 
      cols_width(c(Result,Range) ~ px(120)) |> 
      cols_label_with(fn = ~ "") |> 
      tab_source_note(source_note = footnote_subj_meta_title) |> 
      tab_source_note(source_note = footnote_subj_meta) |> 
      tab_source_note(source_note = notes_comments1) |> 
      tab_source_note(source_note = notes_comments2) |> 
      tab_source_note(source_note = notes_comments3) |> 
      tab_source_note(source_note = notes_comments4) |> 
      tab_source_note(source_note = notes_comments5) |> 
      tab_source_note(source_note = notes_comments6) |> 
      tab_source_note(source_note = gt::md(paste(notes_comments7, '<a href="', "https://thetabiomarkers.com/",'">', "*www.thetabiomarkers.com*", '</a>'))) |> 
      opt_table_font(font = google_font("Roboto Slab")) |> 
      tab_options(table.font.size = 12, table.width = pct(80), table_body.border.bottom.color = "white", table_body.border.top.color = "white", table_body.border.top.width = 3) |> 
      tab_style(style = cell_text(align = "justify"), locations = cells_source_notes())
    } else {
      # logoChoice == 2
      gtbl3 <- tibble(Indices = "", Result = "", Range = "", Distribution = "") |> 
        gt(id = "allNotes") |> 
        tab_header( # ?tab_header
          title = htmltools::tagList(
            htmltools::tags$div(
              style = htmltools::css(
                `display` = "flex", 
                `align-items` = "center" 
              ),
              htmltools::HTML("&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;"), 
              htmltools::tags$div(
                htmltools::tags$strong(notes_comments9)
              )
            )
          )
        ) |> 
        tab_spanner(
          label = htmltools::HTML(
            local_image("./www/theta_biomarkers_logo.png", height = 57)
          ), columns = 1, level = 1, id = "logo1"
        ) |> 
        tab_spanner(
          label = htmltools::HTML(
            local_image("./www/analysis_diagnostic_logo.png", height = 57) 
          ), columns = 2, level = 1, id = "logo2"
        ) |> 
        tab_spanner(
          label = htmltools::HTML("<p style = 'margin-left: 23px'; align = 'left'>", "<br>", paste0("ID: ", selID), "<br><br>", paste0(date_title, format(Sys.Date(),format="%d %B %Y")), "<br><br>", "</p>"), columns = 3:4, level = 1, id = "id_date"
        ) |>
        opt_css(css = "#allNotes .gt_column_spanner {border-bottom-style: hidden;}") |>
        cols_align('center', columns = everything()) |> 
        cols_width(Indices ~ px(150)) |> 
        cols_width(c(Result,Range) ~ px(120)) |> 
        cols_label_with(fn = ~ "") |> 
        tab_source_note(source_note = footnote_subj_meta_title) |> 
        tab_source_note(source_note = footnote_subj_meta) |> 
        tab_source_note(source_note = notes_comments1) |> 
        tab_source_note(source_note = notes_comments2) |> 
        tab_source_note(source_note = notes_comments3) |> 
        tab_source_note(source_note = notes_comments4) |> 
        tab_source_note(source_note = notes_comments5) |> 
        tab_source_note(source_note = notes_comments6) |> 
        tab_source_note(source_note = gt::md(paste(notes_comments7, '<a href="', "https://thetabiomarkers.com/",'">', "*www.thetabiomarkers.com*", '</a>'))) |> 
        opt_table_font(font = google_font("Roboto Slab")) |> 
        tab_options(table.font.size = 12, table.width = pct(80), table_body.border.bottom.color = "white", table_body.border.top.color = "white", table_body.border.top.width = 3) |>  
        tab_style(style = cell_text(align = "justify"), locations = cells_source_notes())
    }
    
    gtblFull <- gt_group(gtbl1, gtbl2, gtbl3)
    
    reportOut$rep <- gtblFull 
    
    gtblFull$gt_tbls$gt_tbl
   
  })

    # message(curl::curl_version()$version) # check curl is installed
    # if (identical(Sys.getenv("R_CONFIG_ACTIVE"), "shinyapps")) {
    #   chromote::set_default_chromote_object(
    #     chromote::Chromote$new(chromote::Chrome$new(
    #       args = c("--disable-gpu", 
    #                "--no-sandbox", 
    #                "--disable-dev-shm-usage", # required bc the target easily crashes
    #                c("--force-color-profile", "srgb"))
    #     ))
    #   )
    # }
    
    output$export <- downloadHandler(
      filename = function() {
        paste0("id_", subject_info$id, "_", format(Sys.Date(),format="%d_%m_%Y"), ".pdf")
        },
      content = function(file) {
        req(reportOut$rep)
        
        gt_save_webshot_rep <- function (data, filename, path = NULL, ..., selector = "table", 
                  zoom = 2, expand = 5) {
          filename <- gt:::gtsave_filename(path = path, filename = filename)
          tempfile_ <- tempfile(fileext = ".html")
          tempfile_ <- gt:::tidy_gsub(tempfile_, "\\\\", "/")
          gt:::gt_save_html(data = data, filename = tempfile_, path = NULL)
          rlang::check_installed("webshot2", "to save gt tables as images.")
          webshot2::webshot(url = paste0("file:///", tempfile_), file = filename, 
                            selector = selector, zoom = zoom, expand = expand, ...)
        }
        
        gtsave_rep <- function (data, filename, path = NULL, ...) {
          gt:::stop_if_not_gt_tbl_or_group(data = data)
          file_ext <- gt:::gtsave_file_ext(filename)
          if (file_ext == "") {
            cli::cli_abort(c("A file extension is required in the provided filename.", 
                             i = "We can use:", `*` = "`.html`, `.htm` (HTML file)", 
                             `*` = "`.png`          (PNG file)", `*` = "`.pdf`          (PDF file)", 
                             `*` = "`.tex`, `.rnw`  (LaTeX file)", `*` = "`.rtf`          (RTF file)", 
                             `*` = "`.docx`         (Word file)"))
          }
          switch(file_ext, htm = , html = gt:::gt_save_html(data = data, 
                                                       filename, path, ...), ltx = , rnw = , tex = gt:::gt_save_latex(data = data, 
                                                                                                                 filename, path, ...), rtf = gt:::gt_save_rtf(data = data, 
                                                                                                                                                         filename, path, ...), png = , pdf = gt_save_webshot_rep(data = data, 
                                                                                                                                                                                                             filename, path, ...), docx = gt:::gt_save_docx(data = data, 
                                                                                                                                                                                                                                                       filename, path, ...), {
                                                                                                                                                                                                                                                         cli::cli_abort(c("The file extension supplied (`.{file_ext}`) cannot be used.", 
                                                                                                                                                                                                                                                                          i = "We can use:", `*` = "`.html`, `.htm` (HTML file)", 
                                                                                                                                                                                                                                                                          `*` = "`.png`          (PNG file)", `*` = "`.pdf`          (PDF file)", 
                                                                                                                                                                                                                                                                          `*` = "`.tex`, `.rnw`  (LaTeX file)", `*` = "`.rtf`          (RTF file)", 
                                                                                                                                                                                                                                                                          `*` = "`.docx`         (Word file)"))
                                                                                                                                                                                                                                                       })
          if (!is.null(path)) {
            filename <- file.path(path, filename)
          }
          invisible(filename)
        }
        
        reportOut$rep |> gtsave_rep(file)
        
        # message(curl::curl_version()$version) # check curl is installed
        # if (identical(Sys.getenv("R_CONFIG_ACTIVE"), "shinyapps")) {
        #   chromote::set_default_chromote_object(
        #     chromote::Chromote$new(chromote::Chrome$new(
        #       args = c("--disable-gpu",
        #                "--no-sandbox",
        #                "--disable-dev-shm-usage", # required bc the target easily crashes
        #                c("--force-color-profile", "srgb"))
        #     ))
        #   )
        # }
        
        # f <- chromote::default_chromote_object()
        # f$close()
        
      }
    )
} 

#---- run app ----
shinyApp(ui = ui, server = server)
