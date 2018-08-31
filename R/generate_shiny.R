#' Default plot function for plsr shiny app
#'
#' @param x A vector of predicted data to plot
#' @param time_steps Number of time steps
#' @param t Which time step to plot
plot_default = function(x,time_steps=10,t){
  time_step_len = length(x)/time_steps

  a = ((t-1)*time_step_len)+1
  b = t*time_step_len

  barplot(x[a:b],ylim = c(min(x),max(x)))
}

#TODO: somewhere bug: works only if plsr_obj in environment
create_shiny = function(plsr_obj,time_steps, app_path=".", plot_func="plot_default",...){
  doc_start = sprintf('time_range <- c(1,%s)
time_start_val <- time_range[2]/2

# Define UI for app
ui <- fluidPage(

  # App title ----
  titlePanel("Emotions"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
sidebarPanel(',time_steps)


#slider generation
slider_string = ""
for (s in 1:ncol(plsr_obj$orig_data$X)){
  new_slider=sprintf('sliderInput(inputId = %s, label = %s, min = -300, max = 300, value = 0),', paste0("'s",s,"'"), paste0("'Slider ", s,"'") )
  slider_string = paste0(slider_string,new_slider)
}

time_slider = 'sliderInput(inputId = "time", label = "time", min = time_range[1], max = time_range[2], value = time_start_val),'
slider_string = paste0(slider_string,time_slider)

  doc_middle = '      fileInput(
        "modelfile", "Upload model",
  accept = c(".Rdata", ".rda", ".Rdat")
  ),
  tags$div(
  id="divActiveModels", checked=NA,
  checkboxGroupInput(
  "activeModels", "Select models for display",
  choices = c("default"),
  selected = c("default")
  )
  ),
  actionButton("play", "Play!")

  ),

  # Main panel for displaying outputs ----
  mainPanel(
  fluidRow(
  #textOutput(outputId="labels"),
  plotOutput(
  outputId = "distPlot"
  )
  )


  )
  )
  )
# Define server logic required to draw the plots
server <- function(input, output, session) {

  # set maximum upload size to 10MB
  options(shiny.maxRequestSize=10*1024^2)

  # init reactive context
  context <- reactiveValues()
  context$running <- FALSE

  models <- reactiveValues()
  models$default <- plsr_obj

  observeEvent(input$play, {
  context$running <- TRUE
  })

  # create non-reactive buffer
  active <- c()

  observeEvent(input$activeModels, {
  active <<- names(models)[names(models) %in% input$activeModels]
  })

  observeEvent(input$modelfile, {
  # event gets triggered when files are uploaded
  # load model
  name <- input$modelfile$name

  models[[name]] <<- get(
  load(input$modelfile$datapath)
  )

  removeUI("#activeModels")

  # and add model to checkbox
  insertUI(
  "#divActiveModels",
  where = "afterEnd",
  ui = checkboxGroupInput(
  "activeModels", "Select models for display",
  choices = names(models),
  selected = c(active, name)
  )
  )
  })

  output$distPlot <- renderPlot({
  if ((input$time < time_range[2]) && context$running){
  # raise time state
  updateSliderInput(
  session,
  inputId = "time",
  label = "time",
  value = {
  # fix this to change moving "speed"
  input$time + 5
  }
  )

  # this causes the plot to re-render every 10 ms
  invalidateLater(10)
  } else if (context$running){
  context$running <- FALSE
  updateSliderInput(
  session,
  inputId = "time",
  label = "time",
  value = 0
  )
  }
'
  #slider readout
  readout_string = 'new_vec = c('
  for (s in 1:ncol(plsr_obj$orig_data$X)){
    if (s==ncol(plsr_obj$orig_data$X)){
      readout_string=paste(readout_string,paste0("input$s",s,")"))
    }else{
      readout_string=paste(readout_string,paste0("input$s",s,","))
    }
  }

  doc_end = sprintf('active <- input$activeModels
    if (!is.null(active)){
  %s(predict(plsr_obj,new_vec),%s)
}

}, height=1000, width=1000
)
}
shiny::shinyApp(ui = ui, server = server)',plot_func,...)

  full_doc = paste0(doc_start,"\n",slider_string,"\n",doc_middle,"\n",readout_string,"\n",doc_end) # complete string of shiny app.R file

  full_path=paste(app_path,"app.R",sep="/")
  fileConn=file(full_path)
  writeLines(full_doc,fileConn)
  close(fileConn)
  shiny::runApp(full_path)
}

#TODO: I feel the passing of additional arguments could be handled in a better way, i.e. no strings
#wrapper function to create shinyapp for face plots
create_shiny_face=function(plsr_obj,tstp,app_path="."){
  create_shiny(plsr_obj,time_steps = tstp,app_path = app_path, plot_func = "plot_frame","single_frame=input$time")
}

create_shiny_default = function(plsr_obj,tstp,app_path="."){
  create_shiny(plsr_obj,time_steps = tstp,app_path = app_path, plot_func = "plot_default","time_steps=time_range[2], t=input$time")
}
