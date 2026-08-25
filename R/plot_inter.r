
#' Interface for plotting of impulse responses and covariance sequences. 
#'
#' @param data     list of matrices containing the data
#' @param type   string; either "IR" for impulse response or "CV" for covariance sequence
#' @param interactive Boolean, indicating, if interactive plot is wanted. Default: FALSE.
#'
#' @return x    invisibly returns the input \code{x} after being coerced to an array.
#'
#' @export 
plot_inter <- function(data = list(),type = "IR",interactive = FALSE){
  
  #s <- dim(data[[1]])[1]
  # put all systems into a vector. 
  if (is.null(names(data)) == FALSE){
    data_v <- list(1)
    data_v[[1]] <- data
  } else {
    data_v <- data
  }

  sout <- dim(data_v[[1]]$C)[1]
  
  if (data_v[[1]]$type == "innov"){
    sin <- sout
  } else {
    sin <- dim(data_v[[1]]$Q)[2]+sout
  }
  
  if (interactive){
    
    # specify the UI 
    ui <- shiny::fluidPage(
      
      # Application title
      shiny::titlePanel(type),
      
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput( 
            "input", 
            "Select input channels:", 
            paste(1:sin), 
            multiple = TRUE,
            selected = paste(1)
          ), 
          shiny::selectInput( 
            "output", 
            "Select output channels:", 
            paste(1:sout), 
            multiple = TRUE,
            selected = paste(1) 
          ), 
          shiny::checkboxInput("all", "Plot all?", value = FALSE
          )
        ),
        
        shiny::mainPanel(
          shiny::plotOutput("pl")
        )
      )
    )
    
    # specify the server 
    server = function(input, output, session){
      
      ## put more server logic here
      
      output$pl <- shiny::renderPlot({
        # get parameters 
        pl_all <- input$all
        inp <- as.numeric(input$input)
        outp <- as.numeric(input$output)
        
        if (pl_all == TRUE){
          # plot a matrix 
          if (type == "IR"){
            p <- plot_IR(syst_vec=data_v, channel_in=1:sin, channel_out=1:sout, L =20)
          } 
          if (type == "CV"){
            p <- plot_cov(syst_vec=data_v, channel_in=1:sin, channel_out=1:sin, L =20)
          } 
          
        } else {
          if (type == "IR"){
            p<-plot_IR(syst_vec=data_v, channel_in=inp, channel_out=outp, L =20)
          } 
          if (type == "CV"){
            p<-plot_cov(syst_vec=data_v, channel_in=inp, channel_out=outp, L =20)
          } 
          
        }
        plot(p)
      })
      
    }
    
    # now start the shiny app
    app <- shiny::shinyApp(ui, server)
    shiny::runApp(app)
  } else { # static plot 
    if (type == "IR"){
      p<-plot_IR(syst_vec=data_v, channel_in=1:sin, channel_out=1:sout, L =20)
    } 
    if (type == "CV"){
      p<-plot_cov(syst_vec=data_v, channel_in=1:sout, channel_out=1:sout, L =20)
    } 
    plot(p)
  }
  
}
