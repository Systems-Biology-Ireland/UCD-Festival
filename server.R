function(input, output) {
  output$distPlot <- renderPlot({
    max_size = 1000
    x_tumor = sim_gompertz(input$obs/10, input$d1, input$d2, max_size)
    res = sim_circle_diam(x_tumor)
    col_cancer = rgb(1, 0.5, 0.5, 0.5)
    plot( res[[1]],res[[2]], col  = col_cancer, pch = 19, xlim = c(-3, 3),
          labels = F, ylim = c(-3,3), xlab = "", ylab = "", tck = 0)
  })
  
  output$statsPlot <- renderPlot({plt_Europe_can_inc(input$Stats)})
  
  output$freqshow <- renderPlot({plt_canfreq(input$cancer_type)})
  
  output$cellshow <- renderImage({ 
    if(input$cancer == "hc"){
      list(src = "./assets/spheroidal_cell.png", alt = "Normal Cell")
    } else if(input$cancer == "hdiv") {
      list(src = "./assets/spheroidal_cell_cluster.png", alt = "Normal Growth")
    } else if (input$cancer == "cdiv") {
      list(src = "./assets/cancer_cells_large.png", alt = "Fast Growth")
    } else {
      list(src = "./assets/cancer_cell.png", alt = "Cancer Cell")
    }
    },deleteFile = FALSE)
  
  clickys <- reactiveValues(x = NULL, y = NULL)
  
  output$netshow <- renderPlot({plt_graph(input$nets, clickys)})
  
  output$celldec <- renderImage({cell_dec(input$nets, clickys)})
  
  observeEvent( input$plot1_click, {
    brush <- input$plot1_click
    if (!is.null(brush)) {
      clickys$x <- c(brush$x, brush$x)
      clickys$y <- c(brush$y, brush$y)
      
    } else {
      clickys$x <- NULL
      clickys$y <- NULL
    }
  })
  
  observeEvent( input$plot1_dblclick,{
    clickys$x <- NULL
    clickys$y <- NULL
  })
  
  

}




