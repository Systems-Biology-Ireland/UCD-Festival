sim_gompertz <- function(t,dose_1, dose_2,max_size){
  # initalize all zeros (every cell is normal)
  alpha = 1
  total_dose = dose_1/10 + dose_2/4 + 5*dose_1*dose_2
  x_tumor = max_size*exp(log(0.1)*exp(-alpha*t)-total_dose*t)
  # rate grows lower 
  return(x_tumor)
}

sim_circle_diam <- function(x_tumor){
  set.seed(1)
  t_x = rnorm(20*round(x_tumor),0,x_tumor/1000)
  t_y =  rnorm(20*round(x_tumor),0,x_tumor/1000)
  return(list(t_x, t_y))
}

plt_Europe_can_inc <- function(location){
  if (location == "europe"){
    #install.packages("sf")
    library(sf)
    #install.packages("dplyr")
    library(dplyr)
    #install.packages("ggplot2")
    library(ggplot2)
    #install.packages("giscoR")
    library(giscoR)
    # Year
    year_ref <- 2016
    europe <- read.delim("./assets/European Countries.txt")[,1]
    # Get countries
    cntries <- gisco_get_countries(year = year_ref,
                                   resolution = 20, country = europe) %>%
      st_transform(3035)
    # Here we are going to read a table with the incidences:
    inc <- read.delim("./assets/Cancer Incidence.txt")
    # choose palette
    pal <- hcl.colors(5, "Inferno", rev = TRUE, alpha = 0.7)
    # reorder inc. 
    inc_ord = matrix(NA,44,1)
    for (i in 1:44){
      w = which(inc[,1] == cntries$NAME_ENGL[i])
      if (length(w)>0){
        inc_ord[i] = inc[w,3]
      }
    }
    
    breaks <- matrix(0,6,1)
    quants = c(0.1,0.25,0.5,0.75,0.9, 1.0)
    for (i in 1:6){
      breaks[i] = quantile(inc_ord, quants[i], na.rm = TRUE)
    }
    
    #create breaks 
    cntries$inc_brc<- cut(inc_ord,
                          breaks = breaks,
                          dig.lab = 5)
    # reorder breaks to match the cntries object.
    # Focus on European Countries. 
    labs <- breaks[,1]
    labs_plot <- paste0("(", labs[1:5], "-", labs[2:6], "]")
    # Plot
    choro_adv  <- ggplot() +
      geom_sf(data = cntries,
              aes(fill = inc_brc), color = NA) +
      # Labs
      labs(title = "Cancer Incidence",
           caption = "Source: World Cancer Research Fund International",
           fill = "ASR/100k Inhabitants") +
      # Custom palette
      scale_fill_manual(values = pal,
                        drop = FALSE,
                        na.value = "grey80",
                        label = labs_plot,
                        # Legend
                        guide = guide_legend(direction = "horizontal",
                                             nrow = 1,
                                             label.position = "bottom")) +
      # Theme
      theme_void() +
      theme(plot.caption = element_text(size = 7, face = "italic"),
            legend.position = "bottom") +
      # Set limits
      xlim(c(2200000, 7150000)) +
      ylim(c(1380000, 5500000))
    
    return(choro_adv)}
  else{
    library(grid)
    library(tidyverse)
    library(shadowtext)
    names <- c(
      "Skin", "Prostate", "Breast", "Bowel", 
      "Lung"
    )
    
    # Name is an ordered factor. We do this to ensure the bars are sorted.
    data <- data.frame(
      count = c(13311,3890,3704,2819,2753), 
      name = factor(names, levels = names),
      y = seq(length(names)) * 0.9
    )
    BLUE <- "#076fa2"
    RED <- "#E3120B"
    BLACK <- "#202020"
    GREY <- "grey50"
    
    plt  <- ggplot(data) +
      geom_col(aes(count, name), fill = BLUE, width = 0.6) +
      
      # The vertical axis only extends upwards 
      #scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
      theme(
        # Set background color to white
        panel.background = element_rect(fill = "white"),
        # Set the color and the width of the grid lines for the horizontal axis
        panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
        # Remove tick marks by setting their length to 0
        axis.ticks.length = unit(0, "mm"),
        # Remove the title for both axes
        axis.title = element_blank(),
        # Only left line of the vertical axis is painted in black
        axis.line.y.left = element_line(color = "black"),
        # Remove labels from the vertical axis
        axis.text.y = element_blank(),
        # But customize labels for the horizontal axis
        axis.text.x = element_text(family = "Econ Sans Cnd", size = 16)
      ) + 
      geom_text(
        data = data,
        aes(0, y = name, label = name),
        hjust = 0,
        nudge_x = 0.3,
        colour = "white",
        family = "Econ Sans Cnd",
        size = 7
      ) + 
      labs(
        title = "Cancer Cases in Ireland 2020",
        subtitle = "Source: Irish Cancer Society"
      ) + 
      theme(
        plot.title = element_text(
          family = "Econ Sans Cnd", 
          face = "bold",
          size = 22
        ),
        plot.subtitle = element_text(
          family = "Econ Sans Cnd",
          size = 20
        )
      )
    
    return(plt)
    
    
  }
}

plt_canfreq <- function(can_type){
  
  library(grid)
  library(tidyverse)
  library(shadowtext)
  
  # read data.
  freqs = as.matrix(read.delim("./assets/cancer_mutation_frequencies.txt", header = F, sep = " "))
  gens = as.matrix(read.delim("./assets/genes.txt", header = F,sep = " "))
  # read can_type
  if (can_type == "brc"){
    freqs_i = freqs[1,]
  } else if (can_type == "bc"){
    freqs_i = freqs[2,]
  } else if (can_type == "lc"){
    freqs_i = freqs[3,]
  } else if (can_type == "pc"){
    freqs_i = freqs[4,]
  } else if (can_type == "sc"){
    freqs_i = freqs[5,]
  }
  # Get top 10.
  O = order(-freqs_i)[1:10]
  
  # Name is an ordered factor. We do this to ensure the bars are sorted.
  data <- data.frame(
    count = freqs_i[O], 
    name = factor(gens[O], levels = gens[O]),
    y = seq(10) * 0.9
  )
  BLUE <- "#076fa2"
  RED <- "#E3120B"
  BLACK <- "#202020"
  GREY <- "grey50"
  
  plt  <- ggplot(data) +
    geom_col(aes(count, name), fill = BLUE, width = 0.6) +
    
    # The vertical axis only extends upwards 
    #scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
    theme(
      # Set background color to white
      panel.background = element_rect(fill = "white"),
      # Set the color and the width of the grid lines for the horizontal axis
      panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
      # Remove tick marks by setting their length to 0
      axis.ticks.length = unit(0, "mm"),
      # Remove the title for both axes
      axis.title = element_blank(),
      # Only left line of the vertical axis is painted in black
      axis.line.y.left = element_line(color = "black"),
      # Remove labels from the vertical axis
      axis.text.y = element_blank(),
      # But customize labels for the horizontal axis
      axis.text.x = element_text(family = "Econ Sans Cnd", size = 16)
    ) + 
    geom_text(
      data = data,
      aes(0, y = name, label = name),
      hjust = 0,
      nudge_x = 0.0,
      colour = "white",
      family = "Econ Sans Cnd",
      size = 4
    ) + 
    labs(
      title = "Top Mutated Genes",
      subtitle = "Source: The Cancer Genome Atlas"
    ) + 
    theme(
      plot.title = element_text(
        family = "Econ Sans Cnd", 
        face = "bold",
        size = 22
      ),
      plot.subtitle = element_text(
        family = "Econ Sans Cnd",
        size = 20
      )
    )
  
  return(plt)
  
  
}

plt_graph<- function(net,click){
  BLUE <- "#076fa2"
  RED <- "#E3120B"
  BLACK <- "#202020"
  GREY <- "grey50"
  if (net == "net1") {
    adj = matrix(0,5,5)
    adj[1,2] = 1
    adj[2,1] = 1
    adj[2,3] = 1
    adj[2,4] = 1
    adj[3,2] = 1
    adj[4,2] = 1
    adj[2,5] = 1
    adj[5,2] = 1
    cols = matrix(BLUE,5)
    cols[2] = RED
  } else if (net == "net2"){
    adj = matrix(0,5,5)
    adj[1,2] = 1
    adj[2,1] = 1
    adj[2,3] = 1
    adj[3,2] = 1
    adj[3,4] = 1
    adj[4,3] = 1
    adj[4,5] = 1
    adj[5,4] = 1
    adj[1,5] = 1
    adj[5,1] = 1
    cols = matrix(BLUE,5)
    cols[5] = RED
  } else {
    adj = as.matrix(read.delim("./assets/adj_symm.txt", sep = " ", header = F))
    # read _freqs.
    # read data.
    freqs = as.matrix(read.delim("./assets/cancer_mutation_frequencies.txt", header = F, sep = " "))
    if (net == "net3"){
      freq_i = freqs[1,]
    } else  if (net == "net4"){
      freq_i = freqs[2,]
    } else  if (net == "net5"){
      freq_i = freqs[3,]
    } else  if (net == "net6"){
      freq_i = freqs[4,]
    } else  if (net == "net7"){
      freq_i = freqs[5,]
    }
    cols = matrix("i",196)
    quant = quantile(freq_i, 0.90)
    cols[freq_i>quant] = RED
    cols[freq_i<=quant] = BLUE
  }
  # load igraph.
  library(igraph)
  Net = graph_from_adjacency_matrix(adj)
  # read freqs.
  coords = layout_with_kk(Net)
  # load packages
  if (net == "net1" || (net == "net2")){
    plot(Net, layout=coords, edge.color = GREY, vertex.color = cols,
         vertex.size = 30, edge.arrow.size = 0.1 , labels.font = 0, vertex.label = NA,
         rescale = F, xlim = c(min(coords[,1])-1, max(coords[,1])+1), ylim = c(min(coords[,2]), max(coords[,2])))
    
  }
  else {
    plot(Net, layout=coords, edge.color = GREY, vertex.color = cols,
         vertex.size = 50, edge.arrow.size = 0.1 , labels.font = 0, vertex.label = NA, 
         rescale = F, xlim = c(min(coords[,1]), max(coords[,1])), ylim = c(min(coords[,2]), max(coords[,2])))
    
  }
  # delete nodes with clicker.
  coord_x = click$x
  coord_y = click$y
  # match to a node
  if (!is.null(coord_x)){
    dist = ((coord_x - coords[,1])^2 + (coord_y-coords[,2])^2)^1/2
    del_id = which.min(dist) 
    adj[del_id,] = 0
    adj[,del_id] = 0
    Net = graph_from_adjacency_matrix(adj)
    if (net == "net1" || (net == "net2")){
      plot(Net, layout=coords, edge.color = GREY, vertex.color = cols,
           vertex.size = 30, edge.arrow.size = 0.1 , labels.font = 0, vertex.label = NA,
           rescale = F, xlim = c(min(coords[,1])-1, max(coords[,1]))+1, ylim = c(min(coords[,2]), max(coords[,2])))
      
    } 
    else {
      plot(Net, layout=coords, edge.color = GREY, vertex.color = cols,
           vertex.size = 50, edge.arrow.size = 0.1 , labels.font = 0, vertex.label = NA, 
           rescale = F, xlim = c(min(coords[,1]), max(coords[,1])), ylim = c(min(coords[,2]), max(coords[,2])))
    }
  }
  
  
}

cell_dec <- function(net,click){
  if (net == "net1") {
    adj = matrix(0,5,5)
    adj[1,2] = 1
    adj[2,1] = 1
    adj[2,3] = 1
    adj[2,4] = 1
    adj[3,2] = 1
    adj[4,2] = 1
    adj[2,5] = 1
    adj[5,2] = 1
    mutant = 2
  } else if (net == "net2"){
    adj = matrix(0,5,5)
    adj[1,2] = 1
    adj[2,1] = 1
    adj[2,3] = 1
    adj[3,2] = 1
    adj[3,4] = 1
    adj[4,3] = 1
    adj[4,5] = 1
    adj[5,4] = 1
    adj[1,5] = 1
    adj[5,1] = 1
    mutant = 5
    
  } else {
    adj = as.matrix(read.delim("./assets/adj_symm.txt", sep = " ", header = F))
    # read _freqs.
    # read data.
    freqs = as.matrix(read.delim("./assets/cancer_mutation_frequencies.txt", header = F, sep = " "))
    if (net == "net3"){
      freq_i = freqs[1,]
    } else  if (net == "net4"){
      freq_i = freqs[2,]
    } else  if (net == "net5"){
      freq_i = freqs[3,]
    } else  if (net == "net6"){
      freq_i = freqs[4,]
    } else  if (net == "net7"){
      freq_i = freqs[5,]
    }
    quant = quantile(freq_i, 0.90)
    mutant =  which(freq_i>quant)
  }
  mut_deg = sum(adj[mutant,])
  # load igraph.
  library(igraph)
  Net = graph_from_adjacency_matrix(adj)
  # read freqs.
  coords = layout_with_kk(Net)
  # delete nodes with clicker.
  coord_x = click$x
  coord_y = click$y
  # match to a node
  if (!is.null(coord_x)){
    dist = ((coord_x - coords[,1])^2 + (coord_y-coords[,2])^2)^1/2
    del_id = which.min(dist) 
    adj[del_id,] = 0
    adj[,del_id] = 0
  }
  
  change_centrality = sum(adj[mutant,])
  
  if (change_centrality<mut_deg){
    return(list(src = "./assets/cancer_cells_apoptosis.png", alt = "Cancer Cell"))
  } else {
    return(list(src = "./assets/cancer_cell.png", alt = "Cancer Cell"))
  }
  
}