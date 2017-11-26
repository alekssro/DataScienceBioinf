
nature_theme <- function(base_size=7, base_family="") {
   library(ggplot2)
   library(grid)
   library(ggthemes)
   tmp = theme_foundation(base_size=base_size, base_family=base_family) + 
         theme(plot.title    = element_text(family=base_family, face="plain", size = base_size), 
               text             = element_text(family=base_family, face="plain", size = base_size), 
               panel.background = element_rect(colour = NA),
               plot.background  = element_rect(colour = NA),
               panel.border     = element_rect(colour = NA),
               axis.title       = element_text(family=base_family, face="plain", size = base_size), 
               axis.title.y     = element_text(angle=90,vjust =2),
               axis.title.x     = element_text(vjust = -0.2),
               axis.text        = element_text(family=base_family, face="plain", size = base_size), 
               axis.line        = element_line(size = 0.3, colour="#000000"), # 0.3 mm is ~1 pt
               axis.ticks       = element_line(size = 0.3, colour="#000000"), # 0.3 mm is ~1 pt
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key       = element_rect(colour = NA),
               legend.key.size  = unit(2, "mm"),
               legend.position  = "right",
               legend.direction = "vertical",
               legend.title     = element_text(family=base_family, face="plain", size = base_size),
               plot.margin      = unit(c(3,3,3,3),"mm"),
               strip.background = element_blank(),
               strip.text       = element_text(family=base_family, color="#000000",face="plain", size = base_size, margin = margin(1,1,1,1,unit="mm")),
               complete = TRUE )
   return(tmp)
}

discrete_colors = c("#386cb0","#ef3b2c","#fdb462","#7fc97f","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
parula_colors = c('#352a87', '#0f5cdd', '#1481d6', '#06a4ca', '#2eb7a4', '#87bf77', '#d1bb59', '#fec832', '#f9fb0e')
