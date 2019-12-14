theme_Publication <- function(base_size=16, base_family="Helvetica", legend_direction="horizontal", legend_position='bottom') {
   library(grid)
   library(ggthemes)
   (theme_foundation(base_size=base_size, base_family=base_family)
      + theme(plot.title = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              text = element_text(),
              panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              #panel.border = element_rect(colour = NA),
              panel.border = element_rect(colour = "black"),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(), 
              axis.line = element_line(colour="black"),
              axis.ticks = element_line(),
              legend.key = element_rect(colour = NA),
              legend.position = legend_position,
              legend.direction = legend_direction,
              legend.key.size= unit(0.5, "cm"),
              legend.key.height = unit(.5,"cm"),
              legend.margin = unit(0, "cm"),
              legend.title = element_text(face="italic"),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_blank(),
              strip.text = element_text(face="bold"),
              strip.placement = "outside"
      ))
   
}

scale_fill_Publication <- function(...){
   library(scales)
   discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
   
}

scale_colour_Publication <- function(...){
   library(scales)
   discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
   
}