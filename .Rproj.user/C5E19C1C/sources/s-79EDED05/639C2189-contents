
data <- get(load("../data/data.rda"))

X <- t(data$data)
latlong2 <- data$latlong

max_n = NULL
numfac = 5
covars = NULL
numfac_init = NULL
override=FALSE
NUM_ITER = 2
cl = NULL

out <- lgamfac(X, numfac = 4, NUM_ITER = 3)



library(ggplot2)
library(maps)
library(mapdata)
library(mapplots)

world_map <- map_data("world")
world_map <- world_map[world_map$region != "Antarctica",] # intercourse antarctica

world_map <- world_map[world_map$long > 90 & world_map$long < 160, ]
world_map <- world_map[world_map$lat > -18 & world_map$lat < 20, ]


p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")

#Add map to base plot
base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), colour="light green", fill="light green")

cleanup <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), legend.position="none",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

base_world <- base_world_messy + cleanup

base_world



PlotAssemblageIdx <- function(vec){
  dat <- cbind.data.frame(latlong2, vec)
  colnames(dat) <- c("Latitude", "Longitude", "Value")
  map_data_coloured <- 
    base_world +
    geom_point(data=dat, 
               aes(x=Latitude, y=Longitude, colour=Value), size=0.5) +
    scale_colour_gradient(low = "white", high = "black") 
  
  map_data_coloured
}

PlotAssemblageIdx(out$F[,1])
PlotAssemblageIdx(out$F[,2])
PlotAssemblageIdx(out$F[,3])
