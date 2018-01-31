
###############   Test lgamfac ()  model   #########################

data <- get(load("../data/data.rda"))

X <- t(data$data)
latlong2 <- data$latlong

d <- 8
m <- nrow(X)
n <- ncol(X)

# check for d validity
if(d != as.integer(d)){
  stop("d should be integer")
} else if(d < 1){
  stop("d should be at least 1")
} else if(d == 1){
  return(matrix(1, n, 1))
} else if(d >1){
  d <- d #for the svd stuff
}

NA_IND <- is.na(X)

#center the matrix...
mean_X <- rowMeans(X, na.rm=TRUE)
norm_X <- X - mean_X

# ...then 'impute'
norm_X[NA_IND] <- 0
override <- FALSE
adjust <- 15

d_init = 30
mysvd <- trunc.svd(norm_X, d=d_init, adjust=adjust, tol=1e-13, override=override)

rm(norm_X)
D <- diag(mysvd$d, d_init, d_init)
U <- mysvd$u
V <- mysvd$v
rm(mysvd)

# form projection
z <- U %*% D %*% t(V)

z <- z + mean_X

#zmin <- apply(z, 1, min)
#zmax <- apply(z, 1, max)
#ind  <- (zmax<(1-2/n)) & (zmin>(2/n))

#z <- z[ind,]
z1 <- z
zmin <- min(z1)
zmax <- max(z1)
z <- t(apply(z, 1, function(x) return(2/n + (1 - 4/n)*((x - zmin)/(zmax - zmin)))))


z <- log(z/(1-z))

norm_z <- centerscale(z)

v <- trunc.svd(norm_z, d=d, adjust=adjust, tol=1e-13, override=override)$v
v <- cbind(v,1)


onerun <- function(xj, argl){
  argl$y <- xj
  fit <- do.call(cv.gamlr,argl)
  out <- c(fit$gamlr$beta[,fit$seg.min], fit$gamlr$alpha[fit$seg.min])
  return(out)
}

cl <- makeCluster(10,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))
print(cl)

NUM_ITER <- 5
for(iter in 1:NUM_ITER){
      verb <- 0
      argl <- list()
      argl$family <- "binomial"
      argl$nlamda = NULL
      argl$x = v[,1:4]
      argl$cv = TRUE
      if(is.null(argl$nlambda))
        argl$nlambda <- formals(gamlr)$nlambda
      argl$verb <- max(verb-1,0)
      
      
      C <- dim(X)[1]
      chunks <- lapply(1:C, 
                       function(i) X[i,])
      
      cat("distributed run.\n") 
      print(cl)
      mods <- parLapply(cl,chunks,onerun,argl=argl) 
      h <- do.call(rbind, mods)
      
      
      D <- dim(X)[2]
      chunks2 <- lapply(1:D, 
                        function(i) X[,i])
      argl2 <- argl
      argl2$x <- h[,1:4]
      mods <- parLapply(cl,chunks2,onerun,argl=argl2) 
      v <- do.call(rbind, mods)
      
      cat("we are at iteration,", iter, "\n")
}

stopCluster(cl)

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

PlotAssemblageIdx(v[,1])
PlotAssemblageIdx(v[,2])
PlotAssemblageIdx(v[,3])
PlotAssemblageIdx(v[,4])
PlotAssemblageIdx(v[,5])

PlotAssemblageIdx(fits2[,1])


