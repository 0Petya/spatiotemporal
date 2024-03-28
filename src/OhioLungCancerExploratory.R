
####  install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

####  Code taken from
####  https://groups.google.com/g/r-inla-discussion-group/c/tPoKFpK5vqs

##  setwd("C:/Users/srigd/Dropbox/MyCourses/BST_5620_Spring_2023")
setwd("C:/Users/serig/Dropbox/MyCourses/BST_5620_Spring_2023")
library(maptools)
library(sp)
library(INLA)
library(tigris)
library(SpatialEpi)
library(tmap)
library(spdep)
library(nimble)
library(dplyr)

##  url <- "https://sites.google.com/a/r-inla.org/stbook/Ohio_data.zip"
##  temp <- tempfile()
##  download.file(url, temp)
##  unzip(temp, "OhioRespMort1.csv")
OH.data <- read.csv("OhioRespMort.csv")

##  unzip(temp, exdir = tempdir())
##  shapefile <- paste0(tempdir(),"\\Ohio_data\\tl_2010_39_county00")
##  ohio <- readShapeSpatial(shapefile)
##  ohio = rgdal::readOGR(tl_2010_39_county00)

ohio <- counties(state = 'OH')    ##  Get shape file from tigris
View( ohio )                      ##  Let's investigate the object

class( ohio )
ohio.nb = poly2nb( ohio )

ohio.lw = nb2listw( ohio.nb )

k = nrow( ohio )
num = rep(0,k)
for (i in 1:k) num[i] = length( ohio.lw$neighbours[[i]] )
adj = c()
for (i in 1:k) adj = c(adj,ohio.lw$neighbours[[i]] )
L = length(adj)

plotLC.Ohio = function( year )
{
  ohio.year = OH.data[ OH.data$year == year , ]
  P = ohio.year$n
  y = ohio.year$y
  LCrate = y/P
  ohio.year$LCrate = LCrate
  ohio.temp = left_join( ohio , ohio.year , by="NAME" )  
  brks = seq( 0 , 0.0014 , 0.0002 )
  tm_shape( ohio.temp ) + 
  tm_polygons( col = "LCrate" , breaks=brks , legend.show = FALSE ) + 
  tm_layout( title = paste0( round(1967+year,0) ), 
             title.position = c('left', 'top') ,
             title.size = 0.8 )  
}

p1 = plotLC.Ohio( 1 )
p2 = plotLC.Ohio( 2 )
p3 = plotLC.Ohio( 3 )
p4 = plotLC.Ohio( 4 )
p5 = plotLC.Ohio( 5 )
p6 = plotLC.Ohio( 6 )
p7 = plotLC.Ohio( 7 )
p8 = plotLC.Ohio( 8 )
p9 = plotLC.Ohio( 9 )
p10 = plotLC.Ohio( 10 )
p11 = plotLC.Ohio( 11 )
p12 = plotLC.Ohio( 12 )
p13 = plotLC.Ohio( 13 )
p14 = plotLC.Ohio( 14 )
p15 = plotLC.Ohio( 15 )
p16 = plotLC.Ohio( 16 )
p17 = plotLC.Ohio( 17 )
p18 = plotLC.Ohio( 18 )
p19 = plotLC.Ohio( 19 )
p20 = plotLC.Ohio( 20 )
p21 = plotLC.Ohio( 21 )

windows( 15 , 9 )
tmap_arrange( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
              p13,p14,p15,p16,p17,p18,p19,p20,p21,
              ncol = 7,nrow = 3 )

LC.matrix = matrix( 0 , nrow=88 , ncol=21 )    ####  rows are counties; columns are years
for ( year in 1:21 )
{
  ohio.year = OH.data[ OH.data$year == year , ]
  P = ohio.year$n
  y = ohio.year$y
  LCrate = y/P
  LC.matrix[,year] = LCrate
}

windows( 9 , 7 )
plot( 1:21 , LC.matrix[1,] , col="gray" , ylim=c(0,0.0014) )
for (cnty in 2:88 ) lines( 1:21 , LC.matrix[cnty,] , col="gray" )

OH.data$LCrate = OH.data$y / OH.data$n






