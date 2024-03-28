################################################################################
####
####  Lab 2-3:  Exploratory Data Analysis & Descriptive Statistics
####
####  from Spatio-Temporal Statistics with R by WZC
####
################################################################################

library( CCA )
library( dplyr )
library( ggplot2 )
library( gstat )
library( sp )
library( spacetime )
library( STRbook )
library( tidyr )
library( fields )    ####  Add this one

set.seed(1)

data( "NOAA_df_1990" , package="STRbook" )

locs = read.table(system.file("extdata", "Stationinfo.dat",
                              package = "STRbook"),
                  col.names = c("id", "lat", "lon"))

####  Select only certain rows of this data set (those rows for which proc=Tmax
####  and month is 5, 6, 7, 8, or 9)

####  Select only certain rows of this data set (those rows for which proc=Tmax)

Tmax = filter( NOAA_df_1990 ,
               proc=="Tmax" & month %in% 5:9 & year == 1993 )
head( Tmax )
Tmax$t = Tmax$julian - 728049     ## New time variable

Tmax_1 = subset( Tmax , t %in% c(1,15,30) )
Tmax_2 = subset( Tmax , t %in% c(14,15,16) )

NOAA_plot = ggplot( Tmax_2 ) +
  geom_point( aes( x=lon , y=lat , color=z , size="3" ) ) +
  col_scale( name = "degF" ) +
  xlab( "Longitude (deg)" ) +
  ylab( "Latitude (deg)" ) +
  geom_path( data=map_data("state") , aes( x=long , y=lat , group=group ) ) +
  facet_grid( ~ date ) + 
  coord_fixed( xlim=c(-105,-75) , ylim=c(25,50) ) +
  theme_bw()

windows( 12 , 5 );  NOAA_plot

unique( locs$id )
unique( Tmax$id )

####  Empirical Spatial Means

spat_av = group_by( Tmax , lat , lon )  %>%
          summarize( mu_emp = mean(z) )

spat_av_plot = ggplot( spat_av ) +
  geom_point( aes( x=lon , y=lat , color=mu_emp , size=1 ) ) +
  col_scale( name = "degF" ) +
  xlab( "Longitude (deg)" ) +
  ylab( "Latitude (deg)" ) +
  geom_path( data=map_data("state") , aes( x=long , y=lat , group=group ) ) +
  coord_fixed( xlim=c(-105,-75) , ylim=c(25,50) ) +
  theme_bw()

windows(9,7);  spat_av_plot

lat_means = ggplot( spat_av ) +
  geom_point( aes( x=lat , y=mu_emp )) +
  xlab( "Latitude (deg)" ) +
  ylab( "Average Maximum Temperature (degF)" )
windows( 9 , 7 );  lat_means

lon_means  = ggplot( spat_av ) +
  geom_point( aes( x=lon , y=mu_emp )) +
  xlab( "Longitude (deg)" ) +
  ylab( "Average Maximum Temperature (degF)" )
windows( 9 , 7 );  lon_means

####  Empirical Temporal Means

Tmax_av = group_by( Tmax , date )  %>%  
  summarize( meanTmax = mean(z) )

gTmaxav = ggplot() +
  geom_line( data = Tmax ,
             aes( x=date , y=z , group=id ) ,
             color="blue" , alpha=0.04 ) +
  geom_line( data=Tmax_av , 
             aes( x=date , y=meanTmax ) ,
             size=1.1 ) +
  xlab( "Month" ) +
  ylab( "Maximum Temperature (degF)" ) +
  theme_bw()
windows( 9 , 7 );  gTmaxav

####  Empirical Covariances
####  Start here

lm1 = lm( z ~ lat + t + I(t^2) ,
          data = Tmax )
Tmax$residuals = residuals( lm1 )

spat_df = filter( Tmax , t==1 )  %>%
  select( lon , lat )  %>%
  arrange( lon , lat )
m = nrow( spat_av )

X = select( Tmax , lon , lat , residuals , t )  %>%    ## Select desired columns
  spread( t , residuals )  %>%                         ## Make time-wide
  select( -lon , -lat )  %>%                           ## Drop lon & lat
  t()                                                  ## Transpose (make space-wide)
?cov

Lag0_cov = cov( X , use='complete.obs' )
Lag0_cov = cov( X , X , use='complete.obs' )
Lag1_cov = cov( X[-1, ] , X[-nrow(X), ] , use='complete.obs' )
####  strips out 1st row of X  &  strips out last row of X

Lag0_cor = cor( X , use='complete.obs' )
Lag1_cor = cor( X[-1, ] , X[-nrow(X), ] , use='complete.obs' )

####  Too many correlations to make sense of
####  Add some STRUCTURE to reduce complexity of model
####  Split the domain into four longitudinal strips

spat_df$n = 1:nrow( spat_df )        ##  Each station gets an index
lim_lon = range( spat_df$lon )       ##  Get the range of all the long. strips
lon_strips = seq( lim_lon[1] , lim_lon[2] , length=5 )
                                     ##  5 cutoffs make four strips
spat_df$lon_strip = cut(  spat_df$lon , lon_strips , 
                          labels=FALSE , include.lowest = TRUE )

plot_cov_strips( Lag0_cov , spat_df )  
plot_cov_strips( Lag1_cov , spat_df )

data( "STObj3" , package="STRbook" )        ##  Object of class STFDF

STObj4 = STObj3[ , "1993-07-01::1993-07-31" ]
vv = variogram( object = z ~  + lat ,       ##  fixed effect component
                data = STObj4 ,             ##  use only July data
                width = 80 ,                ##  spatial bin = 80 km
                cutoff = 1000 ,             ##  consider points < 1000 km apart
                tlags = 0.01:6.01           ##  time lags of 0 to 6 days
               )
plot( vv )


####
####  from https://rdrr.io/github/andrewzm/STRbook/src/R/Cemp.R
####
plot_cov_strips <- function(C,spat_df) {  
  # for each longitudinal strip 
  require(fields)   # load fields for plotting
  for(i in seq_along(unique(spat_df$lon_strip))){                        
    spat_strip <- spat_df %>%         # take spat_df
      filter(lon_strip == i)  %>%   # extract the ith strip
      arrange(lat)                  # sort by latitude
    idx <- spat_strip$n               # extract indices of locations
    jitter <- seq(0,0.0001,           # add jitter for locations that
                  length=length(idx)) # have same latitude component
    image.plot(spat_strip$lat+jitter, # plot the matrix using fields
               spat_strip$lat+jitter,
               C[idx,idx],            # subset and permute C
               xlab="latitude",
               ylab="latitude",
               zlim=c(-15,85),
               col=tim.colors(10),
               cex=200)
  }
}


