################################################################################
####
####  Lab 2-2:  Data Visualization
####
####  from Spatio-Temporal Statistics with R by WZC
####
################################################################################

library( dplyr )
library( animation )
library( ggplot2 )
library( gstat )
library( maps )
library( STRbook )

set.seed(1)

data( "NOAA_df_1990" , package="STRbook" )

locs = read.table(system.file("extdata", "Stationinfo.dat",
                              package = "STRbook"),
                  col.names = c("id", "lat", "lon"))

####  Select only certain rows of this data set (those rows for which proc=Tmax)
Tmax = filter( NOAA_df_1990 ,
               proc=="Tmax" & month %in% 5:9 & year == 1993 )
head( Tmax )
Tmax$t = Tmax$julian - 728049     ## New time variable

Tmax_1 = subset( Tmax , t %in% c(1,15,30) )

NOAA_plot = ggplot( Tmax_1 ) +
  geom_point( aes( x=lon , y=lat , color=z , size=3 ) ) +
  col_scale( name = "degF" ) +
  xlab( "Longitude (deg)" ) +
  ylab( "Latitude (deg)" ) +
  geom_path( data=map_data("state") , aes( x=long , y=lat , group=group ) ) +
  facet_grid( ~ date ) + 
  coord_fixed( xlim=c(-105,-75) , ylim=c(25,50) ) +
  theme_bw()

# windows( 12 , 7 );  
NOAA_plot
####  This is Figure 2.1, p. 19

####  Skip p. 63 and top of p. 64

####  Time Series Plots

UIDs = unique( Tmax$id )
UIDs_sub = sample( UIDs , 10 )
Tmax_sub = filter( Tmax , id %in% UIDs_sub )

Tmax_TS = ggplot( Tmax_sub ) +
  geom_line( aes( x=t , y=z ) ) +
  facet_wrap( ~id , ncol=5 ) +
  xlab( "Day number (days)" ) +
  ylab( "Max Temperature (degF)" ) + 
  theme_bw() +
  theme( panel.spacing=unit(1,"lines") )  ## Remove last line & see what happens

# windows( 12 , 8 );
Tmax_TS

####  Hovmoller Plots

##    Figure 2.11, p. 30
##    Latitudinal Hovmoller Plot
##    Create a 25 by 100 grid   (spatial by temporal)

lim_lat = range( Tmax$lat )
lim_t   = range( Tmax$t )
lat_axis = seq( lim_lat[1] , lim_lat[2] , length=7 )   ## Change to 25
t_axis   = seq( lim_t[1] , lim_t[2] , length=100 )
lat_t_grid = expand.grid( lat = lat_axis , t = t_axis )

Tmax_grid = Tmax
##  The  outer  function in R
##  Let's get a handle on  outer  

set.seed(1);   lcs = sample(1:328,5)
lats = locs$lat[lcs]
lons = locs$lon[lcs]

# windows(9,7)
plot( lons , lats , pch=19 , ylim=c(32,48))
  abline( h=lat_axis )
  for (i in 1:5) text( lons[i] , lats[i]+1 , as.character(round(i)) )
  for (j in 1:7) text( -84 , lat_axis[j]+0.5 , 
                       as.character(round(lat_axis[j], 2)) , col="red" )
dists = outer( lats , lat_axis , "-" )
dists = abs( dists )

apply( dists , 1 , which.min)
lat_axis[ apply( dists , 1 , which.min) ]

##  Scale up!

lim_lat = range( Tmax$lat )
lim_t   = range( Tmax$t )
lat_axis = seq( lim_lat[1] , lim_lat[2] , length=25 )   
t_axis   = seq( lim_t[1] , lim_t[2] , length=100 )
lat_t_grid = expand.grid( lat = lat_axis , t = t_axis )                      

dists = abs( outer( Tmax$lat , lat_axis , "-" ) )
Tmax_grid$lat = lat_axis[ apply( dists , 1 , which.min) ]
unique( Tmax_grid$lat )

Tmax_lat_Hov = group_by( Tmax_grid , lat , t ) %>%
               summarize( z = mean(z) )

Hovmoller_lat = ggplot( Tmax_lat_Hov ) +
  geom_tile( aes( x=lat , y=t , fill=z ) ) +
  fill_scale( name = "degF" ) +
  scale_y_reverse() +                ## Optionally, leave this out
  ylab( "Day number (days)" ) +
  xlab( "Latitude (degrees)" ) +
  theme_bw()
  
# windows( 9 , 7 );
Hovmoller_lat

Hovmoller_lat_1 = ggplot( Tmax_lat_Hov ) +
  geom_tile( aes( x=t , y=lat , fill=z ) ) +
  fill_scale( name = "degF" ) +
  xlab( "Day number (days)" ) +
  ylab( "Latitude (degrees)" ) +
  theme_bw()

# windows( 9 , 7 );
Hovmoller_lat_1

####  Animations

Tmax_t = function( tau ) 
{
  Tmax_sub = filter (Tmax , t==tau )
  ggplot( Tmax_sub ) +
    geom_point( aes( x=lon , y=lat , size=3 , color=z ) ) +
    col_scale( name="z" , limits=c(40,110) ) +
    theme_bw()
}

# windows( 9 , 7 );
Tmax_t(21)


gen_anim = function()
{
  for( t in lim_t[1]:lim_t[2] )
  {
    plot( Tmax_t(t) )
  }
}


gen_anim <- function() {
  for(t in lim_t[1]:lim_t[2]){  # for each time point
    plot(Tmax_t(t))            # plot data at this time point
  }
}

library( animation )
ani.options( interval=0.2 )     ##  time interval in sec between frames

####  saveHTML() saves the animation to a file, given here as the last argument

dir.create("./figures/NOAA_anim/")
saveHTML( gen_anim(),            ##  run function that draws graph
          autoplay = FALSE,      ##  don't automatically play on load
          loop = FALSE,          ##  don't loop automatically
          verbose = FALSE,       ##  suppress lots of output
          outdir = ".",          ##  save to current directory
          single.opts = "'controls': ['first', 'previous',
                                      'play', 'next', 'last',
                                      'loop', 'speed'],
                        'delayMin': 0",
         htmlfile = "./figures/NOAA_anim/NOAA_anim.html")  

