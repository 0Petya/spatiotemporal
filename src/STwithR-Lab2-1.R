################################################################################
####
####  Lab 2-1:  Data Wrangling
####
####  from Spatio-Temporal Statistics with R by WZC
####
################################################################################

library( dplyr )
library( tidyr )
library( STRbook )

locs = read.table(system.file("extdata", "Stationinfo.dat",
                              package = "STRbook"),
                  col.names = c("id", "lat", "lon"))
Times = read.table(system.file("extdata", "Times_1990.dat",
                               package = "STRbook"),
                   col.names = c("julian", "year", "month", "day"))
Tmax = read.table(system.file("extdata", "Tmax_1990.dat",
                              package = "STRbook"))

plot( locs$lon , locs$lat )
names(Tmax)
names(Tmax) = locs$id

Tmax = cbind( Times , Tmax )    ####  Add time stamps to Tmax data
head( names(Tmax) , 20 )

####  Note: Data are in space-wide format
####  gather() is a tidyr function
####  R documentation suggests switching to  pivot_longer

Tmax_long = gather( Tmax ,       ## the space-wide object to gather
                    id ,         ## create an id variable consisting of all
                                 ##   columns of Tmax except see below
                    z ,          ## create a column for the measurement Tmax
                    -julian ,    ## except the julian column
                    -year ,      ##        the year column
                    -month ,     ##        the month column
                    -day ,       ## and    the day column
                  )  
View( Tmax_long )
Tmax_long$id = as.integer( Tmax_long$id )

####
####  Let's filter out missing values (i.e. those with Tmax = -9999 )
####
nrow( Tmax_long )
Tmax_long = filter( Tmax_long , z > -9998 )  ## Notice -9998 and different from WZC
nrow( Tmax_long )

####  Note that more than half of the days had a missing value

Tmax_long = mutate( Tmax_long , proc="Tmax" )  ##  dplyr cmnd to create new column

data( Tmin_long   , package="STRbook" )
data( TDP_long    , package="STRbook" )
data( Precip_long , package="STRbook" )

####  NOAA = National Oceanic Atmospheric Administration
####  https://www.noaa.gov/
####  Noah Rigdon = my grandson

NOAA_df_1990 = rbind( Tmax_long , Tmin_long , TDP_long , Precip_long )

summ = group_by( NOAA_df_1990 , year , proc )  %>% summarize( mean_proc = mean(z) )
View( summ )

####  How many days did it not rain at each station in June of every year?

NOAA_precip = filter( NOAA_df_1990 , proc=="Precip" & month==6 )
summ = group_by( NOAA_precip , year , id ) %>%
       summarize( days_no_precip = sum( z==0 ) )
View( summ )

####  How dplyr's select() works
df1 = select( NOAA_df_1990, julian, z)
df2 = select( NOAA_df_1990, -julian)

####  Pull locs (locations) into NOAA_df_1990
NOAA_df_1990 = left_join( NOAA_df_1990 , locs , by="id" )

####  Select three columns in Tmax_long
Tmax_long_sel = select( Tmax_long , julian , id , z )

####  Rearrange NOAA_df
NOAA_df_sorted = arrange( NOAA_df_1990 , julian , id )


####  spread() does the opposite of gather

Tmax_wide = spread( Tmax_long_sel , ##  Long object 
                    id ,            ##  Column in above obj used to name columns
                    z               ##  Which variable to use to populate data frame
                  )

M = select( Tmax_wide , -julian )   ##  Select all columns of Tmax_wide except julian
str( M )
M = as.matrix(M)                    ##  Convert M to a matrix
str( M )

M = select( Tmax_wide , -julian ) %>% as.matrix()    ##  Book uses the pipe operator
####
####  Recall that the pipe operator is defined as
####  x %>% f(y) means f(x,y)
####  When used repeatedly
####  x %>% f(y) %>% g(z)  means  g(f(x,y),z)
####

library( sp )
library( spacetime )

####
####  Constructing STIDF Object
####  Spatio Temporal Irregular Data Frame
####

NOAA_df_1990$date = with( NOAA_df_1990 , paste( year , month , day , sep="-" ) )
NOAA_df_1990$date = as.Date( NOAA_df_1990$date )
str( NOAA_df_1990 )
class( NOAA_df_1990$date )

####  Now construct an  STIDF  object using  stConstruct()

Tmax_long2 = filter( NOAA_df_1990 , proc=="Tmax" )
STObj = stConstruct( x = Tmax_long2 ,              ##  Data frame
                     space=c("lon","lat") ,        ##  Spatial fields
                     time="date"                   ##  Time field
                   )
class( STObj )
str( STObj )

####  Now construct an  STIDF  object using  STIDF()
####  The Spatial field must be a  Spatial  or  SpatialPoints  object

spat_part = SpatialPoints( coords = Tmax_long2[ , c( "lon" , "lat" ) ] )
temp_part = Tmax_long2$date                         ##  already in Date format
STObj2 = STIDF( sp = spat_part ,                   ##  SpatialPoints data frame
                time = temp_part ,                 ##  Temporal part
                data = select( Tmax_long2 , -date , -lon , -lat )  ##  Data frame
)

str( STObj2 )
class( STObj2 )


####  Now construct an  STIDF  object using  STIDF()
####  The Spatial field must be a  Spatial  or  SpatialPoints  object
####  Spatial points must be fixed in time

spat_part = SpatialPoints( coords = locs[ , c("lon","lat") ] )
temp_part = with( Times , paste( year , month , day , sep="-" ) )
temp_part = as.Date( temp_part )
Tmax_long3 = gather( Tmax , id , z , -julian , -year , -month , -day )

Tmax_long3$id = as.integer( Tmax_long3$id )
Tmax_long3 = arrange( Tmax_long3 , julian , id )  ##  When sending info to
                                                  ##  STFDF the spatial index
                                                  ##  must be moving faster 
                                                  ##  than the temporal index.
?all
all( unique(Tmax_long3$id) == locs$id )

STObj3 = STFDF( sp = spat_part ,
                time = temp_part ,
                data = Tmax_long3 )
proj4string(STObj3) = CRS("+proj=longlat +ellps=WGS84")

STObj3$z[ STObj3$z < -9998 ] = NA
View(STObj3@data)
