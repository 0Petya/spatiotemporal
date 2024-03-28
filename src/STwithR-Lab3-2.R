################################################################################
####
####  Lab 3-2:  Trend Prediction
####
####  from Spatio-Temporal Statistics with R by WZC
####
################################################################################

library( leaps );           library( lmtest );         library( nlme )
library( spacetime );       library( STRbook );        library( ape )
library( broom );           library( FRK );            library( purrr )
library( lattice );         library( ggplot2 );        library( RColorBrewer )
library( dplyr );           library( gstat );          
library( sp );              library( tidyr )

data( "NOAA_df_1990" , package="STRbook" )
Tmax = filter( NOAA_df_1990 ,            ## subset the data
               proc == "Tmax" &          ## only max temperature
               month == 7 &              ## July
               year == 1993)             ## year of 1993
?auto_basis
G = auto_basis( data = Tmax[,c("lon","lat")] %>% SpatialPoints() ,   
                nres = 1 ,
                type = "Gaussian" )      ## Basis functions shown in Figure 3.6
class(G)
str(G)
##  G is a complicated object!

S = eval_basis( basis=G , s=Tmax[,c("lon","lat")] %>% as.matrix() ) %>% as.matrix() 
S = eval_basis( basis=G , s=as.matrix( Tmax[,c("lon","lat")] ) ) %>% as.matrix() 
S = as.matrix( eval_basis( basis=G , s=as.matrix( Tmax[,c("lon","lat")] ) ) )
colnames(S) = paste0("B",1:ncol(S))
print( paste0("B",1:ncol(S)) )

Tmax2 = cbind( Tmax , S )  %>%  
        select( -year , -month , -proc , -julian , -date )
##  cbind Tmax and S, but remove variables that are not in the model

Tmax_no_14 = filter( Tmax2 , !(day==14) )

####  Tangent
x1 = c( 2 , 4 , 6, 8, 10 )
x2 = c( 20, 15 , 25 , 14 , 23 )
y =  c( 100 , 105 , 166 , 157 , 171 )
lm1 = lm( y ~ (x1 + x2)^2 )     ## includes interaction (cross product) terms
                                ## but not squared terms
summary( lm1 )
####  End of Tangent

Tmax_July_lm = lm( z ~ (lon + lat + day)^2 + . ,    ##  "dot" indicates all other vars
                   data = select(Tmax_no_14,-id) )
summary( Tmax_July_lm )

#-------------------------------------------------------------------------------
#-
#-  Assume some *structure*
#-  1. covariance has Gaussian structure
#-  2. isotropic
#-  3. range parameter of 0.5
#-
#-------------------------------------------------------------------------------

?gls
##  gls is a function in the nlme (NonLinear Mixed Effects) package
##  Uses *generalized* least squares

Tmax_July_gls = gls( z ~ (lon + lat + day)^2 + . , 
                     data = select(Tmax_no_14,-id) ,
                     correlation=corGaus( value = 0.5 ,
                                          form = ~ lon + lat + day ,
                                          fixed = TRUE ) ) 
summary( Tmax_July_gls )

pred_grid = expand.grid( lon = seq(-100, -80, length = 20),    ## 20 longitude
                         lat = seq(32, 46, length = 20),       ## 20 latitide
                         day = seq(1, 31, length = 31) )       ## 31 length

