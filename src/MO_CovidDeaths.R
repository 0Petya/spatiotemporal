# setwd("C:/Users/serig/Dropbox/MyCourses/BST_5620_Spring_2023")
MO.data = read.csv("./data/MO_COVID_Deaths.csv")
library( dplyr )

MO.cum.deaths = MO.data[,5:ncol(MO.data)]
MO.cum.deaths.mat = as.matrix( MO.cum.deaths )
I = nrow( MO.cum.deaths.mat )
J = ncol( MO.cum.deaths.mat )
K = J/7

MO.daily.deaths.mat = MO.cum.deaths.mat[,2:J] - MO.cum.deaths.mat[,1:(J-1)]

start = 41
plot( MO.cum.deaths.mat[start,] , type="l" , col="gray" , ylim=c(0,3000) )
for( i in start:(start+19) ) lines( MO.cum.deaths.mat[i,] , col="gray" )
lines( MO.cum.deaths.mat[48,] , col="red" )

MO.daily.deaths.mat = MO.daily.deaths.mat[,6:(J-1)]
## Above guarantees that the number of points in time series is a multiple of 7
## Note:  J = 1091 which is not a multiple of 7.  J = 1085 = 155*7
str( MO.daily.deaths.mat )

MO.weekly.deaths.mat = matrix( -999 , nrow=I , ncol=K )
for (i in 1:I)
{
  for (k in 1:K)
  {
    MO.weekly.deaths.mat[i,k] = sum( MO.daily.deaths.mat[i,(7*(k-1)+1):(7*k)] )
  }
}
write.csv( MO.weekly.deaths.mat , "./data/MO.weekly.deaths.csv" )

MO.weekly.deaths.mat1 = MO.weekly.deaths.mat

for (i in 1:I)
{
  for (k in 2:(K-1))
  {
    y = MO.weekly.deaths.mat1[i,k]
    y.left = MO.weekly.deaths.mat1[i,k-1]
    y.right = MO.weekly.deaths.mat1[i,k+1]
    if( y < 0 & y.left > abs(y)/2 & y.right > abs(y)/2 ) 
    {
      MO.weekly.deaths.mat1[i,k-1] = y.left - floor(0.6+abs(y/2))
      MO.weekly.deaths.mat1[i,k]   = 0
      MO.weekly.deaths.mat1[i,k+1] = y.right - floor(abs(y/2))
    }
  }
}

matxMax <- function(mtx)
{
  colmn <- which(mtx == max(mtx)) %/% nrow(mtx) + 1
  row <- which(mtx == max(mtx)) %% nrow(mtx)
  return( matrix(c(row, colmn), 1))
}
matxMax( abs(MO.weekly.deaths.mat1) )
##  Column 146 is very wierd

MO.weekly.deaths.mat1 = MO.weekly.deaths.mat1[, c(1:145,147:K)]
str( MO.weekly.deaths.mat1 )

sum( MO.weekly.deaths.mat1 < 0 )






windows( 16 , 14 )
par( mfrow=c(5,3) )
start=1
for( i in start:(start+14) )
ts.plot( MO.weekly.deaths.mat[i,] , main=paste0("County ",round(i)) )

start = 341
for (j in start:(start+19) )
{
  if( sum( MO.daily.deaths.mat[,j] < 0 ) )  print( c( j , which( MO.daily.deaths.mat[,j] < 0 ) ) )
}
  








write.csv( MO.daily.deaths.mat , "MO_COVID_Deaths_Daily.csv")



