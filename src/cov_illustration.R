personal = data.frame(
  age = c( 45, 22, 56, 41, 66, 33, 39, 59, 60, 48 ) ,
  height = c( 70, 68, 72, 75, 69, 71, 66, 73, 70, 68 ),
  weight = c( 185, 150, 210, 268, 192, 183, 265, 221, 165, 138 )
)

company = data.frame(
  years = c(12, 1, 28, 15, 40, 8, 12, 4, 5, 21 ),
  salary = c(55, 41, 53, 29, 67, 55, 70, 42, 64, 84)
)

cov( personal )
cor( personal )

cov( company )
cor( company )

cov( personal , company )
cor( personal , company )

####
####  Lag 1 correlation
####
set.seed(1)
x = rnorm( 20 , 100 , 1 )
x0 = x[-1]
x1 = x[-length(x)]
cor( x0 , x1 )

cbind( x0 , x1 )
