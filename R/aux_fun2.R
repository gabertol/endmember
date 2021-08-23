#' auxilary code
#' @export

r2=function(x,y)

{
  p<-ncol(x)
  R2<-vector()

  for(l in c(1:p)){

    R2[l]<-(sd(x[,l])^2-sd(x[,l]-y[,l])^2)/sd(x[,l])^2

  }
  list(R2=R2)

}
