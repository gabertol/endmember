#' T-matrix calculation see Weltje, 1997 Appendix C, k=matrix M, c=weighting Factor
#' @export

tcalc=function(k,c)
{
  if(min(k)>=0){print("STOP, no negative entries in M")}
  if(min(k)>=0)break
  #c5<-as.numeric(readline(paste("please enter weighting exponent c5: "))) #should be activated when used as single function

  c5<-c

  n<-nrow(k)
  p<-ncol(k)

  T<-diag(rep(1,each=p))

  #calculate expansion factor; see Eq.38 (Weltje, 1997)

  prop=vector()

  for(s in c(1:p)){

    prop[s]<-length(k[,s][k[,s]<0])/length(k[,s])}

  mk=vector()

  for(s in c(1:p)){
    mk[s]<-mean(k[,s][k[,s]<0])}

  for(s in c(1:p)){

    if(is.nan(mk[s])==TRUE){mk[s]<-0}}


  pol_exp<-mk*prop^c5           #pol_exp=vector()

  print("mode of polytope expansion of the end member")
  print(pol_exp)

  #----------------------------------------------------------------------

  for(j in c(1:p)){       #calculate T-matrix (Appendix C in Weltje 1997)

    z=vector()
    y<-matrix(0,nrow=n,ncol=p)
    Y=vector

    for(i in c(1:n)){
      if(k[i,j]<0){y[i,]<-k[i,]}}				#scan columns for negative values

    Y<-colSums(y)/length(k[,j][k[,j]<0])

    if(length(k[,j][k[,j]<0])==0){Y[1:p]<-0}
    if(length(Y[Y==0])!=p) z<-which.max(Y)
    if(length(Y[Y==0])!=p){T[z,j]<-pol_exp[j]}	#replace entries in T-Matrix; Weltje 1997 (Appendix C)

    #print("gescannte Spalte von A:")
    #print(j)
    #print("current correction vector Y")
    #print(Y)
    #print("current T")
    #print(T)
  }

  for(s in c(1:p)){T[s,s]<-T[s,s]-sum(T[s,])+1}

  list(T=T)
}
