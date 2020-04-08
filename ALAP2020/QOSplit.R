## Michalski AI, Grigoriev P, Gorlishchev VP  ##

## MPIDR Technical Report ##

## 22.11.2017 ##

## Splitting Fx (age-specific fertility rates) 
## using quadratic optimization (QO) method  (HFD data) ##


## Arguments of 'QOSplit' function ## 
## Fx - vector of aggregated ASFRs to be split into single ages ##
## L - vector of ages  (lower limits. e.g 10-15 is 10, 20-24 is 20, and so on) ##
## AgeInt - vector of age intervals ##
QOSplit<-function(Fx,L,AgeInt){
  
  ## (0) Closing left (-99) and right (99) open intervals, if any ##
  n<-length(AgeInt)
  if (AgeInt[1]==-99) {L[1]<-L[1]-AgeInt[2]+1;AgeInt[1]<-AgeInt[2] } 
  if (AgeInt[n]==99) {AgeInt[n]<-AgeInt[n-1]}
  
  h<-AgeInt			# number of 1-years age groups in the 5-years age groups
  N<-sum(h)			# numberof columns  in G
  
  ## (1) Getting G matrix
  G<-matrix(0,nrow=n,ncol=N)
  G[1,1:h[1]]<-1/h[1]			
  j1<-h[1]+1				
   for (i in 2:n){ 
    j2<-j1+h[i]-1	
    G[i,j1:j2]<-1/h[i]			
    j1<-j2+1  } 
  
  # (2) Quadratic minimization
    
  # Matrix for the 2-nd order derivatives
  F<-matrix(0,nrow=N-2,ncol=N)
  for (i in 1:(N-2)){  F[i,i]<-1; F[i,i+1]<--2; F[i,i+2]<-1;}
  G<-G[,-N]; G<-G[,-1];
  F<-F[,-N]; F<-F[,-1];
  Dmat=t(F)%*%F			# matrix for Quadratic form	
  dvec<-rep(0,N-2)
  Amat<-t(rbind(G,diag(N-2)))	# constrains matrix
  meq<-n					# number of equations
  
  bvec<-c(Fx,dvec)		# vector of birth
  aux<-solve.QP(Dmat, dvec, Amat, bvec, meq=meq)
  aux<-c(0,aux$solution,0)	#  the first and the last values are assumed to have no fertility!
  aux[aux<0]<-0   ## extremely small negative values can be neglected and treated as '0's 
  
  solutionN<-round(aux,8)
  #OUT<-rbind(OUT,solutionN)
  Age<-L[1]+c(1:sum(AgeInt))-1
  out<-as.data.frame(cbind(Age,solutionN))
  names(out)<-c("Age","ASFR")
  return(out)  #uncomment this to see the output 
  #plot(Age ,solutionN, type='l',lwd=2,col='red',xlab='Age',ylab="ASFR",main="") 
}

######################################## END ################################################

