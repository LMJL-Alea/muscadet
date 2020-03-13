mfunbessel=function(r,alpha){2*besselJ(2*r/alpha,1)/(2*r/alpha)}
mfunspecbessel=function(r,alpha){ifelse(pi^2*alpha^2*r^2<1,pi*alpha^2,0)}
gbessel=function(r,alpha,nu=10){ 1-mfunbessel(r,alpha)^2}
gbessel12=function(r,tau,alpha12,nu=10){1-tau^2*mfunbessel(r,alpha12)^2}
Kspecbessel=function(r,rho1=100,rho2=100,alpha1=0.0044,alpha2=0.0044,alpha12=0.0044,tau=0){
  matrix(c(rho1*mfunspecbessel(r,alpha1),tau*sqrt(rho1*rho2)*mfunspecbessel(r,alpha12),tau*sqrt(rho1*rho2)*mfunspecbessel(r,alpha12),rho2*mfunspecbessel(r,alpha2)),2,2)
}
testtaubessel=function(tau,rho1,rho2,alpha1,alpha2,alpha12,nu1=10,nu2=10,nu12=10){
  bound1<-alpha1*alpha2/alpha12^2
  bound2<-sqrt(4*(alpha1*alpha2/alpha12^2)^2*(1/(rho1*pi*alpha1^2)-1)*(1/(rho2*pi*alpha2^2)-1))
  if(tau>bound1 || tau>bound2 || rho1*pi*alpha1^2>1 || rho2*pi*alpha2^2>1 || alpha12<alpha1 || alpha12<alpha2 ){return(FALSE)}
  else{return(TRUE)}}
validtaubessel=function(rho1=100,rho2=100,alpha1=0.03,alpha2=0.03,alpha12=0.03){
  bound1<-alpha1*alpha2/alpha12^2
  bound2<-sqrt(4*(alpha1*alpha2/alpha12^2)^2*(1/(rho1*pi*alpha1^2)-1)*(1/(rho2*pi*alpha2^2)-1))
  return(paste("tau must be smaller than ", bound1," and ",bound2 ))
}
