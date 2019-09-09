ggauss=function(r,alpha,nu=10){ 1-exp(-2*r^2/alpha^2)}
ggauss12=function(r,tau,alpha12,nu=10){1-tau^2*exp(-2*r^2/alpha12^2)}
Kspecgauss=function(r,rho1=100,rho2=100,alpha1=0.03,alpha2=0.03,alpha12=0.03,tau=0){
  matrix(c(rho1*alpha1^2*pi*exp(-pi^2*alpha1^2*r^2),tau*sqrt(rho1*rho2)*alpha12^2*pi*exp(-pi^2*alpha12^2*r^2),tau*sqrt(rho1*rho2)*alpha12^2*pi*exp(-pi^2*alpha12^2*r^2),rho2*alpha2^2*pi*exp(-pi^2*alpha2^2*r^2)),2,2)
}
validtaugauss=function(rho1=100,rho2=100,alpha1=0.03,alpha2=0.03,alpha12=0.03){
  bound1<-alpha1*alpha2/alpha12^2
  bound2<-sqrt(4*(alpha1*alpha2/alpha12^2)^2*(1/(pi*rho1*alpha1^2)-1)*(1/(pi*rho2*alpha2^2)-1))
  return(paste("tau must be smaller than ", bound1," and ",bound2 ))
}
testtaugauss=function(tau,rho1,rho2,alpha1,alpha2,alpha12){
  bound1<-alpha1*alpha2/alpha12^2
  bound2<-sqrt(4*(alpha1*alpha2/alpha12^2)^2*(1/(pi*rho1*alpha1^2)-1)*(1/(pi*rho2*alpha2^2)-1))
  if(tau>bound1 || tau>bound2 || pi*rho1*alpha1^2>1 || pi*rho2*alpha2^2>1 || alpha12^2<(alpha1^2+alpha2^2)/2){return(FALSE)}
  else{return(TRUE)}}
