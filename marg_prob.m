function [f,f1,f2]=marg_prob(y)
sigmax=.1;mu=0;h=1;b=1;
f=quadgk(@(x) S1(y,x,h,b,mu,sigmax),-inf,inf);
sigma=1;
f1=quadgk(@(x) S2(y,x,h,sigma,sigmax),0,inf);
sigmax=.1;sigma1=1;sigma2=.1; mu1=1;mu2=-1;
f2=quadgk(@(x) S3(y,x,h,sigma1,sigma2,mu1,mu2,sigmax),-inf,inf);

function f=S1(y,x,h,b,mu,sigmax)
f=1/2/b*exp(-abs(y-h*x-mu)/b).*exp(-x.^2/2/sigmax^2);

function f=S2(y,x,h,sigma,sigmax)
%x<0 f=0;
f=1/2/sigma^2*exp(-(y-h*x).^2/2/sigma^2).*exp(-x.^2/2/sigmax^2);
        
function f=S3(y,x,h,sigma1,sigma2,mu1,mu2,sigmax)   
epsilon=.1;
f=(1-epsilon)*1/sqrt(2*pi)/sigma1/sigmax*exp(-(y-h*x-mu1).^2/2/sigma1^2).*exp(-x.^2/2/sigmax^2)...
        +epsilon/sqrt(2*pi)/sigmax/sigma2*exp(-(y-h*x-mu2).^2/2/sigma2^2).*exp(-x.^2/2/sigmax^2);
 