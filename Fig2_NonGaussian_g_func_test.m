% g=d/dy(log(marg_prob(y))
clear all
y=[-5:.01:5];

for i=1:length(y)
 [f1(i),f2(i),f3(i)]=marg_prob(y(i));
end
plot(y,f1);
hold on
plot(y,f2,'r');
plot(y,f3,'c');

figure
g=diff(log(f1));
g1=diff(log(f2));
g2=diff(log(f3));

y=y(2:end);
plot(y,g);
hold on
plot(y,g1,'r');
plot(y,g2,'c');
legend('Laplacian','Rayleigh','Impulsive')


% distributions
y=[-5:.01:5];
sigmax=1;mu=0;h=1;b=1;x=0;
f1=1/2/b*exp(-abs(y-h*x-mu)/b);
sigma=1;

%x<0 f=0;
f2=1/2/sigma^2*exp(-(y-h*x).^2/2/sigma^2);      
mid=floor(length(f2)/2);
f2(1:mid)=f2(1:mid)*0;

sigma1=1;sigma2=.1; mu1=1;mu2=-1;
epsilon=.1;
f3=(1-epsilon)*1/sqrt(2*pi)/sigma1*exp(-(y-h*x-mu1).^2/2/sigma1^2)...
        +epsilon*1/sqrt(2*pi)/sigma2*exp(-(y-h*x-mu2).^2/2/sigma2^2);
figure;
plot(y,f1);
hold on
plot(y,f2,'r');
plot(y,f3,'c');
legend('Laplacian','Rayleigh','Impulsive')
title('Distribution')


c=0;
sigma_x=.1;mu=0;h=1;
sigma_1=1;sigma_2=.1; mu1=1;mu2=-1;
for y=[-5:.01:5]
    c=c+1;
lambda1=sqrt(sigma_1^2+sigma_x^2*h^2);
kappa1=(1-h^2*sigma_x^2/lambda1^2);
a1(c)=1/(lambda1)*exp(-kappa1*(y-mu1).^2/sigma_1^2/2);
da1=-(y-mu1)*kappa1/lambda1/sigma_1^2.*exp(-kappa1*(y-mu1).^2/sigma_1^2/2);
lambda2=sqrt(sigma_2^2+sigma_x^2*h^2);
kappa2=(1-h^2*sigma_x^2/lambda2);
a2(c)=1/(lambda2)*exp(-kappa2*(y-mu2).^2/sigma_2^2/2);
da2=-(y-mu2)*kappa2/lambda2/sigma_2^2.*exp(-kappa2*(y-mu2).^2/sigma_2^2/2);
f(c)=(da1+da2)/(a1(c)+a2(c));
f1(c)=a1(c)+a2(c);
end
figure
plot([-5:.01:5],f,'k')
figure
plot([-5:.01:5],f1,'k')

