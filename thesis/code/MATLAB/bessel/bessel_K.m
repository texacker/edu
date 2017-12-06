function K = bessel_K(x,n)
%
% function for Bessel function
%
if (angle(x)>-pi)&(angle(x)<=(pi/2))
   H=bessel_H1(x*exp(pi*(1i)/2),n);
   K=pi*(1i)/2*exp(n*pi*(1i)/2)*H;
else
   H=bessel_H2(x*exp(-pi*(1i)/2),n);
   K=-pi*(1i)/2*exp(-n*pi*(1i)/2)*H;
end
return