function I = bessel_I(x,n)
%
% function for Bessel function
%
if (angle(x)>-pi)&(angle(x)<=(pi/2))
   J=bessel_J(x*exp(pi*(1i)/2),n);
   I=exp(-n*pi*(1i)/2)*J;
else
   J=bessel_J(x*exp(-3*pi*(1i)/2),n);
   I=exp(3*n*pi*(1i)/2)*J;
end
return