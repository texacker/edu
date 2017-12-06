function H1 = bessel_H1(x,n)
%
% function for Bessel function
%
H1=bessel_J(x,n)+bessel_Y(x,n)*(1i);
return