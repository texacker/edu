function Y = bessel_Y(x,n)
%
% function for Bessel function
%
J=bessel_J(x,n);
Y=(cos(n*pi)-(-1)^n)*J/sin(n*pi);
return