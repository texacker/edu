function H2 = bessel_H2(x,n)
%
% function for Bessel function
%
H2=bessel_J(x,n)-bessel_Y(x,n)*(1i);
%H2=conj(bessel_H1(x,n));
return