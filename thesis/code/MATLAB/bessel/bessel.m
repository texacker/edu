function J = bessel(x,n)
%
% function for Bessel function
%
J=zeros(1,4);
J(1,1)=quad8('fun_bessel',-pi,pi,[],[],n,x)/2/pi;
J(1,2)=(cos(n*pi)-(-1)^n)*J(1,1)/sin(n*pi);
if (angle(x)>-pi)&(angle(x)<=(pi/2))
   J(1,3)=0;
   J(1,4)=0;
else
   J(1,3)=1;
   J(1,4)=1;
end
return