function miles = miles_fun_s(a,b,h,n,mu)

miles=zeros(1,size(mu,2));

for i=1:size(miles,2)
   d = zeros(2);
   m = mu(1,i);
   k = sqrt(1-m^2)*n*pi/h;
   l = sqrt(m^2-1)*n*pi/h;
   d(1,1) = l*a*besselj(0,l*a)+(m-1)*besselj(1,l*a);
   d(1,2) = l*a*bessely(0,l*a)+(m-1)*bessely(1,l*a);
   d(2,1) = m^2*l*b*besselj(0,l*b)+(4-5*m^2+m^3)*besselj(1,l*b);
   d(2,2) = m^2*l*b*bessely(0,l*b)+(4-5*m^2+m^3)*bessely(1,l*b);
   miles(i) = det(d);
end

return