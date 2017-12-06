function miles = miles_fun_v(a,b,h,n,mu)

miles=zeros(1,size(mu,2));

for i=1:size(miles,2)
   d = zeros(2);
   m = mu(1,i);
   k = sqrt(1-m^2)*n*pi/h;
   l = sqrt(m^2-1)*n*pi/h;
   d(1,1) =  k*a*besseli(0,k*a)+(m-1)*besseli(1,k*a);
   d(1,2) = -k*a*besselk(0,k*a)+(m-1)*besselk(1,k*a);
   d(2,1) =  m^2*k*b*besseli(0,k*b)+(4-5*m^2+m^3)*besseli(1,k*b);
   d(2,2) = -m^2*k*b*besselk(0,k*b)+(4-5*m^2+m^3)*besselk(1,k*b);
   miles(i) = det(d);
end

return