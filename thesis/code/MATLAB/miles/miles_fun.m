function [miles,dd] = miles_fun(a,b,h,m,n,mu)

miles=zeros(1,size(mu,2));
d = zeros(2);

if m==1
   for i=1:size(miles,2)
      t = mu(1,i);
      if t^2 > 1%abs(t) > 1
         l = sqrt(t^2-1)*n*pi/h;
         d(1,1)=l*a*besselj(0,l*a)+(t-1)*besselj(1,l*a);
         d(1,2)=l*a*bessely(0,l*a)+(t-1)*bessely(1,l*a);
         d(2,1)=t^2*l*b*besselj(0,l*b)+(4-5*t^2+t^3)*besselj(1,l*b);
         d(2,2)=t^2*l*b*bessely(0,l*b)+(4-5*t^2+t^3)*bessely(1,l*b);
      else
         k = sqrt(1-t^2)*n*pi/h;
         d(1,1)= k*a*besseli(0,k*a)+(t-1)*besseli(1,k*a);
         d(1,2)=-k*a*besselk(0,k*a)+(t-1)*besselk(1,k*a);
         d(2,1)= t^2*k*b*besseli(0,k*b)+(4-5*t^2+t^3)*besseli(1,k*b);
         d(2,2)=-t^2*k*b*besselk(0,k*b)+(4-5*t^2+t^3)*besselk(1,k*b);
      end
      miles(i) = det(d);
   end
else
   for i=1:size(miles,2)
      t = mu(1,i);
      if t^2 > 1%abs(t) > 1
         l = sqrt(t^2-1)*n*pi/h;
         d(1,1)=l*a*besselj(0,l*a)+(t-1)*besselj(1,l*a);
         d(1,2)=l*a*bessely(0,l*a)+(t-1)*bessely(1,l*a);
         d(2,1)=t^2*l*b*besselj(0,l*b)+(4-5*t^2+t^3)*besselj(1,l*b);
         d(2,2)=t^2*l*b*bessely(0,l*b)+(4-5*t^2+t^3)*bessely(1,l*b);
      else
         k = sqrt(1-t^2)*n*pi/h;
         d(1,1)= k*a*besseli(0,k*a)+(t-1)*besseli(1,k*a);
         d(1,2)=-k*a*besselk(0,k*a)+(t-1)*besselk(1,k*a);
         d(2,1)= t^2*k*b*besseli(0,k*b)+(4-5*t^2+t^3)*besseli(1,k*b);
         d(2,2)=-t^2*k*b*besselk(0,k*b)+(4-5*t^2+t^3)*besselk(1,k*b);
      end
      miles(i) = i;%det(d);
   end
end

if nargout==2
   dd = d;
end

return