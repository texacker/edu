function gen_ac(q)

thesis_data;

Ae = zeros(3,3);
Ce = zeros(3,3);
Ae1 = zeros(3,3);
Ae2 = zeros(3,3);

r = zeros(2,size(Gauss,2));
z1 = zeros(2,size(Gauss,2));
z2 = zeros(2,size(Gauss,2));
zz = [1 3];

for j=1:2
   i=zz(j);
   r(j,:)=(Gauss*(sx(2,q)-sx(i,q))+(sx(i,q)+sx(2,q)))/2;
   z1(j,:)=(Gauss+1)*(sy(3,q)-sy(1,q))*(sx(2,q)-sx(i,q))/(sx(3,q)-sx(1,q))/2+sy(i,q);
   z2(j,:)=(Gauss*(sy(2,q)-sy(i,q))+(sy(i,q)+sy(2,q)))/2;
end

for i=1:3
   for j=1:3
      Ce(i,j)  = tc(i,q)*tc(j,q)*xc(q)/4/ts(q);
      Ae1(i,j) = tb(i,q)*tb(j,q)*xc(q)/4/ts(q);
      
      tmp = 0;
      for k=1:2
         for m=1:size(Gauss,2)
            eA=(ta(i,q)+tb(i,q)*r(k,m))*(ta(j,q)+tb(j,q)*r(k,m));
            eB=((ta(i,q)+tb(i,q)*r(k,m))*tc(j,q)+(ta(j,q)+tb(j,q)*r(k,m))*tc(i,q))/2;
            eC=(tc(i,q)*tc(j,q))/3;
            if r(k,m)>Prec
               tmp=tmp+(eA*abs(z2(k,m)-z1(k,m))+...
                  eB*abs(z2(k,m)^2-z1(k,m)^2)+...
                  eC*abs(z2(k,m)^3-z1(k,m)^3))*...
                  abs(sx(2,q)-sx(zz(k),q))/r(k,m);
            end
         end
      end
      Ae2(i,j) = tmp*M0^2/(8*ts(q)^2);
   end
end
Ae = Ae1 + Ae2 + Ce;