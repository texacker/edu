function gen_d(q)

thesis_data;

De = zeros(3,3);

r = zeros(size(Gauss));
z = zeros(size(Gauss));
tmp = zeros(size(Gauss));

for i=1:3
   if e_type(i,q)==1
      for j=1:3
         for k=1:3
            r = (Gauss*(nx(i,q)-tx(i,q))+(nx(i,q)+tx(i,q)))/2;
            z = (Gauss*(ny(i,q)-ty(i,q))+(ny(i,q)+ty(i,q)))/2;
            if all(r)
               tmp=(ta(j,q)+tb(j,q)*r+tc(j,q)*z).*...
                  (ta(k,q)+tb(k,q)*r+tc(k,q)*z)./...
                  sqrt(1+(G0(ones(size(Gauss)))/...
                  (r*w0^2*Fea_Len)).^2);
               De(j,k)=De(j,k)+sum(tmp)*len(i,q)/...
                  (8*ts(q)^2)*d_init_matrix(j,k,i);
            end
         end
      end
   end
end