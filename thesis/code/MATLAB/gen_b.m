function gen_b(q)

thesis_data;

Be = zeros(3,3);

for i=1:3
   if boundary(i,q)==1
      Be = Be + b_init_matrix(:,:,i) / 6 * nr(i,q) * len(i,q);
   end
end