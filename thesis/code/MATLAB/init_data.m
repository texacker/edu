function init_data
thesis_data

tx = zeros(3,size(t_t, 2));
ty = zeros(3,size(t_t, 2));
nx = zeros(3,size(t_t, 2));
ny = zeros(3,size(t_t, 2));
dx = zeros(3,size(t_t, 2));
dy = zeros(3,size(t_t, 2));
sx = zeros(3,size(t_t, 2));
sy = zeros(3,size(t_t, 2));
ta = zeros(3,size(t_t, 2));
tb = zeros(3,size(t_t, 2));
tc = zeros(3,size(t_t, 2));
len= zeros(3,size(t_t, 2));
nr = zeros(3,size(t_t, 2));
e_type = zeros(3,size(t_t, 2));
boundary = zeros(3,size(t_t, 2));
free_point = zeros(1,size(t_p,2));
free_edge = zeros(1,size(t_e,2));
distance = zeros(1,size(t_p,2));

ts = zeros(1,size(t_t, 2));
xc = zeros(1,size(t_t, 2));
yc = zeros(1,size(t_t, 2));

if(SpheroidDomain)
   distance = abs(t_p(1,:)-alpha);
else
   if abs(r2-r1)>abs(h2-h1)
      distance = abs(((t_p(1,:).^2-r1^2)*w0*w0/2/G0+h1)-t_p(2,:));
   else
      distance = abs(sqrt((t_p(2,:)-h1)*2*G0/w0/w0+r1^2)-t_p(1,:));
   end
end

free_point(find(distance < Prec)) = 1;

sort_flag = 1;
while sort_flag == 1
   sort_flag = 0;
   for i=1:size(free_point, 2)-1
      if free_point(i) < free_point(i+1)
         p_tmp = t_p(:,i);
         t_p(:,i) = t_p(:,i+1);
         t_p(:,i+1) = p_tmp;
         for j=1:2
            e_tmp1=find(t_e(j,:)==i);
            e_tmp2=find(t_e(j,:)==i+1);
            t_e(j,e_tmp1) = (i+1)*ones(size(e_tmp1));
            t_e(j,e_tmp2) = i*ones(size(e_tmp2));
         end
         for j=1:3
            t_tmp1=find(t_t(j,:)==i);
            t_tmp2=find(t_t(j,:)==i+1);
            t_t(j,t_tmp1) = (i+1)*ones(size(t_tmp1));
            t_t(j,t_tmp2) = i*ones(size(t_tmp2));
         end
         tmp = free_point(i);
         free_point(i) = free_point(i+1);
         free_point(i+1)= tmp;
         sort_flag = 1;
      end
   end
end

for i=1:size(t_t, 2)
   for j=1:3
      tx(j,i) = t_p(1,t_t(j,i));
      ty(j,i) = t_p(2,t_t(j,i));
   end
end

nx(1,:) = tx(2,:);
nx(2,:) = tx(3,:);
nx(3,:) = tx(1,:);

ny(1,:) = ty(2,:);
ny(2,:) = ty(3,:);
ny(3,:) = ty(1,:);

dx = nx-tx;
dy = ny-ty;

sx = tx;
sy = ty;

len = sqrt(dx.^2+dy.^2);
nr = dy./len;

xc = mean(tx,1);
yc = mean(ty,1);

ta(1,:) = tx(2,:).*ty(3,:)-tx(3,:).*ty(2,:);
ta(2,:) = tx(3,:).*ty(1,:)-tx(1,:).*ty(3,:);
ta(3,:) = tx(1,:).*ty(2,:)-tx(2,:).*ty(1,:);

tb(1,:) = ty(2,:)-ty(3,:);
tb(2,:) = ty(3,:)-ty(1,:);
tb(3,:) = ty(1,:)-ty(2,:);

tc(1,:) = tx(3,:)-tx(2,:);
tc(2,:) = tx(1,:)-tx(3,:);
tc(3,:) = tx(2,:)-tx(1,:);

ts = sum(ta,1)/2;

for i=1:size(t_t, 2)
   for j=1:3
      for k=1:2
         if sx(k,i) > sx(k+1,i)
            swap = sx(k,i);
            sx(k,i) = sx(k+1,i);
            sx(k+1,i) = swap;
            
            swap = sy(k,i);
            sy(k,i) = sy(k+1,i);
            sy(k+1,i) = swap;
         end
      end
   end
end

if(SpheroidDomain)
   distance = abs(t_p(1,:)-alpha);
else
   if abs(r2-r1)>abs(h2-h1)
      distance = abs(((t_p(1,:).^2-r1^2)*w0*w0/2/G0+h1)-t_p(2,:));
   else
      distance = abs(sqrt((t_p(2,:)-h1)*2*G0/w0/w0+r1^2)-t_p(1,:));
   end
end

free_point(find(distance < Prec)) = 1;
free_edge(find((free_point(t_e(1,:))+free_point(t_e(2,:))) == 2)) = 1;

for i=1:size(t_t,2)
   for j=1:3
      start_point = find(t_e(1,:)==t_t(j,i));
      end_point = find(t_e(2,:)==t_t(rem(j,3)+1,i));
      if ~(isempty(start_point)|isempty(end_point))
         if start_point == end_point
            e_type(j,i)=free_edge(start_point);
            boundary(j,i) = 1;
         end
      end
      
      start_point = find(t_e(2,:)==t_t(j,i));
      end_point = find(t_e(1,:)==t_t(rem(j,3)+1,i));
      if ~(isempty(start_point)|isempty(end_point))
         if start_point == end_point
            e_type(j,i)=free_edge(start_point);
            boundary(j,i) = 1;
         end
      end
      
   end
end

disp(['Total Points...... ' num2str(size(t_p,2))]);
disp(['Total Borders..... ' num2str(size(t_e,2))]);
disp(['Total Triangles... ' num2str(size(t_t,2))]);

return