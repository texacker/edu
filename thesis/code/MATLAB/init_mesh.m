function init_mesh

thesis_data;

W0 = sqrt(2*G0*H/R/R);

if(SpheroidDomain)
   
   r1 = acos(alpha);
   h1 = 0;
   r2 = 0;
   h2 = 0;
   
   t_g = zeros(12,nodes);
   t_x = R*cos(linspace(-r1,r1,nodes));
   t_y = H*(1+sin(linspace(-r1,r1,nodes)));
   
   t_g(1,:) = 2; t_g(6,:) = 1; t_g(7,:) = 0;
   t_g(2,:) = t_x; t_g(4,:) = t_y;
   t_g(3,1:nodes-1) = t_x(2:nodes); t_g(3,nodes) = t_x(1);
   t_g(5,1:nodes-1) = t_y(2:nodes); t_g(5,nodes) = t_y(1);
   %plot(t_x,t_y);pause;

else
   
if (W0 > w0)
   a0 = (w0*R).^2/(4*G0*H);
   disp('W0 > w0');

   if (alpha < a0)
      disp('alpha < a0');
      
      r1 = sqrt(R*R-sqrt(4*alpha*R*R*H*G0/w0/w0));
      h1 = 0;
      r2 = R;
      h2 = sqrt(alpha*w0*w0*R*R*H/G0);
      
      t_g = zeros(12, nodes+2);
      t_x = linspace(r2,r1,nodes+1); t_x(nodes+2)=r2;
      t_y = (t_x.^2-r1^2)*w0*w0/2/G0; t_y(nodes+2) = h1;
      t_g(1,:) = 2; t_g(6,:) = 1; t_g(7,:) = 0;
      t_g(2,:) = t_x; t_g(4,:) = t_y;
      t_g(3,1:nodes+1) = t_x(2:nodes+2); t_g(3,nodes+2) = t_x(1);
      t_g(5,1:nodes+1) = t_y(2:nodes+2); t_g(5,nodes+2) = t_y(1);
      
      %plot(t_x, t_y);
      
   elseif (alpha > (1-a0))
      disp('1 - a0 < alpha');
      
      r1 = 0;
      h2 = H;
      r2 = sqrt(4*G0*R*R*H*(1-alpha)/w0/w0);
      h1 = H-w0*w0*r2/2/G0;
      r2 = sqrt(r2);
      
      t_g = zeros(12, nodes+4);
      t_x = linspace(r2,r1,nodes+1); t_x(nodes+2)=r1;t_x(nodes+3)=R;t_x(nodes+4)=R;
      t_y = t_x.^2*w0*w0/2/G0+h1; t_y(nodes+2) = 0;t_y(nodes+3) = 0;t_y(nodes+4) = h2;
      t_g(1,:) = 2; t_g(6,:) = 1; t_g(7,:) = 0;
      t_g(2,:) = t_x; t_g(4,:) = t_y;
      t_g(3,1:nodes+3) = t_x(2:nodes+4); t_g(3,nodes+4) = t_x(1);
      t_g(5,1:nodes+3) = t_y(2:nodes+4); t_g(5,nodes+4) = t_y(1);
      
      %plot(t_x, t_y);
      
   else
      disp('a0 < alpha < 1 - a0');
      r1 = 0;
      h1 = alpha*H-w0*w0*R*R/4/G0;
      r2 = R;
      h2 = h1+w0*w0*R*R/2/G0;
      
      t_g = zeros(12, nodes+3);
      t_x = linspace(r2,r1,nodes+1); t_x(nodes+2)=r1;t_x(nodes+3)=r2;
      t_y = t_x.^2*w0*w0/2/G0+h1; t_y(nodes+2) = 0;t_y(nodes+3) = 0;
      t_g(1,:) = 2; t_g(6,:) = 1; t_g(7,:) = 0;
      t_g(2,:) = t_x; t_g(4,:) = t_y;
      t_g(3,1:nodes+2) = t_x(2:nodes+3); t_g(3,nodes+3) = t_x(1);
      t_g(5,1:nodes+2) = t_y(2:nodes+3); t_g(5,nodes+3) = t_y(1);
      
      %plot(t_x, t_y);
      
   end
   
else
   
   a0 = G0*H/(w0*R).^2;
   disp('W0 <= w0');
   
   if (alpha < a0)
      disp('alpha < a0');
      
      r1 = sqrt(R*R-sqrt(4*alpha*R*R*H*G0/w0/w0));
      h1 = 0;
      r2 = R;
      h2 = sqrt(alpha*w0*w0*R*R*H/G0);
      
      t_g = zeros(12, nodes+2);
      t_x = linspace(r2,r1,nodes+1); t_x(nodes+2)=r2;
      t_y = (t_x.^2-r1^2)*w0*w0/2/G0; t_y(nodes+2) = h1;
      
      t_g(1,:) = 2; t_g(6,:) = 1; t_g(7,:) = 0;
      t_g(2,:) = t_x; t_g(4,:) = t_y;
      t_g(3,1:nodes+1) = t_x(2:nodes+2); t_g(3,nodes+2) = t_x(1);
      t_g(5,1:nodes+1) = t_y(2:nodes+2); t_g(5,nodes+2) = t_y(1);
      
      %plot(t_x, t_y);
      
   elseif (alpha > (1-a0))
      disp('1 - a0 < alpha');
      
      r1 = 0;
      h2 = H;
      r2 = sqrt(4*G0*R*R*H*(1-alpha)/w0/w0);
      h1 = H-w0*w0*r2/2/G0;
      r2 = sqrt(r2);
      
      t_g = zeros(12, nodes+4);
      t_x = linspace(r2,r1,nodes+1); t_x(nodes+2)=r1;t_x(nodes+3)=R;t_x(nodes+4)=R;
      t_y = t_x.^2*w0*w0/2/G0+h1; t_y(nodes+2) = 0;t_y(nodes+3) = 0;t_y(nodes+4) = h2;
      t_g(1,:) = 2; t_g(6,:) = 1; t_g(7,:) = 0;
      t_g(2,:) = t_x; t_g(4,:) = t_y;
      t_g(3,1:nodes+3) = t_x(2:nodes+4); t_g(3,nodes+4) = t_x(1);
      t_g(5,1:nodes+3) = t_y(2:nodes+4); t_g(5,nodes+4) = t_y(1);
      
      %plot(t_x, t_y);
      
   else
      disp('a0 < alpha < 1 - a0');
      r1 = R*R*(1-alpha)-G0*H/w0/w0;
      r2 = sqrt(r1+2*G0*H/w0/w0);
      r1 = sqrt(r1);
      h1 = 0;
      h2 = H;
      
      t_g = zeros(12, nodes+3);
      t_y = linspace(h2,h1,nodes+1);t_y(nodes+2)=h1;t_y(nodes+3)=h2;
      t_x = sqrt(t_y*2*G0/w0/w0+r1^2);t_x(nodes+2)=R;t_x(nodes+3)=R;
                  
      t_g(1,:) = 2; t_g(6,:) = 1; t_g(7,:) = 0;
      t_g(2,:) = t_x; t_g(4,:) = t_y;
      t_g(3,1:nodes+2) = t_x(2:nodes+3); t_g(3,nodes+3) = t_x(1);
      t_g(5,1:nodes+2) = t_y(2:nodes+3); t_g(5,nodes+3) = t_y(1);
      
      %plot(t_x, t_y);
   end
end
end

if(RegularMesh)
   %r1 = R*R*(1-alpha)-G0*H/w0/w0;
   %r2 = sqrt(r1+2*G0*H/w0/w0);
   %r1 = sqrt(r1);
   %h1 = 0;
   %h2 = H;
   startx = R*sqrt(1-alpha);
   inc_x = (R-startx)/(nodes_x-1);
   inc_y = H/(nodes_y-1);
   t_p = zeros(2,nodes_x*nodes_y);
   t_e = zeros(7,2*((nodes_x-1)+(nodes_y-1)));
   t_t = zeros(4,2*(nodes_x-1)*(nodes_y-1));
   
   for I = 1:nodes_x,
      for J = 1:nodes_y,
         ct = (I-1)*nodes_y+J;
         t_p(1,ct) = (I-1)*inc_x + startx;
         t_p(2,ct) = (J-1)*inc_y;
      end
   end
   
   for J = 1:nodes_y-1,
      t_e(1,J) = nodes_y-J+1;
      t_e(2,J) = nodes_y-J;
      t_e(5,J) = 1;
   end
   for I = 1:nodes_x-1,
      ct = (nodes_y-1)+I;
      t_e(1,ct) = (I-1)*nodes_y+1;
      t_e(2,ct) = I*nodes_y+1;
      t_e(5,ct) = 2;
   end
   for J = 1:nodes_y-1,
      ct = (nodes_x-1)+(nodes_y-1)+J;
      t_e(1,ct) = (nodes_x-1)*nodes_y+J;
      t_e(2,ct) = (nodes_x-1)*nodes_y+J+1;
      t_e(5,ct) = 3;
   end
   for I = 1:nodes_x-1,
      ct = (nodes_x-1)+(nodes_y-1)+(nodes_y-1)+I;
      t_e(1,ct) = (nodes_x-I+1)*nodes_y;
      t_e(2,ct) = (nodes_x-I)*nodes_y;
      t_e(5,ct) = 4;
   end
   t_e(3,:) = 0;
   t_e(4,:) = 1;
   t_e(6,:) = 0;
   t_e(7,:) = 1;
   
   for I = 1:nodes_x-1,
      for J = 1:floor((nodes_y-1)/2),
         ct = 2*((nodes_y-1)*(I-1)+(J-1))+1;
         t_t(1,ct) = (I-1)*nodes_y+J+1;
         t_t(2,ct) = (I-1)*nodes_y+J;
         t_t(3,ct) = I*nodes_y+J+1;
         
         ct = ct+1;
         t_t(1,ct) = I*nodes_y+J+1;
         t_t(2,ct) = (I-1)*nodes_y+J;
         t_t(3,ct) = I*nodes_y+J;
      end
      for J = floor((nodes_y-1)/2)+1:nodes_y-1,
         ct = 2*((nodes_y-1)*(I-1)+(J-1))+1;
         t_t(1,ct) = (I-1)*nodes_y+J+1;
         t_t(2,ct) = I*nodes_y+J;
         t_t(3,ct) = I*nodes_y+J+1;
         
         ct = ct+1;
         t_t(1,ct) = (I-1)*nodes_y+J+1;
         t_t(2,ct) = (I-1)*nodes_y+J;
         t_t(3,ct) = I*nodes_y+J;
      end
   end
   t_t(4,:) = 1;
else
   [t_p,t_e,t_t]=initmesh(t_g,'hmax',2);
   %[t_p,t_e,t_t]=refinemesh(t_g,t_p,t_e,t_t);
   %[t_p,t_e,t_t]=refinemesh(t_g,t_p,t_e,t_t);
   %[t_p,t_e,t_t]=refinemesh(t_g,t_p,t_e,t_t);
   %[t_p,t_e,t_t]=refinemesh(t_g,t_p,t_e,t_t);
end
%plot(t_x,t_y);pause;