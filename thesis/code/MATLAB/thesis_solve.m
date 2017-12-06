function thesis_solve
thesis_data;

SpheroidDomain = 0;
RegularMesh= 1;
Fea_Len = 1;
nodes = 19;
w0 = 1.2*pi;
G0 = 0;
R = 1;
M0 = 1;

nodes_x = 6;
nodes_y = 11;%2n+1!
alpha = 0.75;
H = R/1;%(1-sqrt(1-alpha));

start_time = cputime;
VORTEX = 0;

%--------------------------------------------------------------------------
switch 3
%0: alpha
%1: Mesh
%2: H
%3: Direct-Vortex Method
%4: Subplot
%5: Sphere Cavity
%6: Mesh for Sphere Cavity
%--------------------------------------------------------------------------
case 0,
%--------------------------------------------------------------------------
nodes_x = 6;
nodes_y = 11;
for H=[4 1 0.25],
hold on;
spec_start = 0.05;
spec_step = 0.05;
for j=0:18
   alpha = spec_start+j*spec_step;
   init_mesh;%pdemesh(t_p,t_e,t_t); axis equal; pause;
   init_data;
   init_solve;
   
   plot_x = init_spectrum;
   spectrum = [plot_x plot_x];
   for i=1:size(spectrum,1)
      plot(spectrum(i,:),[0 1]*spec_step+alpha);
   end
   disp(['finished ' num2str(alpha)]);
end
title(['固有频谱(' init_string ',有限元网格' num2str(nodes_x) 'x' num2str(nodes_y) ',腔体高度' num2str(H) ')']);
xlabel('固有频率');
ylabel('充液比');
hold off;
pause;
close all;
end
disp('r u still alive ?'); pause;
%--------------------------------------------------------------------------
case 1,
%--------------------------------------------------------------------------
psnum=0;
for alpha=[0.75 0.5 0.25],
for H=[4 1 0.25],
psnum = psnum + 1;   
hold on;
for j=5:18,
   jj = j-4;
   nodes_x = fix(j/2)+1;%nodes_x = ceil(j/2);
   nodes_y = ceil(j/2)*2-1;%nodes_y = floor(j/2)*2+1;
   alpha, H, nodes_x, nodes_y
   init_mesh;%pdemesh(t_p,t_e,t_t); axis equal; pause;
   init_data;
   init_solve;
   
   plot_x = init_spectrum;
   spectrum = [plot_x plot_x];
   for i=1:size(spectrum,1),
      if abs(spectrum(i,1)) < 100
         plot(spectrum(i,:),[0 1]+jj);
      end
   end
   disp(['finished ' num2str(nodes_x) ' x ' num2str(nodes_y)]);
end
title(['固有频谱(' init_string ',腔体高度' num2str(H) ',充液比' num2str(alpha) ')']);
xlabel('固有频率');
ylabel('网格划分');
hold off;
print('-dps2',['PS_1_' num2str(psnum) '.ps']);
%print('-dcdeskjet',['DJ_1_' num2str(psnum) '.dj']);
%pause;
close all;
end
end
disp('r u still alive ?'); pause;
%--------------------------------------------------------------------------
case 2,
%--------------------------------------------------------------------------
H_step = 1;
nodes_x = 6;
nodes_y = 11;
for alpha=[0.5],%[0.75 0.25 0.5],
hold on;
for jj=1:50,
   j=jj*H_step;
   H = R/j;%H = (1-sqrt(1-alpha))/j;
   init_mesh;%pdemesh(t_p,t_e,t_t); axis equal; pause;
   init_data;
   init_solve;
   
   plot_x = init_spectrum;
   spectrum = [plot_x plot_x];
   for i=1:size(spectrum,1)
      if abs(spectrum(i,1)) < 5,
         plot(spectrum(i,:),[0 1]*H_step-j);
      end
   end
   disp(['finished ' num2str(j)]);
end

for jj=1:50,
   j=jj*H_step;
   H = R*j;%H = (1-sqrt(1-alpha))*j;
   init_mesh;%pdemesh(t_p,t_e,t_t); axis equal; pause;
   init_data;
   init_solve;
   
   plot_x = init_spectrum;
   spectrum = [plot_x plot_x];
   for i=1:size(spectrum,1)
      if abs(spectrum(i,1)) < 5,
         plot(spectrum(i,:),[-1 0]*H_step+j);
      end
   end
   disp(['finished ' num2str(j)]);
end
title(['固有频谱(' init_string ',有限元网格' num2str(nodes_x) 'x' num2str(nodes_y) ',充液比' num2str(alpha) ')']);
xlabel('固有频率');
%ylabel('圆柱腔半径(圆柱腔高度1)       圆柱腔高度(圆柱腔半径1)');
ylabel('');
hold off;
pause;
close all;
end
disp('r u still alive ?'); pause;
%--------------------------------------------------------------------------
case 3,
%--------------------------------------------------------------------------
nodes_x = 6;
nodes_y = 11;
for alpha=[0.75],%[0.75 0.5 0.25],
   for H=[4],%[4 1 0.25],
      fid = fopen(['3_' num2str(alpha) '_' num2str(H) '.DAT'],'a');
      init_mesh;%pdemesh(t_p,t_e,t_t); axis equal; pause;
      init_data;
      init_solve;
      
      VORTEX = 1; plot_x = (init_spectrum).';
      spectrum = [plot_x; plot_x; plot_x];
      spectrum(2,:) = 0;
      spectrum(3,:) = 1;
      fprintf(fid,'%12.8f %2d %2d\n',spectrum);
      
      VORTEX = 0; plot_x = (init_spectrum).';
      spectrum = [plot_x; plot_x; plot_x];
      spectrum(2,:) = 1;
      spectrum(3,:) = 2;
      fprintf(fid,'%12.8f %2d %2d\n',spectrum);
      
      fclose(fid);
   end
end
disp('r u still alive ?'); pause;
%--------------------------------------------------------------------------
case 4,
%--------------------------------------------------------------------------
nodes_x = 6;
nodes_y = 11;
for alpha=[0.75 0.5 0.25],
   for H=[4 1 0.25],
      init_mesh;%pdemesh(t_p,t_e,t_t); axis equal; pause;
      init_data;
      init_solve;
      
      subplot(2,1,1);
      VORTEX = 0; plot_x = init_spectrum;
      plot_y = ones(size(plot_x));
      bar(plot_x, plot_y, 0.01);
      title('固有频谱(直接方法)');
      subplot(2,1,2);
      VORTEX = 1; plot_x = init_spectrum;
      plot_y = ones(size(plot_x));
      bar(plot_x, plot_y, 0.01);
      title('固有频谱(涡旋分解方法)');
      pause;
      subplot(1,1,1);
   end
end

%--------------------------------------------------------------------------
case 5,
%--------------------------------------------------------------------------
SpheroidDomain = 1;
RegularMesh= 0;
R = 1;
H=1;
for nodes = [19 29],
for alpha=[0.75 0.5 0.25],
   hold on;
   init_mesh;%pdemesh(t_p,t_e,t_t); axis equal; pause;
   init_data;
   init_solve;
      
   VORTEX = 1; plot_x = init_spectrum;
   spectrum = [plot_x plot_x];
   for i=1:size(spectrum,1)
      plot(spectrum(i,:),[0 1]);
   end
   VORTEX = 0; plot_x = init_spectrum;
   spectrum = [plot_x plot_x];
   for i=1:size(spectrum,1)
      plot(spectrum(i,:),[1 2]);
   end
   title(['固有频谱(椭球腔:长轴:' num2str(R) ',短轴:' num2str(H) ',有限元网格:节点数:' num2str(size(t_p,2)) ',单元数:' num2str(size(t_t,2)) ')']);
   xlabel('固有频率');
   ylabel('涡旋分解方法    直接方法');
   hold off;
   pause;
   close all;
end
end
%--------------------------------------------------------------------------
case 6,
%--------------------------------------------------------------------------
SpheroidDomain = 1;
RegularMesh= 0;
R = 1;
H=1;
for nodes = [19 29],
for alpha=[0.75 0.5 0.25],
   alpha_2 = (1-alpha*alpha)*sqrt(1-alpha*alpha);
   init_mesh; pdemesh(t_p,t_e,t_t); axis equal;
   %init_data;
   %init_solve;
   title(['长轴:' num2str(R) ',短轴:' num2str(H) ' 充液比:' num2str(alpha_2)]);
   xlabel(['有限元网格(椭球腔):节点数:' num2str(size(t_p,2)) ',单元数:' num2str(size(t_t,2))]);
   ylabel('');
   pause;
   close all;
end
end
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------

for j=5:18,
   jj = j-4;
   nodes_x = fix(j/2)+1%nodes_x = ceil(j/2);
   nodes_y = ceil(j/2)*2-1%nodes_y = floor(j/2)*2+1;
   
   init_mesh;
   pdemesh(t_p,t_e,t_t); axis equal;
   title(['有限元网格' num2str(nodes_x) ' x ' num2str(nodes_y) ',节点数:' num2str(size(t_p,2)) ',单元数:' num2str(size(t_t,2))]);
   xlabel(['半径:' num2str(R) ',高度:' num2str(H) ' 充液比:' num2str(alpha)]);
   disp(['finished ' num2str(nodes_x) ' x ' num2str(nodes_y)]);
   pause;   
end

   init_mesh; pdemesh(t_p,t_e,t_t); axis equal; pause;
   init_data;
   init_solve;
   
   subplot(2,1,1);
   VORTEX = 0; plot_x = init_spectrum;
   plot_y = ones(size(plot_x));
   bar(plot_x, plot_y, 0.01);
   title('固有频谱(直接方法)');
   subplot(2,1,2);
   VORTEX = 1; plot_x = init_spectrum;
   plot_y = ones(size(plot_x));
   bar(plot_x, plot_y, 0.01);
   title('固有频谱(涡旋分解方法)');
   pause;
   subplot(1,1,1);
pdeplot(t_p,t_e,t_t...
   ,'xydata',zeros(size(t_p,2),1)...
   ,'zdata',zeros(size(t_p,2),1)...
   ,'mesh','on'...
   ,'title',['有限元网格' num2str(nodes_x) 'x' num2str(nodes_y)]...
   ,'colorbar','off'...
   ,'colormap',[1 1 1]);
pause;
   
for i=1:size(s_mode,2)
   pdeplot(t_p,t_e,t_t...
      ,'xydata',s_mode(:,i)...
      ,'zdata',s_mode(:,i)...
      ,'mesh','on'...
      ,'title',['固有频率: ' num2str(s_freq(i))]...
      ,'colorbar','off'...
      ,'colormap','default');%[1 1 1]);
   disp(['S Mode:' num2str(i) '=' num2str(s_freq(i))]);
   pause;
end

for i=1:size(e_vector,2)
   pdeplot(t_p,t_e,t_t...
      ,'xydata',e_vector(1:size(t_p,2),i)...
      ,'zdata',e_vector(1:size(t_p,2),i)...
      ,'mesh','on'...
      ,'title',['固有频率: ' num2str(e_value(i,i))]...
      ,'colorbar','off'...
      ,'colormap','default');%[1 1 1]);
   disp(['Mode:' num2str(i) '=' num2str(e_value(i,i))]);
   pause
end
x1 = 1:15;
y1 = [1.0863 1.3358 1.5777 1.8273 2.0877 2.3671 2.6639 2.9828 3.3173 3.6651 4.0086 4.3308 4.5969 4.7778 4.8401];
y2 = [1.4318 1.5991 1.7825 1.9883 2.2155 2.4693 2.7465 3.0503 3.3732 3.7123 4.0491 4.3665 4.6293 4.8083 4.8700];
y3 = [1.2462 1.4652 1.6817 1.9112 2.1557 2.4226 2.7096 3.0210 3.3496 3.6930 4.0331 4.3528 4.6172 4.7971 4.8591];
z2 = (y2-y1)./y1;
z3 = (y3-y1)./y1;
for i=1:15
   disp([num2str(y1(i)) ' ' num2str(y2(i)) ' ' num2str(y3(i)) ' ' num2str(z2(i)) ' ' num2str(z3(i))]);
end
plot(x1,y1,x1,y2,'x-',x1,y3,'+-');
title('绝对误差');
xlabel('频率阶数');
pause;
plot(x1,z2,x1,z3,'+-');
title('相对误差');
xlabel('频率阶数');
pause;

x1 = [3*5 4*5 4*7 5*7 5*9 6*9 6*11 7*11 7*13 8*13 8*15 9*15 9*17 10*17];
y1 = [1.4781 1.3557 1.2577 1.2082 1.1747 1.1510 1.1348 1.1209 1.1127 1.1041 1.0991 1.0934 1.0902 1.0863];
y2 = [1.5424 1.4597 1.3618 1.3351 1.3022 1.2905 1.2756 1.2695 1.2615 1.2579 1.2532 1.2508 1.2478 1.2462];
y3 = [1.6685 1.5969 1.5224 1.4992 1.4756 1.4653 1.4549 1.4494 1.4440 1.4407 1.4375 1.4354 1.4333 1.4318];
plot(x1,y1,x1,y2,'x-',x1,y3,'+-');
pause;

for j=5:18
%   nodes_x = fix(j/2)+1;
%   nodes_y = ceil(j/2)*2-1;
   nodes_x = ceil(j/2);
   nodes_y = floor(j/2)*2+1;
   disp([num2str(nodes_x) 'x' num2str(nodes_y)]);
   init_mesh;
   init_data;
   disp(' ');
end
pause;

%plot_x = 1:size(G,1);
%plot_y = sort(diag(e_value));
%plot(plot_x,(plot_y)');

%axis([-8 8 0 1]);

%pdemesh(p,e,t); axis equal;
%pdeplot(p,e,t);
%pdesurf(p,t,e_vector(1:size(p,2),1));
%pdecont(p,t,e_vector(1:size(p,2),1));

%plot_x = 1:size(p,2);
%while 1

%plot(1:size(t_p,2),e_vector(1:size(t_p,2),1)...
%   ,1:size(t_p,2),e_vector(1:size(t_p,2),4));
%   ,1:size(t_p,2),s_mode(:,1));

%subplot(3, ceil(size(s_mode,2)/3), i);


%plot_x = 1:size(G,1);
%plot_y = zeros(1,size(G,1));
%for j=1:size(G,1)
%   for i=1:size(G,1)
      %k=e_value(i,i);
%      plot_y(i) = norm(abs(e_vector(:,j)) - abs(e_vector(:,i)));
      %plot_y(i) = norm(e_vector(:,j) - e_vector(:,i));
%   end
%   plot(plot_x,plot_y);
%   disp(['Mode:' num2str(j) '=' num2str(e_value(j,j))]);
%   pause;
%end

%for i=1:size(u_mode,2)
%   pdeplot(t_p,t_e,t_t,'xydata',u_mode(:,i),'zdata',u_mode(:,i)),'mesh','on','title','Sloshing Modes');
%   disp(['Mode:' num2str(i) '=' num2str(u_freq(i))]);
%   pause;
%end
%   ,'colormap',[1 1 1]);
%   ,'colormap','gray');

return