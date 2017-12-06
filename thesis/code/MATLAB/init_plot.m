function init_plot

thesis_data;

H_spectrum = figure('NumberTitle','off','Name','自振频谱');%,'CreateFcn','fcn_c_s','DeleteFcn','fcn_d_s');
hold on;
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
hold off;
title(['固有频谱(有限元网格' num2str(nodes_x) 'x' num2str(nodes_y) ',充液比' num2str(alpha) ',腔体高度' num2str(H) ')']);
xlabel('固有频率');
ylabel('涡旋分解方法-直接方法');
zoom(H_spectrum, 'on');

H_mesh = figure('NumberTitle','off','Name','有限元网格');%,'CreateFcn','fcn_c_m','DeleteFcn','fcn_d_m');
pdemesh(t_p,t_e,t_t);
axis equal;
if SpheroidDomain == 1
   alpha = alpha;%(1-alpha*alpha)*sqrt(1-alpha*alpha);
end   
title(['长轴:' num2str(R) ',短轴:' num2str(H) ' 充液比:' num2str(alpha)]);
xlabel(['有限元网格(椭球腔):节点数:' num2str(size(t_p,2)) ',单元数:' num2str(size(t_t,2))]);
ylabel('');
zoom(H_mesh, 'on');

return