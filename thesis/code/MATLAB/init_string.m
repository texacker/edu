function plot_string = init_string
thesis_data;

if VORTEX
   plot_string = '涡旋分解方法';
else
   plot_string = '直接方法';
end

return