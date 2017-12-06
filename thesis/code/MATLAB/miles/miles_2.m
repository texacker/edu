echo off; clear; clc;

miles_data;

m=1;
n=1;
a=1;
b=0.5;
h=0.8;

mile=100;
s_mile=ones(1,mile);
v_mile=-1*ones(1,mile);

for k=2:mile
   n=k;
   s_mile(k)=1/fzero('miles_zero',[ 0.01  0.99]);
   v_mile(k)=1/fzero('miles_zero',[-0.99 -0.01]);
   
   %disp(['Mode' num2str(xx) ':' num2str(fzero('miles_zero',[0.01 0.99]))]);
   %pause;
end

plot_x = sort([v_mile s_mile]);
plot_y = ones(size(plot_x));
bar(plot_x, plot_y);
pause;