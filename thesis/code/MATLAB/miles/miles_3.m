echo off; clear; clc;

n=2;
x=linspace(1,1.1,100);
for k=1:100
   y=miles_fun(1,0.5,0.8,1,k,x);
   disp(['Mode:' num2str(k)]);
   plot(x,y);grid;
   pause;
end