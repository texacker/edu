%
echo off; clear; clc;

k=30;
a=1;
b=0.5;
h=0.7;
n=3;
%for n=1:30
x=0.9196:0.00001:0.9198;
y=miles_fun(1,0.6,1.2,1,n,x);
disp(['Mode:' num2str(n)]);
plot(x,abs(y));
pause;
%end

%x=(-1:0.001:-0.98);
%y = abs(miles_fun_v(a,b,h,20,x));
%plot(x,y);
%pause;

%x=(-1:0.01:1);
%y = abs(miles_fun_v(a,b,h,20,x));
%plot(x,y);
%pause;

x=(-5:0.1:0);
%x=(-1:0.025:1);
%x=(-1:0.001:-0.98);

for i=0:k
   y = abs(miles_fun_s(a,b,h,i,x))*0.9;
   z = abs(miles_fun(a,b,h,1,i,x));
   %z = abs(miles_fun(a,b,h,1,i,x));
   %y = miles_fun_s(a,b,h,n,x);
   %y = abs(miles_fun_v(a,b,h,n,x+0.075));
   plot(x,y,'b',x,z,'r');
   disp(['Mode:' num2str(i)]);
   pause;
end