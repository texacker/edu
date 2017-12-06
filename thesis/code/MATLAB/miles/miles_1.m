echo off; clear; clc;

A=1;
B=1;
nodes=80;

m=1;
n=1;
a=1;
b=0.6;
h=1*2*n*b;
t=1.451;
   
l=sqrt(t^2-1)*n*pi/h;
r=linspace(b,a,nodes);
lr=l*r;
f=A*besselj(m,lr)+B*bessely(m,lr);
f1=(A*(besselj(m-1,lr)-besselj(m+1,lr))+B*(bessely(m-1,lr)-bessely(m+1,lr)))*l/2;
f2=(A*(besselj(m-2,lr)+besselj(m+2,lr)-2*besselj(m,lr))+B*(bessely(m-2,lr)+bessely(m+2,lr)-2*bessely(m,lr)))*l*l/4;

u=f1+t*m*f./r;
v=f1*t+m*f./r;
w=(t^2-1)*n*pi/h*f;
tmp=f2+f1./r-((1-t^2)*(n*pi/h)^2+(m./r).^2).*f;
%a*f1(nodes)+m*t*f(nodes)
%t^2*b*f1(1)+(4*(1-t^2)+m*t^3)*f(1)
plot(r,tmp);pause;
%abs(miles_fun_s(a,b,h,1,t))
plot(r,f,'b',r,f1,'r',r,f2,'k');
%plot(r,u,'b',r,v,'r',r,w,'k');
%tmp=(f(2:nodes)-f(1:nodes-1))./(r(2:nodes)-r(1:nodes-1))-(f1(1:nodes-1)+f1(2:nodes))/2;
%plot(r(1:nodes-1),tmp);
%pause;

%x=(-1:0.01:1);
%y = abs(miles_fun_v(a,b,h,20,x));
%plot(x,y);
%pause;