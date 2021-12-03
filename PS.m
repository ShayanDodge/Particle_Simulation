% A Fast Algorithm for Particle Simulations
%% Initialize
clc
clear all
N=input('N=');
n=floor(log(N)./log(4));
x=input('x=');
y=input('y=');
x_min=min(x);
y_min=min(y);
x_max=max(x);
y_max=max(y);
x_domain=linspace(x_min,x_max,sqrt(4^n));
y_domain=linspace(y_min,y_max,sqrt(4^n));
Phi=zeros(sqrt(4^n),sqrt(4^n));
psi=zeros(sqrt(4^n),sqrt(4^n));
psi_p=zeros(sqrt(4^n),sqrt(4^n));
%% step 1



