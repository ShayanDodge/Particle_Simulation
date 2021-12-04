%===================================
%    Shayan Dodge
%    MSc in Physics
% M: +98 935 741 67 63 
% E: dodgeshayan@gmail.com 
% A: Tehran, Iran 
%===================================
% A Fast Algorithm for Particle Simulations
%===================================
clc
clear all
%===================================
%% ==============Initialize====================
N=input('N=');%<N=number of particles><n=level>
n=floor(log(N)./log(4));
q=input('q=');
Q_total=sum(q);
x=input('x=');
y=input('y=');
x_min=min(x);
y_min=min(y);
x_max=max(x);
y_max=max(y);
x_domain=linspace(x_min,x_max,sqrt(4^n));
y_domain=linspace(y_min,y_max,sqrt(4^n));
z=zeros(sqrt(4^n),sqrt(4^n),n+1);% center at the mesh
a=zeros(sqrt(4^n),sqrt(4^n),n+1);
Q_cell=zeros(sqrt(4^n),sqrt(4^n),n+1);
for k=1:n+1
for i=1:k
for j=1:k   
y_length=y_min+(y_max-y_min)./(2^(k-1));
x_length=x_min+(x_max-x_min)./(2^(k-1));
z(i,j,k)= i.*(x_length./2) + j.*(y_length./2);
for charge=1:N
if (i-1).*x_length<=x(charge) & x(charge)<=i.*x_length & (j-1).*y_length<=y(charge) & y(charge)<=j.*y_length
a(i,j,k)=a(i,j,k)-q(charge).*sqrt(x(charge).^2+y(charge)^2);
Q_cell(i,j,k)=Q_cell(i,j,k)+q(charge);
end
end
end
end
end
phi=zeros(sqrt(4^n),sqrt(4^n),n+1);
psi=zeros(sqrt(4^n),sqrt(4^n),n+1);
psi_p=zeros(sqrt(4^n),sqrt(4^n),n+1);
%% ==============step 1======================= 
% From multipole expansions of potential field
% due to particles in each box about the box 
% center at the finest mesh level
phi = Q_cell.*log(z)+(a./z);
%% ==============step 2=======================

