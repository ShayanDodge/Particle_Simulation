%==========================================================================
%    Shayan Dodge
%    MSc in Physics
% M: +98 935 741 67 63 
% E: dodgeshayan@gmail.com 
% A: Tehran, Iran 
%==========================================================================
% A Fast Algorithm for Particle Simulations
%==========================================================================
clc
clear all
%==========================================================================
%% ==============Initialize================================================
N=input('N=');% number of particles
n=floor(log(N)./log(4));% levels of refinement
q=input('q='); % charge of particles
Q_total=sum(q); % total charge
x=input('x='); % location of particles (x,y)
y=input('y='); 
if size(x)~=N | size(y)~=N | size(q)~=N
    disp('Error')
else
% Determining the simulation domain
x_min=min(x);
y_min=min(y);
x_max=max(x);
y_max=max(y);
x_domain=linspace(x_min,x_max,sqrt(4^n));
y_domain=linspace(y_min,y_max,sqrt(4^n));
% Defining some parameters
z=zeros(sqrt(4^n),sqrt(4^n),n+1);% location of each box center
a_1=zeros(sqrt(4^n),sqrt(4^n),n+1);% the coefficient in the multipole expansion
Q_cell=zeros(sqrt(4^n),sqrt(4^n),n+1);% total charge inside of each box
N_cell=zeros(sqrt(4^n),sqrt(4^n),n+1);% number of particles inside of each box
% the p-term multipole expansion (about the box center) of the potential field 
% created by the particles contained inside box i at level l.
phi=zeros(sqrt(4^n),sqrt(4^n),n+1);
% the p-term expansion about the center of box i at level l, describing 
% the potential field due to all particles outside the box and its nearest neighbors.
psi=zeros(sqrt(4^n),sqrt(4^n),n+1);
% the p-term local expansion about the center of box i at level l, describing
% the potential field due to all particles outside i’s parent box and the parent box’s
% nearest neighbors
psi_p=zeros(sqrt(4^n),sqrt(4^n),n+1);
% This loop detects that which particles are inside of each box.
% a_1 and Q_cell and N_cell are calculated in this loop
for k=1:n+1
 for i=1:2^(k-1)
  for j=1:2^(k-1)   
y_length=y_min+(y_max-y_min)./(2^(k-1));
x_length=x_min+(x_max-x_min)./(2^(k-1));
z(i,j,k)= i.*(x_length./2) + j.*(y_length./2);
   for charge=1:N
    if (i-1).*x_length<=x(charge) &...
       x(charge)<=i.*x_length     &...
       (j-1).*y_length<=y(charge) &...
       y(charge)<=j.*y_length
a_1(i,j,k)=a_1(i,j,k)-q(charge).*sqrt(x(charge).^2+y(charge)^2);
Q_cell(i,j,k)=Q_cell(i,j,k)+q(charge);
N_cell(i,j,k)=N_cell(i,j,k)+1;
    end
   end
  end
 end
end
% 
for k=1:n
  cells(k).z_i=zeros(N_cell(k),1);
  cells(k).phi_z_i=zeros(N_cell(k),1);
end

for k=1:n+1
 for i=1:k
  for j=1:k   
y_length=y_min+(y_max-y_min)./(2^(k-1));
x_length=x_min+(x_max-x_min)./(2^(k-1));
z(i,j,k)= i.*(x_length./2) + j.*(y_length./2);
number=1;
   for charge=1:N
    if (i-1).*x_length<=x(charge) &...
       x(charge)<=i.*x_length     &...
       (j-1).*y_length<=y(charge) &...
       y(charge)<=j.*y_length
   
cells(k).z_i(number)=sqrt(x(charge).^2+y(charge)^2);
number=number+1;
    end
   end
  end
 end
end
%% ==============step 1======================= 
% From multipole expansions of potential field
% due to particles in each box about the box 
% center at the finest mesh level
phi(:,:,n+1) = Q_cell(:,:,n+1).*log(z(:,:,n+1))+(a_1(:,:,n+1)./z(:,:,n+1));
%% ==============step 2=======================
% Form multipole expansions about the centers
% of all boxes at all coarser mesh levels, each
% expansion representing the potential field
% due to all particles contained in one box.
for k=1:n
    for i=1:(sqrt(4^(k-1)))
        for j=1:(sqrt(4^(k-1)))
            for ci=i.*sqrt(4^(n-(k-1))):(i-1).*sqrt(4^(n-(k-1)))+1
                for cj=j.*sqrt(4^(n-(k-1))):(j-1).*sqrt(4^(n-(k-1)))+1
                phi(i,j,k)=phi(i,j,k)+...
                Q_cell(ci,cj,k).*log(z(ci,cj,k))+Q_cell(ci,cj,k).*z(ci,cj,k);
                end
            end
        end
    end
end
%% ==============step 3=======================
% Form a local expansion about the center of each...
% box at each mesh level l<=n-1. This local expansion...
% describes the field due to all particles in the...
% system that are not contained in the current box...
% or its nearest neighbors. Once the local expansion...
% is obtained for a given box, it is shifted,...
% in the second inner loop to the centers of the...
% box’s children, forming the initial expansion for...
% the boxes at the next level.
for k=1:n-1
    for i=1:(sqrt(4^(k-1)))
        for j=1:(sqrt(4^(k-1)))
            for ci=1:sqrt(4^(k-1))
                for cj=1:sqrt(4^(k-1))
                   if cj==j & ci==i 
                   else
                psi(i,j,k)=psi(i,j,k)+...
                -a_1(ci,cj,k)./z(ci,cj,k)+Q(ci,cj,k).*log(z(ci,cj,k))+(a_1(ci,cj,k)./z(ci,cj,k).^2-Q./z(ci,cj,k)).*z(ci,cj,k);
                   end
                end
            end
        end
    end
 psi_p(i,j,k)=psi(i,j,k);
 kp1=k+1;
     for ci=1:(sqrt(4^(kp1-1)))
        for cj=1:(sqrt(4^(kp1-1)))
                psi_p(ci,cj,kp1)=psi_p(ci,cj,kp1)+...
                    psi(ceil(ci./sqrt(4^(n-(kp1-1)))),ceil(cj./sqrt(4^(n-(kp1-1)))),k).*(Q-a_1(ci,cj,kp1).*z(ci,cj,kp1)).*z(ci,cj,kp1);
        end
     end 
end
%% ==============step 4=======================
% Compute interactions at finest mesh level.
    for i=1:(sqrt(4^(n)))
        for j=1:(sqrt(4^(n)))
            for ci=1:sqrt(4^(n))
                for cj=1:sqrt(4^(n))
                   if cj==j & ci==i 
                   else
                psi(i,j,k)=psi(i,j,k)+...
                -a_1(ci,cj,k)./z(ci,cj,k)+Q(ci,cj,k).*log(z(ci,cj,k))+(a_1(ci,cj,k)./z(ci,cj,k).^2-Q./z(ci,cj,k)).*z(ci,cj,k);
                   end
                end
            end
        end
    end
%% ==============step 5=======================
% Evaluate local expansions at particle positions.
for k=1:4^n
    for num_i=1:N_cell(k)
        for num_j=1:N_cell(k)
cells(k).phi_z_i(num_i)=-log(cells(k).z_i(num_j)-cells(k).z_i(num))
        end
    end
end
%% ==============step 6=======================





end
