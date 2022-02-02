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
%% test input==============================================================
N=64;% l=0,1,2,3
q=ones(1,64);
a_x=linspace(1,2,8);
x=[a_x(1),a_x(1),a_x(1),a_x(1),a_x(1),a_x(1),a_x(1),a_x(1),...
    a_x(2),a_x(2),a_x(2),a_x(2),a_x(2),a_x(2),a_x(2),a_x(2),...
    a_x(3),a_x(3),a_x(3),a_x(3),a_x(3),a_x(3),a_x(3),a_x(3),...
    a_x(4),a_x(4),a_x(4),a_x(4),a_x(4),a_x(4),a_x(4),a_x(4),...
    a_x(5),a_x(5),a_x(5),a_x(5),a_x(5),a_x(5),a_x(5),a_x(5),...
    a_x(6),a_x(6),a_x(6),a_x(6),a_x(6),a_x(6),a_x(6),a_x(6),...
    a_x(7),a_x(7),a_x(7),a_x(7),a_x(7),a_x(7),a_x(7),a_x(7),...
    a_x(8),a_x(8),a_x(8),a_x(8),a_x(8),a_x(8),a_x(8),a_x(8)];
y=[a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8)];
%==========================================================================
%  N=input('N=');% number of particles
% q=input('q='); % charge of particles
% x=input('x='); % location of particles (x,y)
% y=input('y='); 
n=floor(log(N)./log(4));% levels of refinement
Q_total=sum(q); % total charge
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
z_cg=zeros(sqrt(4^n),sqrt(4^n),n+1);% location of each box center
a_1=zeros(sqrt(4^n),sqrt(4^n),n+1);% the coefficient in the multipole expansion
a_2=zeros(sqrt(4^n),sqrt(4^n),n+1);
a_3=zeros(sqrt(4^n),sqrt(4^n),n+1);
a_4=zeros(sqrt(4^n),sqrt(4^n),n+1);
a_5=zeros(sqrt(4^n),sqrt(4^n),n+1);
Q_cell=zeros(sqrt(4^n),sqrt(4^n),n+1);% total charge inside of each box
N_cell=zeros(sqrt(4^n),sqrt(4^n));% number of particles inside of each
% the p-term multipole expansion (about the box center) of the potential field 
% created by the particles contained inside box i at level l.
phi=zeros(sqrt(4^n),sqrt(4^n),n+1);
% the p-term expansion about the center of box i at level l, describing 
% the potential field due to all particles inside the interaction list of ibox
psi=zeros(sqrt(4^n),sqrt(4^n),n+1);
% the p-term local expansion about the center of box i at level l, describing
% the potential field due to all particles inside the interaction list of
% ibox and ibox's parent
% nearest neighbors
psi_p=zeros(sqrt(4^n),sqrt(4^n),n+1);
%% ==============step 1====================================================
% From multipole expansions due to particles in each box
% about the box center at the finest mesh level

% This loop detects which particles are inside each box.
% z_cg, a_1, Q_cell and N_cell are calculated in this loop
for k=1:n+1
 for i=1:2^(k-1)
  for j=1:2^(k-1)   
y_length=(y_max-y_min)./(2^(k-1));
x_length=(x_max-x_min)./(2^(k-1));
z_cg(i,j,k)= sqrt((x_min+(x_length./2)+(i-1).*(x_length)).^2 +...
               (y_min+(y_length./2)+(j-1).*(y_length)).^2);
           
   for charge=1:N
    if x_min+(i-1).*x_length<=x(charge) &...
       x(charge)<=x_min+(i).*x_length   &...
       y_min+(j-1).*y_length<=y(charge) &...
       y(charge)<=y_min+(j).*y_length

a_1(i,j,k)=a_1(i,j,k)-q(charge).*abs(sqrt(x(charge).^2+y(charge)^2)-z_cg(i,j,k));
a_2(i,j,k)=a_2(i,j,k)-q(charge).*abs((sqrt(x(charge).^2+y(charge)^2)-z_cg(i,j,k)).^2)./2;
% a_3(i,j)=a_3(i,j)-q(charge).*((sqrt(x(charge).^2+y(charge)^2)-z_cg(i,j,k)).^3)./3;
% a_4(i,j)=a_4(i,j)-q(charge).*((sqrt(x(charge).^2+y(charge)^2)-z_cg(i,j,k)).^4)./4;
% a_5(i,j)=a_5(i,j)-q(charge).*((sqrt(x(charge).^2+y(charge)^2)-z_cg(i,j,k)).^5)./5;
Q_cell(i,j,k)=Q_cell(i,j,k)+q(charge);
if k==n+1
N_cell(i,j)=N_cell(i,j)+1;
end
    end
   end
  end
 end
end

%% ==============step 2====================================================
% Form multipole expansions about the centers of all boxes at all coarser 
% mesh levels, each expansion representing the potential field due to all
% particles contained in one box.
for k=n:-1:1
    for i=1:(sqrt(4^(k-1)))% i and j = location of parent
        for j=1:(sqrt(4^(k-1)))
            for ci=(i-1).*sqrt(4^(n-(k-1)))+1:(i).*sqrt(4^(n-(k-1)))% ci & cj = location of children
                for cj=(j-1).*sqrt(4^(n-(k-1)))+1:(j).*sqrt(4^(n-(k-1)))
                phi(i,j,k)=phi(i,j,k)+...
                Q_cell(ci,cj,k+1).*log10(abs(z_cg(i,j,k)-z_cg(ci,cj,k+1)))+...
                (a_1(ci,cj,n+1)-Q_cell(ci,cj,k+1).*z_cg(ci,cj,k+1))./z_cg(i,j,k)+...
                (a_2(ci,cj,n+1)+a_1(ci,cj,n+1).*z_cg(ci,cj,k+1)-(Q_cell(ci,cj,k+1).*z_cg(ci,cj,k+1).^2)./2)./z_cg(i,j,k)^2;
                end
            end
        end
    end
end
%% ==============step 3====================================================
% Form a local expansion about the center of each box at each mesh level 
% l<=n-1. This local expansion describes the field due to all particles in
% the system that are not contained in the current box or its nearest 
% neighbors. Once the local expansion is obtained for a given box, it is 
% shifted, in the second inner loop to the centers of the box’s children, 
% forming the initial expansion for the boxes at the next level.
for k=1:n
    for i=1:(sqrt(4^(k-1)))
        for j=1:(sqrt(4^(k-1)))
            for ci=1:sqrt(4^(k-1))
                for cj=1:sqrt(4^(k-1))
                   if ceil(ci./4)==ceil(i./4)-1 & ceil(cj./4)==ceil(j./4)  |...
                      ceil(ci./4)==ceil(i./4)   & ceil(cj./4)==ceil(j./4)  |...
                      ceil(ci./4)==ceil(i./4)+1 & ceil(cj./4)==ceil(j./4  )|...
                      ceil(ci./4)==ceil(i./4)-1 & ceil(cj./4)==ceil(j./4)+1|...
                      ceil(ci./4)==ceil(i./4)   & ceil(cj./4)==ceil(j./4)+1|...
                      ceil(ci./4)==ceil(i./4)+1 & ceil(cj./4)==ceil(j./4)+1|...
                      ceil(ci./4)==ceil(i./4)-1 & ceil(cj./4)==ceil(j./4)-1|...
                      ceil(ci./4)==ceil(i./4)   & ceil(cj./4)==ceil(j./4)-1|...
                      ceil(ci./4)==ceil(i./4)+1 & ceil(cj./4)==ceil(j./4)-1
                   if ci==i   & cj==j     |...
                      ci==i   & cj==j-1   |...
                      ci==i   & cj==j+1   |...
                      ci==i+1 & cj==j     |...
                      ci==i+1 & cj==j-1   |... 
                      ci==i+1 & cj==j+1   |...
                      ci==i-1 & cj==j     |...
                      ci==i-1 & cj==j-1   |... 
                      ci==i-1 & cj==j+1   |...
                      z_cg(i,j,k)==z_cg(ci,cj,k)
                   else
              
psi(i,j,k)=psi(i,j,k)+log10(abs((z_cg(i,j,k)-z_cg(ci,cj,k))))...
-a_1(ci,cj,k)./(z_cg(i,j,k)-z_cg(ci,cj,k))+Q_cell(ci,cj,k).*log10(abs(-(z_cg(i,j,k)-z_cg(ci,cj,k))))+...
(-a_1(ci,cj,k)./(z_cg(i,j,k)-z_cg(ci,cj,k)).^2-Q_cell(ci,cj,k)./(z_cg(i,j,k)-z_cg(ci,cj,k))).*(z_cg(i,j,k)-z_cg(ci,cj,k));
                   end
                   end
                end
            end
        end
    end
psi_p(:,:,k)=psi_p(:,:,k)+psi(:,:,k);
psi_p(1,:,:)=0;
psi_p(end,:,:)=0;
psi_p(:,1,:)=0;
psi_p(:,end,:)=0;
kp1=k+1;
     for ci=1:(sqrt(4^(kp1-1)))
        for cj=1:(sqrt(4^(kp1-1)))
psi_p(ci,cj,kp1)=psi_p(ci,cj,kp1)+psi(ceil(ci./sqrt(4^(n-(kp1-1)))),...
ceil(cj./sqrt(4^(n-(kp1-1)))),k).*(a_1(ceil(ci./sqrt(4^(n-(kp1-1)))),...
ceil(cj./sqrt(4^(n-(kp1-1)))),k).*(-z_cg(ceil(ci./sqrt(4^(n-(kp1-1)))),...
ceil(cj./sqrt(4^(n-(kp1-1)))),k)+z_cg(ci,cj,kp1)));
psi_p(1,:,:)=0;
psi_p(end,:,:)=0;
psi_p(:,1,:)=0;
psi_p(:,end,:)=0;
        end
     end 
end
%% ==============step 4====================================================
% Compute interactions at finest mesh level.
% particles that there are in the interaction list will be defined
    for i=1:(sqrt(4^(n)))
        for j=1:(sqrt(4^(n)))
            for ci=1:sqrt(4^(n))
                for cj=1:sqrt(4^(n))
                   if ceil(ci./4)==ceil(i./4)-1 & ceil(cj./4)==ceil(j./4)  |...
                      ceil(ci./4)==ceil(i./4)   & ceil(cj./4)==ceil(j./4)  |...
                      ceil(ci./4)==ceil(i./4)+1 & ceil(cj./4)==ceil(j./4  )|...
                      ceil(ci./4)==ceil(i./4)-1 & ceil(cj./4)==ceil(j./4)+1|...
                      ceil(ci./4)==ceil(i./4)   & ceil(cj./4)==ceil(j./4)+1|...
                      ceil(ci./4)==ceil(i./4)+1 & ceil(cj./4)==ceil(j./4)+1|...
                      ceil(ci./4)==ceil(i./4)-1 & ceil(cj./4)==ceil(j./4)-1|...
                      ceil(ci./4)==ceil(i./4)   & ceil(cj./4)==ceil(j./4)-1|...
                      ceil(ci./4)==ceil(i./4)+1 & ceil(cj./4)==ceil(j./4)-1
                   if ci==i   & cj==j     |...
                      ci==i   & cj==j-1   |...
                      ci==i   & cj==j+1   |...
                      ci==i+1 & cj==j     |...
                      ci==i+1 & cj==j-1   |... 
                      ci==i+1 & cj==j+1   |...
                      ci==i-1 & cj==j     |...
                      ci==i-1 & cj==j-1   |... 
                      ci==i-1 & cj==j+1   |...
                      z_cg(i,j,k)==z_cg(ci,cj,k)
                   else
psi(i,j,n+1)=psi(i,j,n+1)+log10(abs((z_cg(i,j,n+1)-z_cg(ci,cj,n+1))))...
-a_1(ci,cj,n+1)./(z_cg(i,j,n+1)-z_cg(ci,cj,n+1))+Q_cell(ci,cj,n+1).*log10(abs(-(z_cg(i,j,n+1)-z_cg(ci,cj,n+1))))+...
(-a_1(ci,cj,n+1)./(z_cg(i,j,n+1)-z_cg(ci,cj,n+1)).^2-Q_cell(ci,cj,n+1)./(z_cg(i,j,n+1)-z_cg(ci,cj,n+1))).*(z_cg(i,j,n+1)-z_cg(ci,cj,n+1));
psi_p(i,j,n+1)=psi_p(i,j,n+1)+psi(i,j,n+1);
                       end
                   end
                end
            end
        end
    end
%% ==============step 5====================================================
% Evaluate local expansions at particle positions.
for k=1:4^(n)
  cells(k).z_i=zeros(N_cell(k),1);% location of each particles inside each box
  cells(k).phi_z_i=zeros(N_cell(k),1);% the multiple expansion term at each particles location
  cells(k).dir_z_i=zeros(N_cell(k),1);% the direct term at each particles location
  cells(k).psi_z_i=zeros(N_cell(k),1);% the far field term at each particles location
  cells(k).total_phi_z_i=zeros(N_cell(k),1);% the total potential at each particles location
end
% This loop detects location of each particles inside each box (z=x^2+y^2).
cell_number=1;
 for i=1:2^(n)
  for j=1:2^(n) 
y_length=(y_max-y_min)./(2^(n));
x_length=(x_max-x_min)./(2^(n));
charge_number=1;
   for charge=1:N
    if x_min+(i-1).*x_length<=x(charge) &...
       x(charge)<=x_min+i.*x_length     &...
       y_min+(j-1).*y_length<=y(charge) &...
       y(charge)<=y_min+j.*y_length
cells(cell_number).z_i(charge_number)=sqrt(x(charge).^2+y(charge)^2);
charge_number=charge_number+1; 
    end
   end
cell_number=cell_number+1;
  end
 end

for i=1:2^n
 for j=1:2^n
     k=create_linear_index_list(z_cg,i,j,n+1)-(n.*(2^n).^2);
   for num_i=1:N_cell(k)
cells(k).psi_z_i(num_i)=psi_p(i,j,n+1);
cells(k).phi_z_i(num_i)=phi(1,1,1);
   end
 end
end
%% ==============step 6====================================================
% Compute potential (or force) due to nearest neighbors directly.
for i=2:2^n-1
 for j=2:2^n-1
k=create_linear_index_list(z_cg,i,j,n+1)-(n.*(2^n).^2);
  for ii=-1:1
   for jj=-1:1
    if ii~=0 & jj~=0
     for num_i=1:N_cell(k)
k_n=create_linear_index_list(z_cg,i+ii,j+jj,n+1)-(n.*(2^n).^2);
      for num_j=1:N_cell(k_n)
        if  cells(k).z_i(num_i)==cells(k_n).z_i(num_j)
        else
cells(k).dir_z_i(num_i)=cells(k).dir_z_i(num_i)-...
log10(abs(cells(k).z_i(num_i)-cells(k_n).z_i(num_j)));
      end
      end
     end
    end
   end
  end  
 end
end
%% ==============step 7====================================================
% Adding direct and farfield terms together
for k=1:4^n
   for num_i=1:N_cell(k)
cells(k).total_phi_z_i(num_i)=cells(k).phi_z_i(num_i)+cells(k).dir_z_i(num_i)+...
    cells(k).phi_z_i(num_i);
   end
end
end
