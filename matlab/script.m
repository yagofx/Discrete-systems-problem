%-------------------------------
% Problema sistemas discretos
%-------------------------------
    
young = 1.1e10; % N/m^2
gamma = 24010;  % N/m^3
rmax = 2;       % m
rmin = 1;       % m
l = 3;          % m
p = 20e3;       % N

n_inicial = 3;
n_final = 25;
max_u =1:n_final-n_inicial+1;

for i=n_inicial:n_final
%Discretization
n = i;          % n of beams
nnodes = n+1;   % n of nodes
y_nodes = linspace(0,l,nnodes); % y coordinates at the nodes 
y_mid = (l/n)/2:l/n:l-(l/n)/2;  % y coordinates at the mid sections

%weight of section G's
%radius at the nodes
r_nodes =  (rmin-rmax)/l*y_nodes+rmax;  % radius at the nodes
r_mid = (rmin-rmax)/l*y_mid+rmax;       % radius at the mid section

%Distributed load

fg = zeros(n,1);
G = zeros(n,1);

% G in N/m
for i=1:n
    G(i) = gamma*(l/n)/3*(r_nodes(i)^2+r_nodes(i+1)^2+r_nodes(i)*r_nodes(i+1));
end

%fg in N
fg = G*(l/n)/2;

%Stiffnes matrix

k_element = young*pi*r_mid.^2/(l/n); % k of each beam
k = zeros(nnodes); % global stiffnes matrix

for i=1:n % building global stiffnes matrix
    k(i:i+1,i:i+1) = [1,-1;-1,1].*k_element(i)+k(i:i+1,i:i+1);
end

% %Global force vectors
% 
% f = zeros(nnodes,1); % global force vector
% 
% for i=1:nnodes  % global force vector
%     f(i)=1;
% end

%Get mini stiffnes matrix and mini force vector

mini_k = k(2:nnodes,2:nnodes);

%Creating mini_f

mini_f = zeros(n,1);

for i=1:n
    if i<n
    mini_f(i)=-fg(i)-fg(i+1);
    else
    mini_f(i)= -fg(i)-p;
    end
end

%Find displacement matrix

m_u = mini_k\mini_f;

%Get global displacement vector

u = zeros(nnodes,1);

u(2:end)=m_u;

%Get global force vector

r_f = k*u;

%Deformation of each beam

e = zeros(n,1);

for i=1:n
    e(i)= (u(i+1)-u(i))/(l/n);
end


%Strain stress

sigma = young*e;

%Strain

N = zeros(n,1);

for i=1:n
    N(i)= sigma(i)*pi*r_mid(i)^2;
end

%Verification

if sum(r_f) < 1e-6
    disp('Verification OK!')
else
    disp('Math error!')
    ends
end

max_u(i-2)= u(end);

end

x=n_inicial:n_final;
y=max_u;
plot(x,y,'-o');

title('Problema sistemas discretos')
xlabel('N (number of beams)')
ylabel('u, (max displacement)')
