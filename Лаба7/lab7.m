clc
clear all
g = 9.81;
s = 0.000144;
dens = 7900;
TS = 1000000;
dt  = 1/TS;
l = 0.1;
nodes = 11;
els = 10;
I = 5.027555555555557e-9;
E = 2e11;
x = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
asc = [[1,2];[2,3];[3,4];[4,5];[5,6];[6,7];[7,8];[8,9];[9,10];[10,11]];
ftmp = l/2*-g*dens*s*[1;l/6;1; -l/6];
F = zeros(2*nodes,1);
% for i=1:els
%     F(2*i-1:2*i+2,1)=F(2*i-1:2*i+2,1) + ftmp;
% end
F(1)=0;
F(2)=0;
%  fe = [f_e_tmp(1);f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(1);f_e_tmp(2)];

for i=2:nodes-1
    F(2*i-1)=ftmp(1)+ftmp(3);
    F(2*i)=ftmp(2)-ftmp(4);
end
F(2*nodes-1)=ftmp(3);
F(2*nodes)=ftmp(4);
% F(11)=-dens*g*s;
matrix_for_k_e = zeros(4);
matrix_for_k_e(1,1)=12;
matrix_for_k_e(1,2)=6*l;
matrix_for_k_e(1,3)=-12;
matrix_for_k_e(1,4)=6*l;
matrix_for_k_e(2,1)=6*l;
matrix_for_k_e(2,2)=4*l*l;
matrix_for_k_e(2,3)=-6*l;
matrix_for_k_e(2,4)=2*l*l;
matrix_for_k_e(3,1)=-12;
matrix_for_k_e(3,2)=-6*l;
matrix_for_k_e(3,3)=12;
matrix_for_k_e(3,4)=-6*l;
matrix_for_k_e(4,1)=6*l;
matrix_for_k_e(4,2)=2*l*l;
matrix_for_k_e(4,3)=-6*l;
matrix_for_k_e(4,4)=4*l*l;

me = [156, 12 * l, 54, -13 * l; 22 * l, 4 * l * l, 13 * l, -3 * l * l; 54, 13 * l, 156, -22 * l; -13 * l, -3 * l * l, -22 * l, 4 * l * l;];
me = s * dens * l / 420 * me;
ke= E*I/(0.1^3)*matrix_for_k_e;

K = zeros(2*nodes);
m = zeros(2*nodes);
for i=1:els
    j = asc(i,1);
    v = asc(i,2);
    a_v = 2*j-1;
    a_t = 2*j;
    b_v = 2*v-1;
    b_t = 2*v;
    A = zeros(4,2*nodes);
    A(1,a_v)=1;
    A(2,a_t)=1;
    A(3,b_v)=1;
    A(4,b_t)=1;
    mtmp = A'*me*A;
    ktmp = A'*ke*A;
    K = K+ktmp;
    m = m+mtmp;
end
K(1,:)=0;
K(:,1)=0;
K(1,1)=1;

K(2,:)=0;
K(:,2)=0;
K(2,2)=1;

r = linsolve(K,F);
m(1,:)=0;
m(:,1)=0;
m(1,1)=1;

m(2,:)=0;
m(:,2)=0;
m(2,2)=1;

u = zeros(22);
du = zeros(22);
ddu = zeros(22);
result_u = zeros;
result_du = zeros;
result_ddu = zeros;

dt =0.01;
%  ddu(:, 1) = [0;0;-2.221E-07;0;-2.06545E-07;0;-1.99988E-07;0;-1.97584E-07;0;-1.96706E-07;0;-1.96385E-07;0;-1.96267E-07;0;-1.96223E-07;0;-1.96206E-07;0;-1.96178E-07;0]
ddu(:, 1) = linsolve(m,zeros(2*nodes,1));
preu = zeros(2*nodes,1)
prev = zeros(2*nodes,1)
for k=1:99
    fnext = zeros(2*nodes,1);
    tmp = k-1;
   if(k+1<=50)
       fnext= F*2*(tmp+1)/100;
   end 
   preu = u(:,k)+du(:,k)*dt + 0.25*ddu(:,k)*dt^2;
   prev = du(:,k) + 0.5*ddu(:,k)*dt;
   ddu(:,k+1) = inv(m+0.25*K*dt^2)*(fnext - K*preu);
   u(:,k+1)=preu+ddu(:,k+1)*0.25*dt^2;
   du(:,k+1)=prev+ddu(:,k+1)*0.5*dt;
   for j = 1:11
       result_ddu(j, k+1) = ddu(2 * j - 1, k+1);
       result_du(j, k+1) = du(2 * j - 1, k+1);
       result_u(j, k+1) = u(2 * j - 1, k+1);
   end
end

% for i=1:TS
%    Fs = zeros(2*nodes,1);
%    if(i<=TS/2)
%        Fs = F*((i-1)/TS/2);
%    end
%    ddu(:, i) = m \ (Fs - K * u(:, i));
%    du(:, i + 1) = du(:, i) + dt * ddu(:, i);
%    u(:, i + 1) = u(:, i) + dt * du(:, i + 1);
%    for j = 1:11
%        result_ddu(j, i) = ddu(2 * j - 1, i);
%        result_du(j, i) = du(2 * j - 1, i);
%        result_u(j, i) = u(2 * j - 1, i);
%    end
% end