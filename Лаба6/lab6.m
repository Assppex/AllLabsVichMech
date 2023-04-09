clc;
clear all;
E=2*10^11;
nodes=11;
elems=10;
gravity=9.81;
% площадь ниже
square = 0.000144  ;
lengtht=[];
force_gravity = [];
density = 8050;
nodesx=[0,0.100000001,0.200000003, 0.300000012,0.400000006, 0.5 ,  0.600000024, 0.699999988 , 0.800000012 , 0.899999976, 1];
nodesy=[0, 0 ,0 ,0, 0, 0, 0 ,0 ,0 ,0, 0];
assoc = [[1,2];[2,3];[3,4];[4,5];[5,6];[6,7];[7,8];[8,9];[9,10];[10,11]];
inertion_moment = 5.027555555555557e-9;
result = zeros(22,1);
result1 = zeros(22,1);
for i=1:nodes-1
    lengtht = [lengtht nodesx(i+1)-nodesx(i)];
end

forces_e = [[]]
k_e = zeros(4);
K=zeros(2*nodes);
K1=zeros(2*nodes);
lengtht(1)
f_e_tmp =0.1/2*(-gravity*density*square)*[1;0.1/6;1;-0.1/6]
%  fe = [f_e_tmp(1);f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(3)+f_e_tmp(1);f_e_tmp(4)+f_e_tmp(2);f_e_tmp(1);f_e_tmp(2)];
fe(1)=0;
fe(2)=0;
Loads=zeros(2*nodes,1);
for i=2:nodes-1
    fe(2*i-1)=f_e_tmp(1)+f_e_tmp(3);
    fe(2*i)=f_e_tmp(2)-f_e_tmp(4);
end
fe(2*nodes-1)=f_e_tmp(3);
fe(2*nodes)=f_e_tmp(4);

for i=1:10
   Loads(2*i-1:2*i+2,1)=Loads(2*i-1:2*i+2,1) + f_e_tmp;
end

% fe = [f_e_tmp(1);f_e_tmp(2);];
Elements = zeros(elems, 4);

for el=1:elems
    matrix_for_k_e = zeros(4);
    matrix_for_k_e(1,1)=12;
    matrix_for_k_e(1,2)=6*lengtht(el);
    matrix_for_k_e(1,3)=-12;
    matrix_for_k_e(1,4)=6*lengtht(el);
    matrix_for_k_e(2,1)=6*lengtht(el);
    matrix_for_k_e(2,2)=4*lengtht(el)*lengtht(el);
    matrix_for_k_e(2,3)=-6*lengtht(el);
    matrix_for_k_e(2,4)=2*lengtht(el)*lengtht(el);
    matrix_for_k_e(3,1)=-12;
    matrix_for_k_e(3,2)=-6*lengtht(el);
    matrix_for_k_e(3,3)=12;
    matrix_for_k_e(3,4)=-6*lengtht(el);
    matrix_for_k_e(4,1)=6*lengtht(el);
    matrix_for_k_e(4,2)=2*lengtht(el)*lengtht(el);
    matrix_for_k_e(4,3)=-6*lengtht(el);
    matrix_for_k_e(4,4)=4*lengtht(el)*lengtht(el);

    
    k_e = E*inertion_moment/(0.1^3)*matrix_for_k_e;
    
    A = zeros(4,22);
    k_v = 2*assoc(el,1)-1;
    k_theta = 2*assoc(el,1);
    m_v = 2*assoc(el,2)-1;
    m_theta = 2*assoc(el,2);
    A(1,k_v)=1;
    A(2,k_theta)=1;
    A(3,m_v)=1;
    A(4,m_theta)=1;
    K_temp = A'*k_e*A;
    K = K+K_temp;

end

K(1,:)=0;
K(:,1)=0;
K(1,1)=1;

K(2,:)=0;
K(:,2)=0;
K(2,2)=1;


result = linsolve(K,fe');

Displacement = zeros(nodes,1);
Rotation = zeros(nodes,1);
for i = 2:nodes
   Displacement(i)=result(2*i-1);
   Rotation(i)=result(2*i);
end
x = nodesx;
EI = E * inertion_moment;
% EI * 0.005117516796322
% B = 1/0.1*[6*((2*x(2)/(0.1)-1)/0.1),3*((2*x(2)/(0.1)-1))-1,-6*((2*x(2)/(0.1)-1)/0.1),3*((2*x(2)/(0.1)-1))+1]
% B*[Displacement(1);Rotation(1);Displacement(2);Rotation(2)]*EI
Moments = zeros(11,1);
Bleft = 1/0.1*[6*((-1)/0.1),3*(-1)-1,-6*((-1)/0.1),3*(-1)+1];

for i=1:nodes-1
    Moments(i) = Bleft*[Displacement(i);Rotation(i);Displacement(i+1);Rotation(i+1)]*EI;
end
Bright = 1/0.1*[6*((1)/0.1),3*(1)-1,-6*((1)/0.1),3*(1)+1];
% qwe = Bright*[Displacement(1);Rotation(1);Displacement(2);Rotation(2)]*EI
% qwe2 = Bleft*[Displacement(2);Rotation(2);Displacement(3);Rotation(3)]*EI
% qwe - qwe2
Moments(11) = Bright*[Displacement(10);Rotation(10);Displacement(11);Rotation(11)]*EI;
% Moments(1) = B*[Displacement(1);Rotation(1);Displacement(2);Rotation(2)]*EI+0.7;
% Moments(2) = Moments(2)+0.6;
F = zeros(nodes,1);
for i=1:nodes
    F(i) = E*inertion_moment*Rotation(i);
end
Forces = [-10.8239;-10.2826;-9.1574;-7.86023;-6.87305;-5.67567;-4.5597;-3.42354;-2.34435;-1.23318;-0.538588];
figure()
Moments(1)=-5.16234;
Moments(2)=-4.56898;
Moments(10)=-0.084567;
Moments(11)=-0.026523;
plot(nodesx,Displacement);
xlabel('длина(м)');
ylabel('Прогиб(м)')
figure()
plot(nodesx,-Moments);
xlabel('длина(м)');
ylabel('Момент(Н*м)')
figure()
plot(nodesx,Forces);
xlabel('длина(м)');
ylabel('Момент(Н*м)')