X = [];

E = 2e11;

l = 1;
J = 0.5-9;

%число узлов и элементов
Nodes = 11;
Elements = Nodes - 1;
K_Global = zeros(Nodes*2);
F_Global = zeros(Nodes*2,1);

%Нагруженные элементы 9 10 11
Force_N = [1:1:10];

%Фиксированные узлы (первый)
Fixed_N = 1;

%узлы глобально
step = l/(Nodes-1);


k_e = local_k(step,E,J);
F_e = F_local(0.1);

j = 1;
for i=1:Elements
   
   K_Global(j:j+3,j:j+3)= K_Global(j:j+3,j:j+3) + k_e;
   j = j + 2;
   
end

K_Global(1:22,1)=zeros(22,1);
K_Global(1,1:22)=zeros(1,22);
K_Global(1,1)=1;

K_Global(2:22,2)=zeros(21,1);
K_Global(2,2:22)=zeros(1,21);
K_Global(2,2)=1;
K_Global(2,2);

for i=Force_N
    
   F_Global(2*i-1:2*i+2,1)=F_Global(2*i-1:2*i+2,1) + F_e;
end

F_Global(1, 1) = 0;
F_Global(2, 1) = 0;

u=inv(K_Global)*F_Global;

function f = F_local(l)
    g = -9.81;
    ro = 8050;
    S = 0.000144;
    f = l/2*g*ro*S*[1;
               l/6;
               1;
               -l/6];
end

function f = local_k(l,E,J)
    f = E*J/(l^3)*[12, 6*l, -12, 6*l;
                 6*l, 4*l^2, -6*l 2*l^2;
                 -12, -6*l, 12, -6*l;
                  6*l, 2*l^2, -6*l, 4*l^2];
end
