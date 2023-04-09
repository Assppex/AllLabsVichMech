l = 0.1;
E = 2e11;
J = 0.000001633226782466667;
x = [0,0.1 , 0.2 , 0.3 ,0.4 ,0.5 , 0.6 , 0.7 , 0.8, 0.9 , 1];
el = 10;
nds = 11;
ke = E*J/(l)^3*[12,6*l,-12,6*l;6*l 4* l^2, - 6*l,2*l^2;-12, -6*l, 12 -6*l;6*l,2*l^2,-6*l,4*l^2];
K = zeros(2*nds,2*nds);
F = [0;-10000;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]
for i=1:el
    i_v =2*i-1;
    i_t = 2*i;
    
    j_v = 2*(i+1)-1;
    j_t = 2*(i+1);
    
    A = zeros(4,2*nds);
    A(1,i_v) = 1;
    A(2,i_t) = 1;
    A(3,j_v) = 1;
    A(4,j_t) = 1;
    
    T = A'*ke*A;
    K=K+T;
end
% √”
K(22,:) = 0;
K(:,22) = 0;
K(22,22) = 1;

K(21,:) = 0;
K(:,21) = 0;
K(21,21) = 1;

U = linsolve(K,F);
progibs = zeros(nds,1);
for i=1:nds
progibs(i) = U(2*i-1); 
end
rotations = zeros(nds,1);
for i=1:nds
rotations(i) = U(2*i); 
end
B_l = 1/l*[-6/l,-4,6/l,-2];
moment = zeros(nds,1)
for i=1:el
    moment(i)= -E* J*B_l*[progibs(i);rotations(i);progibs(i+1);rotations(i+1)];
end
B_r = 1/l*[6/l,2,-6/l,4];
moment(11)= -E* J*B_l*[progibs(10);rotations(10);progibs(11);rotations(11)];

force=[-4.13854E-10;
-3.5983E-10;
-5.17032E-10;
-1.00599E-09;
-1.28089E-09;
-2.31819E-10;
-4.77229E-10;
-6.69985E-10;
-6.69985E-10;
-7.47992E-10;
-7.67412E-10;
-7.60863E-10];
