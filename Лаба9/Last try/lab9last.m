clc;
clear all;

E1=55e6;
E2=35e6;
E3=23e9;

v1=0.22;
v2=0.22;
v3=0.25;

d1=2300;
d2=2300;
denwat=1000;

water_level=290;

xnds=[];
ynds=[];
g=9.8;

Temp=readtable('xydata.txt');
xnds=Temp(:,2);
ynds=Temp(:,3);
xnds=table2array(xnds);
ynds=table2array(ynds);
Temp=readtable('assoc.txt'); 
a=Temp(:,2);
a=table2array(a);
b=Temp(:,3);
b=table2array(b);
c=Temp(:,4);
c=table2array(c);


nds=size(xnds,1);
els=size(Temp,1);
assoc=[];

for i=1:els
    temp=[a(i),b(i),c(i)];
    assoc(i,:)=temp;
end

[ige1,ige2,base,kol,kol1,kol2]=get_nodes(xnds,ynds,nds);
[ige1m,ige2m,basem,kol3,kol4,kol5]=get_elements(xnds, ynds, assoc, els);

% K1 = km1(ige1m, assoc, nds, xnds, ynds, E1,v1);
% K2 = km1(ige2m, assoc, nds, xnds, ynds, E2,v2);
% K3 = km1(base, assoc, nds, xnds, ynds, E3,v3);

K= getK(ige1m,ige2m,basem, assoc, els, xnds, ynds, E1,v1,E2,v2,E3,v3,nds);

fe = zeros(2*nds,1);
fe = gravity_force(ige1m, ige2m, xnds, ynds, assoc, els, fe);
% fe(8)=-100000;

fe = fe + NPUf(xnds,ynds,nds,ige1,base);

bc_left=[];
for j=1:nds
    x_c1 = xnds(j);
%     ищем узлы заполняющие рассматриваемую поверхность
    if((x_c1==-1.028400020000000e+03))
        bc_left=[bc_left j];
    end
end
bc_right=[];
for j=1:nds
    x_c1 = xnds(j);
%     ищем узлы заполняющие рассматриваемую поверхность
    if(x_c1==2.056800050000000e+03)
        bc_right=[bc_right j];
    end
end
bc_bottom=[];
for j=1:nds
    y_c1 = ynds(j);
%     ищем узлы заполняющие рассматриваемую поверхность
    if((y_c1==-293))
        bc_bottom=[bc_bottom j];
    end
end
% учитываем гу в матрице К
for j=1:size(bc_bottom,2)
    K(2*bc_bottom(j),:)=0;
    K(:,2*bc_bottom(j))=0;
    K(2*bc_bottom(j),2*bc_bottom(j))=1;
end

for j=1:size(bc_left,2)
    K(2*bc_left(j)-1,:)=0;
    K(:,2*bc_left(j)-1)=0;
    K(2*bc_left(j)-1,2*bc_left(j)-1)=1;
    

end

for j=1:size(bc_right,2)
    K(2*bc_right(j)-1,:)=0;
    K(:,2*bc_right(j)-1)=0;
    K(2*bc_right(j)-1,2*bc_right(j)-1)=1;
    
end
% c=0;
% uz=[];
% for i=1:2*nds
%    if(K(i,i)==0)
%     c=c+1;
%     K(i,i)=1;
%     if(mod(i,2)==0)
%     uz = [uz i/2]
%     end
%    end
% end
u=linsolve(K,fe);
ux=zeros(nds,1);
uy=zeros(nds,1);

for i=1:nds
   ux(i)=u(2*i-1);
   uy(i)=u(2*i);
end
[strain, stress] = get_strain_and_stress(ux,uy, els, xnds, ynds, assoc, ige1m, ige2m,basem,E1,v1,E2,v2,E3,v3);

maxx = max(ux)

maxy =max(uy)