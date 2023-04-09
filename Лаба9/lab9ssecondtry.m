clc;
clear all;
format short G;
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
% матрица K
K=zeros(2*nds,2*nds);
for i=1:els
    xi=xnds(assoc(i,1));
    yi=ynds(assoc(i,1));
    
    xj=xnds(assoc(i,2));
    yj=ynds(assoc(i,2));
    
    xk=xnds(assoc(i,3));
    yk=ynds(assoc(i,3));
%     какому материалу принадлежит узел
    if(((yi>=0) || (yj>=0) || (yk>=0)) && ((xi<=9/5*(yi+71/3)&& xi>=0) || (xj<=9/5*(yj+71/3)&& xj>=0) || (xk<=9/5*(yk+71/3)&& xk>=0)))
%         ige-1
       D=E1*(1-v1)/((1+v1)*(1-2*v1))*[[1,v1/(1-v1),0];[v1/(1-v1),1,0];[0,0,(1-2*v1)/(2*(1-v1))]];
    elseif((yi>=0) || (yj>=0) || (yk>=0)) && ((xi>9/5*(yi+71/3)&& xi<=1028.4) || (xj>9/5*(yj+71/3)&& xj<=1028.4) || (xk>9/5*(yk+71/3)&& xk<=1028.4))
%         ige-2
       D=E2*(1-v2)/((1+v2)*(1-2*v2))*[[1,v2/(1-v2),0];[v2/(1-v2),1,0];[0,0,(1-2*v2)/(2*(1-v2))]];
    else
%         ige-3
       D=E3*(1-v3)/((1+v3)*(1-2*v3))*[[1,v3/(1-v3),0];[v3/(1-v3),1,0];[0,0,(1-2*v3)/(2*(1-v3))]];
    end
    
    der = [[-1,0,1];[-1,1,0]];
    J = der*[[xi,yi];[xj,yj];[xk,yk]];
    Ji=inv(J);
%     считаем В
    B=zeros(3,6);
    tmp = zeros(2,1);
    tmp = Ji*der(:,1);
    B(1,1)=tmp(1);
    B(3,1)=tmp(2);
    B(2,2)=tmp(2);
    B(3,2)=tmp(1);
    
    tmp = Ji*der(:,2);
    B(1,3)=tmp(1);
    B(3,3)=tmp(2);
    B(2,4)=tmp(2);
    B(3,4)=tmp(1);
    
    tmp = Ji*der(:,3);
    B(1,5)=tmp(1);
    B(3,5)=tmp(2);
    
    B(2,6)=tmp(2);
    B(3,6)=tmp(1);
    
    Tmp=zeros(3);
    
    Tmp(:,3)=1;
    
    Tmp(1,1)=xi;
    Tmp(1,2)=yi;
    Tmp(2,1)=xj;
    Tmp(2,2)=yj;
    Tmp(3,1)=xk;
    Tmp(3,2)=yk;
    
    S=1/2*abs(det(Tmp));
    
    ke=B'*D*B*S;
    
    A=zeros(6,2*nds);
    
    it=assoc(i,1);
    jt=assoc(i,2);
    kt=assoc(i,3);
    
    A(1,2*it-1)=1;
    A(2,2*it)=1;
    A(3,2*jt-1)=1;
    A(4,2*jt)=1;
    A(5,2*kt-1)=1;
    A(6,2*kt)=1;
    
    K=K+A'*ke*A;
    
end


fe=zeros(2*nds,1);
grf1=[];
grf2=[];
grf3=[];

for i=1:nds
    xt=xnds(i);
    yt=ynds(i);
    xt
    (yt+71/3)*9/5
    if(xt<=(yt+71/3)*9/5 && yt>=0 && xt>=0)
        grf1=[grf1 i];
    end
    if(xt>(yt+71/3)*9/5 && yt>=0 && xt<=1028.4)
        grf2=[grf2 i];
    end
end


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

u=linsolve(K,fe);
ux=zeros(nds,1);
uy=zeros(nds,1);

for i=1:nds
   ux(i)=u(2*i-1);
   uy(i)=u(2*i);
end