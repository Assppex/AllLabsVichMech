clc;
clear all;
k1 = 1.5;
k2 = 1.75;
Te = zeros(1,3);
n_els = 581;
n_nds = 337;
derivatives = [[-1,0,1];[-1,1,0]];
x=[];
y=[];
Temp = readtable('D:\Учеба\3 курс\Вычислительная механика\Лаба8\Fixed_dimensions\nodes.txt');
x= Temp(:,2);
x=table2array(x);

y= Temp(:,3);
y=table2array(y);

Temp = readtable('D:\Учеба\3 курс\Вычислительная механика\Лаба8\Fixed_dimensions\data.txt');
a=[];
b=[];
c=[];
a = Temp(:,2);
b = Temp(:,3);
c = Temp(:,4);

a=table2array(a);
b=table2array(b);
c=table2array(c);
assoc=zeros(294,3);

for i=1:n_els
    temp = [a(i),b(i),c(i)];
    assoc(i,1)=a(i);
    assoc(i,2)=b(i);
    assoc(i,3)=c(i);
end
Kc=zeros(n_nds);
% составляем матрицу Kc
for m=1:n_els
    Kce=zeros(3);
    
    i=assoc(m,1);
    j=assoc(m,2);
    k=assoc(m,3);

    xi=x(i);
    yi=y(i);

    xj=x(j);
    yj=y(j);

    xk=x(k);
    yk=y(k);

    % матрица якоби и якобиан
    J = derivatives*[[xi,yi];[xj,yj];[xk,yk]];
    detJ = det(J);
    % Матрица Be
    Be=zeros(2,3);
    for n=1:3
        Be(:,n)=inv(J)*derivatives(:,n);
    end
    Kce = 1/2*Be'*Be*detJ;
    %3 узла, в каждом по 1 температуре , тогда матрица А из прошлыш лаб перейдет в:     
    A = zeros(3,n_nds);
    
    A(1,i)=1;
    A(2,j)=1;
    A(3,k)=1;
        %     проверяем принадлжеит ли элемент материалу ИГЭ-1 или ИГЭ-2 (сетка построена так, что если хоть один узел лежит внутри зоны одного из материалов, то и весь элемент лежит там)
    if((xi>42.6 && yi<=-71/3+5/9*xi && xi<1028.4) || (xj>42.6 && yj<=-71/3+5/9*xi && xj<1028.4 ) || (xk>42.6 && yk<=-71/3+5/9*xk && xk<1028.4))
        Kce=k2*Kce;
    else
        Kce=k1*Kce;
    end
    Kc =Kc+A'*Kce*A;
end
% теперь выставляем ГУ
% вода
bcw=[];
% воздух
bca=[];
for m=1:n_els
    
    i=assoc(m,1);
    j=assoc(m,2);
    k=assoc(m,3);

    xi=x(i);
    yi=y(i);

    xj=x(j);
    yj=y(j);

    xk=x(k);
    yk=y(k);
    %     проверяем принадлжеит узел участку где вода
    if(yi<=290 && xi<=462 && sum(bcw(:) == i) == 0)
        if((yi==275 && xi<=424.5 && xi>=412.5) || (abs(abs(xi*2/3)-abs(yi))<10^(-2)) || (abs(abs(-424.5+xi*2/5)-abs(yi))<10^(-2)))
            bcw=[bcw i];
        end
    end
    
     if(yj<=290 && xj<=462 && sum(bcw(:) == j) == 0)
        if((yj==275 && xj<=424.5 && xj>=412.5) || (abs(abs(xj*2/3)-abs(yj))<10^(-2)) || (abs(abs(-424.5+xj*2/5)-abs(yj))<10^(-2)))
            bcw=[bcw j];
        end
     end
    
      if(yk<=290 && xk<=462 && sum(bcw(:) == k) == 0)
        if((yk==275 && xk<=424.5 && xk>=412.5) || (abs(abs(xk*2/3)-abs(yk))<10^(-2)) || (abs(abs(-424.5+xk*2/5)-abs(yk))<10^(-2)))
            bcw=[bcw k];
        end
      end
      
%       && yj>=290 && yj<293
%        && yk>=290 && yk<293
%       (abs(-424.5+xi*2/5)-abs(yi)<10^(-8))||
    %     проверяем принадлжеит узел участку где водздух
      if(xi>462 && yi~=0 && sum(bca(:) == i) == 0)
        if((abs(-424.5+xi*2/5)-abs(yi)<=10^(-8) && yi>=290 && yi<293 && xi<=469.5)||(yi==293)||(abs(-1/2*xi+535.25)-abs(yi)<10^(-5))||(yi==270.5)||(abs(-5/9*xi+1714/3)-abs(yi)<10^(-5)))
            bca=[bca i];
        end
      end
    
      if(xj>462 && yj~=0 && sum(bca(:) == j) == 0)
        if((abs(-424.5+xj*2/5)-abs(yj)<=10^(-8) && yj>=290 && yj<293 && xj<=469.5)||(yj==293)||(abs(-1/2*xj+535.25)-abs(yj)<10^(-5) )||(yj==270.5)||(abs(-5/9*xj+1714/3)-abs(yj)<10^(-5)))
            bca=[bca j];
        end
      end
     
      if(xk>462 && yk~=0 && sum(bca(:) == k) == 0)
        if((abs(-424.5+xk*2/5)-abs(yk)<=10^(-8) && yk>=290 && yk<293 && xk<=469.5)||(yk==293)||(abs(-1/2*xk+535.25)-abs(yk)<10^(-5))||(yk==270.5)||(abs(abs(-5/9*xk+1714/3)-abs(yk))<10^(-5)))
            bca=[bca k];
        end
      end
      

end
R=zeros(n_nds,1);

T_water = 5;
for i=1:length(bcw)
   Kc(bcw(i),:)= 0;
   Kc(bcw(i),bcw(i))=1;
   R(bcw(i))=T_water;
end

T_air = 25;
for i=1:length(bca)
   Kc(bca(i),:)= 0;
   Kc(bca(i),bca(i))=1;
   R(bca(i))=T_air;
end
T_res = linsolve(Kc,R);


