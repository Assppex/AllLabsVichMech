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
Temp = readtable('D:\Учеба\3 курс\Вычислительная механика\Лаба8\nodes.txt');
x= Temp(:,2);
x=table2array(x);

y= Temp(:,3);
y=table2array(y);

Temp = readtable('D:\Учеба\3 курс\Вычислительная механика\Лаба8\data.txt');
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
%     elseif((xi>=3202.5 && xi<=4402.5 && yi<=270.5) ||(xj>=3202.5 && xj<=4402.5 && yj<=270.5) || (xk>=3202.5 && xk<=4402.5 && yk<=270.5))
%         Kce=k2*Kce;
%     elseif ((xi>2715.6 && yi==-2715.6+1.8*xi && xi<3202.5) || (xj>2715.6 && yj==-2715.6+1.8*xj && xj<3202.5 ) || (xk>2715.6 && yk==-2715.6+1.8*xk && xk<3202.5))
%         Kce=k2*Kce;
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
        if((yi==275 && xi<=424.5 && xi>=412.5) || (xi*2/3==yi) || (-424.5+xi*2/5==yi))
            bcw=[bcw i];
        end
    end
    
     if(yj<=290 && xj<=462 && sum(bcw(:) == j) == 0)
        if((yj==275 && xj<=424.5 && xj>=412.5) || (xj*2/3==yj) || (-424.5+xj*2/5==yj))
            bcw=[bcw j];
        end
     end
    
      if(yk<=290 && xk<=462 && sum(bcw(:) == k) == 0)
        if((yk==275 && xk<=424.5 && xk>=412.5) || (xk*2/3==yj) || (-424.5+xk*2/5==yk))
            bcw=[bcw k];
        end
      end
      
      
      
    %     проверяем принадлжеит узел участку где водздух
      if(xi>=462 && yi~=0 && sum(bca(:) == i) == 0)
        if((yi==293 && xi<=484.5 && xi>=469.5) || (yi==270.5 && xi<=541.5 && xi>=529.5) || (-5/9*xi+1871.75==yi)|| (-0.5*xi+535.25==yi))
            bca=[bca i];
        end
      end
    
      if(xj>=462 && yj~=0 && sum(bca(:) == j) == 0)
        if((yj==293 && xj<=484.5 && xj>=469.5) || (yj==270.5 && xj<=541.5 && xj>=529.5) || (-5/9*xj+1871.75==yj)|| (-0.5*xj+535.25==yj))
            bca=[bca j];
        end
      end
     
      if(xk>=462 && yk~=0 && sum(bca(:) == k) == 0)
        if((yk==293 && xk<=484.5 && xk>=469.5) || (yk==270.5 && xk<=541.5 && xk>=529.5) || (-5/9*xk+1871.75==yk)|| (-0.5*xk+535.25==yk))
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


