clc;
clear all;
% коэффициенты Пуассона
v1=0.14;
v2=0.14;
v3=0.22;
% Плотности игэ-1 и игэ-2 (
% основание не подвержено нагрузки собственного веса)
den1=2500;
den2=2500;
denwat=1000;
% Модули Юнга
E1=27500e6;
E2=32500e6;
E3=17e6;
% уровен воды
water_level=280;
% Получаем узлы
x=[];
y=[];

nds=293;
els=511;

assoc=[];

Temp=readtable('xydata.txt');
x=Temp(:,2);
y=Temp(:,3);
x=table2array(x);
y=table2array(y);
% матрица в которой номер строки - номер элемента, а значения в строке -
% номера узлов, которых состоит элемент с этим номером
Temp=readtable('assoc.txt')

a=Temp(:,2);
a=table2array(a);
b=Temp(:,3);
b=table2array(b);
c=Temp(:,4);
c=table2array(c);

for i=1:els
    temp=[a(i),b(i),c(i)];
    assoc(i,:)=temp;
end

% Приступаем к решению
%  1) Считаем матрицу жесткости элемента
%       1.1) Ищем матрицу градиентов, воспользуемся функциями форм
%       изопарамтерического элемента и якобианом
B_iso_derivatives=[[-1,0,1];[-1,1,0]];
% Идет пересчет по элементам,тк матрица градиентов для элемента
nums_of_nodes_in_element=zeros(1,3);
 K=zeros(2*nds);
 kol=0;
 kol1=0;
for i=1:els
    nums_of_nodes_in_element=assoc(i,:);
    xi=x(nums_of_nodes_in_element(1));
    yi=y(nums_of_nodes_in_element(1));
    
    xj=x(nums_of_nodes_in_element(2));
    yj=y(nums_of_nodes_in_element(2));
    
    xk=x(nums_of_nodes_in_element(3));
    yk=y(nums_of_nodes_in_element(3));
    
    coordinates_of_elem_nds=[[xi,yi];[xj,yj];[xk,yk]];
    Jacobi=B_iso_derivatives*coordinates_of_elem_nds;
%     далее считаем столбцы матрицы градиентов уже не в изопараметрической,
%     а в актуальной СК
    inv_Jacobi = inv(Jacobi);
%     пересчитываем первые два столбца (точнее ненулевые элемениы из них)
% матрица градиентов элемента в актульной ск:
    tmp=zeros(2,3);
    B_real=zeros(3,6);
    for j=1:3
        tmp(:,j) = inv_Jacobi*B_iso_derivatives(:,j);
    end
        B_real(:,1)=[tmp(1,1);0;tmp(2,1)];
        B_real(:,2)=[0;tmp(2,1);tmp(1,1)];
        
        B_real(:,3)=[tmp(1,2);0;tmp(2,2)];
        B_real(:,4)=[0;tmp(2,2);tmp(1,2)];
        
        B_real(:,5)=[tmp(1,3);0;tmp(2,3)];
        B_real(:,6)=[0;tmp(2,3);tmp(1,3)];
%      матрица B для элемента найдена,теперь ищем матрицу упругих
%      характеристик, тут надо понять, к какому материалу отсносится
%      материал: игэ-1, игэ-2, основание

        
%      проверяем, принадлежит ли элемент бетону

        if(yi<0 || yj<0 || yk<0)
            kol1=kol1+1;
            D=E3*(1-v3)/((1+v3)*(1-2*v3))*[[1,v3/(1-v3),0];[v3/(1-v3),1,0];[0,0,(1-2*v3)/(2*(1-v3))]];
            
%             проверяем принадлежит ли к ИГЭ-1
        elseif((xi>=0 && xi<=12.6) || (xj>=0 && xj<=12.6) || (xk>=0 && xk<=12.6))
            D=E2*(1-v2)/((1+v2)*(1-2*v2))*[[1,v2/(1-v2),0];[v2/(1-v2),1,0];[0,0,(1-2*v2)/(2*(1-v2))]];
            
            kol=kol+1;
%             проверяем принадлежит ли к ИГЭ-2
        else
            D=E1*(1-v1)/((1+v1)*(1-2*v1))*[[1,v1/(1-v1),0];[v1/(1-v1),1,0];[0,0,(1-2*v1)/(2*(1-v1))]];  
            kol=kol+1;
            
                
        end
        
%         теперь можем считать матрицу ke
    
%         ищем площадь элемента по формуле определения площади через
%         матрицу


% через Герона
%         l_a = sqrt((xj-xi)^2+(yj-yi)^2);
%         l_b = sqrt((xj-xk)^2+(yj-yk)^2);
%         l_c = sqrt((xk-xi)^2+(yk-yi)^2);
%         p=(l_a+l_b+l_c)/2;
%         Square_of_element=sqrt(p*(p-l_a)*(p-l_b)*(p-l_c));

%         ищем площадь элемента по формуле определения площади через
%   
        Tmp=zeros(3);
        Tmp(:,3)=1;
          
        Tmp(1,1)=xi;
        Tmp(2,1)=xj;
        Tmp(3,1)=xk;
          
        Tmp(1,2)=yi;
        Tmp(2,2)=yj;
        Tmp(3,2)=yk;
          
        Square_of_element=abs(det(Tmp)/2);

        ke=B_real'*D*B_real*Square_of_element;
        
        A = zeros(6,2*nds);
%         степени свободы по х - имеют нечетные номера (2*номер узла -1), по
%         y - четные номера (2*номер узла)
        A(1,2*nums_of_nodes_in_element(1)-1)=1;
        A(2,2*nums_of_nodes_in_element(1))=1;
        
        A(3,2*nums_of_nodes_in_element(2)-1)=1;
        A(4,2*nums_of_nodes_in_element(2))=1;
        
        A(5,2*nums_of_nodes_in_element(3)-1)=1;
        A(6,2*nums_of_nodes_in_element(3))=1;
        K_tmp=A'*ke*A;
        K=K+K_tmp;
end

fe=zeros(2*nds,1);
% гидростатическое давление
% y=0
elems_under_force=[];
for i=1:nds
   xt=x(i);
   yt=y(i);
   if(yt==0 && xt<=0)
      elems_under_force = [elems_under_force i]; 
   end
end
% сортировка
for i=1:size(elems_under_force,2)-1
    for j=1:size(elems_under_force,2)-i
        if(x(elems_under_force(j))>x(elems_under_force(j+1)))
            tmp =  elems_under_force(j+1);
            elems_under_force(j+1) = elems_under_force(j);
            elems_under_force(j)=tmp;  
        end
    end
end
% подсчет силы
for i=1:size(elems_under_force,2)-1
   l = sqrt((x(elems_under_force(i))-x(elems_under_force(i+1)))^2+(y(elems_under_force(i))-y(elems_under_force(i+1)))^2)
   fe(2*elems_under_force(i))=fe(2*elems_under_force(size(elems_under_force,2)))-denwat*9.8*(water_level)*l;
end
% fe(2*elems_under_force(size(elems_under_force,2)))=fe(2*elems_under_force(size(elems_under_force,2)))-denwat*9.8*(water_level)*l;


% y=2/3*x
elems_under_force=[];
for i=1:nds
   xt=x(i);
   yt=y(i);
   if(yt<=280 && yt>=0 && xt==0)
      elems_under_force = [elems_under_force i]; 
   end
end
% сортировка
for i=1:size(elems_under_force,2)-1
    for j=1:size(elems_under_force,2)-i
        if(x(elems_under_force(j))>x(elems_under_force(j+1)))
            tmp =  elems_under_force(j+1);
            elems_under_force(j+1) = elems_under_force(j);
            elems_under_force(j)=tmp;  
        end
    end
end
% подсчет силы
for i=1:size(elems_under_force,2)-1
   l = sqrt((x(elems_under_force(i))-x(elems_under_force(i+1)))^2+(y(elems_under_force(i))-y(elems_under_force(i+1)))^2)
   fe(2*elems_under_force(i)-1)=fe(2*elems_under_force(size(elems_under_force,2)))-denwat*9.8*(water_level-y(elems_under_force(i)))*l;
  
end
%  fe(2*elems_under_force(size(elems_under_force,2)))=fe(2*elems_under_force(size(elems_under_force,2)))-denwat*9.8*(water_level-y(elems_under_force(size(elems_under_force,2))))*l*sin(pi/2-atan(2/3));
%  fe(2*elems_under_force(size(elems_under_force,2)))=fe(2*elems_under_force(size(elems_under_force,2)))-denwat*9.8*(water_level-y(elems_under_force(size(elems_under_force,2))))*l*cos(pi/2-atan(2/3));
 
 




% считаем вектор внешних нагрузок
%  сначала учитываем силу тяжести

kol3=0;

for i=1:nds
%    надо нйти все элементы, которые окружают эелмент, ясно, 
% что все эти элементы содержат этот элемент, однако собственный вес нужно
% учитывать только для самой плотины, не для основания
x_c = x(i);
y_c = y(i);
%     гарантированно имеем дело не с основанием, потому что основание
%     находится ниже нуля
    if(y_c>0)
       elements_around=[];
       for j=1:els
           i
           assoc(j,1)
           assoc(j,2)
           assoc(j,3)
           if(i==assoc(j,1) || i==assoc(j,2) || i==assoc(j,3))
               elements_around=[elements_around j];
           end
       end
       V=0;
       for j=1:size(elements_around,2)
           xi=x(assoc(j,1));
           xj=x(assoc(j,2));
           xk=x(assoc(j,3));
           
           yi=y(assoc(j,1));
           yj=y(assoc(j,2));
           yk=y(assoc(j,3));
        Tmp=zeros(3);
        Tmp(:,3)=1;
          
        Tmp(1,1)=xi;
        Tmp(2,1)=xj;
        Tmp(3,1)=xk;
          
        Tmp(1,2)=yi;
        Tmp(2,2)=yj;
        Tmp(3,2)=yk;
          
        V=V+abs(det(Tmp)/2);
       end
       if(x_c>=0 && y_c>=0 && x_c<=12.6)
           cden=den2;
       elseif(x_c>=0 && y_c>=0 && x_c>=12.6)
           cden=den1;
       end
       fe(2*i)=fe(2*i)-cden*9.8*1/3*V;
    end
end

% for i=1:els
%   a=assoc(i,1);
%   b=assoc(i,2);
%   c=assoc(i,3);
%     
%   x1=x(a)
%   y1=y(a);
%   
%   x2=x(a)
%   y2=y(a);
%   
%   x3=x(a)
%   y3=y(a);
% 
%   
%   
%   if(y1>=0 || y2>=0 || y3>=0)
%       
%       Tmp=zeros(3);
%       Tmp(:,3)=1;
%           
%       Tmp(1,1)=x1;
%       Tmp(2,1)=x2;
%       Tmp(3,1)=x3;
%           
%       Tmp(1,2)=y1;
%       Tmp(2,2)=y2;
%       Tmp(3,2)=y3;
%             
%       S=abs(det(Tmp)/2);
%         if(xi>=0 && yi>=0 && xi<=12.6)
%            c_den=den2;
%        elseif(xi>=0 && yi>=0 && xi>=12.6)
%            c_den=den1;
%        end
%       fe(2*a)=fe(2*a)-c_den*9.8*S;
%       fe(2*b)=fe(2*b)-c_den*9.8*S;
%       fe(2*b)=fe(2*b)-c_den*9.8*S;
%   end
% end



% теперь считаем силу гидростатического давления, заметим, что она давит
% только на самую левую грань + на левую часть основания, тут у силы будет
% 2 компаненты по х и у, надо домножать на угол и считать проекции из
% геометрии получим угол между силой давления и наклонными поверхностями:
% beta = pi/2 - alpha, где alpha = arctg(1/m),m - дано в чертеже
nums_of_surround_nds1=[];
%    сначала рассматриваем горизонатальные участки:один - это основание на глубине с у=0, 
% другое - горизонталь левой части поверхности  с у = 275


% горизонталь левой части поверхности  с у = 275
%     ищем узлы заполняющие рассматриваемую поверхность

 









% Ищем номера узлов, которые подлежат границчным условиям
% левый край основая имеет координату x=-1028,4
bc_left=[];
for j=1:nds
    x_c1 = x(j);
%     ищем узлы заполняющие рассматриваемую поверхность
    if((x_c1==-28.512500800000000))
        bc_left=[bc_left j];
    end
end
bc_right=[];
for j=1:nds
    x_c1 = x(j);
%     ищем узлы заполняющие рассматриваемую поверхность
    if(x_c1==65.487503100000000)
        bc_right=[bc_right j];
    end
end
bc_bottom=[];
for j=1:nds
    y_c1 = y(j);
%     ищем узлы заполняющие рассматриваемую поверхность
    if((y_c1==-47))
        bc_bottom=[bc_bottom j];
    end
end
% учитываем гу в матрице К
for j=1:size(bc_bottom,2)
%     K(2*bc_bottom(j)-1,:)=0;
%     K(:,2*bc_bottom(j)-1)=0;
%     K(2*bc_bottom(j)-1,2*bc_bottom(j)-1)=1;
    
    K(2*bc_bottom(j),:)=0;
    K(:,2*bc_bottom(j))=0;
    K(2*bc_bottom(j),2*bc_bottom(j))=1;
end

for j=1:size(bc_left,2)
    K(2*bc_left(j)-1,:)=0;
    K(:,2*bc_left(j)-1)=0;
    K(2*bc_left(j)-1,2*bc_left(j)-1)=1;
    
%     K(2*bc_left(j),:)=0;
%     K(:,2*bc_left(j))=0;
%     K(2*bc_left(j),2*bc_left(j))=1;
end

for j=1:size(bc_right,2)
    K(2*bc_right(j)-1,:)=0;
    K(:,2*bc_right(j)-1)=0;
    K(2*bc_right(j)-1,2*bc_right(j)-1)=1;
    
%     K(2*bc_right(j),:)=0;
%     K(:,2*bc_right(j))=0;
%     K(2*bc_right(j),2*bc_right(j))=1;
end
u=linsolve(K,fe);
ux=zeros(nds,1);
uy=zeros(nds,1);
for i=1:nds
   ux(i)=u(2*i-1);
   uy(i)=u(2*i);
end