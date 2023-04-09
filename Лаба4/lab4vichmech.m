%конечный элемент - линейный гексаэдр
X=zeros(8);
x_nde=[];
y_nde=[];
z_nde=[];
%1
x_nde(end+1)=-1;
y_nde(end+1)=-1;
z_nde(end+1)=-1;
%2
x_nde(end+1)=-1;
y_nde(end+1)=-1;
z_nde(end+1)=1;
%3
x_nde(end+1)=-1;
y_nde(end+1)=1;
z_nde(end+1)=-1;
%4
x_nde(end+1)=-1;
y_nde(end+1)=1;
z_nde(end+1)=1;
%5
x_nde(end+1)=1;
y_nde(end+1)=-1;
z_nde(end+1)=-1;
%6
x_nde(end+1)=1;
y_nde(end+1)=-1;
z_nde(end+1)=1;
%7
x_nde(end+1)=1;
y_nde(end+1)=1;
z_nde(end+1)=-1;
%8
x_nde(end+1)=1;
y_nde(end+1)=1;
z_nde(end+1)=1;

X(:,1)=1;
for i=1:8
    X(i,2)=x_nde(i);
    X(i,3)=y_nde(i);
    X(i,4)=z_nde(i);
    X(i,5)=x_nde(i)*y_nde(i);
    X(i,6)=x_nde(i)*z_nde(i);
    X(i,7)=y_nde(i)*z_nde(i);
    X(i,8)=x_nde(i)*y_nde(i)*z_nde(i);
end
A=zeros(8);
A=inv(X);
syms x y z;
P=[1;x;y;z;x*y;x*z;y*z;x*y*z];
P=P.';
N=P*A

f(1)={@(x,y,z) (x.*y)/8 - y./8 - z./8 - x./8 + (x.*z)./8 + (y.*z)./8 - (x.*y.*z)./8 + 1/8};
f(2)={@(x,y,z) z./8 - y./8 - x./8 + (x.*y)./8 - (x.*z)./8 - (y.*z)./8 + (x.*y.*z)./8 + 1/8};
f(3)={@(x,y,z) y./8 - x./8 - z./8 - (x.*y)./8 + (x.*z)./8 - (y.*z)./8 + (x.*y.*z)./8 + 1/8};
f(4)={@(x,y,z) y./8 - x./8 + z./8 - (x.*y)./8 - (x.*z)./8 + (y.*z)./8 - (x.*y.*z)./8 + 1/8};
f(5)={@(x,y,z) x./8 - y./8 - z./8 - (x.*y)./8 - (x.*z)./8 + (y.*z)./8 + (x.*y.*z)./8 + 1/8};
f(6)={@(x,y,z) x./8 - y./8 + z./8 - (x.*y)./8 + (x.*z)./8 - (y.*z)./8 - (x.*y.*z)./8 + 1/8};
f(7)={@(x,y,z) x./8 + y./8 - z./8 + (x.*y)./8 - (x.*z)./8 - (y.*z)./8 - (x.*y.*z)./8 + 1/8};
f(8)={@(x,y,z) x./8 + y./8 + z./8 + (x.*y)./8 + (x.*z)./8 + (y.*z)./8 + (x.*y.*z)./8 + 1/8};
%  1 1 -1 - 7 узел; -1 1 -1 - 3 узел; -1 1 1 - 4 узел; 1 1 1 - 8 узел;
% 1 -1 -1 - 5 узел; 1 -1 1 - 6 узел; -1 -1 1 - 2 узел; -1 -1 -1 - 1 узел
for i=1:8
    figure()
        hold on
        % fill3 ([координаты Х вершин грани],[координаты Y вершин грани],[координаты Z вершин грани],[цвет определяется значениями iой функции формы в узлах (те в верщинах граней), те задавая так цвета мы получаем проекцию значений функции формы на исходный конечный элемент (линейный гексаэдр)])
        fill3([1 1 -1 -1],[1 1 1 1],[-1 1 1 -1],[f{i}(x_nde(7),y_nde(7),z_nde(7)),f{i}(x_nde(8),y_nde(8),z_nde(8)),f{i}(x_nde(4),y_nde(4),z_nde(4)),f{i}(x_nde(3),y_nde(3),z_nde(3))]);
        fill3([-1 -1 -1 -1],[1 1 -1 -1],[-1 1 1 -1],[f{i}(x_nde(3),y_nde(3),z_nde(3)),f{i}(x_nde(4),y_nde(4),z_nde(4)),f{i}(x_nde(2),y_nde(2),z_nde(2)),f{i}(x_nde(1),y_nde(1),z_nde(1))]);
        fill3([-1 -1 1 1],[-1 -1 -1 -1],[-1 1 1 -1],[f{i}(x_nde(1),y_nde(1),z_nde(1)),f{i}(x_nde(2),y_nde(2),z_nde(2)),f{i}(x_nde(6),y_nde(6),z_nde(6)),f{i}(x_nde(5),y_nde(5),z_nde(5))]);
        fill3([1 1 1 1],[-1 -1 1 1],[-1 1 1 -1],[f{i}(x_nde(5),y_nde(5),z_nde(5)),f{i}(x_nde(6),y_nde(6),z_nde(6)),f{i}(x_nde(8),y_nde(8),z_nde(8)),f{i}(x_nde(7),y_nde(7),z_nde(7))]);
        fill3([-1 -1 1 1],[1 -1 -1 1],[1 1 1 1],[f{i}(x_nde(4),y_nde(4),z_nde(4)),f{i}(x_nde(2),y_nde(2),z_nde(2)),f{i}(x_nde(6),y_nde(6),z_nde(6)),f{i}(x_nde(8),y_nde(8),z_nde(8))]);
        fill3([-1 -1 1 1],[1 -1 -1 1],[-1 -1 -1 -1],[f{i}(x_nde(3),y_nde(3),z_nde(3)),f{i}(x_nde(1),y_nde(1),z_nde(1)),f{i}(x_nde(5),y_nde(5),z_nde(5)),f{i}(x_nde(7),y_nde(7),z_nde(7))]);
        colorbar;
        xlabel('x');
        ylabel('y');
        zlabel('z');
        view(125,21);
        title(sprintf('Проекция значений %d-ой ([%d,%d,%d]) функции формы на конечный элемент',i,x_nde(i),y_nde(i),z_nde(i)));
end