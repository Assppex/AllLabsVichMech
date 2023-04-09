t=[0:0.01:0.5];
x=[0:0.1:1];
dt=0.01;
h=0.1;
n=(1/0.1)+1;
k=0.5/0.01+1;
u=zeros(k,n);
u1=zeros(k,n);
u_dop=zeros(n,1);
%Заполняю ГУ
for i=1:n
     u(k,i)=((i-1)*h)^2*cos(pi*(i-1)*h);
     u1(k,i)=((i-1)*h)^2*cos(pi*(i-1)*h);
%    u(k,i)=2*((i-1)*h)*((i-1)*h+1)+0.3;
%    u1(k,i)=2*((i-1)*h)*((i-1)*h+1)+0.3;
end
for i=k:-1:1
     u(i,1)=0.5*((k-i)*dt);
     u1(i,1)=0.5*((k-i)*dt);
%    u(i,1)=0.3;
%    u1(i,1)=0.3;
end
for i=k:-1:1
   u(i,n)=(k-i)*dt-1;
   u1(i,n)=(k-i)*dt-1;
%    u(i,n)=4.3+(i-1)*dt;
%    u1(i,n)=4.3+(i-1)*dt;
end
% Первая строка перемещений
for i=2:n-1
    u1(k-1,i)=((i*h)^2*cos(pi*i*h))+((i*h)^2*(i*h+1))*dt+(dt^2/(2*h^2))*((((i-1)*h)^2*cos(pi*(i-1)*h))-2*(((i)*h)^2*cos(pi*i*h))+(((i+1)*h)^2*cos(pi*(i+1)*h)));
    u(k-1,i)=((i*h)^2*cos(pi*i*h))+((i*h)^2*(i*h+1))*dt+(dt^2/(2*h^2))*((((i-1)*h)^2*cos(pi*(i-1)*h))-2*(((i)*h)^2*cos(pi*i*h))+(((i+1)*h)^2*cos(pi*(i+1)*h)));
%    ((i-1)*h)^2*cos(pi*(i-1)*h)+(((i-1)*h)^2*((i-1)*h+1))*dt+(dt^2/(2*h^2))*(((i-1)*h)^2*cos(pi*(i-1)*h)-2*((i)*h)^2*cos(pi*(i)*h)+((i+1)*h)^2*cos(pi*(i+1)*h));
%   u1(k-1,i)=(2*((i-1)*h)*((i-1)*h+1)+0.3)+2*sin((i-1)*h)*dt+(dt^2/(2*h^2))*((2*((i-1)*h)*((i-1)*h+1)+0.3)-2*(2*((i)*h)*((i)*h+1)+0.3)+(2*((i+1)*h)*((i+1)*h+1)+0.3));
%   u(k-1,i)=(2*((i-1)*h)*((i-1)*h+1)+0.3)+2*sin((i-1)*h)*dt+(dt^2/(2*h^2))*((2*((i-1)*h)*((i-1)*h+1)+0.3)-2*(2*((i)*h)*((i)*h+1)+0.3)+(2*((i+1)*h)*((i+1)*h+1)+0.3));
end
% Решение задачи по явному мкр
for j=k-1:-1:2
   for i=2:n-1
       u(j-1,i)=(dt^2/(h^2))*(u(j,i+1)-2*u(j,i)+u(j,i-1))+2*u(j,i)-u(j+1,i);
   end
end
% Решение задачи по неявному мкр
A=-1/(h^2);
B=(h^2+2*dt^2)/((dt^2)*(h^2));
C=-1/(h^2);
% Заполняю матрицу коэффициентов
matrix=zeros(n);
for i=2:n-1
   matrix(i,i)=B;
   matrix(i,i-1)=-A;
   matrix(i,i+1)=-C;
end
matrix(1,1)=1;
matrix(n,n)=1;


for j=k-1:-1:2
    del=zeros(n,1);
    lam=zeros(n,1);
    F=zeros(n,1);
      
      P(1)=-C/B;
% реализуем сам метод прогонки
%  прямой ход
    for i=1:n
       F(i)=(2*u1(j,i))/(dt^2)-u1(j+1,i)/(dt^2); 
    end
    Q(1)=F(1)/B;
    for i=2:n
      P(i)=-C/(B+A*P(i-1));
      Q(i)=(F(i)-A*Q(i-1))/(B+A*P(i-1));
    end
% обратный ход
    for i=n-1:-1:2
       u1(j-1,i)=P(i)*u1(j-1,i+1)+Q(i); 
    end
end
temp=zeros(1,n);
temp1=zeros(1,n);
for i=1:(k-1)/2
    temp=u(i,:);
    u(i,:)=u(k-i+1,:);
    u(k-i+1,:)=temp;
    temp1=u1(i,:);
    u1(i,:)=u1(k-i+1,:);
    u1(k-i+1,:)=temp1;
end
figure()
hold on
grid on
surf(x,t,u);
title('u(x,t) (неявный метод)')
xlabel('x')
ylabel('t')
zlabel('u')
figure()
surf(x,t,u1);
title('u(x,t) (явный метод)')
xlabel('x')
ylabel('t')
zlabel('u')
figure()
hold on
plot(x,u(16,:));
plot(x,u1(16,:));
title('Зависимость u(x) в момент времени 0.15 (оба метода)')
xlabel('x')
ylabel('u')
legend('ЯМКР','НМК')
figure()
hold on
plot(x,u(k,:));
plot(x,u1(k,:));
title('Зависимость u(x) в момент времени 0.5 (оба метода)')
xlabel('x')
ylabel('u')
legend('ЯМКР','НМКP')
figure()
hold on
plot(x,u(26,:));
plot(x,u1(26,:));
title('Зависимость u(x) в момент времени 0.25 (оба метода)')
xlabel('x')
ylabel('u')
legend('ЯМКР','НМКP')