t=[0:0.001:0.01];
x=[0:0.1:0.6];
dt=0.001;
h=0.1;
n=0.6/0.1+1;
n=floor(n)+1;
k=0.01/0.001+1;
T=zeros(k,n);
T1=zeros(k,n);
%Заполняю ГУ
for i=1:n
   T(k,i)=(((i-1)*0.1-0.2)*((i-1)*0.1+1))+0.2;
   T1(k,i)=(((i-1)*0.1-0.2)*((i-1)*0.1+1))+0.2;
end
for i=k:-1:1
   T(i,1)=6*(k-i)*dt;
   T1(i,1)=6*(k-i)*dt; 
end
for i=1:k
   T(i,n)=0.84;
   T1(i,n)=0.84;
end
% Решение задачи по явному мкр
for j=k:-1:2
   for i=2:n-1
       T(j-1,i)=(dt/(h^2))*(T(j,i+1)-2*T(j,i)+T(j,i-1))+T(j,i);
   end
end
% Решение задачи по неявному мкр
A=1/(h^2);
B=(h^2+2*dt)/(dt*(h^2));
C=1/(h^2);
% Заполняю матрицу коэффициентов
matrix=zeros(n);
for i=2:n-1
   matrix(i,i)=B;
   matrix(i,i-1)=-A;
   matrix(i,i+1)=-C;
end
matrix(1,1)=1;
matrix(n,n)=1;


for j=k:-1:2
    P=zeros(n,1);
    Q=zeros(n,1);
    F=zeros(n,1);
    P(1)=0;
    F(1)=T1(j-1,1);
    F(n)=T1(j-1,n);
    Q(1)=T1(k-1,1);
% реализуем сам метод прогонки
%  прямой ход
    for i=2:n-1
       F(i)=T1(j,i)/dt; 
    end
    for i=2:n
      P(i)=C/(B-A*P(i-1));
      Q(i)=(F(i)+A*Q(i-1))/(B-A*P(i-1));
    end
% обратный ход
    for i=n-1:-1:2
       T1(j-1,i)=P(i)*T1(j-1,i+1)+Q(i); 
    end
end
temp=zeros(1,n);
temp1=zeros(1,n);
for i=1:(k-1)/2
    temp=T(i,:);
    T(i,:)=T(k-i+1,:);
    T(k-i+1,:)=temp;
    temp1=T1(i,:);
    T1(i,:)=T1(k-i+1,:);
    T1(k-i+1,:)=temp1;
end
figure()
hold on
surf(x,t,T);
surf(x,t,T1);
title('T(x,t) (неявный и явный метод)')
xlabel('x')
ylabel('t')
zlabel('T')
figure()
surf(x,t,T1);
title('T(x,t) (явный метод)')
xlabel('x')
ylabel('t')
zlabel('T')
figure()
hold on
plot(x,T(2,:));
plot(x,T1(2,:));
title('Зависимость T(x) в момент времени 0.001 (оба метода)')
xlabel('x')
ylabel('t')
legend('НМК','ЯМКР')
figure()
hold on
plot(x,T(k-7,:));
plot(x,T1(k-7,:));
title('Зависимость T(x) в момент времени 0.003 (оба метода)')
xlabel('x')
ylabel('t')
legend('НМК','ЯМКР')
figure()
hold on
plot(x,T(k-6,:));
plot(x,T1(k-6,:));
title('Зависимость T(x) в момент времени 0.004 (оба метода)')
xlabel('x')
ylabel('t')
legend('НМК','ЯМКР')