h=0.2;
x=[0:h:1];
y=[0:h:1];
N=size(x,2);
M=size(y,2);
kol_it=[];
eps=0.01;
%Подбор w и решение
for w=0.1:0.1:1.9
    temp=1; 
    u=zeros(M,N);
    kol=0;
    %Граничные условия
    for i=1:N
        u(M,i)=0;
        u(1,i)=20*(1-((i-1)*h)^2);
    end
    for j=M:-1:1
        u(j,1)=20*((M-j)*h);
        u(j,N)=30*sqrt((M-j)*h)*(1-(M-j)*h);
    end
    while(temp>eps)
        for i=2:N-1
           for j=M-1:-1:2
               u(j,i)=0.25*(u(j,i-1)+u(j,i+1)+u(j-1,i)+u(j+1,i));
           end
        end
        t=zeros(M,N);
        t(1,:)=u(1,:);
        t(M,:)=u(M,:);
        t(:,1)=u(:,1);
        t(:,N)=u(:,N);
        for i=2:N-1
           for j=M-1:-1:2
               t(j,i)=0.25*(u(j,i-1)+u(j,i+1)+u(j-1,i)+u(j+1,i));
           end
        end
        for i=1:N
            for j=M:-1:1
                u1(j,i)=u(j,i)+w*(t(j,i)-u(j,i));
            end
        end
        temp=norm(u1-u,inf);
        u=u1;
        kol=kol+1;
    end
    kol_it(end+1)=kol;
end
figure()
plot(0.1:0.1:1.9,kol_it);
xlabel('w(омега)');
ylabel('количество итераций');
title('Зависимость количесвтва итераций от w');
%Оптимальное омега-1.8 и 1.9 - для него кол-во итераций 10
figure()
surf(x,y,u);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Зависимость u(x,y)');
