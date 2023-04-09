E=2*10^11;
S=0.0001;
nodes=14;
elems=25;
force=10^(3);
forms=[1 -1;-1 1];
nodesx=[12,15,17,15,10,8,10,20,20,5,5,3,0,0];
nodesy=[5.63999987,5.0999999,4.73999977,0.899999976,1.79999995,5.63999987,6,4.19999981,0,0.899999976,5.0999999,4.73999977,0,4.19999981];
F_local=[0, -force,0 ,-force,0,-force,0,0,0,0,0, -force,0 ,-force,0, -force,0,0,0, 0,0 ,-force,0, -force,0, 0,0, -force];
K_global_for_system_temp=zeros(4,2*nodes);
K_global_for_system_on_each_iteration=zeros(2*nodes);
K_global_for_system=zeros(2*nodes);
U=[];
% ������� ����� � ������ ������� �������� (1 ���������� k ���� ������� - ����� ���� � ������� ��-��,2 ���������� k ���� ������� - ����� ���� � ������ ��-��,)
assoc=[[1,2];[2,3];[3,4];[5, 4];[6, 5];[6, 7];[7, 1];[1, 5];[5, 7];[4, 1];[4, 2];[3, 8];[9, 3];[9, 8];[4, 9];[10,5];[10, 11];[12, 11];[13, 12];[13, 10];[10,  6];[12, 10];[11,  6];[14, 12];[13, 14]];
%c��������� ������� ���������� (����������) -� ���������� ��� ��������
for el=1:25
%   ���������� ������ ������ � ����� ���������������� ��������
    i=assoc(el,1);
    j=assoc(el,2);
    length=sqrt((nodesx(j)-nodesx(i))^2+(nodesy(j)-nodesy(i))^2);
    l_ij=(nodesx(j)-nodesx(i))/length;
    m_ij=(nodesy(j)-nodesy(i))/length;
%     �������  ���������� � �� ����������� � ������ ����� �������
    K_local=E*S/length*forms;
    T=[[l_ij m_ij 0 0];[0 0 l_ij m_ij]];
%     ������� ���������� ��� �������� � ���������� ��, � �� ����������� � �������
    K_for_elem_in_global = T'*K_local*T
%     ���������� ������� ��� � ���������� ��� ���� �������
% ����� ���� ���-�� ��������� ����������������� ��������� � ��� ����� 14
% ����� � � ������ �� 2 ������� ������� ---> ������� ���������� ���������
% ����� ����������� 28 �� 28, � ������� ��������� ��������� �����
% ����������� 4 �� 4 , ���� �� �������� � ��������� ��������� �������
% �������������� ������� ��������� � ������� 28 �� 28 ������� �����
% �������� ���������, �� ����� ���� ��� �������������� ����������� ����
% ������� (���� ������� �� ������ ��� ������ ������ ���� ����� �� ������ ��� ������ ����) ��� ����� ������ ������� �� ������� ��������� x ���������� � y
% ���������� ����������� � ����� �������������� ����� ������������ ���
% ������� ��-�� (�.�. ������ ��� ���������� �� 4 ������� �������):
    k_x=2*i-1;
    k_y=2*i;
    m_x=2*j-1;
    m_y=2*j;
    K_global_for_system_temp(1,k_x)=1;
    K_global_for_system_temp(2,k_y)=1;
    K_global_for_system_temp(3,m_x)=1;
    K_global_for_system_temp(4,m_y)=1;
    K_global_for_system_on_each_iteration=K_global_for_system_temp'*K_for_elem_in_global*K_global_for_system_temp;
   
    K_global_for_system=K_global_for_system+K_global_for_system_on_each_iteration;
    
    
    
    K_global_for_system_temp=zeros(4,2*nodes);
    K_global_for_system_on_each_iteration=zeros(2*nodes);
end
%     ����� � ������� ���������� ��� ���� �� ��������� ������� - ��� �����
%     ������� �������������� ������� �������, � �������������� ������������
%     �������� ������� ������� 1 (����� �������� �� ������� ���������� ��� 1 �������� ��� ��� ����� �� ��������� 1, � ��� ��������� ���� - ������� ���������)
    K_global_for_system(25,:) = 0;
    K_global_for_system(:,25) = 0;
    K_global_for_system(25,25) = 1;
    K_global_for_system(26,:) = 0;
    K_global_for_system(:,26) = 0;
    K_global_for_system(26,26) = 1;
    K_global_for_system(17,:) = 0;
    K_global_for_system(:,17) = 0;
    K_global_for_system(17,17) = 1;
    K_global_for_system(18,:) = 0;
    K_global_for_system(:,18) = 0;
    K_global_for_system(18,18) = 1;
    M=inv(K_global_for_system)
    U=linsolve(K_global_for_system,F_local');
    F=zeros(25,1);
%     ������� ���������� ������ � ������ �������
for el=1:elems
    i=assoc(el,1);
    j=assoc(el,2);
    k_x=2*i-1;
    k_y=2*i;
    m_x=2*j-1;
    m_y=2*j;
    length=sqrt((nodesx(j)-nodesx(i))^2+(nodesy(j)-nodesy(i))^2);
    new_length=sqrt((nodesx(j)+U(m_x)-nodesx(i)-U(k_x))^2+(nodesy(j)+U(m_y)-nodesy(i)-U(k_y))^2);
    dL=new_length-length;
    F(el)=E*S*dL/length;
end
Ux=zeros(14,1);
Uy=zeros(14,1);
for i=1:14
   Ux(i,1)=U(2*i-1);
end
for i=1:14
   Uy(i)=U(2*i); 
end
Um=zeros(14,2);
Um(:,1)=Ux;
Um(:,2)=Uy;
new=nodesy'
% vtkwrite('Temp.vtk', 'polydata', 'lines', new_nodesx, new_nodesy, zeros(1,14))