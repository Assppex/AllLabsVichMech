clear all
clc
%Зададим параметры задачи
g = 9.8; %Ускорение свободного падения
rho = 1000; %Плотность жидкости
%Для плотины
%Для 1 матриала
E_1_1 = 5e9; %Модуль Юнга
nu_1_1 = 0.2; %Коэффициент Пуассона
rho_1 = 1500; %Плотность
%Для 2 материала
E_1_2 = 7e9; %Модуль Юнга
nu_1_2 = 0.2; %Коэффициент Пуассона
rho_2 = 1600; %Плотность
%Для основания
E_2 = 17e9; %Модуль Юнга
nu_2 = 0.2; %Коэффициент Пуассона
%Импортируем точки и элементы
[coordinates, N, nodes_1_1, nodes_1_2, nodes_2, BC_X_nodes, BC_Y_nodes, force_nodes, force_elements] = Input();
%Количество узлов
n = size(coordinates, 1);
%Находим элементы, которые относятся к плотине и к основанию
N_1_1 = [];
pointer_1_1 = 1;
N_1_2 = [];
pointer_1_2 = 1;
N_2 = [];
pointer_2 = 1;
for i = 1:size(N, 1)
count_1_1 = 0;
count_1_2 = 0;
for j = 1:size(N, 2)
for k = 1:size(nodes_1_1, 1)
if (N(i, j) == nodes_1_1(k))
count_1_1 = count_1_1 + 1;
break;
end
end
for k = 1:size(nodes_1_2, 1)
if (N(i, j) == nodes_1_2(k))
count_1_2 = count_1_2 + 1;
break;
end
end
end
if(count_1_1 ~= 3 && count_1_2 ~= 3)
N_2(pointer_2, :) = N(i, :);
pointer_2 = pointer_2 + 1;
elseif(count_1_1 ~= 3)
N_1_2(pointer_1_2, :) = N(i, :);
pointer_1_2 = pointer_1_2 + 1;
else
N_1_1(pointer_1_1, :) = N(i, :);
pointer_1_1 = pointer_1_1 + 1;
end
end
%Находим элементы с граничными условиями
BC_X_elements = [];
pointer_1 = 1;
BC_Y_elements = [];
pointer_2 = 1;
for i = 1:size(N, 1)
count_x = 0;
count_y = 0;
for j = 1:size(N, 2)
for k = 1:size(BC_X_nodes, 1)
if (N(i, j) == BC_X_nodes(k))
count_x = count_x + 1;
break;
end
end
for k = 1:size(BC_Y_nodes, 1)
if (N(i, j) == BC_Y_nodes(k))
count_y = count_y + 1;
break;
end
end
end
if(count_x == 2)
BC_X_elements(pointer_1, :) = N(i, :);
pointer_1 = pointer_1 + 1;
end
if(count_y == 2)
BC_Y_elements(pointer_2, :) = N(i, :);
pointer_2 = pointer_2 + 1;
end
end
%Используемые матрицы
K_1_1 = zeros(2 * n, 2 * n); %Матрица жесткости плотины 1 материала
K_1_2 = zeros(2 * n, 2 * n); %Матрица жесткости плотины 2 материала
K_2 = zeros(2 * n, 2 * n); %Матрица жесткости основания
F = zeros(2 * n, 1); %Вектор-столбец нагрузок системы
D_1_1 = Elastic_Matrix(E_1_1, nu_1_1); %Матрица упругих характеристик плотины 1 материала
D_1_2 = Elastic_Matrix(E_1_2, nu_1_2); %Матрица упругих характеристик плотины 2 материала
D_2 = Elastic_Matrix(E_2, nu_2); %Матрица упругих характеристик основания
epsilon = zeros(size(N, 1), 3);
sigma = zeros(size(N, 1), 3);
%Строим матрицу жесткости для плотины 1 материала
for i = 1:size(N_1_1, 1)
%Находим коэффициенты для функций форм
a = [coordinates(N_1_1(i, 2), 1) * coordinates(N_1_1(i, 3), 2) - coordinates(N_1_1(i, 3), 1) * coordinates(N_1_1(i, 2), 2);
coordinates(N_1_1(i, 3), 1) * coordinates(N_1_1(i, 1), 2) - coordinates(N_1_1(i, 1), 1) * coordinates(N_1_1(i, 3), 2);
coordinates(N_1_1(i, 1), 1) * coordinates(N_1_1(i, 2), 2) - coordinates(N_1_1(i, 2), 1) * coordinates(N_1_1(i, 1), 2)];
b = [coordinates(N_1_1(i, 2), 2) - coordinates(N_1_1(i, 3), 2);
coordinates(N_1_1(i, 3), 2) - coordinates(N_1_1(i, 1), 2);
coordinates(N_1_1(i, 1), 2) - coordinates(N_1_1(i, 2), 2)];
c = [coordinates(N_1_1(i, 3), 1) - coordinates(N_1_1(i, 2), 1);
coordinates(N_1_1(i, 1), 1) - coordinates(N_1_1(i, 3), 1);
coordinates(N_1_1(i, 2), 1) - coordinates(N_1_1(i, 1), 1)];
%Находим площадь треугольного элемента
S = zeros(3, 3);
for j = 1:3
S(j, 1) = 1;
for k = 2:3
S(j, k) = coordinates(N_1_1(i, j), k - 1);
end
end
A = det(S) / 2; % Площадь
%Строим матрицу градиентов
B = [b(1), 0, b(2), 0, b(3), 0;
0, c(1), 0, c(2), 0, c(3);
c(1), b(1), c(2), b(2), c(3), b(3)];
B = B / (2 * A);
%Строим матрицу жесткости конечного элемента
k_e = A * B' * D_1_1 * B;
%Строим глобальную матрицу жесткости
for n = 1:3
for j = 1:3
K_1_1(2 * N_1_1(i, n) - 1, 2 * N_1_1(i, j) - 1) = K_1_1(2 * N_1_1(i, n) - 1, 2 * N_1_1(i, j) - 1) + k_e(2 * n - 1, 2 * j - 1);
K_1_1(2 * N_1_1(i, n), 2 * N_1_1(i, j) - 1) = K_1_1(2 * N_1_1(i, n), 2 * N_1_1(i, j) - 1) + k_e(2 * n, 2 * j - 1);
K_1_1(2 * N_1_1(i, n) - 1, 2 * N_1_1(i, j)) = K_1_1(2 * N_1_1(i, n) - 1, 2 * N_1_1(i, j)) + k_e(2 * n - 1, 2 * j);
K_1_1(2 * N_1_1(i, n), 2 * N_1_1(i, j)) = K_1_1(2 * N_1_1(i, n), 2 * N_1_1(i, j)) + k_e(2 * n, 2 * j);
end
end
end
%Строим матрицу жесткости для плотины 2 материала
for i = 1:size(N_1_2, 1)
%Находим коэффициенты для функций форм
a = [coordinates(N_1_2(i, 2), 1) * coordinates(N_1_2(i, 3), 2) - coordinates(N_1_2(i, 3), 1) * coordinates(N_1_2(i, 2), 2);
coordinates(N_1_2(i, 3), 1) * coordinates(N_1_2(i, 1), 2) - coordinates(N_1_2(i, 1), 1) * coordinates(N_1_2(i, 3), 2);
coordinates(N_1_2(i, 1), 1) * coordinates(N_1_2(i, 2), 2) - coordinates(N_1_2(i, 2), 1) * coordinates(N_1_2(i, 1), 2)];
b = [coordinates(N_1_2(i, 2), 2) - coordinates(N_1_2(i, 3), 2);
coordinates(N_1_2(i, 3), 2) - coordinates(N_1_2(i, 1), 2);
coordinates(N_1_2(i, 1), 2) - coordinates(N_1_2(i, 2), 2)];
c = [coordinates(N_1_2(i, 3), 1) - coordinates(N_1_2(i, 2), 1);
coordinates(N_1_2(i, 1), 1) - coordinates(N_1_2(i, 3), 1);
coordinates(N_1_2(i, 2), 1) - coordinates(N_1_2(i, 1), 1)];
%Находим площадь треугольного элемента
S = zeros(3, 3);
for j = 1:3
S(j, 1) = 1;
for k = 2:3
S(j, k) = coordinates(N_1_2(i, j), k - 1);
end
end
A = det(S) / 2; % Площадь
%Строим матрицу градиентов
B = [b(1), 0, b(2), 0, b(3), 0;
0, c(1), 0, c(2), 0, c(3);
c(1), b(1), c(2), b(2), c(3), b(3)];
B = B / (2 * A);
%Строим матрицу жесткости конечного элемента
k_e = A * B' * D_1_2 * B;
%Строим глобальную матрицу жесткости
for n = 1:3
for j = 1:3
K_1_2(2 * N_1_2(i, n) - 1, 2 * N_1_2(i, j) - 1) = K_1_2(2 * N_1_2(i, n) - 1, 2 * N_1_2(i, j) - 1) + k_e(2 * n - 1, 2 * j - 1);
K_1_2(2 * N_1_2(i, n), 2 * N_1_2(i, j) - 1) = K_1_2(2 * N_1_2(i, n), 2 * N_1_2(i, j) - 1) + k_e(2 * n, 2 * j - 1);
K_1_2(2 * N_1_2(i, n) - 1, 2 * N_1_2(i, j)) = K_1_2(2 * N_1_2(i, n) - 1, 2 * N_1_2(i, j)) + k_e(2 * n - 1, 2 * j);
K_1_2(2 * N_1_2(i, n), 2 * N_1_2(i, j)) = K_1_2(2 * N_1_2(i, n), 2 * N_1_2(i, j)) + k_e(2 * n, 2 * j);
end
end
end
%Строим матрицу жесткости для основания
for i = 1:size(N_2, 1)
%Находим коэффициенты для функций форм
a = [coordinates(N_2(i, 2), 1) * coordinates(N_2(i, 3), 2) - coordinates(N_2(i, 3), 1) * coordinates(N_2(i, 2), 2);
coordinates(N_2(i, 3), 1) * coordinates(N_2(i, 1), 2) - coordinates(N_2(i, 1), 1) * coordinates(N_2(i, 3), 2);
coordinates(N_2(i, 1), 1) * coordinates(N_2(i, 2), 2) - coordinates(N_2(i, 2), 1) * coordinates(N_2(i, 1), 2)];
b = [coordinates(N_2(i, 2), 2) - coordinates(N_2(i, 3), 2);
coordinates(N_2(i, 3), 2) - coordinates(N_2(i, 1), 2);
coordinates(N_2(i, 1), 2) - coordinates(N_2(i, 2), 2)];
c = [coordinates(N_2(i, 3), 1) - coordinates(N_2(i, 2), 1);
coordinates(N_2(i, 1), 1) - coordinates(N_2(i, 3), 1);
coordinates(N_2(i, 2), 1) - coordinates(N_2(i, 1), 1)];
%Находим площадь треугольного элемента
S = zeros(3, 3);
for j = 1:3
S(j, 1) = 1;
for k = 2:3
S(j, k) = coordinates(N_2(i, j), k - 1);
end
end
A = det(S) / 2; % Площадь
%Строим матрицу градиентов
B = [b(1), 0, b(2), 0, b(3), 0;
0, c(1), 0, c(2), 0, c(3);
c(1), b(1), c(2), b(2), c(3), b(3)];
B = B / (2 * A);
%Строим матрицу жесткости конечного элемента
k_e = A * B' * D_2 * B;
%Строим глобальную матрицу жесткости
for n = 1:3
for j = 1:3
K_2(2 * N_2(i, n) - 1, 2 * N_2(i, j) - 1) = K_2(2 * N_2(i, n) - 1, 2 * N_2(i, j) - 1) + k_e(2 * n - 1, 2 * j - 1);
K_2(2 * N_2(i, n), 2 * N_2(i, j) - 1) = K_2(2 * N_2(i, n), 2 * N_2(i, j) - 1) + k_e(2 * n, 2 * j - 1);
K_2(2 * N_2(i, n) - 1, 2 * N_2(i, j)) = K_2(2 * N_2(i, n) - 1, 2 * N_2(i, j)) + k_e(2 * n - 1, 2 * j);
K_2(2 * N_2(i, n), 2 * N_2(i, j)) = K_2(2 * N_2(i, n), 2 * N_2(i, j)) + k_e(2 * n, 2 * j);
end
end
end
%Итоговая матрица жесткости
K = K_1_1 + K_1_2 + K_2;
%Зададим вектор-столбец усилия
%Учтём давление жидкости
correct_force = zeros(size(force_nodes, 1), 3);
count = 1;
for i = 1:size(force_elements, 1)
for j = 1:3
for k = 1:size(force_nodes, 1)
if(N(force_elements(i), j) == force_nodes(k))
correct_force(count, 1) = force_nodes(k);
correct_force(count, 2) = coordinates(force_nodes(k), 1);
correct_force(count, 3) = coordinates(force_nodes(k), 2);
count = count + 1;
break;
end
end
end
end %Сохраняем в массив координаты точек, к которым приложена сила
for i = 1:size(correct_force, 1)
for j = (i + 1):size(correct_force, 1)
if(correct_force(i, :) == correct_force(j, :))
correct_force(j, :) = [];
break;
end
end
end %Удаляем повторяющиеся строки
for i = 1:size(correct_force, 1) %Сортировка по возрастанию Y
for j = 1:size(correct_force, 1) - i
if(correct_force(j + 1, 3) < correct_force(j, 3))
temp = correct_force(j + 1, :);
correct_force(j + 1, :) = correct_force(j, :);
correct_force(j, :) = temp;
end
end
end
water_level = correct_force(size(correct_force, 1), 3); %Уровень воды
for i = size(correct_force, 1):-1:2 %Учет давления жидкости
length = sqrt((correct_force(i, 2) - correct_force(i - 1, 2))^2 + (correct_force(i, 3) - correct_force(i - 1, 3))^2);
tg = (correct_force(i, 3) - correct_force(i - 1, 3)) / (correct_force(i, 2) - correct_force(i - 1, 2));
F(2 * correct_force(i - 1, 1) - 1) = rho * g * (water_level - correct_force(i - 1, 3)) * length * cos(pi / 2 - atan(tg));
F(2 * correct_force(i - 1, 1)) = -rho * g * (water_level - correct_force(i - 1, 3)) * length * sin (pi / 2 - atan(tg));
F(2 * correct_force(i, 1) - 1) = rho * g * (water_level - correct_force(i, 3)) * length * cos(pi / 2 - atan(tg));
F(2 * correct_force(i, 1)) = -rho * g * (water_level - correct_force(i, 3)) * length * sin (pi / 2 - atan(tg));
end
%Учтём силю тяжести в плотине 1 материал
for i = 1:size(nodes_1_1, 1)
for j = 1:size(N_1_1, 1)
for k = 1:size(N_1_1, 2)
if(N_1_1(j, k) == nodes_1_1(i))
%Находим площадь треугольного элемента
S = zeros(3, 3);
for n = 1:3
S(n, 1) = 1;
for m = 2:3
S(n, m) = coordinates(N_1_1(j, n), m - 1);
end
end
A = det(S) / 2; % Площадь
F(2 * nodes_1_1(i)) = F(2 * nodes_1_1(i)) - rho_1 * A * g;
end
break;
end
end
end
%Учтём силю тяжести в плотине 2 материал
for i = 1:size(nodes_1_2, 1)
for j = 1:size(N_1_2, 1)
for k = 1:size(N_1_2, 2)
if(N_1_2(j, k) == nodes_1_2(i))
%Находим площадь треугольного элемента
S = zeros(3, 3);
for n = 1:3
S(n, 1) = 1;
for m = 2:3
S(n, m) = coordinates(N_1_2(j, n), m - 1);
end
end
A = det(S) / 2; % Площадь
F(2 * nodes_1_2(i)) = F(2 * nodes_1_2(i)) - rho_2 * A * g;
end
break;
end
end
end
%Используем граничные условия
%Граничные условия по оси Х
for i = 1:size(BC_X_nodes, 1)
F(2 * BC_X_nodes(i) - 1) = 0;
K(2 * BC_X_nodes(i) - 1, :) = 0;
K(:, 2 * BC_X_nodes(i) - 1) = 0;
K(2 * BC_X_nodes(i) - 1, 2 * BC_X_nodes(i) - 1) = 1;
end
%Граничные условия по оси Y
for i = 1:size(BC_Y_nodes, 1)
F(2 * BC_Y_nodes(i)) = 0;
K(2 * BC_Y_nodes(i), :) = 0;
K(:, 2 * BC_Y_nodes(i)) = 0;
K(2 * BC_Y_nodes(i), 2 * BC_Y_nodes(i)) = 1;
end
%Решаем систему
u = K \ F;
u_gl = zeros(size(u, 1) / 2, 3); %Перемещения
u_gl(:, 3) = 0;
for i = 1:(size(u, 1) / 2)
u_gl(i, 1) = u(2 * i - 1);
u_gl(i, 2) = u(2 * i);
end
%Вычисляем напряжения и деформации
for i = 1:size(N, 1)
%Находим коэффициенты для функций форм
a = [coordinates(N(i, 2), 1) * coordinates(N(i, 3), 2) - coordinates(N(i, 3), 1) * coordinates(N(i, 2), 2);
coordinates(N(i, 3), 1) * coordinates(N(i, 1), 2) - coordinates(N(i, 1), 1) * coordinates(N(i, 3), 2);
coordinates(N(i, 1), 1) * coordinates(N(i, 2), 2) - coordinates(N(i, 2), 1) * coordinates(N(i, 1), 2)];
b = [coordinates(N(i, 2), 2) - coordinates(N(i, 3), 2);
coordinates(N(i, 3), 2) - coordinates(N(i, 1), 2);
coordinates(N(i, 1), 2) - coordinates(N(i, 2), 2)];
c = [coordinates(N(i, 3), 1) - coordinates(N(i, 2), 1);
coordinates(N(i, 1), 1) - coordinates(N(i, 3), 1);
coordinates(N(i, 2), 1) - coordinates(N(i, 1), 1)];
%Находим площадь треугольного элемента
S = zeros(3, 3);
for j = 1:3
S(j, 1) = 1;
for k = 2:3
S(j, k) = coordinates(N(i, j), k - 1);
end
end
A = det(S) / 2; % Площадь
%Строим матрицу градиентов
B = [b(1), 0, b(2), 0, b(3), 0;
0, c(1), 0, c(2), 0, c(3);
c(1), b(1), c(2), b(2), c(3), b(3)];
B = B / (2 * A);
%Составляем вектор перемещений узлов элемента
u_3 = zeros(6, 1);
for j = 1:3
u_3(2 * j - 1) = u_gl(N(i, j), 1);
u_3(2 * j) = u_gl(N(i, j), 2);
end
%Вычисляем деформации
epsilon(i, :) = B * u_3;
%Вычисляем напряжения
flag = false;
%Ищем к чему относится конечный элемент(Плотина 1 или 2 материал или
%основание)
for j = 1:size(N_1_2, 1)
if(N(i, :) == N_1_2(j, :))
flag = true;
D = D_1_2;
break;
end
end
if(flag == false)
for j = 1:size(N_1_1, 1)
if(N(i, :) == N_1_1(j, :))
flag = true;
D = D_1_1;
break;
end
end
end
if(flag == false)
for j = 1:size(N_2, 1)
if(N(i, :) == N_2(j, :))
D = D_2;
break;
end
end
end
sigma(i, :) = D * B * u_3;
end