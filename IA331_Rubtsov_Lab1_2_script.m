clear; clc; close all;

% Параметры
I = 13;  % длина вектора-столбца
J = 28;  % длина вектора-строки
a = 0;   % нижняя граница диапазона
b = 10;  % верхняя граница

% Вызов функции
[matrix, avg, var_val] = IA331_Rubtsov_Lab1_2_func(I, J, a, b);

% Вывод результатов
fprintf('Лабораторная работа номер 1.2\n');
fprintf('Размеры: I = %d, J = %d\n', I, J);
fprintf('Диапазон: [%g, %g]\n', a, b);
fprintf('Размер матрицы: %d x %d\n', size(matrix,1), size(matrix,2));
fprintf('Среднее значение: %.4f\n', avg);
fprintf('Дисперсия: %.4f\n', var_val);

% Построение графика матрицы
figure;
imagesc(matrix);
colorbar;
title(['Матрица ', num2str(I), '×', num2str(J), ' (вектор-столбец × вектор-строка)']);
xlabel('Столбцы');
ylabel('Строки');