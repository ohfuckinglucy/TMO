function [matrix, avg, var_val] = IA331_Rubtsov_Lab1_2_func(I, J, a, b)
% Входные параметры
%   I  - длина вектора-столбца
%   J  - длина вектора-строки
%   a, b - диапазон [a, b] для равномерного распределения (a < b)
%
% Выходные параметры:
%   matrix - матрица I x J
%   avg - среднее значение матрицы
%   var_val - дисперсия матрицы

    % 1. Генерируем вектор-столбец длины I
    col_vector = a + (b - a) * rand(I, 1);  % rand(I,1) — I строк, 1 столбец

    % 2. Генерируем вектор-строку длины J
    row_vector = a + (b - a) * rand(1, J);  % rand(1,J) — 1 строка, J столбцов

    % матрица размера I x J
    matrix = col_vector * row_vector;

    % 4. Среднее значение матрицы
    avg = mean(matrix, 'all');

    % 5. Дисперсия матрицы
    var_val = var(matrix, 1, 'all');

end