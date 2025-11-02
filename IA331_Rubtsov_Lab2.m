%% 2cos(2x) [0; pi/4]

clear; clc;

%% 1, 2 - расчет выборочных значений случайной перменной, заданной законом
% Параметры
N_values = [50, 200, 1000];
alpha_levels = [0.1, 0.05, 0.01];
x_limits = [0, pi/4]; % Ограничения по x

% Генерация выборок (метод обратных функций)
sample_50  = 0.5 * asin(rand(N_values(1), 1));
sample_200 = 0.5 * asin(rand(N_values(2), 1));
sample_1000 = 0.5 * asin(rand(N_values(3), 1));

samples = {sample_50, sample_200, sample_1000};

%% 3, 4 - расчет точечных оценок среднего, дисперсии и СКО и интервальных оценок среднего и дисперсии для a

fprintf('=== Доверительные интервалы для среднего ===\n');
fprintf('N\t\tα\t\tLower\t\tUpper\n');
fprintf('--------------------------------------------\n');

for i = 1:length(N_values)
    x = samples{i};
    n = N_values(i);
    mean_x = mean(x);
    std_x = std(x);
    
    for j = 1:length(alpha_levels)
        alpha = alpha_levels(j);
        t_crit = tinv(1 - alpha/2, n-1);      % квантиль t-распределения
        margin = t_crit * std_x / sqrt(n);
        lower_mean = mean_x - margin;
        upper_mean = mean_x + margin;
        
        fprintf('%.4d\t%.2f\t%.6f\t%.6f\n', n, alpha, lower_mean, upper_mean);
    end
end

fprintf('\n=== Доверительные интервалы для дисперсии ===\n');
fprintf('N\t\tα\t\tLower\t\tUpper\n');
fprintf('--------------------------------------------\n');

for i = 1:length(N_values)
    x = samples{i};
    n = N_values(i);
    var_x = var(x);  % несмещённая дисперсия
    
    for j = 1:length(alpha_levels)
        alpha = alpha_levels(j);
        chi2_low = chi2inv(alpha/2, n-1);      % нижняя квантиль
        chi2_high = chi2inv(1 - alpha/2, n-1); % верхняя квантиль
        
        lower_var = (n-1)*var_x / chi2_high;   % нижняя граница
        upper_var = (n-1)*var_x / chi2_low;    % верхняя граница
        
        fprintf('%.4d\t%.2f\t%.6f\t%.6f\n', n, alpha, lower_var, upper_var);
    end
end

%% 6, 7 - построить гистограммы, описывающие закон распределения сл. величины и поверх график теор. пл. распределения вероятности

figure('Position', [100, 100, 900, 600]);
for i = 1:length(N_values)
    N = N_values(i);
    x = samples{i};
    
    k = floor(1 + 3.2 * log(N));
    edges = linspace(x_limits(1), x_limits(2), k+1);
    
    subplot(1, 3, i);
    histogram(x, edges, 'Normalization', 'pdf', 'EdgeColor', 'k', 'LineWidth', 0.5);
    hold on;
    
    % Теоретическая плотность f(x) = 2*cos(2x) — без масштабирования!
    x_grid = linspace(x_limits(1), x_limits(2), 500);
    f_x = 2 * cos(2 * x_grid);
    plot(x_grid, f_x, 'r-', 'LineWidth', 2, 'DisplayName', 'Теоретическая плотность');
    
    title(sprintf('N = %d, k = %d', N, k));
    xlabel('x');
    ylabel('f(x)');
    xlim(x_limits);
    ylim([0, 2.2]);  % чтобы видеть максимум 2
    grid on;
    legend('Location', 'northeast');
    hold off;
end
sgtitle('Гистограммы (PDF) с теоретической плотностью f(x) = 2cos(2x)');

%% 8 - график теор-ой функции распр. вероятности и плотность распр. вероятности

x = linspace(0, pi/4, 1000);

% Теоретические функции
f_x = 2 * cos(2 * x);         % плотность вероятности (PDF)
F_x = sin(2 * x);             % функция распределения (CDF)

figure('Position', [100, 100, 800, 500]);
plot(x, f_x, 'b-', 'LineWidth', 2, 'DisplayName', 'Плотность f(x) = 2cos(2x)');
hold on;
plot(x, F_x, 'r--', 'LineWidth', 2, 'DisplayName', 'Функция распределения F(x) = sin(2x)');

title('Теоретическая плотность и функция распределения');
xlabel('x');
ylabel('Значение');
xlim([0, pi/4]);
ylim([0, 1.1]);  % F(x) до 1, f(x) до 2 — но F(x) важнее для масштаба

grid on;
legend('Location', 'southeast');
hold off;

% ключевые точки
text(pi/4, 0, '\leftarrow x=\pi/4', 'HorizontalAlignment', 'right', 'FontSize', 9);
text(0, 1, '\uparrow F(\pi/4)=1', 'VerticalAlignment', 'bottom', 'FontSize', 9);
text(0, 2, '\uparrow f(0)=2', 'VerticalAlignment', 'bottom', 'FontSize', 9, 'Color', 'blue');

%% 10 P {X = k} = ((0.5^k)/k!)*e^(-0.5); k = 0,1, ... - повторить шаги 1- 9 для дискретно распр. случайно величины

lambda = 0.5;
K_max = 20;

% Шаг 1: Вычисление теоретических вероятностей P(X = k) и CDF
k_vals = 0:K_max;
P_k = exp(-lambda) .* (lambda.^k_vals) ./ factorial(k_vals);  % P(X=k)
F_k = cumsum(P_k);  % F(k) = P(X ≤ k)

sampless = cell(1, length(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    U = rand(N, 1);  % равномерные числа
    
    % Для каждого U находим k: минимальное k, где F_k >= U
    sample = zeros(N, 1);
    for j = 1:N
        idx = find(F_k >= U(j), 1, 'first');
        sample(j) = k_vals(idx);
    end
    sampless{i} = sample;
end

sample_50  = sampless{1};
sample_200 = sampless{2};
sample_1000 = sampless{3};

% Шаг 3: Расчёт точечных оценок
fprintf('\n=== Точечные оценки ===\n');
fprintf('N\tСреднее\t\tДисперсия\tСКО\n');
fprintf('----------------------------------------\n');

for i = 1:length(N_values)
    x = sampless{i};
    mean_x = mean(x);
    var_x = var(x);
    std_x = std(x);
    fprintf('%d\t%.6f\t%.6f\t%.6f\n', N_values(i), mean_x, var_x, std_x);
end

% Шаг 5: Доверительные интервалы (α = 0.1, 0.05, 0.01)
alpha_levels = [0.1, 0.05, 0.01];

fprintf('\n=== Доверительные интервалы для среднего ===\n');
fprintf('N\tα\tLower\t\tUpper\n');
fprintf('--------------------------------------------\n');

for i = 1:length(N_values)
    x = sampless{i};
    n = length(x);
    mean_x = mean(x);
    std_x = std(x);
    
    for j = 1:length(alpha_levels)
        alpha = alpha_levels(j);
        t_crit = tinv(1 - alpha/2, n-1);
        margin = t_crit * std_x / sqrt(n);
        lower_mean = mean_x - margin;
        upper_mean = mean_x + margin;
        fprintf('%d\t%.2f\t%.6f\t%.6f\n', n, alpha, lower_mean, upper_mean);
    end
end

fprintf('\n=== Доверительные интервалы для дисперсии ===\n');
fprintf('N\tα\tLower\t\tUpper\n');
fprintf('--------------------------------------------\n');

for i = 1:length(N_values)
    x = sampless{i};
    n = length(x);
    var_x = var(x);
    
    for j = 1:length(alpha_levels)
        alpha = alpha_levels(j);
        chi2_low = chi2inv(alpha/2, n-1);
        chi2_high = chi2inv(1 - alpha/2, n-1);
        lower_var = (n-1)*var_x / chi2_high;
        upper_var = (n-1)*var_x / chi2_low;
        fprintf('%d\t%.2f\t%.6f\t%.6f\n', n, alpha, lower_var, upper_var);
    end
end

% Шаг 7: Построение гистограммы + теоретическое распределение
figure('Position', [100, 100, 900, 600]);
for i = 1:length(N_values)
    x = sampless{i};
    N = N_values(i);
    
    subplot(1, 3, i);
    
    % Эмпирические частоты: считаем сколько раз встречается каждое k
    k_range = 0:10;  % покажем первые 10 значений — достаточно
    emp_counts = histcounts(x, [k_range, inf]);  % подсчитываем частоты
    emp_probs = emp_counts / N;  % относительные частоты = эмпирическая вероятность
    
    % Теоретические вероятности
    theor_probs = exp(-lambda) .* (lambda.^k_range) ./ factorial(k_range);
    
    % Построение столбцов
    bar(k_range, emp_probs, 'FaceColor', [0.7 0.8 1], 'EdgeColor', 'k');
    hold on;
    plot(k_range, theor_probs, 'ro-', 'LineWidth', 2, 'MarkerSize', 6, ...
         'DisplayName', 'Теоретическая вероятность');
    
    title(sprintf('N = %d', N));
    xlabel('k');
    ylabel('Вероятность');
    xlim([-0.5, 10.5]);
    ylim([0, 0.7]);
    grid on;
    legend('Location', 'northeast');
    hold off;
end
sgtitle('Эмпирические и теоретические распределения Poisson(λ=0.5)');

% Шаг 9: График PMF и CDF (отдельный)
figure('Position', [100, 100, 800, 500]);
k_plot = 0:15;

% PMF (теоретическая)
pmf = exp(-lambda) .* (lambda.^k_plot) ./ factorial(k_plot);

% CDF (накопленная)
cdf = cumsum(pmf);

% Построение
subplot(2,1,1);
stem(k_plot, pmf, 'b-', 'filled', 'LineWidth', 1.5);
title('Распределение Пуассона (\lambda = 0.5)');
ylabel('P(X = k)');
grid on;
xlim([-0.5, 15.5]);

subplot(2,1,2);
stairs(k_plot, cdf, 'r-', 'LineWidth', 1.5);
xlabel('k');
ylabel('F(k) = P(X ≤ k)');
grid on;
xlim([-0.5, 15.5]);

sgtitle('Теоретическая PMF и CDF для Poisson(0.5)');

% Генерация нормальной случайной величины
rng('shuffle');
X_norm = randn(1000, 1);

% Расчёт асимметрии и эксцесса для нормальной выборки
skewness_X_norm = skewness(X_norm);
kurtosis_X_norm = kurtosis(X_norm);  % избыточный эксцесс

% Расчёт асимметрии и эксцесса для (f(x)=2cos(2x))
skewness_X_orig = skewness(sample_1000);
kurtosis_X_orig = kurtosis(sample_1000);

fprintf('\n=== Асимметрия и эксцесс ===\n');
fprintf('Распределение             | Асимметрия (Skewness) | Эксцесс (Kurtosis)\n');
fprintf('------------------------------------------------------------------\n');
fprintf('Нормальное (N(0,1))       | %.6f               | %.6f\n', skewness_X_norm, kurtosis_X_norm);
fprintf('Ваше распределение (2cos) | %.6f               | %.6f\n', skewness_X_orig, kurtosis_X_orig);
fprintf('------------------------------------------------------------------\n');
fprintf('Примечание: Эксцесс — избыточный (для нормального = 0).\n');

% Генерация выборки
sample = 0.5 * asin(rand(1000,1));

% Определяем границы выборки с небольшим запасом
x_min = min(sample);
x_max = max(sample);

% Добавляем запас: расширяем интервал на 1% с каждой стороны
margin = (x_max - x_min) * 0.01;
x_min_extended = max(0, x_min - margin);   % не ниже 0
x_max_extended = min(pi/4, x_max + margin); % не выше pi/4

% Создаём плотную сетку, покрывающую всю выборку + запас
n_points = 2000;  % много точек — для точности
x_grid = linspace(x_min_extended, x_max_extended, n_points);

% Вычисляем теоретическую CDF: F(x) = sin(2x), обрезаем в [0,1]
F_theor_values = sin(2 * x_grid);
F_theor_values = max(0, min(1, F_theor_values));  % защита от артефактов

% Формируем CDF-матрицу: 2 столбца — [x, F(x)]
cdf_matrix = [x_grid', F_theor_values'];

fprintf('Диапазон выборки: [%.6f, %.6f]\n', min(sample), max(sample));
fprintf('Диапазон CDF:     [%.6f, %.6f]\n', min(x_grid), max(x_grid));

% Проводим тест Колмогорова–Смирнова
[h, p, ksstat] = kstest(sample, 'CDF', cdf_matrix, 'Alpha', 0.01);

fprintf('\n=== Тест Колмогорова–Смирнова (α = 0.01) ===\n');
fprintf('Статистика KS: %.6f\n', ksstat);
fprintf('p-значение:   %.6f\n', p);

if h == 0
    result_str = 'НЕ ОТКЛОНЯЕТСЯ';
else
    result_str = 'ОТКЛОНЯЕТСЯ';
end
fprintf('Гипотеза H0 (распределение совпадает): %s\n', result_str);