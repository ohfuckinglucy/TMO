clear; clc; close all;

%% 6

% Параметры из Таблицы 3.2
N = 1000;   % длина реализации
K = 1000;   % количество реализаций
mu_omega = 0;
sigma_omega = 1;

% Генерация случайных блужданий
X_rw = zeros(N, K); % матрица N x K
omega = mu_omega + sigma_omega * randn(N, K); % приращения

% Начальное условие xi[0] = 0 учтено
X_rw(1, :) = omega(1, :);
for n = 2:N
    X_rw(n, :) = X_rw(n-1, :) + omega(n, :);
end

% --- 1. Построение всех реализаций ---
figure('Name', 'Реализации случайного блуждания', 'NumberTitle', 'off');
plot(1:N, X_rw, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
hold on;
plot(1:N, mean(X_rw, 2), 'r-', 'LineWidth', 2, 'DisplayName', 'Среднее по ансамблю');
xlabel('Время n');
ylabel('\xi[n]');
title(sprintf('Случайные блуждания: %d реализаций, N=%d, \\omega ~ N(0,1)', K, N));
legend('Location', 'best'); grid on;

% --- 2. Скаттерограммы: пары (n, n-1) ---
pairs1 = [10, 9; 50, 49; 100, 99; 200, 199];
colors1 = lines(size(pairs1, 1));

figure('Name', 'Диаграммы рассеяния: пары (n, n-1)', 'NumberTitle', 'off');
for p = 1:size(pairs1, 1)
    ni = pairs1(p, 1); nj = pairs1(p, 2);
    scatter(X_rw(nj, :), X_rw(ni, :), 15, colors1(p, :), 'filled');
    hold on;
end
title('\xi[n] vs \xi[n-1]');
xlabel('\xi[n-1]'); ylabel('\xi[n]');
legend(cellstr(arrayfun(@(x,y) sprintf('(%d,%d)',x,y), pairs1(:,1), pairs1(:,2),'UniformOutput',false)));
grid on; axis equal;

% --- 3. Скаттерограммы: пары (n, n-10) ---
pairs2 = [50, 40; 100, 90; 200, 190];
colors2 = lines(size(pairs2, 1));

figure('Name', 'Диаграммы рассеяния: пары (n, n-10)', 'NumberTitle', 'off');
for p = 1:size(pairs2, 1)
    ni = pairs2(p, 1); nj = pairs2(p, 2);
    scatter(X_rw(nj, :), X_rw(ni, :), 15, colors2(p, :), 'filled');
    hold on;
end
title('\xi[n] vs \xi[n-10]');
xlabel('\xi[n-10]'); ylabel('\xi[n]');
legend(cellstr(arrayfun(@(x,y) sprintf('(%d,%d)',x,y), pairs2(:,1), pairs2(:,2),'UniformOutput',false)));
grid on; axis equal;

% --- 4. Теория vs практика ---
fprintf('-------------------------------------------------------------------\n');
fprintf('СРАВНЕНИЕ ТЕОРИИ И ПРАКТИКИ: автокорреляция случайного блуждания\n');

all_pairs = [pairs1; pairs2];
for p = 1:size(all_pairs, 1)
    ni = all_pairs(p, 1);
    nj = all_pairs(p, 2);
    l = ni - nj;
    r_empirical = mean(X_rw(ni, :) .* X_rw(nj, :));
    r_theoretical = nj;
    rho_empirical = corr(X_rw(ni, :)', X_rw(nj, :)');
    rho_theoretical = sqrt(nj / ni);
    fprintf('\nПара (n=%d, m=%d), лаг=%d:\n', ni, nj, l);
    fprintf('  Теория:    r=%d, rho=%.4f\n', r_theoretical, rho_theoretical);
    fprintf('  Практика:  r=%.2f, rho=%.4f\n', r_empirical, rho_empirical);
end
fprintf('-------------------------------------------------------------------\n');

%% 7

n_vals = 2:N;
r_empirical = zeros(size(n_vals));
for i = 1:length(n_vals)
    n = n_vals(i);
    r_empirical(i) = mean(X_rw(n, :) .* X_rw(n-1, :));
end

r_theoretical = n_vals - 1;

figure('Name', 'Автокорреляция: теория vs практика', 'NumberTitle', 'off');
plot(n_vals, r_empirical, 'b-', 'LineWidth', 1.5, 'DisplayName','Эмпирическая');
hold on;
plot(n_vals, r_theoretical, 'r--', 'LineWidth', 2, 'DisplayName','Теоретическая');
xlabel('Время n'); ylabel('r_\xi(n,n-1)');
title('Сравнение r_\xi(n,n-1): теория и практика');
legend('Location','best'); grid on;

fprintf('-------------------------------------------------------------------\n');
fprintf('ЗАДАНИЕ 7: Сравнение r_\\xi(n,n-1)\n\n');
sample_points = [2, 10, 50, 100, 500, 999];
fprintf('%4s | %10s | %10s | %10s\n','n','Теория','Практика','Ошибка');
fprintf('------------------------------------------------\n');
for i = 1:length(sample_points)
    n = sample_points(i);
    if n>=2 && n<=N
        idx = n-1;
        fprintf('%4d | %10.2f | %10.2f | %10.2f\n', n, r_theoretical(idx), r_empirical(idx), abs(r_theoretical(idx)-r_empirical(idx)));
    end
end
fprintf('-------------------------------------------------------------------\n');

fprintf('\nОТВЕТ НА ВОПРОС:\n');
fprintf('Можно ли оценить r_\\xi(n,n-1) по одной реализации?\n');
fprintf('→ НЕТ. Случайное блуждание нестационарно и не эргодично.\n');
fprintf('Для корректной оценки нужна выборка по ансамблю.\n');
fprintf('-------------------------------------------------------------------\n');

%% 8

fprintf('-------------------------------------------------------------------\n');
fprintf('ЗАДАНИЕ 8: Процесс с затуханием\n');
fprintf('xi[n] = 0.9*xi[n-1] + w[n], w ~ N(0,1)\n\n');

fprintf('1. Рекуррентная формула:\n');
fprintf('   sigma_xi^2[n] = 0.81*sigma_xi^2[n-1] + 1\n\n');

fprintf('2. Общая формула:\n');
fprintf('   sigma_xi^2[n] = (1 - 0.81^n)/0.19\n\n');

sigma2_inf = 1/0.19;
fprintf('3. При n→∞: sigma_xi^2[n] → %.4f\n\n', sigma2_inf);

fprintf('4. В стационарном режиме:\n');
fprintf('   r_xi(l) = 0.9^{|l|} * %.4f\n\n', sigma2_inf);

fprintf('5. Процесс нестационарен при малых n, но асимптотически стационарен.\n');
fprintf('-------------------------------------------------------------------\n');
