%% 9
clear; clc; close all;

% Параметры
N = 1000;   % длина реализации
K = 1000;   % количество реализаций
a = 0.9;    % коэффициент затухания
mu_omega = 0;
sigma_omega = 1;

% Генерация AR(1)-процесса
X_ar = zeros(N, K);           % инициализация
omega = mu_omega + sigma_omega * randn(N, K);  % шум

% Рекурсия: xi[n] = a * xi[n-1] + w[n]
for n = 2:N
    X_ar(n, :) = a * X_ar(n-1, :) + omega(n, :);
end
X_ar(1, :) = omega(1, :);  % xi[1] = w[1]

% 1. Построение всех реализаций на одном графике
figure('Name', 'Реализации AR(1)-процесса', 'NumberTitle', 'off');
plot(1:N, X_ar, 'Color', [0.85 0.85 0.9], 'LineWidth', 0.5);  % светло-серые линии
hold on;
plot(1:N, mean(X_ar, 2), 'r-', 'LineWidth', 2, 'DisplayName', 'Среднее по ансамблю');
xlabel('Время n');
ylabel('\xi[n]');
title(sprintf('AR(1)-процесс: %d реализаций, \\xi[n] = %.1f \\cdot \\xi[n-1] + \\omega[n]', K, a));
legend('Location', 'best');
grid on;

% 2. Диаграммы рассеяния — Группа 1: (n, n-1)
pairs1 = [10, 9; 50, 49; 100, 99; 200, 199];
colors1 = lines(size(pairs1, 1));

figure('Name', 'Диаграммы рассеяния: пары (n, n-1) — AR(1)', 'NumberTitle', 'off');
for p = 1:size(pairs1, 1)
    ni = pairs1(p, 1);
    nj = pairs1(p, 2);
    scatter(X_ar(nj, :), X_ar(ni, :), 15, colors1(p, :), 'filled');
    hold on;
end
title('AR(1): \xi[n] vs \xi[n-1]');
xlabel('\xi[n-1]');
ylabel('\xi[n]');
legend(cellstr(arrayfun(@(x,y) sprintf('(%d, %d)', x,y), pairs1(:,1), pairs1(:,2), 'UniformOutput', false)), ...
       'Location', 'best');
grid on; axis equal;

% 3. Диаграммы рассеяния — Группа 2: (n, n-10)
pairs2 = [50, 40; 100, 90; 200, 190];
colors2 = lines(size(pairs2, 1));

figure('Name', 'Диаграммы рассеяния: пары (n, n-10) — AR(1)', 'NumberTitle', 'off');
for p = 1:size(pairs2, 1)
    ni = pairs2(p, 1);
    nj = pairs2(p, 2);
    scatter(X_ar(nj, :), X_ar(ni, :), 15, colors2(p, :), 'filled');
    hold on;
end
title('AR(1): \xi[n] vs \xi[n-10]');
xlabel('\xi[n-10]');
ylabel('\xi[n]');
legend(cellstr(arrayfun(@(x,y) sprintf('(%d, %d)', x,y), pairs2(:,1), pairs2(:,2), 'UniformOutput', false)), ...
       'Location', 'best');
grid on; axis equal;

% 4. Сравнение с теорией (Задание 8)
fprintf('-------------------------------------------------------------------\n');
fprintf('СРАВНЕНИЕ С ТЕОРИЕЙ: AR(1)-процесс\n');

% Теоретическая дисперсия в установившемся режиме
sigma2_theory_inf = 1 / (1 - a^2);  % = 1 / (1 - 0.81) = 1/0.19 ≈ 5.2632

% Эмпирическая дисперсия на последнем шаге
sigma2_empirical_end = mean(X_ar(end, :) .^ 2);

fprintf('Теоретическая дисперсия (n→∞): %.4f\n', sigma2_theory_inf);
fprintf('Эмпирическая дисперсия (n=%d): %.4f\n', N, sigma2_empirical_end);
fprintf('Относительная ошибка: %.2f%%\n', 100 * abs(sigma2_theory_inf - sigma2_empirical_end) / sigma2_theory_inf);

% 5. Сравнение автокорреляции с теорией
fprintf('\nСравнение автокорреляции r(n, n-l):\n');
fprintf('%8s | %8s | %10s | %10s | %10s\n', 'n', 'l', 'Теория', 'Практика', 'Ошибка');
fprintf('--------------------------------------------------------------\n');

all_pairs = [pairs1; pairs2];
for p = 1:size(all_pairs, 1)
    ni = all_pairs(p, 1);
    nj = all_pairs(p, 2);
    l = ni - nj;
    
    % Эмпирическая автокорреляция
    r_emp = mean(X_ar(ni, :) .* X_ar(nj, :));
    
    % Теоретическая автокорреляция (в стационарном режиме)
    % r(l) = a^l * sigma^2
    r_theory = (a^l) * sigma2_theory_inf;
    
    fprintf('%8d | %8d | %10.3f | %10.3f | %10.3f\n', ni, l, r_theory, r_emp, abs(r_theory - r_emp));
end
fprintf('-------------------------------------------------------------------\n');

% 6. Сравнение с Заданием 6 (обычные блуждания) — ВЫВОДЫ
fprintf('\nСРАВНЕНИЕ С ЗАДАНИЕМ 6 (СЛУЧАЙНЫЕ БЛУЖДАНИЯ):\n');
fprintf('1. Визуально:\n');
fprintf('   - В Задании 6 траектории расходятся → дисперсия растёт → НЕСТАЦИОНАРНО.\n');
fprintf('   - Здесь — траектории "колеблются" в ограниченной полосе → дисперсия стабилизируется → СТАЦИОНАРНО.\n\n');

fprintf('2. Диаграммы рассеяния:\n');
fprintf('   - В Задании 6: сильная корреляция, особенно при больших n (ρ→1).\n');
fprintf('   - Здесь: корреляция экспоненциально затухает с ростом l → ρ(l) = 0.9^l.\n\n');

fprintf('3. Автокорреляция:\n');
fprintf('   - В Задании 6: r(n, n-l) = n - l → растёт с n.\n');
fprintf('   - Здесь: r(l) = const * 0.9^l → зависит только от l, не от n (в стационарном режиме).\n');
fprintf('-------------------------------------------------------------------\n');


%% 10

% === ШАГ 0: Генерация данных (если не существуют) ===
% Параметры
N = 1000; K = 1000;
a = 0.9;
mu = 0; sigma = 1;

% Белый шум (если не существует)
if ~exist('X_wn', 'var')
    X_wn = mu + sigma * randn(N, K);
end

% Случайное блуждание (если не существует)
if ~exist('X_rw', 'var')
    X_rw = zeros(N, K);
    omega_rw = mu + sigma * randn(N, K);
    X_rw(1, :) = omega_rw(1, :);
    for n = 2:N
        X_rw(n, :) = X_rw(n-1, :) + omega_rw(n, :);
    end
end

% AR(1) (если не существует)
if ~exist('X_ar', 'var')
    X_ar = zeros(N, K);
    omega_ar = mu + sigma * randn(N, K);
    X_ar(1, :) = omega_ar(1, :);
    for n = 2:N
        X_ar(n, :) = a * X_ar(n-1, :) + omega_ar(n, :);
    end
end

% === ШАГ 1: Выборочная автокорреляция по ансамблю для AR(1) — r(n, n-1) ===
n_vals = 2:N;
r_ar_empirical = zeros(size(n_vals));
for i = 1:length(n_vals)
    n = n_vals(i);
    r_ar_empirical(i) = mean(X_ar(n, :) .* X_ar(n-1, :));
end

% Теоретическая автокорреляция для AR(1) в стационарном режиме: r(1) = a * sigma^2
sigma2_inf = 1 / (1 - a^2);  % ≈ 5.2632
r_ar_theory = a * sigma2_inf;  % const = 0.9 * 5.2632 ≈ 4.7369

% График: эмпирическая vs теоретическая (горизонтальная линия)
figure('Name', 'AR(1): r_\xi(n, n-1) — теория vs практика', 'NumberTitle', 'off');
plot(n_vals, r_ar_empirical, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Эмпирическая \hat{r}_\xi(n, n-1)');
hold on;
plot([2 N], [r_ar_theory r_ar_theory], 'r--', 'LineWidth', 2, 'DisplayName', 'Теория: const = 0.9 / 0.19');
xlabel('Время n');
ylabel('r_\xi(n, n-1)');
title('AR(1): выборочная автокорреляция по ансамблю');
legend('Location', 'best');
grid on;

% === ШАГ 2: Сравнение с Заданием 6 (случайные блуждания) ===
% Рассчитаем r_rw(n, n-1) для сравнения
r_rw_empirical = zeros(size(n_vals));
for i = 1:length(n_vals)
    n = n_vals(i);
    r_rw_empirical(i) = mean(X_rw(n, :) .* X_rw(n-1, :));
end
r_rw_theory = n_vals - 1;  % теория для блужданий

% График сравнения
figure('Name', 'Сравнение r(n, n-1): Блуждания vs AR(1)', 'NumberTitle', 'off');
plot(n_vals, r_rw_empirical, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Блуждания (эмп.)');
hold on;
plot(n_vals, r_rw_theory, 'r--', 'LineWidth', 1, 'DisplayName', 'Блуждания (теория)');
plot(n_vals, r_ar_empirical, 'b-', 'LineWidth', 1.5, 'DisplayName', 'AR(1) (эмп.)');
plot([2 N], [r_ar_theory r_ar_theory], 'b--', 'LineWidth', 1, 'DisplayName', 'AR(1) (теория)');
xlabel('Время n');
ylabel('r_\xi(n, n-1)');
title('Сравнение автокорреляции: случайные блуждания vs AR(1)');
legend('Location', 'best');
grid on;

% === ШАГ 3: Проверка эргодичности — среднее по времени vs ансамблю для AR(1) ===
l1 = 10; l2 = 100;

fprintf('-------------------------------------------------------------------\n');
fprintf('ПРОВЕРКА ЭРГОДИЧНОСТИ ДЛЯ AR(1)-ПРОЦЕССА (две реализации)\n');
fprintf('Лаги: l1=%d, l2=%d\n\n', l1, l2);

% Выберем две произвольные реализации (например, 1 и 2)
k1 = 1; k2 = 2;

for l = [l1, l2]
    % Среднее по ансамблю (теоретически правильное)
    r_ensemble = mean(X_ar(l+1:end, :) .* X_ar(1:end-l, :), 2);  % вектор длины (N-l)
    r_ensemble_mean = mean(r_ensemble);  % усредним по n для сравнения
    
    % Среднее по времени для реализации k1
    xi1 = X_ar(:, k1);
    r_time_k1 = mean(xi1(l+1:end) .* xi1(1:end-l));
    
    % Среднее по времени для реализации k2
    xi2 = X_ar(:, k2);
    r_time_k2 = mean(xi2(l+1:end) .* xi2(1:end-l));
    
    % Теоретическое значение
    r_theory_l = (a^l) * sigma2_inf;
    
    fprintf('Лаг l = %d:\n', l);
    fprintf('  Теория:                   %.4f\n', r_theory_l);
    fprintf('  Среднее по ансамблю:      %.4f\n', r_ensemble_mean);
    fprintf('  Среднее по времени (k=%d): %.4f\n', k1, r_time_k1);
    fprintf('  Среднее по времени (k=%d): %.4f\n', k2, r_time_k2);
    fprintf('  Отклонение (k1):          %.4f\n', abs(r_time_k1 - r_ensemble_mean));
    fprintf('  Отклонение (k2):          %.4f\n', abs(r_time_k2 - r_ensemble_mean));
    fprintf('\n');
end

fprintf('→ При больших N и стационарности процесса, среднее по времени\n');
fprintf('  близко к среднему по ансамблю → процесс эргодичен по автокорреляции.\n');
fprintf('-------------------------------------------------------------------\n');

% === ШАГ 4: Графики НОРМИРОВАННОЙ АКФ (коэффициент корреляции ρ(l)) ===
max_lag = 50;
lags = 0:max_lag;

% Инициализация
acf_wn_norm = zeros(size(lags));
acf_rw_norm = zeros(size(lags));
acf_ar_norm = zeros(size(lags));

% --- Белый шум: стационарен ---
sigma2_wn = mean(X_wn(:).^2);  % ≈ 1
for l = lags
    if l == 0
        acf_wn_norm(l+1) = 1;  % ρ(0) = 1
    else
        cov = mean(mean(X_wn(1:end-l, :) .* X_wn(l+1:end, :)));
        acf_wn_norm(l+1) = cov / sigma2_wn;  % ρ(l) = cov / σ²
    end
end

% --- AR(1): стационарен ---
sigma2_ar = mean(X_ar(:).^2);  % ≈ 5.26
for l = lags
    if l == 0
        acf_ar_norm(l+1) = 1;
    else
        cov = mean(mean(X_ar(1:end-l, :) .* X_ar(l+1:end, :)));
        acf_ar_norm(l+1) = cov / sigma2_ar;
    end
end

% --- Случайные блуждания: НЕСТАЦИОНАРНЫ → берём в момент n = N ---
n_fixed = N;
sigma2_rw_N = mean(X_rw(n_fixed, :).^2);  % Var(ξ[N]) ≈ N * σ² = 1000

for l = lags
    if l == 0
        acf_rw_norm(l+1) = 1;
    else
        if n_fixed - l >= 1
            cov = mean(X_rw(n_fixed, :) .* X_rw(n_fixed - l, :));
            sigma2_rw_N_minus_l = mean(X_rw(n_fixed - l, :).^2);  % Var(ξ[N-l])
            acf_rw_norm(l+1) = cov / sqrt(sigma2_rw_N * sigma2_rw_N_minus_l);
        else
            acf_rw_norm(l+1) = NaN;
        end
    end
end

% --- Теоретические значения ---
acf_wn_theory_norm = [1, zeros(1, max_lag)];  % белый шум — только ρ(0)=1

acf_ar_theory_norm = a.^lags;  % для AR(1): ρ(l) = a^l

% Для блужданий: теоретически
% ρ(N, N-l) = (N - l) / sqrt(N * (N - l)) = sqrt((N - l)/N)
acf_rw_theory_norm = sqrt((n_fixed - lags) / n_fixed);

% --- Построение графика ---
figure('Name', 'Сравнение НОРМИРОВАННЫХ АКФ (ρ(l))', 'NumberTitle', 'off');
hold on;

plot(lags, acf_wn_norm, 'g-o', 'MarkerSize', 4, 'DisplayName', 'Белый шум (эмп.)');
plot(lags, acf_wn_theory_norm, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Белый шум (теория)');

plot(lags, acf_rw_norm, 'r-s', 'MarkerSize', 4, 'DisplayName', 'Случайные блуждания (эмп.)');
plot(lags, acf_rw_theory_norm, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Случайные блуждания (теория)');

plot(lags, acf_ar_norm, 'b-^', 'MarkerSize', 4, 'DisplayName', 'AR(1) (эмп.)');
plot(lags, acf_ar_theory_norm, 'b--', 'LineWidth', 1.5, 'DisplayName', 'AR(1) (теория)');

xlabel('Лаг l');
ylabel('Нормированная автокорреляция \rho(l)');
title('Сравнение нормированных автокорреляционных функций');
legend('Location', 'best');
grid on;
xlim([0 max_lag]);
ylim([-0.1 1.1]);  % чтобы видеть всё в одном масштабе

% === ШАГ 5: Выводы ===
fprintf('\n\n=== ОСНОВНЫЕ ВЫВОДЫ ===\n');
fprintf('1. AR(1)-процесс:\n');
fprintf('   - Выборочная автокорреляция r(n, n-1) стабилизируется со временем →\n');
fprintf('     соответствует теории (константа ≈ 4.737).\n');
fprintf('   - В отличие от случайных блужданий (где r(n, n-1) = n-1 → растёт),\n');
fprintf('     здесь процесс стационарен → автокорреляция зависит только от лага.\n\n');

fprintf('2. Эргодичность:\n');
fprintf('   - Для AR(1) среднее по времени близко к среднему по ансамблю →\n');
fprintf('     процесс эргодичен (при больших N).\n\n');

fprintf('3. Сравнение АКФ:\n');
fprintf('   - Белый шум: корреляция только на нулевом лаге.\n');
fprintf('   - Случайные блуждания: медленно затухающая (почти линейная) АКФ —\n');
fprintf('     признак нестационарности.\n');
fprintf('   - AR(1): экспоненциально затухающая АКФ — признак стационарного\n');
fprintf('     авторегрессионного процесса.\n');
fprintf('-------------------------------------------------------------------\n');