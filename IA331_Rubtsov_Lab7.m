clc; close all;

%% 3 входные параметры

labmda = 5; % интенсивность поступления
mu = 8; % интенсивность обслуживания
N = 1000; % число заявок

%% 4 генерация СМО

% M/M/1
tn_mm1 = exprnd(1/labmda, 1, N);
vn_mm1 = exprnd(1/mu, 1, N);
smo_func('M/M/1', lambda, mu, tn_mm1, vn_mm1, N);

% M/G/1
tn_mg1 = exprnd(1/labmda, 1, N);
vn_mg1 = gamrnd(1/mu, 1, N);
smo_func('M/G/1', lambda, mu, tn_mg1, vn_mg1, N);

% G/M/1
tn_gm1 = gamrnd(1/labmda, 1, N);
vn_gm1 = exprnd(1/mu, 1, N);
smo_func('G/M/1', lambda, mu, tn_gm1, vn_gm1, N);

% G/G/1
tn_gg1 = gamrnd(1/labmda, 1, N);
vn_gg1 = gamrnd(1/labmda, 1, N);
smo_func('G/G/1', lambda, mu, tn_gg1, vn_gg1, N);

%% 1 приложенная функция (написал сам)

function smo_func(name, labmda, mu, tn, vn, N)
    % моменты прихода
    t_arr = zeros(1, length(tn));
    t_arr(1) = tn(1);

    for i = 2:N
        t_arr(i) = t_arr(i - 1) + tn(i);
    end

    % моменты ухода
    v_arr = zeros(1, length(vn));

    for i = 1:N
        if i == 1
            start = t_arr(i);
        else
            start = max(t_arr(i), v_arr(i-1));
        end
        v_arr(i) = start + vn(i);
    end

    sys_time = v_arr - t_arr;

    wait_time = zeros(1, N);

    for i = 1:N
        if i == 1
            wait_time(i) = 0;
        else
            wait_time(i) = max(0, v_arr(i - 1) - t_arr(i));
        end
    end

    all_arr = [t_arr, v_arr]; % все события
    what_arr = [ones(1, N), -ones(1, N)]; % 1 приходы - 1 уходы

    events = [all_arr(:), what_arr(:)];
    events = sortrows(events, 1);
    sorted_time = events(:, 1);
    sorted_type = events(:, 2);

    num_in_sys = zeros(1, 2*N + 1);
    time_line = zeros(1, 2*N + 1);
    cur = 0;

    for i = 1:(2*N)
        cur = cur + sorted_type(i);
        time_line(i + 1) = sorted_time(i);
        num_in_sys(i + 1) = cur;
    end

    % Поступление обслуживание

    figure('Name', name);
    subplot(3, 1, 1);
    plot(t_arr, 0:N-1, 'b', 'DisplayName', 'Поступили');
    hold on;
    plot(v_arr, 0:N-1, 'r', 'DisplayName', 'Обслужены');
    title([name ': поступили/обслужены']);
    xlabel('время');
    ylabel('кол-во заявок');
    legend;
    grid on;

    % Кол-во заявок в системе

    subplot(3, 1, 2);
    stairs(time_line, num_in_sys, 'k');
    title([name ': кол-во заявок в системе']);
    xlabel('время');
    ylabel('кол-во заявок');
    grid on;

    % Распределение заявок

    subplot(3,1,3);
    histogram(num_in_sys, 'Normalization', 'probability');
    title([name ': распределение числа заявок']);
    xlabel('Число заявок'); 
    ylabel('Вероятность');
    grid on;

    % расчет характеристик

    L = mean(num_in_sys); % среднее число заявок в СМО
    tq = mean(wait_time); % среднее время пребывания в очереди
    w = mean(sys_time); % среднее время пребывания в системе
    rho = L; % коэффициент загрузки

    fprintf('\n%s:\n', name);
    fprintf("среднее число заявок в СМО: %.4f\n", L);
    fprintf("среднее время пребывания в очереди: %.4f\n", tq);
    fprintf("среднее время пребывания в системе: %.4f\n", w);
    fprintf("коэффициент загрузки: %.4f\n", rho);
end