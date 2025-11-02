clc; close all;

lambda = 5;    % интенсивность поступления
mu = 8;       % интенсивность обслуживания
N = 1000;      % число заявок

% M/M/1
tn = exprnd(1/lambda, 1, N);
vn = exprnd(1/mu, 1, N);
model_smo(tn, vn, 'M/M/1', lambda, mu);

% M/G/1
tn = exprnd(1/lambda, 1, N);
vn = normrnd(1/mu, 0.3/mu, 1, N);
vn = max(vn, 0.01/mu);  % убираем отрицательные
model_smo(tn, vn, 'M/G/1', lambda, mu);

% G/M/1
tn = unifrnd(0.5/lambda, 1.5/lambda, 1, N);
vn = exprnd(1/mu, 1, N);
model_smo(tn, vn, 'G/M/1', lambda, mu);

% G/G/1
tn = unifrnd(0.5/lambda, 1.5/lambda, 1, N);
vn = normrnd(1/mu, 0.3/mu, 1, N);
vn = max(vn, 0.01/mu);
model_smo(tn, vn, 'G/G/1', lambda, mu);


function model_smo(arr_dt, serv, name, lambda, mu)
    % 1: моменты прихода
    t_arr = zeros(1, length(arr_dt));
    t_arr(1) = arr_dt(1);
    for i = 2:length(arr_dt)
        t_arr(i) = t_arr(i-1) + arr_dt(i);
    end

    % 2: моменты ухода
    t_dep = zeros(1, length(serv));
    for i = 1:length(serv)
        if i == 1
            start = t_arr(i);
        else
            if t_arr(i) > t_dep(i-1)
                start = t_arr(i);
            else
                start = t_dep(i-1);
            end
        end
        t_dep(i) = start + serv(i);
    end

    % 3: время в очереди и в системе
    wait_time = zeros(1, length(serv));
    sys_time = zeros(1, length(serv));
    for i = 1:length(serv)
        if i == 1
            wait_time(i) = 0;
        else
            wait_time(i) = max(0, t_dep(i-1) - t_arr(i));
        end
        sys_time(i) = t_dep(i) - t_arr(i);
    end

    % 4: число заявок в системе во времени
    % все события: приход (+1), уход (-1)
    total_events = 2 * length(t_arr);
    event_time = zeros(1, total_events);
    event_type = zeros(1, total_events);
    k = 1;
    for i = 1:length(t_arr)
        event_time(k) = t_arr(i);   event_type(k) = +1;  k = k + 1;
        event_time(k) = t_dep(i);   event_type(k) = -1;  k = k + 1;
    end
    % Сортируем по времени
    for i = 1:total_events-1
        for j = i+1:total_events
            if event_time(i) > event_time(j)
                tmp = event_time(i); event_time(i) = event_time(j); event_time(j) = tmp;
                tmp = event_type(i); event_type(i) = event_type(j); event_type(j) = tmp;
            end
        end
    end
    % Строим график числа заявок
    time_plot = [0, event_time];
    num_plot = zeros(1, total_events+1);
    current = 0;
    for i = 1:total_events
        current = current + event_type(i);
        num_plot(i+1) = current;
    end

    % 5: графики
    figure('Name', name);

    subplot(3,1,1);
    plot(t_arr, 0:length(t_arr)-1, 'b', t_dep, 0:length(t_dep)-1, 'r');
    title([name ': поступившие и обслуженные']);
    xlabel('Время'); ylabel('Номер заявки');
    legend('Поступили','Обслужены'); grid on;

    subplot(3,1,2);
    stairs(time_plot, num_plot, 'k');
    title([name ': число заявок в системе']);
    xlabel('Время'); ylabel('Число заявок'); grid on;

    subplot(3,1,3);
    histogram(num_plot, 'Normalization', 'probability');
    title([name ': распределение числа заявок']);
    xlabel('Число заявок'); ylabel('Вероятность'); grid on;

    % 6: считаем характеристики
    L = mean(num_plot);               % среднее число в системе
    Wq = mean(wait_time);             % среднее время в очереди
    W = mean(sys_time);               % среднее время в системе
    rho = L;                          % коэффициент загрузки (эмпирически)

    fprintf('\n%s:\n', name);
    fprintf('  Коэффициент загрузки: %.4f\n', rho);
    fprintf('  Среднее число заявок в СМО: %.4f\n', L);
    fprintf('  Среднее время в очереди: %.4f\n', Wq);
    fprintf('  Среднее время в системе: %.4f\n', W);
end