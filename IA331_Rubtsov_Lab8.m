clc; close all;

%% Вычисление характеристик СМО

% M/M/1
[N_q_MM1, N_MM1, W_MM1, T_MM1] = smo_func('M/M/1');

% M/G/1
[N_q_MG1p, N_MG1p, W_MG1p, T_MG1p] = smo_func('M/G/1_p');
[N_q_MG1, N_MG1, W_MG1, T_MG1] = smo_func('M/G/1_Cb');

% M/D/1
[N_q_MD1, N_MD1, W_MD1, T_MD1] = smo_func('M/D/1');

%% Графики

% M/M/1

figure();
subplot(4, 1, 1)
plot(0:0.1:0.9, N_q_MM1);
title("Средняя длина очереди M/M/1")
xlabel("p");
ylabel("N_q");
subplot(4, 1, 2)
plot(0:0.1:0.9, N_MM1);
title("Среднее число заявок в СМО M/M/1")
xlabel("p");
ylabel("N");
subplot(4, 1, 3)
plot(0:0.1:0.9, W_MM1);
title("Среднее время ожидания M/M/1")
xlabel("p");
ylabel("W");
subplot(4, 1, 4)
plot(0:0.1:0.9, T_MM1);
title("Среднее время пребывания требования в системе M/M/1")
xlabel("p");
ylabel("T");

% M/D/1

figure();
subplot(4, 1, 1)
plot(0:0.1:0.9, N_q_MD1);
title("Средняя длина очереди M/D/1")
xlabel("p");
ylabel("N_q");
subplot(4, 1, 2)
plot(0:0.1:0.9, N_MD1);
title("Среднее число заявок в СМО M/D/1")
xlabel("p");
ylabel("N");
subplot(4, 1, 3)
plot(0:0.1:0.9, W_MD1);
title("Среднее время ожидания M/D/1")
xlabel("p");
ylabel("W");
subplot(4, 1, 4)
plot(0:0.1:0.9, T_MD1);
title("Среднее время пребывания требования в системе M/D/1")
xlabel("p");
ylabel("T");

% M/G/1

figure();
subplot(4, 1, 1)
plot(0:10:100, N_q_MG1);
title("Средняя длина очереди M/G/1")
xlabel("Cb");
ylabel("N_q");
subplot(4, 1, 2)
plot(0:10:100, N_MG1);
title("Среднее число заявок в СМО M/G/1")
xlabel("Cb");
ylabel("N");
subplot(4, 1, 3)
plot(0:10:100, W_MG1);
title("Среднее время ожидания M/G/1")
xlabel("Cb");
ylabel("W");
subplot(4, 1, 4)
plot(0:10:100, T_MG1);
title("Среднее время пребывания требования в системе M/G/1")
xlabel("Cb");
ylabel("T");

figure();
subplot(4, 1, 1)
plot(0:0.1:0.9, N_q_MG1p);
title("Средняя длина очереди M/G/1")
xlabel("p");
ylabel("N_q");
subplot(4, 1, 2)
plot(0:0.1:0.9, N_MG1p);
title("Среднее число заявок в СМО M/G/1")
xlabel("p");
ylabel("N");
subplot(4, 1, 3)
plot(0:0.1:0.9, W_MG1p);
title("Среднее время ожидания M/G/1")
xlabel("p");
ylabel("W");
subplot(4, 1, 4)
plot(0:0.1:0.9, T_MG1p);
title("Среднее время пребывания требования в системе M/G/1")
xlabel("p");
ylabel("T");

%% Функция для расчета хар-ик СМО

function [N_q, N, W, T] = smo_func(name)
    x = 1/8; % среднее время обслуживания

    % Коэффициент загрузки и дисперсия
    p = 0:0.1:0.9;
    C_b = 20;
    Cb_vec = 0:10:100;
    p_fix = 0.4;

    switch name
        case 'M/G/1_p'
            N_q = (p.^2 .* (1 + C_b^2)) ./ (2*(1 - p)); % Средняя длина очереди
            N = p + N_q; % Среднее число заявок в СМО
            W = x * N_q ./ p; % Среднее время ожидания
            T = W + x; % Среднее время пребывания требования в системе

        case 'M/G/1_Cb'
            N_q = (p_fix^2 .* (1 + Cb_vec.^2)) ./ (2*(1 - p_fix));
            N = p_fix + N_q;
            W = x * N_q / p_fix;
            T = W + x;

        case 'M/M/1'
            N_q = p.^2 ./ (1 - p);
            N = p ./ (1 - p);
            W = x * p ./ (1 - p);
            T = x ./ (1 - p);

        case 'M/D/1'
            N_q = (p.^2) ./ (2*(1 - p));
            N = p + N_q;
            W = x * p ./ (2*(1 - p));
            T = (x*(1-p))./(2*(1-p));
    end
end
