clear; clc; close all;

%% 1. Функция одной переменной и графики

% Задаём функцию
f = @(x) (x.^sin(x)) ./ cos(x);

% Диапазон x кроме 0 тк там значение будет неопределенно и pi/2, где cos будет 0
x = linspace(0.1, pi/2 - 0.1, 500);

% Считаем значения f(x)
y_f = f(x);

% Численная производная с помощью diff
h = x(2) - x(1);
y_deriv = diff(y_f) / h;           % производная
x_deriv = x(1:end-1) + h/2;        % точки для производной (на 1 меньше)

% Интеграл F(x) от 0 до х_вал
F = @(x_val) integral(f, 0, x_val);
y_integr = arrayfun(F, x);         % считаем интеграл в каждой точке x

% Строим графики
figure;
hold on;
plot(x, y_f, 'b-', 'LineWidth', 2, 'DisplayName', 'f(x)');
plot(x_deriv, y_deriv, 'r--', 'LineWidth', 2, 'DisplayName', 'f''(x)');
plot(x, y_integr, 'g-.', 'LineWidth', 2, 'DisplayName', 'Интеграл F(x)');
title('Графики f(x), f''(x) и F(x)');
xlabel('x');
ylabel('y');
legend('Location', 'best');
grid on;
hold off;

%% 2. Решение уравнения a*x + b = f(x)

% Выбираем коэффициенты
a = 3;
b = 0;

% Графическое решение
figure;
xx = linspace(0.1, pi/2 - 0.1, 500);
plot(xx, f(xx), 'b-', 'LineWidth', 2, 'DisplayName', 'f(x)');
hold on;
plot(xx, a*xx + b, 'k-', 'LineWidth', 2, 'DisplayName', ['Прямая y = ' num2str(a) 'x']);
title('Решение уравнения a*x + b = f(x)');
xlabel('x');
ylabel('y');
legend('Location', 'best');
grid on;

% Численное решение с помощью fsolve
equation = @(x) f(x) - (a*x + b);   % уравнение f(x) = a*x + b
x0 = 0.5;                           % начальное приближение
x_sol = fsolve(equation, x0);       % решаем

% Отмечаем точку пересечения
plot(x_sol, f(x_sol), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
    'DisplayName', ['Решение: x = ' num2str(x_sol, 3)]);
legend;  % обновляем легенду

% Вывод результата
fprintf('Уравнение: %g*x + %g = f(x)\n', a, b);
fprintf('Численное решение: x = %.4f\n', x_sol);
fprintf('Проверка: f(x) = %.4f, a*x + b = %.4f\n', f(x_sol), a*x_sol + b);

hold off;

%% 3. Функция двух переменных

% Задаём функцию F(x,y)
F2 = @(x,y) (x.^2 + y.^2) ./ sqrt(x + y);

% Создаём сетку
x2 = linspace(0.1, 2, 50);
y2 = linspace(0.1, 2, 50);
[X, Y] = meshgrid(x2, y2);

% Считаем Z
Z = F2(X, Y);

% Убираем возможные ошибки (NaN или Inf)
Z(isnan(Z)) = NaN;
Z(isinf(Z)) = NaN;

% 3D-график
figure;
surf(X, Y, Z);
title('F(x,y) = (x^2 + y^2)/sqrt(x + y)');
xlabel('x');
ylabel('y');
zlabel('F(x,y)');
colorbar;
grid on;
shading interp;