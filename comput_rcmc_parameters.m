function [delta_R, delta_tau] = compute_rcmc_parameters(f_azimuth, params)
% =========================================================================
% Вычисление параметров Range Cell Migration Correction
% =========================================================================
% Функция вычисляет величину миграции дальности для каждой
% азимутальной частоты
%
% Вход:
%   f_azimuth - вектор азимутальных частот, Гц
%   params - структура параметров SAR системы
%
% Выход:
%   delta_R - миграция дальности, м
%   delta_tau - миграция во времени, с
%
% Референции:
%   [1] Cumming & Wong, 2005, eq. 6.17, 6.18
%   [2] Sentinel-1 Algorithm, Section 6.3.2
% =========================================================================

% Извлечение параметров
R_0 = params.R_0;           % Минимальная дальность, м
V_r = params.V_r;           % Скорость платформы, м/с
lambda = params.lambda;     % Длина волны, м
theta_sq = params.theta_sq; % Угол squint, рад

% Вычисление миграции дальности [1, eq.6.17]
% ΔR(f_a) = R_0 * (sec(θ_sq) * sqrt(1 - (λ*f_a/(2*V_r))^2) - 1)

cos_theta = cos(theta_sq);
sec_theta = 1/cos_theta;

% Нормированная доплеровская частота
f_norm = lambda * f_azimuth / (2 * V_r);

% Проверка диапазона (|f_norm| должно быть < 1)
if any(abs(f_norm) >= 1)
    warning('Некоторые азимутальные частоты выходят за допустимый диапазон');
    f_norm(abs(f_norm) >= 1) = sign(f_norm(abs(f_norm) >= 1)) * 0.999;
end

% Миграция дальности
delta_R = R_0 * (sec_theta * sqrt(1 - f_norm.^2) - 1);

% Преобразование в миграцию времени [1, eq.6.18]
% Δτ = 2*ΔR/c
c = 3e8;  % Скорость света
delta_tau = 2 * delta_R / c;

% Дополнительная коррекция для больших углов squint [1, Section 6.4]
if abs(theta_sq) > 10*pi/180  % > 10 градусов
    % Учет члена второго порядка
    % ΔR_2 = R_0 * λ² * f_a² / (8*V_r²) * tan²(θ_sq)
    delta_R_second = R_0 * lambda^2 * f_azimuth.^2 / (8*V_r^2) * tan(theta_sq)^2;
    delta_R = delta_R + delta_R_second;
    delta_tau = 2 * delta_R / c;
end

end
