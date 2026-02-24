% =========================================================================
% SAR Range Doppler Algorithm (RDA) Demonstration
% =========================================================================
% Этот скрипт демонстрирует обработку синтетических IQ данных SAR
% с использованием Range Doppler Algorithm (RDA)
%
% Основано на:
% [1] Cumming & Wong, "Digital Processing of SAR Data", 2005, Chapter 6
% [2] Sentinel-1 Level 1 Detailed Algorithm Definition, 2022
% [3] Jansing, "Introduction to SAR", 2021
%
% Автор: Демонстрация для Sentinel-1C
% Дата: 2026
% =========================================================================

clear all; close all; clc;

%% ========================================================================
% ПАРАМЕТРЫ SAR СИСТЕМЫ (близкие к Sentinel-1C IW mode)
% ========================================================================

% Параметры платформы
params.V_r = 7400;              % Скорость платформы, м/с [1, p.39]
params.H = 693e3;               % Высота орбиты, м
params.R_0 = 850e3;             % Минимальная наклонная дальность, м

% Параметры РЛС
params.f_c = 5.405e9;           % Центральная частота, Гц (C-band) [Sentinel-1]
params.lambda = 3e8/params.f_c; % Длина волны, м (~5.55 см)
params.B_r = 56.5e6;            % Полоса ЧМ-импульса, Гц [1, Table 6.1]
params.T_p = 41.74e-6;          % Длительность импульса, с [1, Table 6.1]
params.K_r = params.B_r/params.T_p; % Частотная модуляция (chirp rate), Гц/с

% Параметры антенны
params.L_a = 12.3;              % Длина антенны, м [Sentinel-1]
params.theta_a = params.lambda/params.L_a; % Ширина ДН по азимуту, рад

% Параметры дискретизации
params.F_s = 90e6;              % Частота дискретизации по дальности, Гц
params.PRF = 1700;              % Частота повторения импульсов, Гц [1, Table 6.1]

% Угол косого визирования (squint angle)
params.theta_sq = 5 * pi/180;   % 5 градусов [1, Chapter 6]

% Вычисляемые параметры
params.theta_inc = acos(params.H/params.R_0); % Угол визирования
params.K_a = 2*params.V_r^2/(params.lambda*params.R_0); % FM-rate по азимуту [1, eq.4.27]

% Ширина полосы доплеровских частот [1, eq.4.16]
params.B_a = 2*params.V_r*sin(params.theta_a/2)/params.lambda;

fprintf('=== Параметры SAR системы ===\n');
fprintf('Центральная частота: %.3f GHz\n', params.f_c/1e9);
fprintf('Длина волны: %.2f см\n', params.lambda*100);
fprintf('Скорость платформы: %.0f м/с\n', params.V_r);
fprintf('PRF: %.0f Гц\n', params.PRF);
fprintf('Полоса по дальности: %.2f МГц\n', params.B_r/1e6);
fprintf('Полоса по азимуту: %.2f Гц\n', params.B_a);
fprintf('Угол squint: %.2f град\n', params.theta_sq*180/pi);
fprintf('=============================\n\n');

%% ========================================================================
% ОПРЕДЕЛЕНИЕ ТОЧЕЧНЫХ ЦЕЛЕЙ
% ========================================================================
% Создаем несколько сценариев для демонстрации различных случаев
% соответствующих рисункам из [1, Chapter 6]

scenario = input('Выберите сценарий (1-4, Enter=1): ');
if isempty(scenario)
    scenario = 1;
end

switch scenario
    case 1
        % Одиночная точечная цель в центре сцены [1, Fig.6.3]
        targets = struct('x', 0, 'y', params.R_0, 'rcs', 1);
        scenario_name = 'Одиночная цель';
        
    case 2
        % Три цели на одной дальности [1, Fig.6.4]
        targets = [
            struct('x', -500, 'y', params.R_0, 'rcs', 1),
            struct('x',    0, 'y', params.R_0, 'rcs', 1),
            struct('x',  500, 'y', params.R_0, 'rcs', 1)
        ];
        scenario_name = 'Три цели на одной дальности';
        
    case 3
        % Цели вдоль траектории (азимутальная линия) [1, Fig.6.5]
        targets = [
            struct('x', -800, 'y', params.R_0, 'rcs', 1),
            struct('x', -400, 'y', params.R_0, 'rcs', 1),
            struct('x',    0, 'y', params.R_0, 'rcs', 1),
            struct('x',  400, 'y', params.R_0, 'rcs', 1),
            struct('x',  800, 'y', params.R_0, 'rcs', 1)
        ];
        scenario_name = 'Цели вдоль азимута';
        
    case 4
        % Несколько целей на разных дальностях [1, Fig.6.9]
        targets = [
            struct('x',    0, 'y', params.R_0-1000, 'rcs', 1),
            struct('x', -500, 'y', params.R_0,      'rcs', 1),
            struct('x',  500, 'y', params.R_0,      'rcs', 1),
            struct('x',    0, 'y', params.R_0+1000, 'rcs', 1)
        ];
        scenario_name = 'Цели на разных дальностях';
end

fprintf('Сценарий: %s\n', scenario_name);
fprintf('Количество целей: %d\n\n', length(targets));

%% ========================================================================
% СИНТЕЗ СЫРЫХ IQ ДАННЫХ
% ========================================================================
fprintf('Синтез сырых IQ данных...\n');

% Временные оси
N_range = 4096;                 % Количество отсчетов по дальности
N_azimuth = 2048;               % Количество отсчетов по азимуту

t_range = (0:N_range-1)/params.F_s;          % Быстрое время
t_azimuth = (0:N_azimuth-1)/params.PRF;      % Медленное время
t_azimuth = t_azimuth - mean(t_azimuth);     % Центрирование

[T_az, T_r] = meshgrid(t_azimuth, t_range);

% Инициализация сырых данных
raw_data = zeros(N_range, N_azimuth);

% Генерация сигнала от каждой цели
for i = 1:length(targets)
    x_t = targets(i).x;         % Позиция цели по азимуту, м
    y_t = targets(i).y;         % Позиция цели по дальности, м
    sigma = targets(i).rcs;     % ЭПР
    
    % Траектория платформы с учетом squint
    x_s = params.V_r * T_az;    % Позиция платформы
    
    % Мгновенная наклонная дальность [1, eq.4.11]
    % R(η) = sqrt((x_s - x_t)^2 + y_t^2)
    R_inst = sqrt((x_s - x_t).^2 + y_t^2);
    
    % Двухпроходная задержка
    tau = 2*R_inst/3e8;
    
    % Опорный ЧМ-импульс по дальности [1, eq.6.2]
    % s_r(t) = rect(t/T_p) * exp(j*π*K_r*t^2)
    phase_range = pi * params.K_r * (T_r - tau).^2;
    envelope_range = double(abs(T_r - tau) <= params.T_p/2);
    
    % Фазовая модуляция по азимуту [1, eq.4.19]
    % φ_az(η) = -4π*R(η)/λ
    phase_azimuth = -4*pi*R_inst/params.lambda;
    
    % Амплитудная огибающая по азимуту (ДН антенны) [1, eq.4.15]
    % w_a(η) = sinc^2(0.886 * θ_a * V_r * η / λ)
    azimuth_time_norm = params.V_r * T_az / (params.lambda * params.R_0);
    envelope_azimuth = sinc(0.886 * params.theta_a * azimuth_time_norm).^2;
    
    % Комплексный сигнал
    signal = sigma * envelope_range .* envelope_azimuth .* ...
             exp(1j * (phase_range + phase_azimuth));
    
    raw_data = raw_data + signal;
end

% Добавление шума
SNR_dB = 20;                    % Отношение сигнал/шум, дБ
noise_power = max(abs(raw_data(:)))^2 / (10^(SNR_dB/10));
noise = sqrt(noise_power/2) * (randn(size(raw_data)) + 1j*randn(size(raw_data)));
raw_data = raw_data + noise;

fprintf('Сырые данные синтезированы: %d x %d\n\n', N_range, N_azimuth);

%% ========================================================================
% ВИЗУАЛИЗАЦИЯ СЫРЫХ ДАННЫХ [1, Fig.6.3]
% ========================================================================
figure('Name', 'Сырые IQ данные', 'Position', [100 100 1200 500]);

subplot(1,2,1);
imagesc(t_azimuth, t_range*1e6, abs(raw_data));
xlabel('Азимутальное время (с)');
ylabel('Время дальности (мкс)');
title('Амплитуда сырых данных');
colorbar;
colormap(jet);
axis xy;

subplot(1,2,2);
imagesc(t_azimuth, t_range*1e6, angle(raw_data));
xlabel('Азимутальное время (с)');
ylabel('Время дальности (мкс)');
title('Фаза сырых данных');
colorbar;
colormap(hsv);
axis xy;

sgtitle(sprintf('Сырые IQ данные - %s', scenario_name));

%% ========================================================================
% ШАГИ ОБРАБОТКИ RDA
% ========================================================================

%% -----------------------------------------------------------------------
% ШАГ 1: RANGE COMPRESSION (Сжатие по дальности)
% -----------------------------------------------------------------------
fprintf('ШАГ 1: Range Compression...\n');

% Опорная функция по дальности [1, eq.6.3]
t_ref = t_range - mean(t_range);
ref_range = exp(-1j * pi * params.K_r * t_ref.^2);
ref_range = ref_range .* (abs(t_ref) <= params.T_p/2);

% Согласованная фильтрация в частотной области [2, Section 6.1.2]
% Преобразование в частотную область
S_2D = fft(raw_data, [], 1);
H_range = conj(fft(ref_range'));

% Применение согласованного фильтра
S_rc = S_2D .* repmat(H_range, 1, N_azimuth);

% Обратное преобразование
range_compressed = ifft(S_rc, [], 1);

fprintf('  Разрешение по дальности: %.2f м\n', 3e8/(2*params.B_r));

%% -----------------------------------------------------------------------
% ШАГ 2: AZIMUTH FFT (Переход в Range-Doppler область)
% -----------------------------------------------------------------------
fprintf('ШАГ 2: Azimuth FFT...\n');

% FFT по азимуту [1, Section 6.2]
S_RD = fftshift(fft(range_compressed, [], 2), 2);

% Азимутальная частотная ось
f_azimuth = linspace(-params.PRF/2, params.PRF/2, N_azimuth);

fprintf('  Данные в Range-Doppler области\n');

%% -----------------------------------------------------------------------
% ШАГ 3: DOPPLER CENTROID ESTIMATION (Оценка доплеровского центроида)
% -----------------------------------------------------------------------
fprintf('ШАГ 3: Doppler Centroid Estimation...\n');

% Простая оценка по энергетическому центру спектра [1, Section 4.6]
azimuth_spectrum = sum(abs(S_RD).^2, 1);
f_dc = sum(f_azimuth .* azimuth_spectrum) / sum(azimuth_spectrum);

fprintf('  Оценка Doppler Centroid: %.2f Гц\n', f_dc);

% Коррекция доплеровского центроида (сдвиг в частотной области)
phase_shift = exp(-1j * 2*pi * f_dc * t_azimuth);
S_RD = S_RD .* repmat(phase_shift, N_range, 1);

% Обновление частотной оси
f_azimuth = f_azimuth - f_dc;

%% -----------------------------------------------------------------------
% ШАГ 4: SECONDARY RANGE COMPRESSION (SRC)
% -----------------------------------------------------------------------
fprintf('ШАГ 4: Secondary Range Compression...\n');

% SRC компенсирует азимутально-зависимое расфокусирование по дальности
% возникающее из-за squint angle [2, Section 6.3.1]

% Частотная ось по дальности
f_range = linspace(-params.F_s/2, params.F_s/2, N_range)';

% Вычисление D(f_r, f_a) - SRC фазовой функции [1, eq.6.25]
% Для каждой азимутальной частоты
for k = 1:N_azimuth
    f_a = f_azimuth(k);
    
    % Эффективная скорость [1, eq.6.21]
    V_eff = params.V_r * sqrt(1 - (params.lambda * f_a / (2*params.V_r))^2);
    
    % SRC фазовая функция [1, eq.6.25]
    % H_src = exp(-j*π*λ^2*R_0*f_r^2/(8*V_eff^2))
    H_src = exp(-1j * pi * params.lambda^2 * params.R_0 * f_range.^2 / ...
                (8 * V_eff^2));
    
    % Применение SRC
    S_RD(:, k) = S_RD(:, k) .* H_src;
end

fprintf('  SRC применено\n');

%% -----------------------------------------------------------------------
% ШАГ 5: RANGE CELL MIGRATION CORRECTION (RCMC)
% -----------------------------------------------------------------------
fprintf('ШАГ 5: Range Cell Migration Correction...\n');

% RCMC исправляет кривизну траектории цели в Range-Doppler области
% [1, Section 6.3, 6.4]

% Вычисление RCMC для каждой азимутальной частоты
for k = 1:N_azimuth
    f_a = f_azimuth(k);
    
    % Расчет миграции дальности [1, eq.6.17]
    % ΔR(f_a) = R_0 * (1/cos(θ_sq) * sqrt(1-(λ*f_a/(2*V_r))^2) - 1)
    cos_theta = cos(params.theta_sq);
    delta_R = params.R_0 * (1/cos_theta * ...
              sqrt(1 - (params.lambda*f_a/(2*params.V_r))^2) - 1);
    
    % Преобразование миграции в задержку
    delta_tau = 2 * delta_R / 3e8;
    
    % Интерполяция для коррекции [2, Section 6.3.2]
    % Применяем фазовый сдвиг в частотной области
    H_rcmc = exp(1j * 2*pi * f_range * delta_tau);
    S_RD(:, k) = S_RD(:, k) .* H_rcmc;
end

fprintf('  RCMC применена\n');

%% -----------------------------------------------------------------------
% ШАГ 6: AZIMUTH COMPRESSION (Сжатие по азимуту)
% -----------------------------------------------------------------------
fprintf('ШАГ 6: Azimuth Compression...\n');

% Согласованный фильтр по азимуту [1, Section 6.5]
% Для каждой дальности применяем фазовую коррекцию

for i = 1:N_range
    % Опорная функция по азимуту [1, eq.6.29]
    % H_az(f_a) = exp(j*π*f_a^2/K_a)
    H_azimuth = exp(1j * pi * f_azimuth.^2 / params.K_a);
    
    % Применение согласованного фильтра
    S_RD(i, :) = S_RD(i, :) .* H_azimuth;
end

fprintf('  Разрешение по азимуту: %.2f м\n', params.lambda*params.R_0/(2*params.L_a));

%% -----------------------------------------------------------------------
% ШАГ 7: INVERSE AZIMUTH FFT (Обратное преобразование по азимуту)
% -----------------------------------------------------------------------
fprintf('ШАГ 7: Inverse Azimuth FFT...\n');

% Обратное FFT по азимуту
focused_image = ifft(ifftshift(S_RD, 2), [], 2);

fprintf('  Фокусированное изображение получено\n\n');

%% ========================================================================
% ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ
% ========================================================================

% Преобразование координат для визуализации
range_axis = t_range * 3e8/2;   % Наклонная дальность, м
azimuth_axis = t_azimuth * params.V_r; % Позиция по азимуту, м

%% Рисунок 1: Амплитуда данных после Range Compression [1, Fig.6.3]
figure('Name', 'После Range Compression', 'Position', [100 100 1200 500]);

subplot(1,2,1);
imagesc(t_azimuth, range_axis/1e3, 20*log10(abs(range_compressed)+eps));
xlabel('Азимутальное время (с)');
ylabel('Наклонная дальность (км)');
title('Range Compressed - Амплитуда (дБ)');
colorbar;
caxis([-40 0]);
colormap(jet);
axis xy;

subplot(1,2,2);
imagesc(azimuth_axis, range_axis/1e3, 20*log10(abs(range_compressed)+eps));
xlabel('Позиция по азимуту (м)');
ylabel('Наклонная дальность (км)');
title('Range Compressed - Амплитуда (дБ)');
colorbar;
caxis([-40 0]);
colormap(jet);
axis xy;

sgtitle('После сжатия по дальности (Range Compression)');

%% Рисунок 2: Range-Doppler область [1, Fig.6.4]
figure('Name', 'Range-Doppler Domain', 'Position', [100 100 1200 500]);

subplot(1,2,1);
imagesc(f_azimuth, range_axis/1e3, 20*log10(abs(S_RD)+eps));
xlabel('Доплеровская частота (Гц)');
ylabel('Наклонная дальность (км)');
title('Range-Doppler область - Амплитуда (дБ)');
colorbar;
caxis([-40 0]);
colormap(jet);
axis xy;

subplot(1,2,2);
imagesc(f_azimuth, range_axis/1e3, angle(S_RD));
xlabel('Доплеровская частота (Гц)');
ylabel('Наклонная дальность (км)');
title('Range-Doppler область - Фаза (рад)');
colorbar;
colormap(hsv);
axis xy;

sgtitle('Range-Doppler область (после RCMC и SRC)');

%% Рисунок 3: Фокусированное изображение [1, Fig.6.5, 6.9]
figure('Name', 'Фокусированное изображение', 'Position', [100 100 1200 500]);

subplot(1,2,1);
imagesc(azimuth_axis, range_axis/1e3, 20*log10(abs(focused_image)+eps));
xlabel('Позиция по азимуту (м)');
ylabel('Наклонная дальность (км)');
title('Фокусированное изображение - Амплитуда (дБ)');
colorbar;
caxis([-40 0]);
colormap(jet);
axis xy;
grid on;

subplot(1,2,2);
imagesc(azimuth_axis, range_axis/1e3, angle(focused_image));
xlabel('Позиция по азимуту (м)');
ylabel('Наклонная дальность (км)');
title('Фокусированное изображение - Фаза (рад)');
colorbar;
colormap(hsv);
axis xy;
grid on;

sgtitle(sprintf('Финальное фокусированное изображение - %s', scenario_name));

%% Рисунок 4: Профили точечных целей
if length(targets) > 0
    figure('Name', 'Профили точечных целей', 'Position', [100 100 1200 500]);
    
    % Находим центральную цель для анализа
    center_target_idx = ceil(length(targets)/2);
    x_target = targets(center_target_idx).x;
    y_target = targets(center_target_idx).y;
    
    % Индексы ближайшей позиции цели
    [~, idx_azimuth] = min(abs(azimuth_axis - x_target));
    [~, idx_range] = min(abs(range_axis - y_target));
    
    % Профиль по дальности
    subplot(1,2,1);
    range_profile = abs(focused_image(:, idx_azimuth));
    plot(range_axis/1e3, 20*log10(range_profile/max(range_profile)));
    xlabel('Наклонная дальность (км)');
    ylabel('Нормированная амплитуда (дБ)');
    title('Профиль по дальности');
    grid on;
    ylim([-40 5]);
    
    % Профиль по азимуту
    subplot(1,2,2);
    azimuth_profile = abs(focused_image(idx_range, :));
    plot(azimuth_axis, 20*log10(azimuth_profile/max(azimuth_profile)));
    xlabel('Позиция по азимуту (м)');
    ylabel('Нормированная амплитуда (дБ)');
    title('Профиль по азимуту');
    grid on;
    ylim([-40 5]);
    
    sgtitle('Профили точечной цели');
end

%% ========================================================================
% АНАЛИЗ КАЧЕСТВА ФОКУСИРОВКИ
% ========================================================================
fprintf('=== Анализ качества фокусировки ===\n');

% Теоретическое разрешение
res_range_theory = 3e8/(2*params.B_r);
res_azimuth_theory = params.lambda*params.R_0/(2*params.L_a);

fprintf('Теоретическое разрешение по дальности: %.2f м\n', res_range_theory);
fprintf('Теоретическое разрешение по азимуту: %.2f м\n', res_azimuth_theory);

% Измеренное разрешение (по уровню -3дБ)
if length(targets) > 0
    % Для центральной цели
    [~, idx_azimuth] = min(abs(azimuth_axis - targets(center_target_idx).x));
    [~, idx_range] = min(abs(range_axis - targets(center_target_idx).y));
    
    % Разрешение по дальности
    range_profile = abs(focused_image(:, idx_azimuth));
    [max_val, max_idx] = max(range_profile);
    half_power = max_val/sqrt(2);
    idx_3dB = find(range_profile >= half_power);
    if length(idx_3dB) > 1
        res_range_measured = (range_axis(idx_3dB(end)) - range_axis(idx_3dB(1)));
        fprintf('Измеренное разрешение по дальности: %.2f м\n', res_range_measured);
    end
    
    % Разрешение по азимуту
    azimuth_profile = abs(focused_image(idx_range, :));
    [max_val, max_idx] = max(azimuth_profile);
    half_power = max_val/sqrt(2);
    idx_3dB = find(azimuth_profile >= half_power);
    if length(idx_3dB) > 1
        res_azimuth_measured = abs(azimuth_axis(idx_3dB(end)) - azimuth_axis(idx_3dB(1)));
        fprintf('Измеренное разрешение по азимуту: %.2f м\n', res_azimuth_measured);
    end
end

fprintf('===================================\n');
fprintf('\nОбработка завершена успешно!\n');
fprintf('Демонстрация соответствует Рис. 6.1-6.9 из [1]\n');
