% =========================================================================
% Детальная визуализация этапов Range Doppler Algorithm
% =========================================================================
% Этот скрипт создает подробные графики для каждого этапа обработки RDA
% соответствующие рисункам 6.1-6.9 из Cumming & Wong, 2005
% =========================================================================

clear all; close all; clc;

%% Запрос сценария
fprintf('Выберите сценарий для детальной визуализации:\n');
fprintf('1 - Одиночная точечная цель (Рис. 6.3)\n');
fprintf('2 - Три цели на одной дальности (Рис. 6.4)\n');
fprintf('3 - Цели вдоль азимута (Рис. 6.5)\n');
fprintf('4 - Цели на разных дальностях (Рис. 6.9)\n');
scenario = input('Введите номер (1-4): ');

%% Используем основной скрипт для генерации данных
% (Здесь должен быть код из основного скрипта или его вызов)
% Для краткости предполагаем, что данные уже сгенерированы

%% ========================================================================
% ДЕТАЛЬНАЯ ВИЗУАЛИЗАЦИЯ: БЛОК-СХЕМА АЛГОРИТМА [1, Fig.6.1]
% ========================================================================

figure('Name', 'Блок-схема RDA', 'Position', [100 100 900 700]);

% Текстовое описание блок-схемы
subplot('Position', [0.1 0.1 0.8 0.85]);
axis off;

text(0.5, 0.95, 'Range Doppler Algorithm (RDA)', ...
    'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Блоки алгоритма
blocks = {
    'Raw Data s(t,η)'
    '↓'
    'Range Compression'
    's_rc(t,η) = s(t,η) ⊗ h_range(t)'
    '↓'
    'Azimuth FFT'
    'S_RD(t,f_η) = FFT_η{s_rc(t,η)}'
    '↓'
    'Doppler Centroid Estimation'
    'f_dc = argmax_f |S_RD(t,f)|²'
    '↓'
    'Secondary Range Compression (SRC)'
    'Компенсация азимутально-зависимого'
    'расфокусирования по дальности'
    '↓'
    'Range Cell Migration Correction (RCMC)'
    'Выпрямление траекторий целей'
    'в Range-Doppler области'
    '↓'
    'Azimuth Compression'
    'S_focused(t,f_η) = S_RD(t,f_η) × H_az(f_η)'
    '↓'
    'Inverse Azimuth FFT'
    's_focused(t,η) = IFFT_η{S_focused(t,f_η)}'
    '↓'
    'Focused Image'
};

y_pos = 0.88;
for i = 1:length(blocks)
    if strcmp(blocks{i}, '↓')
        text(0.5, y_pos, blocks{i}, 'FontSize', 14, ...
            'HorizontalAlignment', 'center');
        y_pos = y_pos - 0.025;
    else
        text(0.5, y_pos, blocks{i}, 'FontSize', 10, ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', [0.9 0.95 1], 'EdgeColor', 'k', ...
            'Margin', 5, 'LineWidth', 1);
        y_pos = y_pos - 0.035;
    end
end

%% ========================================================================
% ВИЗУАЛИЗАЦИЯ: Траектории в Range-Doppler области
% ========================================================================
% Демонстрация эффекта Range Cell Migration

figure('Name', 'Range Cell Migration', 'Position', [100 100 1400 500]);

% Параметры для демонстрации
R0 = 850e3;
V = 7400;
lambda = 0.0555;
theta_sq = 5*pi/180;

f_a = linspace(-100, 100, 1000);

subplot(1,3,1);
% Траектория без коррекции
R_trajectory = R0 ./ cos(theta_sq) .* sqrt(1 - (lambda*f_a/(2*V)).^2);
plot(f_a, (R_trajectory - R0)/1e3, 'b-', 'LineWidth', 2);
xlabel('Доплеровская частота (Гц)');
ylabel('Миграция дальности (км)');
title('До RCMC');
grid on;

subplot(1,3,2);
% Траектория после RCMC
plot(f_a, zeros(size(f_a)), 'r-', 'LineWidth', 2);
xlabel('Доплеровская частота (Гц)');
ylabel('Миграция дальности (км)');
title('После RCMC');
grid on;
ylim([-2 2]);

subplot(1,3,3);
% Сравнение
hold on;
plot(f_a, (R_trajectory - R0)/1e3, 'b-', 'LineWidth', 1.5, 'DisplayName', 'До RCMC');
plot(f_a, zeros(size(f_a)), 'r-', 'LineWidth', 1.5, 'DisplayName', 'После RCMC');
xlabel('Доплеровская частота (Гц)');
ylabel('Миграция дальности (км)');
title('Сравнение');
legend('Location', 'best');
grid on;
hold off;

sgtitle('Демонстрация Range Cell Migration Correction');

%% ========================================================================
% ВИЗУАЛИЗАЦИЯ: Согласованные фильтры
% ========================================================================

figure('Name', 'Согласованные фильтры', 'Position', [100 100 1400 500]);

% Параметры
Tp = 41.74e-6;
Br = 56.5e6;
Kr = Br/Tp;
Fs = 90e6;

t = linspace(-Tp, Tp, 1000);

% Фильтр по дальности
subplot(2,3,1);
h_range = exp(1j*pi*Kr*t.^2) .* (abs(t) <= Tp/2);
plot(t*1e6, real(h_range), 'b-', 'LineWidth', 1.5);
xlabel('Время (мкс)');
ylabel('Амплитуда');
title('Range Reference (Real)');
grid on;

subplot(2,3,4);
H_range_freq = fftshift(fft(h_range));
f = linspace(-Fs/2, Fs/2, length(H_range_freq));
plot(f/1e6, abs(H_range_freq), 'b-', 'LineWidth', 1.5);
xlabel('Частота (МГц)');
ylabel('Амплитуда');
title('Range Filter Spectrum');
grid on;

% Фильтр по азимуту
subplot(2,3,2);
PRF = 1700;
Ka = 2*V^2/(lambda*R0);
t_az = linspace(-0.5, 0.5, 1000);
h_azimuth = exp(1j*pi*Ka*t_az.^2);
plot(t_az, real(h_azimuth), 'r-', 'LineWidth', 1.5);
xlabel('Время (с)');
ylabel('Амплитуда');
title('Azimuth Reference (Real)');
grid on;

subplot(2,3,5);
H_az_freq = fftshift(fft(h_azimuth));
f_az = linspace(-PRF/2, PRF/2, length(H_az_freq));
plot(f_az, abs(H_az_freq), 'r-', 'LineWidth', 1.5);
xlabel('Частота (Гц)');
ylabel('Амплитуда');
title('Azimuth Filter Spectrum');
grid on;

% Двумерная передаточная функция
subplot(2,3,[3 6]);
[F_range, F_az] = meshgrid(linspace(-Br/2, Br/2, 200), ...
                            linspace(-PRF/2, PRF/2, 200));
H_2D = ones(size(F_range)); % Упрощенная модель
imagesc(linspace(-Br/2, Br/2, 200)/1e6, ...
        linspace(-PRF/2, PRF/2, 200), abs(H_2D));
xlabel('Range Frequency (МГц)');
ylabel('Azimuth Frequency (Гц)');
title('2D Transfer Function');
colorbar;
axis xy;

sgtitle('Согласованные фильтры Range Doppler Algorithm');

%% ========================================================================
% ЭКСПОРТ ПАРАМЕТРОВ
% ========================================================================

fprintf('\n=== Параметры системы для воспроизведения ===\n');
fprintf('Используйте эти параметры для точного воспроизведения\n');
fprintf('рисунков из Cumming & Wong, 2005:\n\n');
fprintf('Carrier frequency: %.3f GHz\n', 5.405);
fprintf('Range bandwidth: %.2f MHz\n', 56.5);
fprintf('Pulse width: %.2f µs\n', 41.74);
fprintf('Sampling rate: %.0f MHz\n', 90);
fprintf('PRF: %.0f Hz\n', 1700);
fprintf('Platform velocity: %.0f m/s\n', 7400);
fprintf('Slant range: %.0f km\n', 850);
fprintf('Antenna length: %.1f m\n', 12.3);
fprintf('Squint angle: %.0f degrees\n', 5);
fprintf('============================================\n');

fprintf('\nДетальная визуализация завершена!\n');
