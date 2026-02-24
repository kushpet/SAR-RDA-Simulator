function [f_dc, correlation] = estimate_doppler_centroid(data, PRF, method)
% =========================================================================
% Оценка Doppler Centroid
% =========================================================================
% Функция для оценки доплеровского центроида по различным методам
%
% Вход:
%   data - комплексные данные после Range Compression [N_range x N_azimuth]
%   PRF - частота повторения импульсов, Гц
%   method - метод оценки: 'energy', 'correlation', 'clutter'
%
% Выход:
%   f_dc - оценка Doppler centroid, Гц
%   correlation - функция корреляции (для метода 'correlation')
%
% Референции:
%   [1] Cumming & Wong, 2005, Section 4.6
%   [2] Sentinel-1 Algorithm Definition, Section 5.5
% =========================================================================

if nargin < 3
    method = 'energy';
end

[N_range, N_azimuth] = size(data);

switch lower(method)
    case 'energy'
        % Энергетический метод [1, eq.4.62]
        % f_dc находится как центр масс спектра мощности
        
        % FFT по азимуту
        spectrum = fftshift(fft(data, [], 2), 2);
        
        % Усреднение по дальности
        azimuth_spectrum = mean(abs(spectrum).^2, 1);
        
        % Частотная ось
        f_axis = linspace(-PRF/2, PRF/2, N_azimuth);
        
        % Центр масс
        f_dc = sum(f_axis .* azimuth_spectrum) / sum(azimuth_spectrum);
        
        correlation = [];
        
    case 'correlation'
        % Корреляционный метод [1, Section 4.6.2]
        % Использует корреляцию между соседними импульсами
        
        % Корреляция соседних азимутальных линий
        corr_sum = 0;
        for i = 1:N_range
            % Корреляция для каждой дальности
            line1 = data(i, 1:end-1);
            line2 = data(i, 2:end);
            corr_sum = corr_sum + sum(line1 .* conj(line2));
        end
        
        % Фаза корреляции связана с Doppler centroid
        phase = angle(corr_sum);
        f_dc = phase * PRF / (2*pi);
        
        correlation = corr_sum;
        
    case 'clutter'
        % Метод по отражению от подстилающей поверхности
        % [2, Section 5.5.2]
        % Упрощенная реализация
        
        % Выбираем центральную часть по дальности
        center_idx = round(N_range/2);
        range_window = max(1, center_idx-50):min(N_range, center_idx+50);
        
        % Данные центрального участка
        data_center = data(range_window, :);
        
        % FFT по азимуту
        spectrum = fftshift(fft(data_center, [], 2), 2);
        
        % Усреднение
        azimuth_spectrum = mean(abs(spectrum).^2, 1);
        
        % Поиск максимума
        [~, max_idx] = max(azimuth_spectrum);
        f_axis = linspace(-PRF/2, PRF/2, N_azimuth);
        f_dc = f_axis(max_idx);
        
        correlation = [];
        
    otherwise
        error('Неизвестный метод. Используйте: energy, correlation, или clutter');
end

end
