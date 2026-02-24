function S_corrected = apply_src_correction(S_RD, f_range, f_azimuth, params)
% =========================================================================
% Secondary Range Compression (SRC)
% =========================================================================
% Функция применяет вторичное сжатие по дальности для компенсации
% азимутально-зависимого расфокусирования
%
% Вход:
%   S_RD - данные в Range-Doppler области [N_range x N_azimuth]
%   f_range - вектор частот по дальности, Гц
%   f_azimuth - вектор частот по азимуту, Гц
%   params - структура параметров SAR
%
% Выход:
%   S_corrected - данные после SRC [N_range x N_azimuth]
%
% Референции:
%   [1] Cumming & Wong, 2005, eq. 6.25, Section 6.3
%   [2] Sentinel-1 Algorithm, Section 6.3.1
% =========================================================================

[N_range, N_azimuth] = size(S_RD);

% Извлечение параметров
R_0 = params.R_0;
V_r = params.V_r;
lambda = params.lambda;

% Инициализация выходного массива
S_corrected = S_RD;

% Применение SRC для каждой азимутальной частоты
for k = 1:N_azimuth
    f_a = f_azimuth(k);
    
    % Эффективная скорость [1, eq.6.21]
    % V_eff = V_r * sqrt(1 - (λ*f_a/(2*V_r))²)
    f_norm = (lambda * f_a) / (2 * V_r);
    
    % Проверка допустимого диапазона
    if abs(f_norm) >= 1
        continue;
    end
    
    V_eff = V_r * sqrt(1 - f_norm^2);
    
    % SRC фазовая функция [1, eq.6.25]
    % Φ_src(f_r, f_a) = -π * λ² * R_0 * f_r² / (8 * V_eff²)
    phase_src = -pi * lambda^2 * R_0 * f_range.^2 / (8 * V_eff^2);
    
    % Применение фазовой коррекции
    H_src = exp(1j * phase_src);
    S_corrected(:, k) = S_RD(:, k) .* H_src;
end

% Примечание: [2] использует более точную формулу для Sentinel-1:
% Φ_src = -π * λ² * R_ref * f_r² * D(f_a) / (8 * V_eff²)
% где D(f_a) - дополнительный корректирующий член

end
