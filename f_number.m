clc; clear all; close all;
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 1540;               % Speed of sound [m/s]
Fc = 5e6;               % Transducer center frequency [Hz]
Fs = Fc * 4;            % Hardware sampling frequency (e.g., 20 MHz) [Hz]
Lambda = C / Fc;        % Wavelength [m]
Elem_pitch = 298e-6;    % Element pitch [m]
N_RX_element = 128;     % Total number of receiving elements
N_scanline = 128;       % Number of scanlines
depth = 70e-3;          % Imaging depth [m]
N_image_point = 1700;   % Number of image points along depth
Unit_Dis = C / (2 * Fs); % Distance represented by one sample [m]

Inter_factor = 1;       % Interpolation factor for upsampling
Fs_new = Fs * Inter_factor; % New sampling frequency after interpolation [Hz]

% F-number settings
F_numbers = [1.0, 1.5, 2.0]; % Analyze these F-numbers
near_field_limit = 1e-2; % Depth limit for near field (1 cm)

% Storage for results
B_mode_images = cell(1, length(F_numbers));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beamformer Implementation with Dynamic F-number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f_idx = 1:length(F_numbers)
    F_Number_limit_RX = F_numbers(f_idx);
    Sum_out_2D = zeros(N_image_point, N_scanline);

    % Aperture settings
    min_aperture = 32;
    max_aperture = N_RX_element;

    for s = 1:N_scanline
        % Load RF data
        fid = fopen(sprintf('RF_data\\scanline_phantom%03d.bin', s), 'rb');
        rf_data = fread(fid, 'double');
        fclose(fid);

        Nt = length(rf_data) / N_RX_element;
        rf_data_2D = reshape(rf_data, [Nt, N_RX_element]);

        time_orig = (0:Nt-1).' / Fs;
        Nt_up = Nt * Inter_factor;
        time_up = (0:Nt_up-1).' / Fs_new;

        rf_data_2D_up = zeros(Nt_up, N_RX_element);
        for ch = 1:N_RX_element
            rf_data_2D_up(:, ch) = interp1(time_orig, rf_data_2D(:, ch), time_up, 'linear', 0);
        end

        delays = zeros(N_image_point, N_RX_element);
        for i = 1:N_image_point
            imaging_depth = (i - 1) * Unit_Dis;
            if imaging_depth > near_field_limit
                continue; % Only process near field
            end

            % Dynamic aperture configuration
            dynamic_aperture = round(2 * imaging_depth / (F_Number_limit_RX * Elem_pitch));
            bounded_aperture_min = max(min_aperture, dynamic_aperture);
            aperture_size = min(bounded_aperture_min, max_aperture);

            rx_aperture_start = max(1, round(N_RX_element / 2 - aperture_size / 2));
            rx_aperture_end = min(N_RX_element, rx_aperture_start + aperture_size - 1);

            for j = rx_aperture_start:rx_aperture_end
                elem_position = (j - (N_RX_element / 2)) * Elem_pitch;
                lateral_distance = abs(elem_position);

                if lateral_distance <= imaging_depth / F_Number_limit_RX
                    distance_squared = imaging_depth^2 + elem_position^2;
                    delays(i, j) = imaging_depth / C + sqrt(distance_squared) / C;
                else
                    delays(i, j) = Inf;
                end
            end
        end

        Sum_out = zeros(N_image_point, 1);
        for iPoint = 1:N_image_point
            sum_val = 0;
            for m = 1:N_RX_element
                if isinf(delays(iPoint, m))
                    continue;
                end
                delay_samples = delays(iPoint, m) * Fs_new;
                sample_index_low = floor(delay_samples);
                sample_index_high = sample_index_low + 1;
                frac = delay_samples - sample_index_low;

                if sample_index_low >= 1 && sample_index_high <= size(rf_data_2D_up, 1)
                    val_low = rf_data_2D_up(sample_index_low, m);
                    val_high = rf_data_2D_up(sample_index_high, m);
                    interpolated_val = val_low * (1 - frac) + val_high * frac;
                    sum_val = sum_val + interpolated_val;
                end
            end
            Sum_out(iPoint) = sum_val;
        end
        Sum_out_2D(:, s) = Sum_out;
    end

    % Envelope detection and log compression
    Data_Meg = abs(hilbert(Sum_out_2D'));
    Dynamic_Range = 60; % dB
    Ymax = 255;
    Xmax = max(Data_Meg(:));
    Xmin = Xmax * 10^(-Dynamic_Range / 20);
    Log_out = (Data_Meg >= Xmin) .* (Ymax / log10(Xmax / Xmin) * log10(Data_Meg / Xmin));

    % Calculate near field index (maximum depth = 1 cm)
    near_field_idx = round(near_field_limit / Unit_Dis);
    near_field_idx = min(near_field_idx, size(Log_out, 1));

    % Extract near field data for B-mode image
    B_mode_images{f_idx} = Log_out(1:near_field_idx, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization of B-mode Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for f_idx = 1:length(F_numbers)
    subplot(1, length(F_numbers), f_idx);
    imagesc(1:N_scanline, (1:near_field_idx) * Unit_Dis * 1e3, B_mode_images{f_idx});
    xlabel('Scanline Index');
    ylabel('Depth (mm)');
    title(sprintf('F-number = %.1f', F_numbers(f_idx)));
    colormap(gray);
    caxis([-1000,1000]); % Dynamic range
    colorbar;
end
sgtitle('B-mode Images with Different F-numbers (Depth < 1 cm)');
