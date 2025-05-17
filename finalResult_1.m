clc; clear all; close all;
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 1540;               % Speed of sound[m/s]
Fc = 5e6;               % Transducer center frequency [Hz]
Fs = Fc * 4;            % Hardware sampling frequency (e.g., 20 MHz) [Hz]
Lambda = C / Fc;        % Wavelength [m]
Elem_pitch = 298e-6;    % Element pitch (distance between elements) [m]
N_element = 128;        % Total number of elements
N_RX_element = 128;     % Total number of receiving elements
N_scanline = 128;       % Number of scanlines
View_width = Elem_pitch * N_element; % Total view width [m] 
depth = 70e-3;          % Imaging depth (range of RF data) [m]
N_image_point = 1700;   % Number of image points (pixels) along depth
Unit_Dis = C / (2*Fs);  % Distance represented by one sample [m]

Inter_factor = 1;       % Interpolation factor for upsampling
Fs_new = Fs * Inter_factor; % New sampling frequency after interpolation [Hz]
F_Number_limit_RX = 2.0; % F-number limit for receiving aperture


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beamformer Implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sum_out_2D = zeros(N_image_point, N_scanline);

% Aperture settings (minimum and maximum aperture size)
min_aperture = 32;  % Minimum aperture size
max_aperture = N_RX_element;  % Maximum aperture size

for s = 1:N_scanline
    fid = fopen(sprintf('RF_data\\scanline_phantom%03d.bin', s), 'rb');
    rf_data = fread(fid, 'double');
    fclose(fid);

    Nt = length(rf_data) / N_RX_element; % Total number of time samples per element
    % Reshape RF data into a 2D array (time samples x elements)
    rf_data_2D = reshape(rf_data, [length(rf_data) / N_RX_element, N_RX_element]);

    time_orig = (0 : Nt - 1).' / Fs; % Original time axis based on sampling frequency Fs
    Nt_up = Nt * Inter_factor; % Number of time samples after interpolation
    time_up = (0 : Nt_up - 1).' / Fs_new; % Interpolated time axis based on upsampled frequency

    rf_data_2D_up = zeros(Nt_up, N_RX_element); 
    for ch = 1 : N_RX_element
        rf_data_2D_up(:, ch) = interp1( ...
            time_orig, ...
            rf_data_2D(:, ch), ...
            time_up, ...
            'linear', 0 ); % Linear interpolation
    end

    % Pre-compute aperture sizes for each depth
    apertures = zeros(N_image_point, 2);
    for i = 1:N_image_point
        imaging_depth = (i - 1) * Unit_Dis;
        dynamic_aperture = round(2 * imaging_depth / (F_Number_limit_RX * Elem_pitch));
        bounded_aperture_min = max(min_aperture, dynamic_aperture);
        aperture_size = min(bounded_aperture_min, max_aperture);

        rx_aperture_start = max(1, round(N_RX_element / 2 - aperture_size / 2));
        rx_aperture_end = min(N_RX_element, rx_aperture_start + aperture_size - 1);

        apertures(i, :) = [rx_aperture_start, rx_aperture_end];
    end

    % Generate Hann window for apodization
    window = hann(N_RX_element);

    Sum_out = zeros(N_image_point, 1);    % Initialize output for the current scanline
    z_axis = linspace(0, depth, N_image_point).';

    for iPoint = 1:N_image_point
        sum_val = 0;
        imaging_depth = (iPoint - 1) * Unit_Dis;

        % Use pre-computed aperture
        rx_aperture_start = apertures(iPoint, 1);
        rx_aperture_end = apertures(iPoint, 2);

        for m = rx_aperture_start:rx_aperture_end
            elem_position = (m - (N_RX_element / 2)) * Elem_pitch; % Element position
            lateral_distance = abs(elem_position);

            % Apply F-number limit and calculate delays
            if lateral_distance <= imaging_depth / F_Number_limit_RX
                distance_squared = imaging_depth^2 + elem_position^2;
                delay_samples = (imaging_depth / C + sqrt(distance_squared) / C) * Fs_new;
                sample_index_low = floor(delay_samples);
                sample_index_high = sample_index_low + 1;
                frac = delay_samples - sample_index_low;

                if sample_index_low >= 1 && sample_index_high <= size(rf_data_2D_up, 1)
                    val_low = rf_data_2D_up(sample_index_low, m);
                    val_high = rf_data_2D_up(sample_index_high, m);
                    interpolated_val = val_low * (1 - frac) + val_high * frac;
                    sum_val = sum_val + interpolated_val * window(m); % Apply window
                end
            end
        end
        Sum_out(iPoint) = sum_val;
    end

    % Store results in the 2D array for the current scanline
    Sum_out_2D(:, s) = Sum_out;
end  % end for s=1:N_scanline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Final 2D Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save 2D array of beamformed data
fid = fopen('Beamformed_data\\sum_out.bin', 'wb');
count = fwrite(fid, Sum_out_2D, 'double');   % 2D array (size: N_image_point x N_scanline)
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization of Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the beamformed RF data
figure;
imagesc(1:N_scanline, z_axis * 1e3, Sum_out_2D);
xlabel('Scanline Index');
ylabel('Depth (mm)');
title('Beamformed RF Data');
colormap(gray);
colorbar;
clim([-1000,1000]);
grid on;

figure;
plot(z_axis * 1e3, Sum_out_2D(:, 64));
xlabel('Depth (mm)');
ylabel('Amplitude');
title('Beamformed RF (Line #64) with 16x interpolation');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save I/Q Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IQ_data = zeros(size(Sum_out_2D));
for s = 1:N_scanline
    Sum_out_hilbert = hilbert(Sum_out_2D(:, s));
    analytic_signal = abs(Sum_out_hilbert);
    IQ_data(:, s) = analytic_signal;
end

fid_IQ = fopen('IQ_data.bin', 'wb');
fwrite(fid_IQ, IQ_data, 'double');
fclose(fid_IQ);

% Centerline result visualization
centerline_idx = ceil(N_scanline / 2);
figure;
plot(IQ_data(:, centerline_idx));
title('Centerline');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Spectrum Analysis using I + jQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scanline_idx = ceil(N_scanline / 2); % Select center scanline
I = real(Sum_out_2D(:, scanline_idx)); % Extract real part (I)
Q = imag(hilbert(Sum_out_2D(:, scanline_idx))); % Extract imaginary part (Q)

IQ_complex = I + 1i * Q; % Create complex signal

% Compute FFT
freq_axis = linspace(-Fs_new / 2, Fs_new / 2, length(IQ_complex));
IQ_spectrum = fftshift(fft(IQ_complex)); % Shift zero frequency to center

% Plot Frequency Spectrum
figure;
plot(freq_axis / 1e6, abs(IQ_spectrum)); % Convert frequency to MHz
xlabel('Frequency (MHz)');
ylabel('Magnitude');
title('Frequency Spectrum of I + jQ');
grid on;
