% Load the beamformed data
fid = fopen('Beamformed_data\\sum_out.bin', 'rb');
Sum_out_2D = fread(fid, [N_image_point, N_scanline], 'double');
fclose(fid);

% Apply Hilbert Transform to generate I + jQ signal
IQ_data = hilbert(Sum_out_2D); % I + jQ 형태의 복소 신호

% Define analysis depth and lateral position
target_depth = 1.5e-2; % Focus on 1.5 cm depth
[~, depth_idx] = min(abs(z_axis - target_depth)); % Locate depth index
lateral_position = ((1:N_scanline) - N_scanline/2) * Elem_pitch; % In meters

% Define window functions
window_functions = {
    'Rectangular', rectwin(N_RX_element);
    'Hanning', hann(N_RX_element);
    'Hamming', hamming(N_RX_element);
    'Gaussian', gausswin(N_RX_element, 2.5); % Alpha = 2.5
    'Blackman', blackman(N_RX_element);
    'Kaiser', kaiser(N_RX_element, 3); % Beta = 3
};

% Initialize result storage
results = struct();

% Create individual plots for each window function
figure;
for w = 1:size(window_functions, 1)
    window_name = window_functions{w, 1};
    window = window_functions{w, 2};

    % Apply window to the beamformed data
    lateral_profile = abs(IQ_data(depth_idx, :)) .* (window' / max(window));
    lateral_profile = lateral_profile / max(lateral_profile); % Normalize
    lateral_profile_dB = 20 * log10(lateral_profile); % Convert to dB
    
    % Plot each window function result in a separate subplot
    subplot(2, 3, w); % Arrange subplots in 2 rows and 3 columns
    plot(lateral_position * 1e2, lateral_profile_dB); % Convert lateral to cm
    xlabel('Lateral Position (cm)');
    ylabel('Amplitude (dB)');
    title(sprintf('%s Window', window_name));
    xlim([-1.5 1.5]); % Focus on Mainlobe and closest Sidelobes
    ylim([-60 10]); % Adjust dynamic range
    grid on;

    % Mainlobe Width (Full Width at Half Maximum)
    [~, max_idx] = max(lateral_profile); % Find the peak
    half_max = max(lateral_profile) / sqrt(2); % Half maximum (FWHM criterion)
    left_idx = find(lateral_profile(1:max_idx) <= half_max, 1, 'last'); % Left side
    right_idx = find(lateral_profile(max_idx:end) <= half_max, 1, 'first') + max_idx - 1; % Right side
    mainlobe_width = (lateral_position(right_idx) - lateral_position(left_idx)) * 1e2; % In cm

    % Sidelobe Level (Second highest peak)
    sidelobe_profile = lateral_profile;
    sidelobe_profile(left_idx:right_idx) = 0; % Ignore mainlobe region
    sidelobe_level = max(20 * log10(sidelobe_profile(sidelobe_profile > 0))); % Max sidelobe in dB

    % Sidelobe Suppression Ratio (SSR)
    ssr = max(lateral_profile_dB) - sidelobe_level; % Difference between mainlobe and sidelobe

    % Store results
    results(w).Window = window_name;
    results(w).MainlobeWidth = mainlobe_width; % in cm
    results(w).SidelobeLevel = sidelobe_level; % in dB
    results(w).SSR = ssr; % in dB
end

% Add a title for the entire figure
sgtitle(sprintf('Lateral Profile Comparison at Depth = %.1f cm (I + jQ Signal)', target_depth * 1e2));

% Display numerical results
fprintf('Results at Depth = %.1f cm:\n', target_depth * 1e2);
fprintf('%-12s | %-15s | %-15s | %-15s\n', 'Window', 'Mainlobe Width (cm)', 'Sidelobe Level (dB)', 'SSR (dB)');
fprintf('---------------------------------------------------------------\n');
for w = 1:length(results)
    fprintf('%-12s | %-15.3f | %-15.3f | %-15.3f\n', ...
        results(w).Window, results(w).MainlobeWidth, results(w).SidelobeLevel, results(w).SSR);
end

