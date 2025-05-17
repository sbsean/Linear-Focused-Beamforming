clc; clear all; close all;
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 1540; % ���� [m/s]                                             % Speed of sound [m/s]
Fc = 5e6;  % Ʈ�����༭ �߽����ļ� [Hz]                                            % Transducer center frequency [Hz]
Fs = Fc * 4; % �ϵ���� ���ø� ���ļ� (��: 20 MHz)                                          % Sampling frequency [Hz]
Lambda = C / Fc;                                             % lambda
Elem_pitch = 298e-6;                                           % Element Pitch [m]
N_element = 128;                                              % Total # of elements
N_RX_element = 128;                                           % Total # of Rx. elements
N_scanline = 128;                                             % number of scanline    
View_width = Elem_pitch * N_element;
depth = 70e-3;                                                % RF data depth [m]
N_image_point = 1700;
Unit_Dis = C / Fs;

Inter_factor = 4;
Fs_new = Fs * Inter_factor;
F_Number_limit_RX = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �� ���Ӹ� �����Ͻÿ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sum_out_2D = zeros(N_image_point, N_scanline);

for s = 1:N_scanline
    fid = fopen(sprintf('RF_data\\scanline_phantom%03d.bin', s), 'rb');
    rf_data = fread(fid, 'double');
    fclose(fid);

    Nt = length(rf_data) / N_RX_element;
    rf_data_2D = reshape(rf_data, [length(rf_data) / N_RX_element, N_RX_element]);

    % figure; imagesc(rf_data_2D); caxis([-5e2 5e2]); colormap(gray(256)); colorbar;

    time_orig = (0 : Nt - 1).' / Fs;
    Nt_up = Nt * Inter_factor;
    time_up = (0 : Nt_up - 1).' / Fs_new;

    rf_data_2D_up = zeros(Nt_up, N_RX_element);
    for ch = 1 : N_RX_element
        rf_data_2D_up(:, ch) = interp1( ...
            time_orig, ...
            rf_data_2D(:, ch), ...
            time_up, ...
            'linear', 0 );
    end
    
    %--------------------------------------------------------------------------
    % (3) Delay & Sum ������ (���� ��Ŀ��)
    %--------------------------------------------------------------------------
    Sum_out = zeros(N_image_point, 1);    % �� ��ĵ������ ���(1D)
    
    z_axis = linspace(0, depth, N_image_point).';  
    elem_pos = ((1:N_RX_element) - (N_RX_element + 1) / 2) * Elem_pitch;
    
    for iPoint = 1 : N_image_point
        z_val = z_axis(iPoint);

        %----- F-Number �������� Aperture ���� ���� -----
        aperture_half = z_val / (2 * F_Number_limit_RX);
        valid_idx = find(abs(elem_pos) <= aperture_half);

        tmp_sum = 0;
        for m = valid_idx
            d_m = sqrt(z_val^2 + elem_pos(m)^2);
            tau_m = d_m / C;
            sample_index = tau_m * Fs_new;  % �����õ� Fs ����

            % ������ ���
            if sample_index >= 1 && sample_index <= Nt_up
                val = interp1(time_up, rf_data_2D_up(:, m), sample_index, 'linear', 0 );
            else
                val = 0;
            end
            tmp_sum = tmp_sum + val;
        end
        
        Sum_out(iPoint) = tmp_sum;
    end
    
    %--------------------------------------------------------------------------
    % (4) 2D �迭�� ���� (Sum_out_2D�� s��° ����)
    %--------------------------------------------------------------------------
    Sum_out_2D(:, s) = Sum_out;
    
end  % end for s=1:N_scanline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5) ���� 2D ��� ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���� �ڵ尡 'sum_out.bin' ���� �����ϱ� ���Ѵٸ�, ���⼭ 2D ���¸� �� ���� ����
fid = fopen('Beamformed_data\\sum_out.bin', 'wb');
count = fwrite(fid, Sum_out_2D, 'double');   % 2D �迭 (ũ��: N_image_point x N_scanline)
fclose(fid);

% (����) �ʿ��ϴٸ�, Sum_out_2D �� "Sum_out_realigned" �̸����� �Ѱܵ� ��:
% Sum_out_realigned = Sum_out_2D;
% fid = fopen('Beamformed_data\\sum_out.bin','wb');
% fwrite(fid, Sum_out_realigned, 'double');
% fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6) Ȯ�ο� Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����: ��� ����(s=64) ����� ���������� �÷�
figure;
plot(z_axis * 1e3, Sum_out_2D(:, 64));
xlabel('Depth (mm)');
ylabel('Amplitude');
title('Beamformed RF (Line #64) with 16x interpolation');
grid on;

% end

fid = fopen('Beamformed_data\\sum_out.bin', 'wb');
count = fwrite(fid, Sum_out_2D, 'double');   % ���� �����ֵ� RF ������  
fclose(fid);
