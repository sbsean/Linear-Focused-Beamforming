clc;clear all;close all;
format long g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                       = 1540;                                             % Speed of sound [m/s]
Fc                      = 5e6;                                              % Transducer center frequency [Hz]
Fs                      = Fc * 4;                                           % Sampling frequency [Hz]
Lambda                  = C/Fc;                                             % lambda
Elem_pitch              = 298e-6;                                           % Element Pitch [m]
View_lem_pitch          = 298e-6;                                           % Element Pitch [m]
N_element               = 128;                                              % Total # of elements
N_RX_element            = 128;                                              % Total # of Rx. elements
N_scanline              = 128;                                              % number of scanlinewidth              
View_width              = Elem_pitch * N_element;
depth                   = 70e-3;                                            % RF data depth [m]
N_image_point           = 1700;
Dynamic_Range           = 60;                                               % [dB]

Num_pixel_z             = round(depth * 10000);                             % [sample]
Num_pixel_x             = round(View_width * 10000);                        % [sample]


%% Mid & Backend Processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RF data file open
% 빔포밍 결과를 불러오기
fid = fopen('Beamformed_data/sum_out.bin', 'rb');
Sum_out_temp = fread(fid,'double');
fclose all;

Sum_out_temp = reshape(Sum_out_temp, length(Sum_out_temp)/N_scanline, N_scanline);
Sum_out = Sum_out_temp(1:N_image_point,:)';
% Sum_out = Sum_out./max(max(abs(Sum_out)))*2^20;

figure; imagesc(Sum_out'); caxis([-1e3 1e3]); colormap(gray(256)); colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% DC Cancelation
DC_Cancel_out=zeros(N_scanline,N_image_point);
DC_coef = fir1(32,[0.1],'high');
DC_Cancel_out = convn(Sum_out,DC_coef,'same');

disp(sprintf(sprintf('Dc Cancelation\n')));

% Time Gain Compensation
TGC_out=zeros(N_scanline,N_image_point);
TGC_out = DC_Cancel_out;

disp(sprintf('Time Gain Compensation\n'));

% Quadrature Demodulation
Data_I=zeros(N_scanline,N_image_point);
Data_Q=zeros(N_scanline,N_image_point);
Data_I_LPF=zeros(N_scanline,N_image_point);
Data_Q_LPF=zeros(N_scanline,N_image_point);

Cos_t = cos((1:1:N_image_point)*2*pi*Fc/Fs);
Sin_t = sin((1:1:N_image_point)*2*pi*Fc/Fs);

for i=1:1:N_scanline
    Data_I(i,1:1:N_image_point) = TGC_out(i,1:1:N_image_point) .* Cos_t(1,1:1:N_image_point);
    Data_Q(i,1:1:N_image_point) = TGC_out(i,1:1:N_image_point) .* Sin_t(1,1:1:N_image_point);
end

LPF_coef = fir1(48,Fc/Fs);

Data_I_LPF = convn(Data_I,LPF_coef,'same');
Data_Q_LPF = convn(Data_Q,LPF_coef,'same');

disp(sprintf('Quadrature Demodulation\n'));

% Envelope Detection
Data_Meg=zeros(N_scanline,N_RX_element);
Data_Meg = sqrt(Data_I_LPF.^2 + Data_Q_LPF.^2);

disp(sprintf('Envelope Detection\n'));

% Log Compression
Log_out=zeros(N_scanline,N_RX_element);
Ymax = 2^8;
Xmax = max(max(Data_Meg));
Xmin = Xmax*10^(-Dynamic_Range/20);
Log_out = (Data_Meg >= Xmin).* (Ymax/(log10(Xmax/Xmin))*log10(Data_Meg/Xmin));

disp(sprintf('Log Compression\n'));

%% Scan converter
dsc_in=zeros(N_RX_element,N_scanline);
dsc_in_trunc=zeros(N_image_point,N_scanline);
sc_data = zeros(Num_pixel_z,Num_pixel_x);
dsc_in = Log_out';

dsc_in_trunc = dsc_in(1:N_image_point,:);
    
size_data=size(dsc_in_trunc);
Num_data_z = size_data(1);
Num_data_x = size_data(2);

Axis_data_z = ([1:1:Num_data_z] * N_image_point / (Num_data_z + 1)) ;
Axis_data_x = ([1:1:Num_data_x] * View_width / (Num_data_x + 1)) ;
    
Axis_image_z = ([1:1:Num_pixel_z] * N_image_point / (Num_pixel_z + 1)) ;
Axis_image_x = ([1:1:Num_pixel_x] * View_width / (Num_pixel_x + 1)) ;
    
for i = 1:1:Num_pixel_z
    if  (Axis_data_z < Axis_image_z(i))
        cordi_down(i) = Num_data_z + 1;
    else
        cordi_down(i) = find(Axis_data_z>=Axis_image_z(i),1);
    end
end
    
cordi_up  = cordi_down - 1;
 
for i = 1:1:Num_pixel_x
    if  (Axis_data_x < Axis_image_x(i))
        cordi_right(i) = Num_data_x + 1;
    else
        cordi_right(i) = find(Axis_data_x>=Axis_image_x(i),1);
    end
end
    
cordi_left = cordi_right - 1;
    
% boundary treatment

[index_left_0 ]=find(cordi_left == 0);
cordi_left(index_left_0)=1;
index_right_max=find(cordi_right > Num_data_x);
cordi_right(index_right_max)=Num_data_x;

 for z = 1:1:Num_pixel_z
    for x = 1:1:Num_pixel_x
        D = dsc_in_trunc(cordi_up(z),cordi_left(x));
        A = dsc_in_trunc(cordi_down(z),cordi_left(x));    
        a = Axis_data_z(cordi_up(z));
        c = Axis_data_z(cordi_down(z));
        b = Axis_image_z(z);    
        left_value(z,x) = ((D-A)/(a-c))*(b-c)+A; 
        D = dsc_in_trunc(cordi_up(z),cordi_right(x));
        A = dsc_in_trunc(cordi_down(z),cordi_right(x));        
        right_value(z,x) = ((D-A)/(a-c))*(b-c)+A;   
        D = left_value(z,x);
        A = right_value(z,x);    
        a = Axis_data_x(cordi_left(x));
        c = Axis_data_x(cordi_right(x));
        b = Axis_image_x(x);        
        sc_data(z,x) = ((D-A)/(a-c))*(b-c)+A; 
    end
 end

% filling bouandary
for i=1:1:max(index_left_0)
    sc_data(:,i) = sc_data(:,(max(index_left_0)+1));
end
for j=index_right_max(1):1:max(index_right_max)
    sc_data(:,j) = sc_data(:,min(index_right_max)-1);
end

    h = figure(1004);
    set(h,'color',[1 1 1]);
    
    x_axis = View_width*100;
    x_axis_temp = linspace(-1*x_axis/2,x_axis/2,length(sc_data(:,1)));
    z_axis = (N_image_point)/512;
    z_axis_temp = linspace(0,z_axis,length(sc_data(1,:)));
    
    map = [0:255; 0:255; 0:255]'/255;
    imagesc(x_axis_temp,z_axis_temp,sc_data),colormap(map),truesize;
    xlabel('[cm]'); ylabel('[cm]')