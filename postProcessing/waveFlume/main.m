%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read free-surface elevation data and compute error with Stokes theory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% READ DATA
filename = "surfaceElevation_ref4.dat";
addpath('wave_data')

raw_data = table2array(readtable(filename));

idx=find(raw_data(:,1)>0); %Gives indices of positive values of time
data = raw_data(idx,:);
time = data(:,1); % Time vector
surf_elevation = data(:,2:end); % Surface elevation
%plot(time,surf_elevation(:,1))
%% WAVE PARAMETERS
wave_period = 8; % Wave period [s]
height_ratio = 0.05; % Steepness, wave height to wavelength ratio
depth_ratio = 0.4; %Water depth to length ratio
wave_length = 9.81/(2*pi)*wave_period^2*tanh(2*pi*depth_ratio); % Dispersion relation

wave_height = height_ratio*wave_length; % Wave height [m]
water_depth = depth_ratio*wave_length; % Water depth [m]
wave_freq = 2*pi/wave_period; %Wave frequency [rad/s]

%% CONSTRUCT REFERENCE WAVE
n_periods = time(end)/wave_period; % Number of periods in the data
n_samples = length(time)/n_periods; % Mean number of samples per period
ref_time = linspace(0,wave_period,n_samples); % Reference time for cosine wave
gauges_pos = linspace(0,2*wave_length,size(surf_elevation,2)); %Dimensionless position of the gauges (x/lambda)

% Compute second-order phase average
ref_phase_avg = surfaceElevation_2nd(0, ref_time, wave_length, wave_period, wave_height, water_depth);

% Compute the Stokes second-order surface elevation at evry position and time
for i=1:length(time)
    ref_wave_matrix(i,:) = surfaceElevation_2nd(gauges_pos, time(i), wave_length, wave_period, wave_height, water_depth);
end
figure(1)
plot(time,surf_elevation(:,250))
hold on
plot(time,ref_wave_matrix(:,250))
legend('numerical mwl', 'analytical mwl')
xlabel('t/T');
ylabel('mean water level [m]')
title('8 Cells Per Wave Height');

figure
nwt=size(surf_elevation);
nwt_length=linspace(1,nwt(2),500);
plot(nwt_length,surf_elevation(3200,:))
hold on
plot(nwt_length,ref_wave_matrix(3200,:))
%% COMPUTE ERROR ALONG THE DOMAIN
error_amp = zeros(size(surf_elevation,2),1);
error_harm = zeros(size(surf_elevation,2),1);

for i=1:size(surf_elevation,2)
     [phase_avg(:,i)] = averagePhase(surf_elevation(:,i),0,n_samples);
     % nRMSE error, normaslised with wave height
     nrmse(i) = (100/wave_height)*sqrt(mean((phase_avg(:,i)' - ref_phase_avg).^2)); 
     % relative error, normaslised with wave height
     rel_error(:,i) = (100/wave_height)*(surf_elevation(:,i) - ref_wave_matrix(:,i));
     % relative error of the phase average, normaslised with wave height
     avg_rel_error(:,i) = (1/wave_height)*(phase_avg(:,i)' - ref_phase_avg);
end
     
   plot(ref_time,avg_rel_error(:,1))
   hold on
   
[phase_avg1(:,1)] = averagePhase(surf_elevation1(:,1),0,n_samples);
avg_rel_error1(:,1) = (1/wave_height)*(phase_avg1(:,1)' - ref_phase_avg);
plot(ref_time,avg_rel_error1(:,1))
hold on

[phase_avg2(:,1)] = averagePhase(surf_elevation2(:,1),0,n_samples);
avg_rel_error2(:,1) = (1/wave_height)*(phase_avg2(:,1)' - ref_phase_avg);
plot(ref_time,avg_rel_error2(:,1))
hold on

[phase_avg3(:,1)] = averagePhase(surf_elevation3(:,1),0,n_samples);
avg_rel_error3(:,1) = (1/wave_height)*(phase_avg3(:,1)' - ref_phase_avg);
plot(ref_time,avg_rel_error3(:,1))

legend('5 cpwh','10 cpwh','15 cpwh', '20 cpwh')
xlabel('t/T')
ylabel('\deltap/H')
%% amplitude error, 4th peak, wave gauge number 250
    
real_peak=max(ref_wave_matrix(:,250));
stokes_maxima=[];
position_of_maxima=[];
stokes_maxima_new=[];
for i=2:(length(surf_elevation)-1)
            if surf_elevation(i,250)>surf_elevation(i-1,250) && surf_elevation(i,250)>surf_elevation(i+1,250)
        
            stokes_maxima=[stokes_maxima, surf_elevation(i,250)];
            position_of_maxima=[position_of_maxima,i];
        end
end

 ampl_error=((surf_elevation(position_of_maxima(1),250)/ref_wave_matrix(position_of_maxima(1),250)) ...
     + (surf_elevation(position_of_maxima(2),250)/ref_wave_matrix(position_of_maxima(2),250)) ...
     + (surf_elevation(position_of_maxima(3),250)/ref_wave_matrix(position_of_maxima(3),250)) ...
     + (surf_elevation(position_of_maxima(4),250)/ref_wave_matrix(position_of_maxima(4),250)) ...
     + (surf_elevation(position_of_maxima(5),250)/ref_wave_matrix(position_of_maxima(5),250)) ...
     + (surf_elevation(position_of_maxima(6),250)/ref_wave_matrix(position_of_maxima(6),250)))/6 ...
     %+ (surf_elevation(position_of_maxima(7),250)/ref_wave_matrix(position_of_maxima(7),250)) ...
     %+ (surf_elevation(position_of_maxima(8),250)/ref_wave_matrix(position_of_maxima(8),250)))/8;
% position_of_maxima=[1, position_of_maxima];
ampl_rror4th=(surf_elevation(position_of_maxima(4),250)/ref_wave_matrix(position_of_maxima(4),250));
% rel_amplitude_error=[(sqrt((real_peak-surf_elevation(1,1)).^2))/wave_height, (sqrt((real_peak-stokes_maxima).^2))/wave_height];
% relative_time=[1:n_periods];
% figure(2)
% h=semilogy(relative_time,rel_amplitude_error,'-o');
% hold on
% h.LineStyle = '-';
% axis([1 8 10^(-5) 1])
% xlabel('t/T');
% ylabel('log(\deltaa/H)')
% %title('15 Cells Per Wave Height');
%% amplitude error all 4  cpwh (5 10 15 20) together
real_peak=max(ref_wave_matrix(:,1));
stokes_maxima=[];
position_of_maxima=[];
for i=2:length(surf_elevation(:,1)-1)
    if surf_elevation(i,1)>surf_elevation(i-1) && surf_elevation(i,1)>surf_elevation(i+1)
        
        stokes_maxima=[stokes_maxima, surf_elevation(i,1)];
        position_of_maxima=[position_of_maxima,i];
    end
end
position_of_maxima=[1, position_of_maxima];

rel_amplitude_error=[(sqrt((real_peak-surf_elevation(1,1)).^2))/wave_height, (sqrt((real_peak-stokes_maxima).^2))/wave_height];
relative_time=[1:n_periods];
figure()
h=semilogy(relative_time,rel_amplitude_error,'-o');
hold on

real_peak1=max(ref_wave_matrix(:,1));
stokes_maxima1=[];
position_of_maxima1=[];
for i=2:length(surf_elevation1(:,1)-1)
    if surf_elevation1(i,1)>surf_elevation1(i-1) && surf_elevation1(i,1)>surf_elevation1(i+1)
        
        stokes_maxima1=[stokes_maxima1, surf_elevation1(i,1)];
        position_of_maxima1=[position_of_maxima1,i];
    end
end

position_of_maxima1=[1, position_of_maxima1];

rel_amplitude_error1=[(sqrt((real_peak1-surf_elevation1(1,1)).^2))/wave_height, (sqrt((real_peak1-stokes_maxima1).^2))/wave_height];
relative_time=[1:n_periods];
h=semilogy(relative_time,rel_amplitude_error1,'-o');
hold on

real_peak2=max(ref_wave_matrix(:,1));
stokes_maxima2=[];
position_of_maxima2=[];
for i=2:length(surf_elevation2(:,1)-1)
    if surf_elevation2(i,1)>surf_elevation2(i-1) && surf_elevation2(i,1)>surf_elevation2(i+1)
        
        stokes_maxima2=[stokes_maxima2, surf_elevation2(i,1)];
        position_of_maxima2=[position_of_maxima2,i];
    end
end

position_of_maxima2=[1, position_of_maxima2];

rel_amplitude_error2=[(sqrt((real_peak2-surf_elevation2(1,1)).^2))/wave_height, (sqrt((real_peak2-stokes_maxima2).^2))/wave_height];
relative_time=[1:n_periods];
h=semilogy(relative_time,rel_amplitude_error2,'-o');


real_peak3=max(ref_wave_matrix(:,1));
stokes_maxima3=[];
position_of_maxima3=[];
for i=2:length(surf_elevation3(:,1)-1)
    if surf_elevation3(i,1)>surf_elevation3(i-1) && surf_elevation3(i,1)>surf_elevation3(i+1)
        
        stokes_maxima3=[stokes_maxima3, surf_elevation3(i,1)];
        position_of_maxima3=[position_of_maxima3,i];
    end
end

position_of_maxima3=[1, position_of_maxima3];

rel_amplitude_error3=[(sqrt((real_peak3-surf_elevation3(1,1)).^2))/wave_height, (sqrt((real_peak3-stokes_maxima3).^2))/wave_height];
relative_time=[1:n_periods];
h=semilogy(relative_time,rel_amplitude_error3,'-o');


xlabel('t/T');
ylabel('log(\deltaa/H)')
axis([1 8 10^(-5) 1])
legend('5 cpwh','10 cpwh','15 cpwh','20 cpwh')

%% phase error at the peaks of each wave period, it is equal in all 3 cases (maybe wrong)
% Find peaks in numerical and analytical data sets
[~, num_peak_locs] = findpeaks(surf_elevation(:,1), 'MinPeakDistance', round(n_samples/2));
[~, ref_peak_locs] = findpeaks(ref_wave_matrix(:,1), 'MinPeakDistance', round(n_samples/2));

% Calculate phase error at each peak
phase_error = zeros(1, length(num_peak_locs));
for i = 1:length(num_peak_locs)
    num_peak_time = time(num_peak_locs(i));
    ref_peak_time = time(ref_peak_locs(i));
    phase_error(i) = (num_peak_time - ref_peak_time) * wave_freq;
end

% Find maximum phase error and print result
[max_phase_error, max_idx] = max(abs(phase_error));
fprintf("Maximum phase error: %f radians at time t = %f s\n", max_phase_error, time(num_peak_locs(max_idx)));

