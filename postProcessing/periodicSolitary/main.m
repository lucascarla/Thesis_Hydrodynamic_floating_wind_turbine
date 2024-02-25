%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read free-surface elevation data and compute error with Stokes theory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% READ DATA
filename = "surfaceElevation_5.dat";
addpath('wave_data')

raw_data = table2array(readtable(filename));

idx=find(raw_data(:,1)>0); %Gives indices of positive values of time
data = raw_data(idx,:);
time = data(:,1); % Time vector
surf_elevation = data(:,2:end); % Surface elevation
%plot(time,surf_elevation(:,1))
%% READ DATA
filename = "surfaceElevation_10.dat";
addpath('wave_data')

raw_data = table2array(readtable(filename));

idx=find(raw_data(:,1)>0); %Gives indices of positive values of time
data = raw_data(idx,:);
time = data(:,1); % Time vector
surf_elevation1 = data(:,2:end); % Surface elevation
%plot(time,surf_elevation(:,1))
%% READ DATA
filename = "surfaceElevation_15.dat";
addpath('wave_data')

raw_data = table2array(readtable(filename));

idx=find(raw_data(:,1)>0); %Gives indices of positive values of time
data = raw_data(idx,:);
time = data(:,1); % Time vector
surf_elevation2 = data(:,2:end); % Surface elevation
%plot(time,surf_elevation(:,1))
%% READ DATA
filename = "surfaceElevationcph20.dat";
addpath('wave_data')

raw_data = table2array(readtable(filename));

idx=find(raw_data(:,1)>0); %Gives indices of positive values of time
data = raw_data(idx,:);
time = data(:,1); % Time vector
surf_elevation3 = data(:,2:end); % Surface elevation
%plot(time,surf_elevation(:,1))
%% WAVE PARAMETERS
wave_period = 8; % Wave period [s]
height_ratio = 0.005; % Steepness, wave height to wavelength ratio
depth_ratio = 0.4; %Water depth to length ratio
wave_length = 9.81/(2*pi)*wave_period^2*tanh(2*pi*depth_ratio); % Dispersion relation

wave_height = height_ratio*wave_length; % Wave height [m]
water_depth = depth_ratio*wave_length; % Water depth [m]
wave_freq = 2*pi/wave_period; %Wave frequency [rad/s]

%% CONSTRUCT REFERENCE WAVE
n_periods = time(end)/wave_period; % Number of periods in the data
n_samples = length(time)/n_periods; % Mean number of samples per period
ref_time = linspace(0,wave_period,400); % Reference time for cosine wave ref_time = linspace(0,wave_period,n_samples);
gauges_pos = linspace(0,2*wave_length,size(surf_elevation,2)); %Dimensionless position of the gauges (x/lambda)

% Compute second-order phase average
ref_phase_avg = surfaceElevation_2nd(0,ref_time, wave_length, wave_period, wave_height, water_depth); 

% Compute the Stokes second-order surface elevation at evry position and time
for i=1:length(time)
    ref_wave_matrix(i,:) = surfaceElevation_2nd(gauges_pos, time(i), wave_length, wave_period, wave_height, water_depth);
end
figure(1)
plot(time,surf_elevation(:,1))
hold on
plot(time,ref_wave_matrix(:,1))
legend('numerical mwl', 'analytical mwl')
xlabel('t/T');
ylabel('mean water level [m]')
title('15 Cells Per Wave Height');
%% COMPUTE ERROR ALONG THE DOMAIN
error_amp = zeros(size(surf_elevation,2),1);
error_harm = zeros(size(surf_elevation,2),1);
relative_time=[1:n_periods];
for i=1:size(surf_elevation,2)
     [phase_avg(:,i)] = averagePhase(surf_elevation(:,i),0,n_samples);
     % nRMSE error, normaslised with wave height
     [phase_avg1(:,i)] = averagePhase(surf_elevation1(:,i),0,n_samples);
     % nRMSE error, normaslised with wave height
     [phase_avg2(:,i)] = averagePhase(surf_elevation2(:,i),0,n_samples);
     % nRMSE error, normaslised with wave height
     [phase_avg3(:,i)] = averagePhase(surf_elevation3(:,i),0,n_samples);
     % nRMSE error, normaslised with wave height
     nrmse(i) = (100/wave_height)*sqrt(mean((phase_avg(:,i)' - ref_phase_avg).^2)); 
     % relative error, normaslised with wave height
     rel_error(:,i) = (100/wave_height)*(surf_elevation(:,i) - ref_wave_matrix(:,i));
     % relative error of the phase average, normaslised with wave height
     rel_error1(:,i) = (100/wave_height)*(surf_elevation1(:,i) - ref_wave_matrix(:,i));
     % relative error of the phase average, normaslised with wave height
     rel_error2(:,i) = (100/wave_height)*(surf_elevation2(:,i) - ref_wave_matrix(:,i));
     % relative error of the phase average, normaslised with wave height
     rel_error3(:,i) = (100/wave_height)*(surf_elevation3(:,i) - ref_wave_matrix(:,i));
     % relative error of the phase average, normaslised with wave height
     avg_rel_error(:,i) = (1/wave_height).*(phase_avg(:,i)' - ref_phase_avg);
     avg_rel_error1(:,i) = (1/wave_height).*(phase_avg1(:,i)' - ref_phase_avg);
     avg_rel_error2(:,i) = (1/wave_height).*(phase_avg2(:,i)' - ref_phase_avg);
     avg_rel_error3(:,i) = (1/wave_height).*(phase_avg3(:,i)' - ref_phase_avg);
end

plot(ref_time,avg_rel_error(:,1),'Linewidth',1)
hold on
plot(ref_time,avg_rel_error1(:,1),'Linewidth',1)
hold on
plot(ref_time,avg_rel_error2(:,1),'Linewidth',1)
hold on
plot(ref_time,avg_rel_error3(:,1),'Linewidth',1)
xlabel('t/T')
ylabel('\deltap/H')
legend('5 cpwh','10 cpwh','15 cpwh', '20 cpwh')

figure()
   M10=max(avg_rel_error(:,1));
   [phase_avg(:,1)] = averagePhase(surf_elevation(:,1),1,n_samples);
   avg_rel_error(:,1) = (1/wave_height).*(phase_avg(:,1)' - ref_phase_avg);
   M11=max(avg_rel_error(:,1));
   [phase_avg(:,1)] = averagePhase(surf_elevation(:,1),2,n_samples);
   avg_rel_error(:,1) = (1/wave_height).*(phase_avg(:,1)' - ref_phase_avg);
   M12=max(avg_rel_error(:,1));
   [phase_avg(:,1)] = averagePhase(surf_elevation(:,1),3,n_samples);
   avg_rel_error(:,1) = (1/wave_height).*(phase_avg(:,1)' - ref_phase_avg);
   M13=max(avg_rel_error(:,1));
   [phase_avg(:,1)] = averagePhase(surf_elevation(:,1),4,n_samples);
   avg_rel_error(:,1) = (1/wave_height).*(phase_avg(:,1)' - ref_phase_avg);
   M14=max(avg_rel_error(:,1));
   [phase_avg(:,1)] = averagePhase(surf_elevation(:,1),5,n_samples);
   avg_rel_error(:,1) = (1/wave_height).*(phase_avg(:,1)' - ref_phase_avg);
   M15=max(avg_rel_error(:,1));
   M1=[10^(-5),M15,M14,M13,M12,M11,M10];
   h=semilogy([1:7],M1,'-*','Linewidth',1);
   hold on
   h.LineStyle = '-';
   hold on
  
[phase_avg1(:,1)] = averagePhase(surf_elevation1(:,1),0,n_samples);
avg_rel_error1(:,1) = ((1/wave_height)*(phase_avg1(:,1)' - ref_phase_avg));
%plot(ref_time,avg_rel_error1(:,1))
M20=max(avg_rel_error1(:,1));
   [phase_avg1(:,1)] = averagePhase(surf_elevation1(:,1),1,n_samples);
   avg_rel_error1(:,1) = (1/wave_height).*(phase_avg1(:,1)' - ref_phase_avg);
   M21=max(avg_rel_error1(:,1));
   [phase_avg1(:,1)] = averagePhase(surf_elevation1(:,1),2,n_samples);
   avg_rel_error1(:,1) = (1/wave_height).*(phase_avg1(:,1)' - ref_phase_avg);
   M22=max(avg_rel_error1(:,1));
   [phase_avg1(:,1)] = averagePhase(surf_elevation1(:,1),3,n_samples);
   avg_rel_error1(:,1) = (1/wave_height).*(phase_avg1(:,1)' - ref_phase_avg);
   M23=max(avg_rel_error1(:,1));
   [phase_avg1(:,1)] = averagePhase(surf_elevation1(:,1),4,n_samples);
   avg_rel_error1(:,1) = (1/wave_height).*(phase_avg1(:,1)' - ref_phase_avg);
   M24=max(avg_rel_error1(:,1));
   [phase_avg1(:,1)] = averagePhase(surf_elevation1(:,1),5,n_samples);
   avg_rel_error1(:,1) = (1/wave_height).*(phase_avg1(:,1)' - ref_phase_avg);
   M25=max(avg_rel_error1(:,1));
   M2=[10^(-5),M25,M24,M23,M22,M21,M20]
   h=semilogy([1:7],M2,'-*','Linewidth',1);
   hold on
   h.LineStyle = '-';
hold on

[phase_avg2(:,1)] = averagePhase(surf_elevation2(:,1),0,n_samples);
avg_rel_error2(:,1) = ((1/wave_height)*(phase_avg2(:,1)' - ref_phase_avg));
%plot(ref_time,avg_rel_error2(:,1))
M30=max(avg_rel_error2(:,1));
   [phase_avg2(:,1)] = averagePhase(surf_elevation2(:,1),1,n_samples);
   avg_rel_error2(:,1) = (1/wave_height).*(phase_avg2(:,1)' - ref_phase_avg);
   M31=max(avg_rel_error2(:,1));
   [phase_avg2(:,1)] = averagePhase(surf_elevation2(:,1),2,n_samples);
   avg_rel_error2(:,1) = (1/wave_height).*(phase_avg2(:,1)' - ref_phase_avg);
   M32=max(avg_rel_error2(:,1));
   [phase_avg2(:,1)] = averagePhase(surf_elevation2(:,1),3,n_samples);
   avg_rel_error2(:,1) = (1/wave_height).*(phase_avg2(:,1)' - ref_phase_avg);
   M33=max(avg_rel_error2(:,1));
   [phase_avg2(:,1)] = averagePhase(surf_elevation2(:,1),4,n_samples);
   avg_rel_error2(:,1) = (1/wave_height).*(phase_avg2(:,1)' - ref_phase_avg);
   M34=max(avg_rel_error2(:,1));
   [phase_avg2(:,1)] = averagePhase(surf_elevation2(:,1),5,n_samples);
   avg_rel_error2(:,1) = (1/wave_height).*(phase_avg2(:,1)' - ref_phase_avg);
   M35=max(avg_rel_error2(:,1));
   M3=[10^(-5),M35,M34,M33,M32,M31,M30];
   h=semilogy([1:7],M3,'-*','Linewidth',1);
   h.LineStyle = '-';
   hold on
hold on

[phase_avg3(:,1)] = averagePhase(surf_elevation3(:,1),0,n_samples);
avg_rel_error3(:,1) = ((1/wave_height)*(phase_avg3(:,1)' - ref_phase_avg));
%plot(ref_time,avg_rel_error3(:,1))
M40=max(avg_rel_error3(:,1));
   [phase_avg3(:,1)] = averagePhase(surf_elevation3(:,1),1,n_samples);
   avg_rel_error3(:,1) = (1/wave_height).*(phase_avg3(:,1)' - ref_phase_avg);
   M41=max(avg_rel_error3(:,1));
   [phase_avg3(:,1)] = averagePhase(surf_elevation3(:,1),2,n_samples);
   avg_rel_error3(:,1) = (1/wave_height).*(phase_avg3(:,1)' - ref_phase_avg);
   M42=max(avg_rel_error3(:,1));
   [phase_avg3(:,1)] = averagePhase(surf_elevation3(:,1),3,n_samples);
   avg_rel_error3(:,1) = (1/wave_height).*(phase_avg3(:,1)' - ref_phase_avg);
   M43=max(avg_rel_error3(:,1));
   [phase_avg3(:,1)] = averagePhase(surf_elevation3(:,1),4,n_samples);
   avg_rel_error3(:,1) = (1/wave_height).*(phase_avg3(:,1)' - ref_phase_avg);
   M44=max(avg_rel_error3(:,1));
   [phase_avg3(:,1)] = averagePhase(surf_elevation3(:,1),5,n_samples);
   avg_rel_error3(:,1) = (1/wave_height).*(phase_avg3(:,1)' - ref_phase_avg);
   M45=max(avg_rel_error3(:,1));
   M4=[10^(-5),M45,M44,M43,M42,M41,M40]
   h=semilogy([1:7],M4,'-*','Linewidth',1);
   hold on
   h.LineStyle = '-';
legend('5 cpwh','10 cpwh','15 cpwh', '20 cpwh')
ylim([0 10^0])
xlabel('t/T')
ylabel('log(\deltap/H)')
%% amplitude error
    
real_peak=max(ref_wave_matrix(:,1));
stokes_maxima=[];
position_of_maxima=[];
stokes_maxima_new=[];
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
h=semilogy(relative_time,rel_amplitude_error,'-o','Linewidth',1);
hold on
h.LineStyle = '-';
axis([1 8 10^(-5) 1])
xlabel('t/T');
ylabel('log(\deltaa/H)')
%title('15 Cells Per Wave Height');
%% amplitude error all 4  cpwh together
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
h=semilogy(relative_time,rel_amplitude_error,'-o','Linewidth',1);
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
h=semilogy(relative_time,rel_amplitude_error1,'-o','Linewidth',1);
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
h=semilogy(relative_time,rel_amplitude_error2,'-o','Linewidth',1);


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
h=semilogy(relative_time,rel_amplitude_error3,'-o','Linewidth',1);


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

