clear all
close all
clc

%% SYSTEM PARAMETERS

E = 206e9; %[Pa]
rho = 7850; %[kg/m^3]
A = 111e-6; %[m^2]
J = 6370e-12; %[m^4]
L = 0.7; %[m]
Xi_1 = 0.05;
Xi_2 = 0.01;
Xi_3 = Xi_2;

%% READING OF THE DATA ACQUIRED IN THE LABORATORY

load 'Acquired data.mat'
Acc = [Acc_1, -Acc_2]; % Array with the accelerations of both accelerometers; -Acc_2 because the acquisition direction of the 2 accelerometers is opposite (bisogna pareggiare le direzioni di acquisizione)
mean_acc = mean(Acc, 2); % Mean values for the accelerations of the two accelerometers taken at the same instant of time
time = linspace(0, 4, 204800);
time = time';

%% FORCE AND ACELERATIONS PLOTS

% Force plot
figure
plot(time, Force);
xlabel('Time (s)')
ylabel('Force of the hammer (N)')
title('Force')
% Accelerations plots for all acquisition time
figure
subplot(3, 1, 1);
plot(time, Acc_1);
xlabel('Time (s)')
ylabel('Acceleration 1 (m/s^2)')
title('Accelerations of the 1st accelerometer')
xlim([0 4])
subplot(3, 1, 2);
plot(time, Acc_2);
xlabel('Time (s)')
ylabel('Acceleration 2 (m/s^2)')
title('Accelerations of the 2nd accelerometer')
xlim([0 4])
subplot(3, 1, 3);
plot(time, mean_acc);
xlabel('Time (s)')
ylabel('Mean acceleration (m/s^2)')
title('Mean acceleration')
xlim([0 4])
% Accelerations plots only for meaningful time window
figure
subplot(3, 1, 1);
plot(time, Acc_1);
xlabel('Time (s)')
ylabel('Acceleration 1 (m/s^2)')
title('Accelerations of the 1st accelerometer')
xlim([1.23 1.45])
subplot(3, 1, 2);
plot(time, Acc_2);
xlabel('Time (s)')
ylabel('Acceleration 2 (m/s^2)')
title('Accelerations of the 2nd accelerometer')
xlim([1.23 1.45])
subplot(3, 1, 3);
plot(time, mean_acc);
xlabel('Time (s)')
ylabel('Mean acceleration (m/s^2)')
title('Mean acceleration')
xlim([1.23 1.45])

%% EXPERIMENTAL TRANSFER FUNCTION BETWEEN THE ACQUIRED IMPULSE FORCE AND MEAN ACCELERATION FOR ALL TIME

samp_freq = 51200; %[Hz]
[Tf, Fr] = tfestimate(Force, mean_acc, [ ], [ ], [ ], samp_freq); % Tf: transfer function, Fr: frequencies vector
Tm = abs(Tf); % Magnitude of the transfer function
% Abs(Tf) in decibel
%figure 
%plot(Fr, mag2db(Tm));
%xlabel('Frequency (Hz)')
%ylabel('|G(iw)| (dB)')
%xlim([0 2000])

% Abs(Tf) not in decibel, but amplitude; log-log scale
figure 
loglog(Fr, Tm);
xlabel('Frequency (Hz)')
ylabel('|G(iw)| (Abs)')
axis([0 2000 0 2000])

%% EXPERIMENTAL TRANSFER FUNCTION ONLY FOR THE MEANINGFUL TIME WINDOW --> TIME IN [1.2, 1.5]S
Size = 81923-61439;
i=1;
while i < Size
     for j = 62975:74243
      acc_mean(i) = mean_acc(j);
      force(i) = Force(j);
      i=i+1;
     end
end
[Tf, Fr] = tfestimate(force, acc_mean, [ ], [ ], [ ], samp_freq); % Tf: transfer function, Fr: frequencies vector
Tm = abs(Tf); % Magnitude of the transfer function

% Abs(Tf) not in decibel, but amplitude; log-log scale
figure 
loglog(Fr, Tm);
xlabel('Frequency (Hz)')
ylabel('|G(iw)| (Abs)')
axis([0 2000 0 2000])


% Abs(Tf) in decibel
%figure 
%plot(Fr, mag2db(Tm));
%xlabel('Frequency (Hz)')
%ylabel('|G(iw)| (dB)')
%xlim([0 2000])

%% HALF POWER POINTS METHOD
w1 = 124.466;
w2 = 497.864;
w3 = 1120.19;
peaks_Tm = [17.0441, max(Tm)];
peaks_freq = [128.125, 1079.69]; % in [Hz]

% MODE 1
g1 = 0.707*peaks_Tm(1); % constant function of rms value of 1st peak
for i=1:length(Fr)
g1_plot(i) = g1; 
end
g1_plot = g1_plot';
figure 
loglog(Fr, Tm);
hold on
plot(Fr, g1_plot);
hold off
xlabel('Frequency (Hz)')
ylabel('|G(iw)| (Abs)')
axis([0 2000 0 2000])
% prova per trovare punti di intersezione
dy1y2 = g1_plot-Tm; 
idx = find(diff(sign(dy1y2)));  
for k = 1:numel(idx)
    idxrng = max(idx(k)-2,1) : min(numel(Fr),idx(k)+2);        
    xv(k) = interp1(dy1y2(idxrng), Fr(idxrng), 0);             
    yv(k) = interp1(Fr(idxrng), g1_plot(idxrng), xv(k));             
end
figure
loglog(Fr, Tm);
hold on
loglog(Fr, g1_plot);
loglog(xv, yv, 'xr')
hold off
grid
axis([0 2000 0 2000])
text(xv(end-2:end), yv(end-2:end), compose('\\leftarrow (%7.4f, %7.4f)', [xv(end-2:end); yv(end-2:end)].'), 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'Rotation',-15)
xlabel('Frequency (Hz)')
ylabel('|G(iw)| (Abs)')
title('Intersections With Coordinate Values')
int1 = 124.931; int2 = 130.633; % frequencies at intersection
damping_ratio_1 = (int2 - int1)/(2*peaks_freq(1))

% MODE 2
g2 = 0.707*peaks_Tm(2); % constant function of rms value of 1st peak
for i=1:length(Fr)
g2_plot(i) = g2; 
end
g2_plot = g2_plot';
figure 
loglog(Fr, Tm);
hold on
plot(Fr, g2_plot);
hold off
xlabel('Frequency (Hz)')
ylabel('|G(iw)| (Abs)')
axis([0 2000 0 2000])
% prova per trovare punti di intersezione
dy1y2 = g2_plot-Tm; 
idx = find(diff(sign(dy1y2)));  
for k = 1:numel(idx)
    idxrng = max(idx(k)-2,1) : min(numel(Fr),idx(k)+2);        
    xv(k) = interp1(dy1y2(idxrng), Fr(idxrng), 0);             
    yv(k) = interp1(Fr(idxrng), g2_plot(idxrng), xv(k));             
end
figure
loglog(Fr, Tm);
hold on
loglog(Fr, g2_plot);
loglog(xv, yv, 'xr')
hold off
grid
axis([0 2000 0 2000])
text(xv(end-2:end), yv(end-2:end), compose('\\leftarrow (%7.4f, %7.4f)', [xv(end-2:end); yv(end-2:end)].'), 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'Rotation',-15)
xlabel('Frequency (Hz)')
ylabel('|G(iw)| (Abs)')
title('Intersections With Coordinate Values')
int1_2 = 1048.9; int2_2 = 1090.0466; % frequencies at intersection
damping_ratio_2 = (int2_2 - int1_2)/(2*peaks_freq(2))