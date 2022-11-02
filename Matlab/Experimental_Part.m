clear all
close all
clc

%% READING OF THE DATA ACQUIRED IN THE LABORATORY

load 'Acquired data.mat'
Acc = [Acc_1, -Acc_2]; % Array with the accelerations of both accelerometers; -Acc_2 because the acquisition direction of the 2 accelerometers is opposite (bisogna pareggiare le direzioni di acquisizione)
mean_acc = mean(Acc, 2); % Mean values for the accelerations of the two accelerometers taken at the same instant of time
time = linspace(0, 4, 204800);
time = time';

%% EXPERIMENTAL TRANSFER FUNCTION ONLY FOR THE MEANINGFUL TIME WINDOW --> TIME IN [1.2, 1.5]S
samp_freq = 51200; %[Hz]
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

% PROVA IN DECIBEL
%Tm = mag2db(Tm);
%plot(Fr, Tm);
%xlim([0 2000])

%% HALF POWER POINTS METHOD
peaks_Tm = [10.8903, max(Tm), 4.2757]; % in amplitude
%peaks_Tm = [20.7408, 32.8038, 12.6201]; % in dB

% MODE 1
g1 = 0.707*peaks_Tm(1); % constant function of rms value of 1st peak
for r=1:length(Fr)
g1_plot(r) = g1; 
end
g1_plot = g1_plot';
figure 
loglog(Fr, Tm);
hold on
plot(Fr, g1_plot);
hold off
xlabel('Frequency (Hz)')
axis([0 2000 0 2000])
legend('|G(iw)| (Abs)','0.707*Q1')

% MODE 2
g2 = 0.707*peaks_Tm(2); % constant function of rms value of 1st peak
for r=1:length(Fr)
g2_plot(r) = g2; 
end
g2_plot = g2_plot';
figure 
loglog(Fr, Tm);
hold on
plot(Fr, g2_plot);
hold off
xlabel('Frequency (Hz)')
axis([0 2000 0 2000])
legend('|G(iw)| (Abs)','0.707*Q3')

% MODE 3
g3 = 0.707*peaks_Tm(3); % constant function of rms value of 1st peak
for r=1:length(Fr)
g3_plot(r) = g3; 
end
g3_plot = g3_plot';
figure 
loglog(Fr, Tm);
hold on
plot(Fr, g3_plot);
hold off
xlabel('Frequency (Hz)')
axis([0 2000 0 2000])
legend('|G(iw)| (Abs)','0.707*Q2')

%% INTERSECTIONS MODE 1

L11 = [Fr'; Tm']; % first curve
L21 = [Fr'; g1_plot']; % second curve
P = InterX(L11,L21);
w11 = P(1,1); w21 = P(1,2);
damping_ratio_1 = (w21-w11)/(w11+w21);

%% INTERSECTIONS MODE 2

L12 = [Fr'; Tm']; % first curve
L22 = [Fr'; g2_plot']; % second curve
P = InterX(L12,L22);
w12 = P(1,1); w22 = P(1,2);
damping_ratio_2 = (w22-w12)/(w12+w22);

%% INTERSECTIONS MODE 3

L13 = [Fr'; Tm']; % first curve
L23 = [Fr'; g3_plot']; % second curve
P = InterX(L13,L23);
w13 = P(1,3); w23 = P(1,4);
damping_ratio_3 = (w23-w13)/(w13+w23);