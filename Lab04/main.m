clear all; close all;

%% 1. Potentiometer Accelerometer 

% == Part a. ==

weights1 = ([10.9 60.9 110.9 160.9 200.9 215.9]/1000)*9.8 * 3.60; % N -> ozf
e0_weights1 = [619.2 705.7 788 873.4 940.6 965.4]/1000; % mV -> V

figure(1)
plot(weights1, e0_weights1); grid on
xlabel('Weight (oz_{f})'); ylabel('e_{0} (V)')

% slope of the line (V/ozf)
sensitivity = (e0_weights1(end) - e0_weights1(1)) / (weights1(end) - weights1(1));

% == Part b. ==

% the v2p5 folder contains functions to read .tdms files. Make sure to add
% this folder and subfolders to your path when you clone the repo 
filestruct1 = TDMS_readTDMSFile('JPL11.tdms');
e01 = filestruct1.data{5};

% time vector taken from the png plot
t1 = linspace(-0.6, 0.6, length(e01));

figure(2)
plot(t1, e01); grid on

% points taken from plot
overshoot = abs(0.519 - 0.594) / 0.594;
zeta1 = 0.56;

% overshoot is about half of a full oscillation
period1 = t1(54034) * 2;
omega_damped1 = 2*pi/period1;

% natural frequency 
omega_undamped1 = omega_damped1 / sqrt(1 - zeta1^2);

% this is kinda sketchy, not sure if it is correct.


%% 2. Piezoelectric Force Sensor

%2.1 PCB Characteristics

% == Part a. ==

mass = [.2 .4 .6 .8 1 1.2 1.4 1.6] .* 2.2046226218488; %lbf
output = [21.4 43.3 66.1 82.1 106 125 143 163]; %mV

figure(3)
plot(mass, output)
hold on

sens1 = polyfit(mass, output, 1);
output_ideal = mass * sens1(1) + sens1(2);
plot(mass, output_ideal,'-.');
legend('Calibration Curve','Best Fit Line')
title('Calbiration Curve')
ylabel('PCB Force Transducer Output [mV]')
xlabel('Weight [lbf]')

FSerror = (output(3) - output_ideal(3)) / output_ideal(3); %Error in percent Full Scale

% == Part b. ==

PCB = readmatrix('JPL2.13.csv')';

PCB_t = PCB(2,5:end); %s
PCB_v = PCB(1,5:end); %mv

% figure(4)
% plot(PCB_t, PCB_v)

t = 14.25; %seconds - Taken from the plot above

%Cant measure steady state responses because of piezoelectric crystal

% == Part c. ==

PCB_sys = readmatrix('JPL2.14.csv')';
PCB_sys_t = PCB_sys(2,5:end); %s
PCB_sys_v = PCB_sys(1,5:end); %mv

[peak_idx, peak_mag] = peakfinder(PCB_sys_v, 0.01);
peak_idx = peak_idx(16:28);
peak_mag = peak_mag(16:28);
zeta = mean(find_damping_ratios(peak_mag, 0));
[omega_undamped, omega_damped] = find_undamped_natural_frequency(PCB_sys_t, peak_idx, zeta);

% figure(4)
% plot(PCB_sys_t, PCB_sys_v); hold on
% plot(PCB_sys_t(peak_idx), peak_mag, 'o'); 

%2.2 Impulse Loading and Vibration

% == Part a. ==

foam = readmatrix('JPL2.2.csv')';
foam_t = foam(2,49000:end) - foam(2,49000); %s
foam_v = foam(1,49000:end); %mv

[peak_idx1, peak_mag1] = peakfinder(foam_v, 0.01);
%peak_idx = peak_idx(16:28);
%peak_mag = peak_mag(16:28);
zeta1 = mean(find_damping_ratios(peak_mag1, 0));
[omega_undamped2, omega_damped2] = find_undamped_natural_frequency(foam_t, peak_idx1, zeta1);

mass1 = 2 * 32.174; %lbm

figure(4)
plot(foam_t, foam_v); hold on
%plot(out.tout, out.simout)
%legend('Experimental', 'Theoretical')
title('System Response of Impulse Force')


k = (omega_undamped2/(2*pi))^2 * mass1 * (1 / 32.174) * (1 / 12); % [lbf/in.]

% == Part b. ==

B = 2 * zeta1 * omega_undamped2 * mass1 * (1 / 32.174) * (1 / 12); % [lbf-s/in.]

% == Part c. ==

first = 1 / (omega_undamped2^2);
second = 2 * zeta1 / omega_undamped2;
third = 1;
top = 1/k;


%% Functions 
function zeta = find_damping_ratios(peaks, final_value)
    % Calculates the damping ratio between each peak using the log decrement method 
    
    amplitude = peaks - final_value;
    delta = zeros(1, length(peaks)-1);
    zeta = zeros(1, length(peaks)-1);
    
    for n = 1:length(peaks)-1
        delta(n) = (1/n)*log(amplitude(n) / amplitude(n+1));
        zeta(n) = delta(n) / sqrt( (4*(pi^2)) + delta(n)^2 );
    end
end

function [omega_undamped, omega_damped] = find_undamped_natural_frequency(time, peak_indices, zeta)

    wave_length = time(peak_indices(end)) - time(peak_indices(1));
    period = wave_length / (length(peak_indices) - 1);
    omega_damped = 2*pi/period;
    
    omega_undamped = omega_damped / sqrt(1 - zeta^2);
end