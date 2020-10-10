clear all; close all

%% 1 Position and Velocity Measurements 

tach_sensitivity = 0.007 / 0.10472; % V/rpm -> V/(rad./s)
pot_sensitivity = 0.015 / 0.01745; % V/deg. -> V/rad.

Dataset13 = readmatrix('WaveformData_1.3_redo.csv');
t1 = linspace(0, 0.00001408* 100000, 100000);
e_tach1 = Dataset13(6:end, 1);
e_pot1 = Dataset13(6:end, 2);

omega1 = e_tach1 / tach_sensitivity; % rad./s
theta1 = e_pot1 / pot_sensitivity; % rad.

figure(1)
yyaxis left
plot(t1(1:end-5), omega1)
ylabel('\omega (rad./s)')
yyaxis right
plot(t1(1:end-5), theta1)
ylabel('\theta (rad.)')
xlabel('Time (s)')
xlim([0.3 0.75]);

% Theta begins with some initial displacement, and omega is initially 0. As
% the gear oscillates both signal have approximately the same period and damping 
% effects, but omega is lagged behind. 

figure(2)
plot(theta1, omega1)
xlabel('\theta (rad.)'); ylabel('\omega (rad./s)')

theta_int1 = cumtrapz(t1(1:end-5), omega1);

% The signal starts at the outter edge of the spiral and goes inwards as
% time moves forward. This represents the damping response as the magnitude
% of the oscillations decreases gradually. Eventually, the signal dies out
% when omega is equal to 0 at the center of the spiral. 

figure(3)
yyaxis left
plot(t1(1:end-5), theta1) 
ylabel('\theta (rad.)')
yyaxis right
plot(t1(1:end-5), theta_int1)
ylabel('\theta_{int} (rad.)')
xlabel('Time (s)')
xlim([0 0.75]);

% Cumtrapz Effects:
%  - assumes initial condition is 0, which essentially applies a shift
%  - removes noise and smoothes signal 
%  - relative amplitudes of peaks don't exactly match 

%% 2 Second-Order Parameter Identification






