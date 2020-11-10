clear all; close all;

%% 1. Potentiometer Accelerometer 

% == Part a. ==

weights1 = ([10.9 60.9 110.9 160.9 200.9 215.9]/1000)*9.8 * 3.60; % ozf
e0_weights1 = [619.2 705.7 788 873.4 940.6 965.4]/1000; % V

figure(1)
plot(weights1, e0_weights1); grid on
xlabel('Weight (oz_{f})'); ylabel('e_{0} (V)')

sensitivityA = (e0_weights1(end) - e0_weights1(1)) / (weights1(end) - weights1(1)); % V/ozf

% == Part b. ==

% the v2p5 folder contains functions to read .tdms files. Make sure to add
% this folder and subfolders to your path when you clone the repo 
filestruct1 = TDMS_readTDMSFile('JPL11.tdms');
e01 = filestruct1.data{5};

% time vector taken from the png plot
t1 = linspace(-0.6, 0.6, length(e01));

figure(2)
plot(t1, e01); grid on

% effective mass 
eff_mass = ((0.6022 - 0.5250) / sensitivityA); % ozm

% damping ratio
overshoot = abs(0.519 - 0.594) / 0.594;
zeta1 = 0.56;

% natural frequency 
period1 = t1(54034) * 2;
omega_damped1 = 2*pi/period1;
omega_n1 = omega_damped1 / sqrt(1 - zeta1^2); % rad/s

% spring constant
eff_mass_slug = eff_mass * (1/32.174) * (1/16); % lbf-s^2/ft
k = eff_mass_slug * (16) * (1/12) * omega_n1^2; % ozf/in.

% == Part c. ==

sensitivityB = sensitivityA * (1/32.174) * (1/12) * eff_mass; % V/in/s^2

% == Part d. ==

max_accel = 1.01 / sensitivityB; % in/s^2

% == Part e. ==
num = 1;
den = [1/omega_n1^2 2*zeta1/omega_n1 1];

tf = tf(num, den);
bode(tf)









