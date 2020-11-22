%% 2. DC Motor Back EMF Constant
clear all; close all
% == Part A == 
Ktach = 3.0 / 1000 / (2*pi/60); % V/rad/s

% == Part B ==
cali_mut_em = [0.997 1.930 2.870 3.890 4.850 5.920 6.880 7.930 8.900]; % V
cali_tach_eo = [0.560 1.100 1.770 2.420 2.960 3.570 4.200 4.830 5.540]; % V
omega = cali_tach_eo / Ktach; % rad/s

figure(1)
plot(omega, cali_mut_em); grid on
ylabel('e_{m} (V)'); xlabel('\omega (rad/s)')

% Motor voltage constant
% Spec Sheet: Ke = 4.39-5.37 V/krpm -> 0.0466 V/rad/s
Ke = (cali_mut_em(end) - cali_mut_em(1)) / (omega(end) - omega(1)); % V/rad/s

% == Part C ==

% Motor torque constant
% Spec Sheet: Kt = 6.6 +- 10% ozf-in/A; -> see lecture for this conversion
Kt = 141.6*Ke; % ozf-in/A

% == Part D ==
start_volt = 1.4; % V
int_resistance = 4.25; % ohms
start_cur = start_volt / int_resistance; % A
start_torque = Kt * start_cur; % ozf-in

%% 3. Open Loop Response to a Voltage Step Input and Disturbance Torque

% == Part A ==


% == Part B ==
filestruct1 = TDMS_readTDMSFile('3.9Waveform.tdms');
mut_eo = filestruct1.data{5};
t = linspace(-0.600, 0.600, length(mut_eo));
figure(2)
plot(t, mut_eo)


