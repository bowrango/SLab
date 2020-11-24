%% 2. DC Motor Back EMF Constant
clear all;
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
% Spec Sheet: Ke = 4.39-5.37 V/krpm = 0.0466 V/rad/s
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

% Find the transfer function Eo/Ei of Figure 4, ignore the power op-amp,
% Td = 0

% == Part B ==

stall_cur = 6.0/5.0; % A
stall_torque = Kt * stall_cur; % ozf-in

% == Part C ==

% Gain Term K -> See questions 286/298 on piazza
K = (4.2 / 8.0) / Ktach; % (rad/s)/V

% Time Constants for Step and Disturbance 
filestruct37 = TDMS_readTDMSFile('3.7Waveform.tdms'); % Step Response
tach_step_eo = filestruct37.data{5};
t_step = linspace(-0.600, 0.600, length(tach_step_eo));
plot(t_step, tach_step_eo)

t_step_shift = t_step + 0.056;
ef_step = (4.2 - tach_step_eo) ./ (4.2 - -4.1);
tau_step_idx = find_closest_to(0.368, ef_step);
tau_step = t_step_shift(tau_step_idx);

filestruct39 = TDMS_readTDMSFile('3.9Waveform.tdms'); % Disturbance Load
tach_dtb_eo = filestruct39.data{5};
t_dtb = linspace(-0.600, 0.600, length(tach_dtb_eo));

t_dtb_shift = t_dtb + 0.240; 
ef_dtb = (-2.1 - tach_dtb_eo) ./ (-2.1 - -3);
tau_dtb_idx = find_closest_to(0.368, ef_dtb);
tau_dtb = t_dtb_shift(tau_dtb_idx);

% == Part D ==
R = 4.2; % ohms

% See question 276. 
B_rad = ( ((Kt/K) - (Ke*Kt)) / R); % ozf-in/(rad/s)
B_krpm = B_rad * (2*pi/60) * (1000); % ozf-in/krpm

J = (tau_step/R) * ( (Ke*Kt) + (R*B_rad)); % ozf-in-s^2

% == Part E ==
settle_time = 4*tau_dtb;



%% Functions 

function [idx] = find_closest_to(value, array)

    error = abs(array - value);
    best = min(error);
    idx = find(error == best);
    
    idx = idx(1);
end

