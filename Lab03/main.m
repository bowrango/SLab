%% 1. Linear Velocity Transducer 

clear all; close all

e0_lvt = readmatrix('Waveform1.2.csv')';
t1 = linspace( -1.11433336704258, 0.00002403*100000, 100000);

figure(1)
plot(t1, e0_lvt); grid on
xlabel('Time (s)'); ylabel('LVT Voltage (V)')

g = 386; % [in/sec^2]
mass = 0.075/0.45359; % [lbm]

% estimate the slope of the initial response 
fit = polyfit(t1(26371:30632), e0_lvt(26371:30632), 1);
slope = fit(1);

lvt_sensitivity = slope / g; % [V / (in./s)]
core_velocity = e0_lvt  / lvt_sensitivity; % [in./s]
core_position = cumtrapz(t1, core_velocity); % [in.]

figure(2); grid on
xlabel('Time (s)')
yyaxis right
plot(t1, core_velocity); hold on
ylabel('Core Velocity (in./s)')
yyaxis left
plot(t1, core_position); hold on
ylabel('Core Position (in.)')

% can only use the last couple peaks because the core has to be in contact
% with the foam pad for it to be considered second order -> see piazza 
[peak_idx1, peak_mag1] = peakfinder(e0_lvt(42722:end), 0.08, 0.025);
zeta1 = mean(find_damping_ratios(peak_mag1, 0));
[omega_undamped1, omega_damped1] = find_undamped_natural_frequency(t1, peak_idx1+42722, zeta1);

% plot to confirm the peaks look good 
figure(3)
plot(t1, e0_lvt); hold on
plot(t1(peak_idx1+42722), peak_mag1, 'o')

% foam spring constant found from standard form of second order tf
k_foam = omega_undamped1^2 * mass * (1 / 32.174) * (1 / 12); % [lbf/in.]

% damping coeffienct
B = 2 * zeta1 * omega_undamped1 * mass * (1 / 32.174) * (1 / 12); % [lbf-s/in.]

% impact force
v_impact = 53.7; % [in./s] -
F_impact = B*v_impact; % [lbf]

%% 2.1 LDVT Characteristics

e2_1_lvt_above = readmatrix('Waveform2.1.3_Above_Null.csv')';
e2_1_lvt_null = readmatrix('Waveform2.1.3_Null.csv')';
e2_1_lvt_below = readmatrix('Waveform2.1.3_Below_Null.csv')';
t2_1 = linspace(0.33687612, 0.00000015*6668, 6668);

figure(4)
subplot(3,1,1)
title('Voltage Output vs Coil Displacement before Amplitude Demodulation')
hold on
plot(t2_1, e2_1_lvt_above(1,:),'r')
hold on
plot(t2_1, e2_1_lvt_above(2,:),'-or','MarkerIndices',1:300:length(t2_1))
hold on
ylabel('Above Null [volts]')
legend('Coil 1','Coil 2')
grid on
subplot(3,1,2)
plot(t2_1, e2_1_lvt_null(1,:),'b')
hold on
plot(t2_1, e2_1_lvt_null(2,:),'-ob','MarkerIndices',1:300:length(t2_1))
hold on
ylabel('At Null [volts]')
legend('Coil 1','Coil 2')
grid on
subplot(3,1,3)
plot(t2_1, e2_1_lvt_below(1,:),'k')
hold on
plot(t2_1, e2_1_lvt_below(2,:),'-ok','MarkerIndices',1:300:length(t2_1))
hold on
ylabel('Below Null [volts]')
legend('Coil 1','Coil 2')
xlabel('Time [s]')
grid on

e2_2_lvt_above = readmatrix('Waveform2.1.5_Above_Null.csv')';
e2_2_lvt_null = readmatrix('Waveform2.1.5_Null.csv')';
e2_2_lvt_null = e2_2_lvt_null(:,(1:6453));
e2_2_lvt_below = readmatrix('Waveform2.1.5_Below_Null.csv')';
t2_2 = linspace(0.0283703208033347, 0.00000031*6453, 6453);

figure(5)
subplot(3,1,1)
title('Voltage Output vs Coil Displacement after Amplitude Demodulation')
hold on
plot(t2_2, e2_2_lvt_above(1,:),'r')
hold on
plot(t2_2, e2_2_lvt_above(2,:),'-or','MarkerIndices',1:300:length(t2_2))
hold on
ylabel('Above Null [volts]')
legend('Coil 1','Coil 2')
grid on
subplot(3,1,2)
plot(t2_2, e2_2_lvt_null(1,:),'b')
hold on
plot(t2_2, e2_2_lvt_null(2,:),'-ob','MarkerIndices',1:300:length(t2_2))
hold on
ylabel('At Null [volts]')
legend('Coil 1','Coil 2')
grid on
subplot(3,1,3)
plot(t2_2, e2_2_lvt_below(1,:),'k')
hold on
plot(t2_2, e2_2_lvt_below(2,:),'-ok','MarkerIndices',1:300:length(t2_2))
hold on
ylabel('Below Null [volts]')
legend('Coil 1','Coil 2')
xlabel('Time [s]')
grid on

%% 2.2 LVDT System Calibration and Dynamics

%A
weight = ([100 200 300 400 500] + 7.5) / 1000; %Kg
v_weight = [-130.714 -263.856 -391.563 -523.374 -657.848] / 1000; %V

disp = [0.05 0.1 0.15 0.2 0.25]; %in
v_disp = [-130.313 -278.8 -426.208 -574.973 -719.663] / 1000; %V

figure(6)
plot(weight, v_weight,'-o')
grid on
title('Weight vs Output Voltage')
ylabel('Output Voltage [volts]')
xlabel('Weight [Kg]')

figure(7)
plot(disp, v_disp,'-o')
grid on
title('Displacement vs Output Voltage')
ylabel('Output Voltage [volts]')
xlabel('Displacement [in]')

<<<<<<< HEAD
%B
slope_w = (v_weight(1) - v_weight(end)) / (weight(1) - weight(end)); %v/kg
slope_d = (v_disp(1) - v_disp(end)) / (disp(1) - disp(end)); %v/in

slope_w = slope_w * (1/2.2046226218488); %v/lbf

k = slope_d/slope_w; %kg/in

%C
LVDTsens = slope_d; %V/in

%D 
e2_4_lvt = readmatrix('Waveform2.2.4.csv')';
t2_4 = linspace(-3.369857989, 0.00008458*53206, 53206);

[peak_idx2, peak_mag2] = peakfinder(e2_4_lvt, 0.08, 0.0025);
zeta2 = mean(find_damping_ratios(peak_mag2, -.0075));
[omega_undamped2, omega_damped2] = find_undamped_natural_frequency(t2_4, peak_idx2, zeta2);

figure(8)
plot(t2_4, e2_4_lvt)
hold on
plot(t2_4(peak_idx2), peak_mag2,'o')
hold on

eff_m = k/(omega_undamped2^2); %should be lbm, but idk the units of w_undamped

=======
>>>>>>> 1da69901f6e7e4c667c94f635da476ecac276fee
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


