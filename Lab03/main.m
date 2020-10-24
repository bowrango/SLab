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
figure(4)
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

figure(5)
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

figure(6)
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


