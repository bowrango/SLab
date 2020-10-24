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


%% 2.3 


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


