%% 1.1 Step Response of 1st-Order System
clear all; close all

% === +-2 V Square Wave ===
Dataset115 = readmatrix('WaveformData_1.1.5.csv');
t1 = linspace(0, 26002*0.000001, 26002);
V1 = Dataset115(:, 2)';
V_final1 = 2; V_initial1 = -2;

% Time right before system response
response1_idx = find_initial_response(V1, 5, 0.02);
t1o = t1(response1_idx);

% Fit a line where the response behaves linearly 
fit1 = polyfit(t1(response1_idx:response1_idx+800), V1(response1_idx:response1_idx+800), 1);
line1 = @(x) fit1(1)*x + fit1(2);

% Extend the fit line slightly past the final value
slope_time1 = t1(response1_idx:response1_idx+2000);
initial_slope_line1 = line1(slope_time1);

% Find where the initial slope line reaches the final response value (intersection)
tau_idx1_ISM = find_closest_to(V_final1, initial_slope_line1);
tau1_ISM = slope_time1(tau_idx1_ISM) - t1o;

% Find where the error fraction is equal to 0.368
error_fraction1 = (V_final1 - V1) ./ (V_final1 - V_initial1);
tau_idx1_632 = find_closest_to(0.368, error_fraction1);
tau1_632 = t1(tau_idx1_632) - t1o;

R1 = 10.94*10^3;
C1 = tau1_632 / R1;

% Transfer function for ideal 1st order system response
theo_transfer_function1 = tf(1, [R1*C1  1]);

% Compare ideal system response to experimential result
figure(1)
opts1 = stepDataOptions('InputOffset', -2, 'StepAmplitude', 4);
stepplot(theo_transfer_function1, 0.005, opts1, '--'); hold on; ylim([-2 2.1]); grid on
plot(t1(response1_idx:response1_idx+5000)-t1o, V1(response1_idx:response1_idx+5000), 'color', 'r')
legend('Theoretical', 'Actual')

% === -3/+5 V Square Wave ===
Dataset117 = readmatrix('WaveformData_1.1.7.csv');
t2 = linspace(0,  60002* 0.000001,  60002);
V2 = Dataset117(:, 2)';
V_final2 = 5; V_initial2 = -3;

% Time right before system response
response2_idx = find_initial_response(V2(3500:end), 5, 0.02) + 3500;
t2o = t2(response2_idx);

% Fit a line where the response behaves linearly 
fit2 = polyfit(t2(response2_idx:response2_idx+800), V2(response2_idx:response2_idx+800), 1);
line2 = @(x) fit2(1)*x + fit2(2);

% Extend the fit line slightly past the final value
slope_time2 = t2(response2_idx:response2_idx+1200);
initial_slope_line2 = line2(slope_time2);

% Find where the slope line reaches the final response value (intersection)
tau_idx2_ISM = find_closest_to(V_final2, initial_slope_line2);
tau2_ISM = slope_time2(tau_idx2_ISM) - t2o;

% Find where the error fraction is equal to 0.368
error_fraction2 = (V_final2 - V2(3400:end)) ./ (V_final2 - V_initial2);
tau_idx2_632 = find_closest_to(0.368, error_fraction2);
shifted_t2 = t2(3400:end);
tau2_632 = shifted_t2(tau_idx2_632) - t2o;

%% 1.2 Frequency Response of 1st-order System

input_frequencies = [11.25 56.25 112.5 225 337.5 450 900].*6.28; % Hz -> Rad./s
magnitude_ratio = [0.088632296 -0.134667653 -0.939931254 -2.878131514 -4.972838524 -6.52493568 -11.64928871]; % dB
phase_shift = [-2.93625 -13.06125 -25.0695 -44.1 -54.5535 -61.56 -73.872]; % deg.

handle = figure(2);
bode(theo_transfer_function1, '--'); hold on; grid on
children = get(handle, 'Children');
axes(children(3)) % Magnitude Plot
semilogx(input_frequencies, magnitude_ratio, 'color', 'r'); ylim([-40 5]);
axes(children(2))  % Phase Plot
semilogx(input_frequencies, phase_shift, 'color', 'r')
text(150, -45, '{\omega}_{break} = 255 Hz')
legend('Theoretical', 'Actual'); hold off

%% 2.1 Step Response of 2nd-Order System

Dataset214 = readmatrix('WaveformData_2.1.4.csv');
t3 = linspace(0,  28002*0.000001, 28002);
V3 = Dataset214(:, 2)';
V_final3 = 2; V_initial3 = 0;

% Time right before system response
response3_idx = find_initial_response(V3, 30, 0.028);
t3o = t3(response3_idx);

[peak_idx3, peak_mag3] = peakfinder(V3, 0.25);
% Zeta is taken to be the mean over the 4 damping ratios for the 5 peaks 
zeta3 = mean(find_damping_ratios(peak_mag3, V_final3));
omega_ud3 = find_undamped_natural_frequency(t3, peak_idx3, zeta3);

Dataset215 = readmatrix('WaveformData_2.1.5.csv');
t4 = linspace(0,  28002*0.000001, 28002);
V4 = Dataset215(:, 2)';
V_final4 = 5; V_initial4 = -3;

[peak_idx4, peak_mag4] = peakfinder(V4, 0.2);
% Zeta is taken to be the mean over the 4 damping ratios for the 5 peaks 
zeta4 = mean(find_damping_ratios(peak_mag4, V_final4));
omega_ud4 = find_undamped_natural_frequency(t4, peak_idx4, zeta4);

% Capacitance and inductance solved for using the 2nd order d.e. form
R2 = 20.1*10^3;
C2 = sqrt(1 / ((omega_ud3^2)*4*(R2^2)*(zeta3^2)) );
L2 = C2*4*(R2^2)*(zeta3^2);

% Transfer function for ideal 2nd order system response
theo_transfer_function2 = tf(1/(L2*C2), [1 1/(C2*R2) 1/(L2*C2)]);

% Compare ideal system response to experimential result
figure(3)
opts2 = stepDataOptions('InputOffset', 0, 'StepAmplitude', 2);
stepplot(theo_transfer_function2, 0.02, opts2, '--'); hold on; ylim([0 4.5]); grid on
plot(t3-t3o-0.00045, V3, 'color', 'r');
legend('Theoretical', 'Actual')

%% 2.2 Frequency Response of 2nd-order System

input_frequencies_2order = [15 75 150 300 450 600 1200].*6.28; % Hz -> Rad./s
magnitude_ratio_2order = [0.652028655 1.012666838 2.858952064 9.176816756 -2.687605857 -8.74008063 -20.57927392]; % dB
phase_shift_2order = [-0.243 -1.755 -8.1 -86.724 -158.76 -164.16 -172.8]; % deg.

%2.2 Bode Plot theoretical vs experimental
handle = figure(4);
bode(theo_transfer_function2, '--'); hold on; grid on
children = get(handle, 'Children');
axes(children(3)) % Magnitude Plot
semilogx(input_frequencies_2order, magnitude_ratio_2order, 'color', 'r'); ylim([-40 25]);
axes(children(2))  % Phase Plot
semilogx(input_frequencies_2order, phase_shift_2order, 'color', 'r')
text(195, -75, '{\omega}_{break} = 300 Hz')
legend('Theoretical', 'Actual'); hold off

%Solving for undamped natural frequency and damping ratio using
%measurements
sys = tf(theo_transfer_function2);
maxratio = max(magnitude_ratio_2order);
idwnf = find_closest_to(maxratio, magnitude_ratio_2order);
undamp_nat_freq = input_frequencies_2order(idwnf);
damp_ratio = 1 / (2 * maxratio);

%Solving for L and C values with measured undamped natural frequency and
%damping ratio
C2_2 = sqrt(1 / ((undamp_nat_freq^2)*4*(R2^2)*(damp_ratio^2)) );
L2_2 = C2_2*4*(R2^2)*(damp_ratio^2);

%% 3.1 Frequency Response of a First Order System Using LabView

%Importing the excel sheets of data
importdata1 = xlsread('Porter_Lab1_3.1_3.2.xlsx',1);
importdata2 = xlsread('Porter_Lab1_3.1_3.2.xlsx',2);
importdata3 = xlsread('Porter_Lab1_3.1_3.2.xlsx',3);
importdata4 = xlsread('Porter_Lab1_3.1_3.2.xlsx',4);

firstorderphase = importdata3((8:1010),:);
firstordermag = importdata4((8:1010),:);
secondorderphase = importdata1((8:1010),:);
secondordermag = importdata2((8:1010),:);

firstorderphase(:,1) = firstorderphase(:,1) .* 6.28; %Hz ---> Db
firstordermag(:,1) = firstordermag(:,1) .* 6.28; %Hz ---> Db
secondorderphase(:,1) = secondorderphase(:,1) .* 6.28; %Hz ---> Db
secondordermag(:,1) = secondordermag(:,1) .* 6.28; %Hz ---> Db

handle = figure(5);
bode(theo_transfer_function1, '--'); hold on; grid on
children = get(handle, 'Children');
axes(children(3)) % Magnitude Plot
semilogx(firstordermag(:,1), firstordermag(:,2), 'color', 'r');
axes(children(2))  % Phase Plot
semilogx(firstorderphase(:,1), firstorderphase(:,2), 'color', 'r');
%text(50, -45, '{\omega}_{break} = 255 Hz')
legend('Theoretical', 'Actual'); hold off

wbreak_firstorder = 240;
wbreak_secondorder = 355;

tau_firstorder = 1 / (2* pi * wbreak_firstorder);

C3_1 = tau_firstorder / R1;

%% 3.2 Frequency Response of a Second Order System Using LabView

handle = figure(6);
bode(theo_transfer_function2, '--'); hold on; grid on
children = get(handle, 'Children');
axes(children(3)) % Magnitude Plot
semilogx(secondordermag(:,1), secondordermag(:,2), 'color', 'r');
axes(children(2))  % Phase Plot
semilogx(secondorderphase(:,1), secondorderphase(:,2), 'color', 'r');
%text(50, -45, '{\omega}_{break} = 255 Hz')
legend('Theoretical', 'Actual'); hold off

maxratio_labview = max(secondordermag(:,2));
idwnf_labview = find_closest_to(maxratio_labview, secondordermag(:,2));
undamp_nat_freq_labview = secondordermag((idwnf_labview),1);
damp_ratio_labview = 1 / (2 * maxratio_labview);

%Solving for L and C values with labview undamped natural frequency and
%damping ratio
C3_2 = sqrt(1 / ((undamp_nat_freq_labview^2)*4*(R2^2)*(damp_ratio_labview^2)) );
L3_2 = C3_1*4*(R2^2)*(damp_ratio_labview^2);


%% Functions 
function [idx] = find_initial_response(voltage, window, thres)

    idx = 1;
    change = 0;
    while change < thres
        change = abs(voltage(idx+window) - voltage(idx));
        idx = idx + 1;
    end
    
    return
end

function [idx] = find_closest_to(value, array)

    error = abs(array - value);
    best = min(error);
    idx = find(error == best);
    
    idx = idx(1);
end

function omega_natural = find_undamped_natural_frequency(time, peak_indices, zeta)

    wave_length = time(peak_indices(end)) - time(peak_indices(1));
    period = wave_length / (length(peak_indices) - 1);
    omega_damped = 2*pi/period;
    
    omega_natural = omega_damped / sqrt(1 - zeta^2);
end

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
