clear all; 
close all;

%% 1 Position and Velocity Measurements 

tach_sensitivity = 0.007 / 0.10472; % V/rpm -> V/(rad./s)
pot_sensitivity = 0.015 / 0.01745; % V/deg. -> V/rad.

Dataset13 = readmatrix('WaveformData_1.3_redo.csv');
t1 = linspace(0, 0.00001408* 100000, 100000);
V_tach1 = Dataset13(6:end, 1);
V_pot1 = Dataset13(6:end, 2);

omega1 = V_tach1 / tach_sensitivity; % rad./s
theta1 = V_pot1 / pot_sensitivity; % rad.

%1a
figure(1)
yyaxis left
plot(t1(1:end-5), omega1)
ylabel('\omega (rad./s)')
yyaxis right
plot(t1(1:end-5), theta1)
ylabel('\theta (rad.)')
xlabel('Time (s)')
xlim([0.3 0.75]);

figure(2)
plot(theta1, omega1)
xlabel('\theta (rad.)'); ylabel('\omega (rad./s)')

%1b
% Theta begins with some initial displacement, and omega is initially 0. As
% the gear oscillates both signal have approximately the same period and damping 
% effects, but omega is lagged behind. 

%1c
% The signal starts at the outter edge of the spiral and goes inwards as
% time moves forward. This represents the damping response as the magnitude
% of the oscillations decreases gradually. Eventually, the signal dies out
% when omega is equal to 0 at the center of the spiral. 

%1d
theta_int1 = cumtrapz(t1(1:end-5), omega1);

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

%2a 
% This is given in the module 12 velocity sensors video (~29mins)

%2b
wire_length = 118.07 * (0.00328); % mm -> ft
wire_dia = 6.81 * (0.00328); % mm -> ft
k_exp = (0.253) * ( (5/12)/wire_length );

%2c
% I confirmed that using the theta and omega responses produce very similar
% zeta and omega_n parameters. I'm choosing the use the omega response
% because peakfinder is able to get more peaks. 
[peak_idx1, peak_mag1] = peakfinder(omega1, 2.5);
zeta1 = mean(find_damping_ratios(peak_mag1, 0));
omega_natural1 = find_undamped_natural_frequency(t1, peak_idx1, zeta1);

% Plot to confirm the peaks look good 
% figure(4)
% plot(t1(1:end-5), omega1); hold on
% plot(t1(peak_idx1), peak_mag1, 'o')

%2d 
J_sys_exp = k_exp / omega_natural1^2; % lbf-ft-s^2
B_sys_exp = 2*zeta1*k_exp / omega_natural1; 


% 2.e
J_pot = 2.44*10^-3 / (16*12); % ozf-in-sec^2 -> lbf-ft-s^2
J_tach = 1.32*10^-4 / (16*12); % ozf-in-sec^2 > lbf-ft-s^2
J_gear = 0.5 * (0.268/32.174) * (2.5/2)^2 * (1/12) * (1/12); % lbf-ft-s^2
J_coup_pt = 0.5 * (0.0249/32.174) * (0.744/2)^2 * (1/12) * (1/12); % lbf-ft-s^2
J_coup_ps = 0.5 * (0.0253/32.174) * (0.75/2)^2 * (1/12) * (1/12); % lbf-ft-s^2

shaft_length = 114.8 * (0.00328); % mm -> ft: 
shaft_dia = 11.44 * (0.00328); % mm -> ft: 
shaft_volume = shaft_length * pi*(shaft_dia/2)^2; % ft^3
shaft_mass = shaft_volume * (501.297); % lbm
J_shaft = 0.5 * (shaft_mass/32.174) * (shaft_dia/2)^2; % lbf-ft-s^2

J_sys_theo = J_pot + J_tach + J_gear + J_coup_pt + J_coup_ps + J_shaft; % lbf-ft-s^2

J_wire = (pi/32)*(wire_dia)^4; % ft^4
G_wire = wire_length*k_exp/J_wire;
k_theo = J_wire*G_wire / wire_length;

%% 3 Integration of Velocity and Signal Drift

%A Ideal Transfer Function
R1 = 10000; %Ohms
C1 = 10 * (10^-6); %Farads

num = 1;
den = [-1/(C1*R1) 0];

Tf = tf(num,den);

bode(Tf)

%Ideal OP Amp has the 2 other terms that give rise to a drift. Integrating
%the error and integrating the input current
% Module 6 Part 2 end of the video

%B Angular Displacements

Datasetiii = readmatrix('WaveformData(iii).csv');

E_int = Datasetiii(6:end, 1);
E_pot = Datasetiii(6:end, 2); 
E_time = Datasetiii(6:end, 3);

theta_int = E_int / tach_sensitivity; % rad./s
gain = 1/(R1*C1);
theta_int = theta_int / gain; %rad

theta_pot = E_pot / pot_sensitivity; %rad

figure(4)
plot(E_time, theta_int, 'r')
hold on
plot(E_time, theta_pot, 'b')
hold on
legend('Theta Integrator', 'Theta Potentiometer')
ylabel('Theta [rad]')
xlabel('Time');
title('Integrator Angle vs Potentiometer Theta')

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

function omega_natural = find_undamped_natural_frequency(time, peak_indices, zeta)

    wave_length = time(peak_indices(end)) - time(peak_indices(1));
    period = wave_length / (length(peak_indices) - 1);
    omega_damped = 2*pi/period;
    
    omega_natural = omega_damped / sqrt(1 - zeta^2);
end
