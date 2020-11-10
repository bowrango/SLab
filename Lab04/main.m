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

%% 3. Vibration Analysis 

data=readmatrix('3.3');
t=data(6:end, 3);
acc = data(6:end, 1);
lvt = data(6:end, 2);
LVTslope=(lvt(60474)-lvt(56360))/(t(60474)-t(56360));
mass = 75.6/28.35; % [ozm]
g=32.17405*12; %in/sec^2
lvt_sensitivity = LVTslope / g; % [mV / (in./s)]
core_velocity = lvt  / lvt_sensitivity; % [in./s]
vdispl=(602.2-525); %mV
forcesense=g*mass;


accsens=-g/acc(44585); %mV/in/sec^2

% == Part a. ==

%  [peak_idx, peak_mag] = peakfinder(lvt, .25);
% [peak_idx1, peak_mag1] = peakfinder(lvt, .25,500,-1);
% [peak_idx2, peak_mag2] = peakfinder(acc, .5,100);
% [peak_idx3, peak_mag3] = peakfinder(acc, .1,500,-1);

 figure(3)
yyaxis left
plot(t,lvt)
hold on
plot(t(53914),lvt(53914),'o')
plot(t(45400),lvt(45400),'o')

ylabel('LVT (Velocity in/sec)')
hold off
yyaxis right
plot(t, acc)
 hold on
%plot(t(),acc(),'o')
%plot(t(),acc(),'o')
ylabel('Accelerometer in/sec^2)')
xlabel('Time (s)')

% == Part b. == 

new_acc=acc*accsens;
new_acc=detrend(new_acc);
v=cumtrapz(t(41274:end),new_acc(41274:end));

figure(4)
plot(t(41274:end),v,'b'); hold on
plot(t(41274:end),core_velocity(41274:end),'r')
xlabel('Time (s)'); ylabel('in/sec.')
legend('Integrated Signal', 'Core Velocity')

% == Part d. == 

ac=abs(core_velocity);
b=max(ac);

% == Part e. == 

data1=readmatrix('3.4');

t1=data1(6:end, 3);
force = data1(6:end, 2);
lvt1 = data1(6:end, 1);
 
figure(5)
yyaxis left
plot(t1,lvt1)
hold on

ylabel('LVT (mV)')
hold off
yyaxis right
plot(t1, force)
hold on
ylabel('Force Sensor mV)')
xlabel('Time (s)')

figure(6)
yyaxis left
plot(t1,lvt1/lvt_sensitivity)
hold on

ylabel('LVT (in/sec)')
hold off
yyaxis right
plot(t1, force/500)
hold on
ylabel('Force Sensor lbf)')
xlabel('Time (s)')

core_velocity1 = lvt1  / lvt_sensitivity; % [in./s]
ac1=abs(core_velocity1);
b=max(ac1);
%@index 40716 v=max
forcemaxv=force(40716)/500; %lbf
maxforce=max(abs(force))/500; %lbf
ssf=force(83856)/500; %lbf
massofcore=ssf/32.174049;









