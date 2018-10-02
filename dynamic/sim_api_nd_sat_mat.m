% Antithetic integrator
% Activation version with saturation and maturation


%% Initialize

clear all
close all
clc

% Load fitted parameters
load theta_20180410.mat

% Unpack parameters
h0 = x(1);
h1 = x(2);
a0 = x(3);
a1 = x(4);
k3 = x(5);
Vmax = x(6);
K = x(7);
d = log(2)/25;
k5 = 0.5*d;
k6 = x(8);
k7 = x(9);
C = x(10);

% Set disturbance
w = d;

% Set HSL and ARA
hsl = 10;
ara = 0.15;

% Calculate k1 and k2
k1 = h0 + hsl*h1;
k2 = a0 + ara*a1;

% Set initial conditions
x0 = [0 0 0 0];

% Set accuracy parameters
% Define accuracy threshold
thr = 1e-6;
% Define max final time
Tmax = 100000;
% Keep track of maximum simulation error
maxerr = 0;


%% Solve

% Simulate with disturbance
TI = 0; TF = 100;
% Solution in first time interval
solw = ode15s (@api_nd_sat_mat, [TI TF], [x0 k1 k2 k3 Vmax K k5 d w k6 k7]);
% Check the achieved norm
while (norm(solw.y(1:4,end) - solw.y(1:4,end-1), Inf) > thr && TF < Tmax)
    % If accuracy not met, extend simulation time
    TF = TF + 100;
    % Extend the solution
    solw = odextend (solw, [], TF);
end

% Simulate with disturbance
TI = 0; TF = 100;
% Solution in first time interval
solnw = ode15s (@api_nd_sat_mat, [TI TF], [x0 k1 k2 k3 Vmax K k5 d 0 k6 k7]);
% Check the achieved norm
while (norm(solnw.y(1:4,end) - solnw.y(1:4,end-1), Inf) > thr && TF < Tmax)
    % If accuracy not met, extend simulation time
    TF = TF + 100;
    % Extend the solution
    solnw = odextend (solnw, [], TF);
end

% Calculate steady state error
ss1 = solnw.y(4,end);
ss2 = solw.y(4,end);

disp ('Steady state error')
disp (abs(ss1-ss2)/ss1);


%% Plots


% Model fitting
figure;
load 20180406_steps.mat
data = outm - outm(:,1);
plot(time, data, solnw.x, C*solnw.y(4,:))
xlabel('Time (min)')
ylabel('Fluorescence (a.u.)')
title('Model fitting')

figure;
plot(solnw.x, solnw.y(3:4,:))
legend('AraC', 'GFP', 'Location', 'best')
xlabel('Time (min)')
ylabel('Concentration (nM)')
title([num2str(hsl) ' nM HSL - ' num2str(ara) '% ARA'])

figure;
plot(solnw.x, solnw.y(1:2,:))
legend('SigW', 'RsiW', 'Location', 'best')
xlabel('Time (min)')
ylabel('Concentration (nM)')
title([num2str(hsl) ' nM HSL - ' num2str(ara) '% ARA'])


figure;
plot(solnw.x,solnw.y(4,:), solw.x, solw.y (4,:), 'LineWidth', 3)
legend ('GFP no disturbance', 'GFP w/ disturbance', 'Location', 'best')
xlabel ('Time (min)')
ylabel ('Concentration (nM)')
title([num2str(hsl) ' nM HSL - ' num2str(ara) '% ARA'])
