% Antithetic integrator
% Activation version


%% Initialize

clear all
close all
clc

% NOMINAL parameters
% Production rate of controller species z1
k1 = 1;
% Induction of controller species z2 by Y
k2 = 1;
% Annihilation rate
k3 = 1;
% Production of output Y by controller species z1
Vmax = 1;
K = 1;
% Native degradation of Y
k5 = 1;
% Dilution rate
d = 1*log(2)/25;
% Disturbance
w = 5*d;

% Initial conditions
x0 = [0 0 0];

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
solw = ode15s (@api_nd_growth, [TI TF], [x0 k1 k2 k3 Vmax K k5 d w]);
% Check the achieved norm
while (norm(solw.y(1:3,end) - solw.y(1:3,end-1), Inf) > thr && TF < Tmax)
    % If accuracy not met, extend simulation time
    TF = TF + 100;
    % Extend the solution
    solw = odextend (solw, [], TF);
end

% Simulate with disturbance
TI = 0; TF = 100;
% Solution in first time interval
solnw = ode15s (@api_nd_growth, [TI TF], [x0 k1 k2 k3 Vmax K k5 d 0]);
% Check the achieved norm
while (norm(solnw.y(1:3,end) - solnw.y(1:3,end-1), Inf) > thr && TF < Tmax)
    % If accuracy not met, extend simulation time
    TF = TF + 100;
    % Extend the solution
    solnw = odextend (solnw, [], TF);
end

% Calculate steady state error
ss1 = solnw.y(3,end);
ss2 = solw.y(3,end);

disp ('Steady state error')
disp (abs(ss1-ss2)/ss1);


%% Plot

figure;
plot(solnw.x,solnw.y(3,:), solw.x, solw.y (3,:), 'LineWidth', 3)
legend ('no disturbance', 'disturbance', 'Location', 'best')
xlabel ('time (min)')
ylabel ('protein concentration (nM)')
