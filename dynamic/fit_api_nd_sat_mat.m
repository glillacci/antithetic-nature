% Least-squares fitting of antithetic PI model
% 20180409
% Gabriele Lillacci


%% Initialize

clear all
close all
clc

% Load data
load 20180406_steps.mat
% Select first step and remove background
data = outm - outm(:,1);


%% Initial guesses for the parameters

% Coefficients of linear map between HSL and k1
h0 = 0.1; h1 = 32;
% Coefficients of linear map between ARA and k2
a0 = -0.03; a1 = 0.316;
% Production rate of controller species z1
k1 = 10^(1.5);
% Induction of controller species z2 by Y
k2 = 10^(-1.5);
% Annihilation rate
k3 = 0.02;
% Production of output Y by controller species z1
Vmax = 5000;
K = Vmax/(2*1);
% Maturation parameters
k6 = log(2)/30;
k7 = k6;
% Output scaling factor
C = 0.0005;

% Gather all initial guesses
x0 = [h0 h1 a0 a1 k3 Vmax K k6 k7 C];

% Set upper and lower bounds for each parameter
ub = [Inf  Inf Inf  Inf 0.05 10000 Inf log(2)/15 Inf Inf];
lb = [-Inf eps -Inf eps eps  eps   eps log(2)/50 eps eps];

%% Minimize the cost function

% Randomize the initial guesses
%x0 = x0 + x0.*randn(size(x0));

% Set up optimization function
fun = @(theta) cf_api_nd_sat_mat (theta, @api_nd_sat_mat, data, time);

% Set up optimization options
oo = optimset ('Display', 'iter', 'MaxIter', 1000, 'MaxFunEval', 1e6);

% Run the least-squares fitting
x = fmincon (fun, x0, [], [], [], [], lb, ub, [], oo);


%% Save the results

% save theta_nosat.mat x


%% Plot

[F,sol] = fun(x);
plot(time,data,time,x(10)*sol)
