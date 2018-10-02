% Least-squares fitting of antithetic PI model
% 20180630
% Gabriele Lillacci


%% Initialize

clear all
close all
clc

% Load data
load 20180704_steps.mat
data = (outm-outm(:,1))/1000;


%% Initial guesses for the parameters

% Coefficients of linear map between HSL and k1
% h0 = -978; h1 = 207;
% Coefficients of linear map between ARA and k2
% a0 = 0.005; a1 = 0.16;
% Production rate of controller species z1
k1 = 10^(1.5);
% Induction of controller species z2 by Y
k2 = 10^(-1.5);
% Annihilation rate
k3 = 0.05;
% Production of output Y by controller species z1
Vmax = 10000;
K = Vmax/(2*5);
% Maturation parameters
k6 = log(2)/20;
k7 = k6;
% Output scaling factor
C = 0.0005;


% Compute matrix A
A = zeros(4,4);
for ii=[1 2 3 4]
    A(ii,:) = [C*1 C*cond(1,ii) -mean(data(ii,7:end)) -mean(data(ii,7:end))*cond(2,ii)];
end
% Compute a rank-2 approximation of A
[U,S,V] = svd(A);
S2 = S; S2(3,3) = 0; S2(4,4) = 0;
A2 = U*S2*V';
% Compute the null space of A2
N = null(A2);
n = 10*N(:,1) - 10*N(:,2);
% Use it to set the initial conditions
h0=n(1); h1=n(2); a0=n(3); a1=n(4);

% Gather all initial guesses
x0 = [h0 h1 a0 a1 k3 Vmax K k6 k7 C];

% Set upper and lower bounds for each parameter
ub = [ Inf  Inf  Inf  Inf 1    20000 Inf log(2)/15 Inf 1e-3];
lb = [-Inf  eps -Inf  eps eps  eps   eps log(2)/50 eps eps];

% Set the linear constraints for HSL and ARA
Ain = [-1 -5  0  0    0 0 0 0 0 0;
      0  0 -1 -0.15 0 0 0 0 0 0 ];
bin = [0; 0];


%% Minimize the cost function

% Randomize the initial guesses
% x0 = x0 + 0.05*x0.*randn(size(x0));

% Set up optimization function
fun = @(theta) cf_api_nd_sat_mat_3 (theta, @api_nd_sat_mat, data, cond, time);

% Set up optimization options
oo = optimset ('Display', 'iter', 'MaxIter', 1000, 'MaxFunEval', 1e6);

% Run the least-squares fitting
x = fmincon (fun, x0, Ain, bin, [], [], lb, ub, [], oo);


%% Save the results

% save theta_nosat.mat x


%% Plot

figure;
[F,sol,solnw] = fun(x0);
plot(time,data,time,x0(10)*sol)
title('Initial conditions')


figure;
[F,sol,solnw,sols] = fun(x);
hold on
markers = ['b', 'r', 'g', 'k'];
for ii=1:4
    plot(time, data(ii,:), 'o', 'MarkerEdgeColor', markers(ii), 'MarkerSize', 8)
end
for ii=1:4
    plot(sols(ii).x, sols(ii).y(4,:)*x(10), markers(ii), 'LineWidth', 2)
end
title('Fit')
xlabel('Time (min)')
ylabel('sfGFP (a.u. x 1000)')
legend('5 nM HSL 0.15% ARA', '10 nM HSL 0.15% ARA', '7.5 nM HSL 0.2% ARA', '8.5 nM HSL 0.2% ARA', 'Location', 'best')


figure;
plot(sols(4).x, sols(4).y(1:2,:))
legend('z1', 'z2', 'Location', 'best')

figure;
plot(sols(4).x, sols(4).y(3:4,:))
legend('arac', 'gfp', 'Location', 'best')

