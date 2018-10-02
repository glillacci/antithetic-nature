% Least-squares fitting of antithetic PI model
% 20180409
% Gabriele Lillacci


%% Initialize

clear all
close all
clc

% Load data
load 20180504_steps.mat


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
k3 = 0.2;
% Production of output Y by controller species z1
Vmax = 5000;
K = Vmax/(2*1);
% Maturation parameters
k6 = log(2)/30;
k7 = k6;
% Output scaling factor
C = 0.0005;


% Compute matrix A
A = zeros(4,4);
for ii=[1 2 3 4]
    A(ii,:) = [C*1 C*cond(1,ii) -mean(data(ii,15:end)) -mean(data(ii,15:end))*cond(2,ii)];
end
% Compute a rank-2 approximation of A
[U,S,V] = svd(A);
S2 = S; S2(3,3) = 0; S2(4,4) = 0;
A2 = U*S2*V';
% Compute the null space of A2
N = null(A2);
n = -50*N(:,1) -250*N(:,2);
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
[F,sol,solnw] = fun(x);
plot(time,data,time,x(10)*sol)

figure;
plot(solnw.x, solnw.y(1:2,:))
legend('z1', 'z2', 'Location', 'best')

[F,sol,solw] = cf_api_nd_sat_mat_2 (x, @api_nd_sat_mat, data, time);

figure;
plot(solnw.x, solnw.y(3:4,:))
legend('arac', 'gfp', 'Location', 'best')


figure;
[F,sol,solnw] = fun(x0);
plot(time,data,time,x0(10)*sol)