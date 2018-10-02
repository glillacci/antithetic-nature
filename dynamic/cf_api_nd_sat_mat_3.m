% Antithetic integrator
% Cost function for fitting step responses

function [F, sol, solnw, sols] = cf_api_nd_sat_mat_3 (theta, model, data, cond, time)

% Extract parameters from theta vector
h0 = theta(1);
h1 = theta(2);
a0 = theta(3);
a1 = theta(4);
k3 = theta(5);
Vmax = theta(6);
K = theta(7);
k6 = theta(8);
k7 = theta(9);
C = theta(10);
d = log(2)/25;
k5 = 0.5*d;

opts = odeset ('Jacobian', @api_nd_sat_mat_jac, 'AbsTol', 1e-10);

% Set the initial conditions for the model
x0 = [0 0 0 0];

% Initialize cost function
F = 0;

% Initialize model solution
sol = zeros(2, length(time));

% Solve the model for all the conditions
for ii = 1:size(cond,2)
    % Solve model (without disturbance)
    solnw = ode15s(model, [time(1) time(end)], [x0 h0+h1*cond(1,ii)...
        a0+a1*cond(2,ii) k3 Vmax K k5 d cond(3,ii) k6 k7], opts);
    
    % Compare model solution and data
    model_sol = deval(solnw, time);
    error = data(ii,:) - C*model_sol(4,:);
    F = F + sum(error.^2);
    
    % Output model solution
    sol(ii,:) = model_sol(4,:);
    sols(ii) = solnw;
end

end
