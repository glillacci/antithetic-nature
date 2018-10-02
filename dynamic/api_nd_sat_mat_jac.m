% Jacobian for Model of antithetic integrator
% Activation version with saturation and maturation
% Gabriele Lillacci
%
% NULL --(k1*Y)--> Z1
% NULL --(k2)--> Z2
% Z1 + Z2 --(k3)--> NULL
% NULL --(k4*Z1)--> Y
% Y --(k6)--> X


%%
function J = api_nd_sat_mat_jac (t,y)

% State variables
Z1 = y(1);
Z2 = y(2);
Y = y(3);
X = y(4);

% Estimated parameters
k1 = y(5);
k2 = y(6);
k3 = y(7);
Vmax = y(8);
K = y(9);
k5 = y(10);
d = y(11);
w = y(12);
k6 = y(13);
k7 = y(14);

% Initialize the Jacobian
J = zeros(14);

% Fill in the non-zero elements
J(1,:) = [-d-Z2*k3, -Z1*k3, 0, 0, 1, 0, -Z1*Z2, 0, 0, 0, -Z1, 0, 0, 0];
J(2,:) = [-Z2*k3, -d-Z1*k3, k2, 0, 0, Y, -Z1*Z2, 0, 0, 0, -Z2, 0, 0, 0];
J(3,:) = [(K*Vmax)/(K + Z1)^2, 0, -d-k5-w, 0, 0, 0, 0, Z1/(K + Z1), -(Vmax*Z1)/(K + Z1)^2, -Y, -Y, -Y, 0, 0];
J(4,:) = [0, 0, k6, -k7, 0, 0, 0, 0, 0, 0, 0, 0, Y, -X];

end