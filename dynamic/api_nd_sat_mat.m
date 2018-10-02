% Model of antithetic integrator
% Activation version with saturation and maturation
% Gabriele Lillacci
%
% NULL --(k1*Y)--> Z1
% NULL --(k2)--> Z2
% Z1 + Z2 --(k3)--> NULL
% NULL --(k4*Z1)--> Y
% Y --(k6)--> X


%%
function yp = api_nd_sat_mat (t,y)

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

% Dynamics
Z1p = k1 - k3*Z1*Z2 - d*Z1;
Z2p = k2*Y - k3*Z1*Z2 - d*Z2;
Yp = Vmax*Z1./(K+Z1) - (k5+w+d)*Y;
Xp = k6*Y - k7*X;

yp = [Z1p; Z2p; Yp; Xp; zeros(10,1)];

end