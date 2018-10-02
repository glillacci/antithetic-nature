% Model of antithetic integrator
% Activation version with saturation, maturation, and explicit LuxR
% Gabriele Lillacci
%
% NULL --(k1*Y)--> Z1
% NULL --(k2)--> Z2
% Z1 + Z2 --(k3)--> NULL
% NULL --(k4*Z1)--> Y
% Y --(k6)--> X


%%
function yp = api_nd_sat_mat_LuxR (t,y)

% State variables
Z1 = y(1);
Z2 = y(2);
Y = y(3);
X = y(4);
L = y(5);

% Estimated parameters
k1 = y(6);
k2 = y(7);
k3 = y(8);
Vmax = y(9);
K = y(10);
k5 = y(11);
d = y(12);
w = y(13);
k6 = y(14);
k7 = y(15);

% Dynamics
Lp = 0.1 - d*L;
Z1p = k1*L - k3*Z1*Z2 - d*Z1;
Z2p = k2*Y - k3*Z1*Z2 - d*Z2;
Yp = Vmax*Z1./(K+Z1) - (k5+w+d)*Y;
Xp = k6*Y - k7*X;

yp = [Z1p; Z2p; Yp; Xp; Lp; zeros(10,1)];

end