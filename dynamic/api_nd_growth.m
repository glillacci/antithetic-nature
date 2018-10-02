% Model of antithetic integrator
% Activation version with saturation
% Gabriele Lillacci
%
% NULL --(k1*Y)--> Z1
% NULL --(k2)--> Z2
% Z1 + Z2 --(k3)--> NULL
% NULL --(k4*Z1)--> Y


%%
function yp = api_nd_growth (t,y)

% State variables
Z1 = y(1);
Z2 = y(2);
Y = y(3);

% Estimated parameters
k1 = y(4);
k2 = y(5);
k3 = y(6);
Vmax = y(7);
K = y(8);
k5 = y(9);
d = y(10);
w = y(11);

% Dynamics
Z1p = k1 - k3*Z1*Z2 - d*Y*Z1;
Z2p = k2*Y - k3*Z1*Z2 - d*Y*Z2;
Yp = Vmax*Z1./(K+Z1) - (k5+w+d*Y)*Y;

yp = [Z1p; Z2p; Yp; zeros(8,1)];

end