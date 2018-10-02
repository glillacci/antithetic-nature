% Model of antithetic integrator
% "Open-loop" topology suggested by Science review
% Gabriele Lillacci
%


%%
function yp = api_nd_sat_science_review (t,y)

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

% Dynamics
Z1p = k1 - k3*Z1*Z2 - d*Z1;
Z2p = k2*Y - k3*Z1*Z2 - d*Z2;
Yp = Vmax*0.1 - (k5+w+d)*Y;
Xp = Vmax*Z1./(K+Z1) - (k5+w+d)*X;

yp = [Z1p; Z2p; Yp; Xp; zeros(8,1)];

end