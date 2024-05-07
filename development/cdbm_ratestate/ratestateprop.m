format long
%% Parameters Rate-and-state
% frictional properties
L = 0.02; 
a = 0.008; 
b = 0.012;
fo = 0.6; 
Vini = 1e-12; 
Vo = 1e-6;
tau_init = 70e6; 
sigma_init = 120e6;
% material properties
nu = 0.25;
mu = 3.204e10;
%% Calculate theta_init
theta_init = ( L / Vo ) * exp( ( a * log( 2 * sinh( tau_init / ( a * sigma_init ))) - fo - a * log(Vini / Vo) ) / b );
disp(theta_init)
%% Calculate nucleation size L_nuc
L_nuc = ( 2 * mu / ( 1 - nu ) * L * b ) / ( pi * sigma_init * ( b - a ) ^ 2 );
%% Calculate process zone size Lb
Lb = mu / ( 1 - nu ) * L / ( sigma_init * b );
%% Estimate dx
dx = Lb / 8;
%% Output Parameters
fprintf("mu (GPa): %.3f \n", mu);
fprintf("nu (-): %.3f \n", nu);
fprintf("Lb (m): %.3f \n", Lb);
fprintf("dx_max (m): %.3f \n", dx);