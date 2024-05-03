format long
%% Calculate theta_init
L = 0.02; a = 0.008; b = 0.012;
fo = 0.6; Vini = 1e-12; Vo = 1e-6;
tau_init = 70e6; sigma_init = 120e6;
%
theta_init = ( L / Vo ) * exp( ( a * log( 2 * sinh( tau_init / ( a * sigma_init ))) - fo - a * log(Vini / Vo) ) / b );
disp(theta_init)