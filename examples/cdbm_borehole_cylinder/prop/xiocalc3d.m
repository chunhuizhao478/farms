function [xi_out] = xiocalc3d(lambda_in, mu_in, phi_in)
%calculate xi_o in 3d
phi_in = phi_in * pi / 180;

q = sin(phi_in) / ( 1 - sin(phi_in) / 3 );

xi_out = -sqrt(3) / sqrt( 2 * q ^ 2 * ( lambda_in / mu_in + 2 / 3 ) ^ 2 + 1 );

end

