function [xi_out] = xiocalc2d(lambda_in, mu_in, phi_in)
%calculate xi_o in 2d
phi_in = phi_in * pi / 180;
xi_out = -sqrt(2) / sqrt(1+(lambda_in/mu_in+1)^2*sin(phi_in)*sin(phi_in));
end

