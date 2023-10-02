function [var] = comp_sigma(param,var)
%Compute stress (solid state + granular state)
%var.sigma_s = ( var.lambda_ - var.gamma_ / var.xit ) * var.I1t * eye(3) ...
%            + ( 2 * var.mu_ - var.gamma_ * var.xit ) * var.eps; %total strain
var.sigma_s = ( var.lambda_ - var.gamma_ / var.xi ) * var.I1 * eye(3) ...
            + ( 2 * var.mu_ - var.gamma_ * var.xi ) * var.eps_e; %total strain
var.sigma_b = ( 2 * param.a2 + param.a1 / var.xi + 3 * param.a3 * var.xi ) * var.I1 * eye(3) ...
            + ( 2 * param.a0 + param.a1 * var.xi - param.a3 * var.xi ^ 3 ) * var.eps_e; %elastic strain
var.sigma   = ( 1 - var.B ) * var.sigma_s + var.B * var.sigma_b;
end