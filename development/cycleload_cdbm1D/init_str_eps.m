function [var] = init_str_eps(param,var)
%% Initialize sigma_s, sigma_b, sigma
var.sigma_s(1,1) = -param.P_c;
var.sigma_s(2,2) = -param.P_c;
var.sigma_s(3,3) = -param.P_c;
var.sigma = (1 - var.B) * var.sigma_s + var.B * var.sigma_b;
%% Initialize eps, eps_e, eps_p
%%total strain tensor (under 3D compaction)%%
var.eps(1,1) = 1 / var.E * ( -param.P_c + 2 * var.nu * param.P_c );
var.eps(2,2) = 1 / var.E * ( -param.P_c + 2 * var.nu * param.P_c );
var.eps(3,3) = 1 / var.E * ( -param.P_c + 2 * var.nu * param.P_c );
var.eps_e = var.eps - var.eps_p;
end