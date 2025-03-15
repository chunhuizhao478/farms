function [var] = comp_eps_p(y,)
%Compute plastic strain
%eps_p_pre: plastic strain at the last time step (not changing within this
%time step)
% sigma_d = var.sigma - 1/3 * sum(diag(var.sigma)) * eye(3);
% var.eps_p = eps_p_pre + var.dt * param.C_g * var.B ^ param.m1 * sigma_d ^ param.m2;

end