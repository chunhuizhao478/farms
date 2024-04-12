function [var] = comp_eps_sys_loop(param,var)
%Compute total strain / elastic strain / plastic strain
%% Variables @ previous time step
eps_t_pre = var.eps;
eps_e_pre = var.eps_e;
eps_p_pre = var.eps_p;
eps_t_inc11 = var.depsx;

eps_t_iter11 = eps_t_pre(1,1) + eps_t_inc11;
eps_t_iter22 = eps_t_pre(2,2);
eps_t_iter33 = eps_t_pre(3,3);

eps_p_iter11 = eps_p_pre(1,1);
eps_p_iter22 = eps_p_pre(2,2);
eps_p_iter33 = eps_p_pre(3,3);

eps_e_iter11 = eps_t_iter11 - eps_p_iter11;
eps_e_iter22 = eps_t_iter22 - eps_p_iter22;
eps_e_iter33 = eps_t_iter33 - eps_p_iter33;

I1_iter = eps_e_iter11 + eps_e_iter22 + eps_e_iter33;
I2_iter = eps_e_iter11 ^ 2 + eps_e_iter22 ^ 2 + eps_e_iter33 ^ 2;
xi_iter = I1_iter / sqrt(I2_iter);


end