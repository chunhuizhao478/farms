function [var] = comp_eps_sys(param, var)
%Compute total strain/plastic strain/elastic strain
%% 1.Get variables from var struct
eps_pre = var.eps;
%% 2.Update strain in the x-dir
eps11 = eps_pre(1,1) + var.depsx;
%% 3.Define symbolic unknowns eps_22 eps_33
syms eps22 eps33 
%I1 = eps11 + eps22 + eps33;
%I2 = eps11 ^ 2 + eps22 ^ 2 + eps33 ^ 2;
%xi = I1 / sqrt(I2);
sigma11 = var.E / ( ( 1 + var.nu ) * ( 1 - 2 * var.nu ) ) * ( ( 1 - var.nu ) * eps11 + var.nu * ( eps22 + eps33 ) ); 
%sigma11 = ( var.lambda_ - var.gamma_ / xi ) * I1 + ( 2 * var.mu_ - var.gamma_ * xi ) * eps11;
%% 4.Setup equations
eqn2 = eps22 - 1 / var.E * ( - param.P_c - var.nu * ( sigma11 - param.P_c ));
eqn3 = eps33 - 1 / var.E * ( - param.P_c - var.nu * ( sigma11 - param.P_c ));
%eqn2 = eps22 - 1 / ( 2 * mu_in - gamma_in * xi) * ( - param.P_c - (lambda_in - gamma_in/xi ) * I1 );
%eqn3 = eps33 - 1 / ( 2 * mu_in - gamma_in * xi) * ( - param.P_c - (lambda_in - gamma_in/xi ) * I1 );
%% 5.Solve equations
aout = vpasolve(eqn2==0,eqn3==0,eps22,eps33);
eps22_new = eval(aout.eps22);
eps33_new = eval(aout.eps33);
%% 6.Update variables
eps11_new = eps11;
%% 7.Obtain eps (total strain)
var.eps(1,1) = eps11_new;
var.eps(2,2) = eps22_new;
var.eps(3,3) = eps33_new;
%% 8.while loop for computing elastic strain
[var] = loop_eps_e(param,var);
end