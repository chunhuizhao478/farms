function [F] = eps_func_epst(x,eps11_inc,eps_p_pre,param,var)
%Construct system of equations to solve
%{
%%%input%%%
%%eps11t: total strain in vertical direction
eps11_inc: strain incremenet in vertical direction
eps_p_pre: plastic strain from previous time step
x(1):eps22t_inc
x(2):eps33t_inc
%%%output%%%
F(1),F(2)
%}
%% Define total strain
eps11t = var.eps(1,1) + eps11_inc;
eps22t = var.eps(2,2) + x(1);
eps33t = var.eps(3,3) + x(2);
%% Define plastic strain (assume using previous iteration)
eps11p = eps_p_pre(1,1); 
eps22p = eps_p_pre(2,2); 
eps33p = eps_p_pre(3,3); 
%% Define elastic strain
eps11 = eps11t - eps11p;
eps22 = eps22t - eps22p;
eps33 = eps33t - eps33p;
%% Represent I1,I2,xi in terms of unknown variables
%total strain
%I1t = eps11t + eps22t + eps33t;
%I2t = eps11t ^ 2 + eps22t ^ 2 + eps33t ^ 2;
%xit = I1t / sqrt(I2t);
%elastic strain
I1 = eps11 + eps22 + eps33;
I2 = eps11 ^ 2 + eps22 ^ 2 + eps33 ^ 2;
xi = I1 / sqrt(I2);
%% Represent sigma (solid(s) + granular(b))
%sigma11_s = ( var.lambda_ - var.gamma_ / xit ) * I1t ...
%          + ( 2 * var.mu_ - var.gamma_ * xit ) * eps11t;
sigma11_s = ( var.lambda_ - var.gamma_ / xi ) * I1 ...
          + ( 2 * var.mu_ - var.gamma_ * xi ) * eps11;
sigma11_b = ( 2 * param.a2 + param.a1 / xi + 3 * param.a3 * xi ) * I1 ...
          + ( 2 * param.a0 + param.a1 * xi - param.a3 * xi ^ 3 ) * eps11;
%sigma22_s = ( var.lambda_ - var.gamma_ / xit ) * I1t ...
%          + ( 2 * var.mu_ - var.gamma_ * xit ) * eps22t;
sigma22_s = ( var.lambda_ - var.gamma_ / xi ) * I1 ...
          + ( 2 * var.mu_ - var.gamma_ * xi ) * eps22;
sigma22_b = ( 2 * param.a2 + param.a1 / xi + 3 * param.a3 * xi ) * I1 ...
          + ( 2 * param.a0 + param.a1 * xi - param.a3 * xi ^ 3 ) * eps22;
%sigma33_s = ( var.lambda_ - var.gamma_ / xit ) * I1t ...
%          + ( 2 * var.mu_ - var.gamma_ * xit ) * eps33t;
sigma33_s = ( var.lambda_ - var.gamma_ / xi ) * I1 ...
          + ( 2 * var.mu_ - var.gamma_ * xi ) * eps33;
sigma33_b = ( 2 * param.a2 + param.a1 / xi + 3 * param.a3 * xi ) * I1 ...
          + ( 2 * param.a0 + param.a1 * xi - param.a3 * xi ^ 3 ) * eps33;
%% Setup five equations, five unknowns
F(1) = param.P_c + (1 - var.B) * sigma22_s + var.B * sigma22_b;
F(2) = param.P_c + (1 - var.B) * sigma33_s + var.B * sigma33_b;
end