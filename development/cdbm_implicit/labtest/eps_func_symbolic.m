function [eps22t_inc_out,eps33t_inc_out,eps11p_inc_out,eps22p_inc_out,eps33p_inc_out] = eps_func_symbolic(eps11_inc,eps_p_pre,param,var)
%Construct system of equations to solve
%{
%%%input%%%
%%eps11t: total strain in vertical direction
eps11_inc: strain incremenet in vertical direction
eps_p_pre: plastic strain from previous time step
x(1):eps22t_inc
x(2):eps33t_inc
x(3):eps11p_inc
x(4):eps22p_inc
x(5):eps33p_inc
%%%output%%%
F(1),F(2),F(3),F(4),F(5)
%}
syms eps22t_inc eps33t_inc eps11p_inc eps22p_inc eps33p_inc
%% Define total strain
eps11t = var.eps(1,1) + eps11_inc;
eps22t = var.eps(2,2) + eps22t_inc;
eps33t = var.eps(3,3) + eps33t_inc;
%% Define plastic strain
eps11p = eps_p_pre(1,1) + eps11p_inc;
eps22p = eps_p_pre(2,2) + eps22p_inc;
eps33p = eps_p_pre(3,3) + eps33p_inc;
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
sigma11_s = ( var.lambda_ - var.gamma_ / xi ) * I1 ...
          + ( 2 * var.mu_ - var.gamma_ * xi ) * eps11;
sigma11_b = ( 2 * param.a2 + param.a1 / xi + 3 * param.a3 * xi ) * I1 ...
          + ( 2 * param.a0 + param.a1 * xi - param.a3 * xi ^ 3 ) * eps11;
sigma22_s = ( var.lambda_ - var.gamma_ / xi ) * I1 ...
          + ( 2 * var.mu_ - var.gamma_ * xi ) * eps22;
sigma22_b = ( 2 * param.a2 + param.a1 / xi + 3 * param.a3 * xi ) * I1 ...
          + ( 2 * param.a0 + param.a1 * xi - param.a3 * xi ^ 3 ) * eps22;
sigma33_s = ( var.lambda_ - var.gamma_ / xi ) * I1 ...
          + ( 2 * var.mu_ - var.gamma_ * xi ) * eps33;
sigma33_b = ( 2 * param.a2 + param.a1 / xi + 3 * param.a3 * xi ) * I1 ...
          + ( 2 * param.a0 + param.a1 * xi - param.a3 * xi ^ 3 ) * eps33;
%% Represent deviatoric stress
sigma11_t = (1 - var.B) * sigma11_s + var.B * sigma11_b;
sigma22_t = (1 - var.B) * sigma22_s + var.B * sigma22_b;
sigma33_t = (1 - var.B) * sigma33_s + var.B * sigma33_b;

sigma_d11 = sigma11_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
sigma_d22 = sigma22_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
sigma_d33 = sigma33_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
%% Setup five equations, five unknowns
eqn1 = param.P_c + (1 - var.B) * sigma22_s + var.B * sigma22_b;
eqn2 = param.P_c + (1 - var.B) * sigma33_s + var.B * sigma33_b;
eqn3 = eps11p_inc - var.dt * param.C_g * var.B ^ param.m1 * sigma_d11 ^ param.m2;
eqn4 = eps22p_inc - var.dt * param.C_g * var.B ^ param.m1 * sigma_d22 ^ param.m2;
eqn5 = eps33p_inc - var.dt * param.C_g * var.B ^ param.m1 * sigma_d33 ^ param.m2;

aout = vpasolve(eqn1==0,eqn2==0,eqn3==0,eqn4==0,eqn5==0,eps22t_inc,eps33t_inc,eps11p_inc,eps22p_inc,eps33p_inc);

eps22t_inc_out = eval(aout.eps22t_inc);
eps33t_inc_out = eval(aout.eps33t_inc);
eps11p_inc_out = eval(aout.eps11p_inc);
eps22p_inc_out = eval(aout.eps22p_inc);
eps33p_inc_out = eval(aout.eps33p_inc);

end