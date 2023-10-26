function[a0_out,a1_out,a2_out,a3_out] = fb_coeff(param, xi_given)
%input: param (struct); output: None (change in place "param")
%Solve for parameters of free energy in granular state 
%clear all; close all; clc
%xi_given = 0;
%% Intermediate parameters %%
mu_cr = param.mu_0 + comp_alpha_cr(xi_given,param) * param.xi_0 * param.gamma_r;
gamma_cr = comp_alpha_cr(xi_given,param) * param.gamma_r;

mu_ = param.mu_0 + 1 * param.xi_0 * param.gamma_r;
gamma_ = 1 * param.gamma_r;
lambda_ = param.lambda_0;
%% Equations Setup %%
syms a0 a1 a2 a3 
eqn1 = 2 * a2 + a1 / param.xi_1 + 3 * a3 * param.xi_1;
eqn2 = 2 * a0 + a1 * param.xi_1 - a3 * param.xi_1 ^ 3;
%eqn3 = a0 - param.chi * mu_cr;
eqn3 = a0 + a1 * (xi_given) + a2 * (xi_given) ^ 2 + a3 * (xi_given) ^ 3 - ( param.chi * ( mu_cr - gamma_cr * (xi_given) + lambda_/2 * (xi_given) ^ 2 ) );
eqn4 = a0 + a1 * param.xi_d + a2 * param.xi_d ^ 2 + a3 * param.xi_d ^ 3 - ( mu_ - gamma_ * param.xi_d + lambda_ / 2 * param.xi_d ^ 2 );
aout = solve(eqn1==0,eqn2==0,eqn3==0,eqn4==0,a0,a1,a2,a3);
a0_out = aout.a0;
a1_out = aout.a1;
a2_out = aout.a2;
a3_out = aout.a3;
end
