function [var] = comp_modulus(param,var)
%Compute lambda/mu/gamma
%input: param,var output: update in var
%only depends on damage variable
var.lambda_ = param.lambda_0;
var.mu_ = param.mu_0 + var.alpha_ * param.xi_0 * param.gamma_r;
var.gamma_ = var.alpha_ * param.gamma_r;
end