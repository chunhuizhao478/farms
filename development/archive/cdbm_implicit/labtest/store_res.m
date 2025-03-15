function [var] = store_res(var)
%Store results
var.time_list    = [var.time_list; var.t];
var.alpha_list   = [var.alpha_list; var.alpha_];
var.xi_list      = [var.xi_list; var.xi];
var.eps_x_list   = [var.eps_x_list; var.eps(1,1)];
var.eps_y_list   = [var.eps_y_list; var.eps(2,2)];
var.sigma_x_list = [var.sigma_x_list; var.sigma(1,1)];
var.sigma_y_list = [var.sigma_y_list; var.sigma(2,2)];
var.sigmas_x_list = [var.sigmas_x_list; var.sigma_s(1,1)];
var.sigmas_y_list = [var.sigmas_y_list; var.sigma_s(2,2)];
var.sigmab_x_list = [var.sigmab_x_list; var.sigma_b(1,1)];
var.sigmab_y_list = [var.sigmab_y_list; var.sigma_b(2,2)];
var.epsp_x_list  = [var.epsp_x_list; var.eps_p(1,1)];
var.epsp_y_list  = [var.epsp_y_list; var.eps_p(2,2)];
var.epse_x_list  = [var.epse_x_list; var.eps_e(1,1)];
var.epse_y_list  = [var.epse_y_list; var.eps_e(2,2)];
var.B_list = [var.B_list var.B];

%
meansts = 1/3 * ( var.sigma(1,1) + var.sigma(2,2) + var.sigma(3,3) );

sigma_d_11 = var.sigma(1,1) - meansts;
sigma_d_22 = var.sigma(2,2) - meansts;
sigma_d_33 = var.sigma(3,3) - meansts;

devsts2ndI = sigma_d_11 * sigma_d_22 + sigma_d_22 * sigma_d_33 + sigma_d_33 * sigma_d_11;

%
var.meansts_list = [ var.meansts_list meansts];
var.devsts2ndInvariant = [ var.devsts2ndInvariant devsts2ndI ];
end