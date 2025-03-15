%% import parameter %%
[param] = struct_param();
%% import variable %%
[var] = struct_var();
%% initialize material variable %%
[var] = comp_modulus(param,var);
[var] = comp_young_nu(param,var,"init");

alpha_val = 0.8;
eqn = @(xi) ( (var.mu_ - alpha_val * 0.8 * param.gamma_r) - param.a0) - (param.a1 + alpha_val*(param.gamma_r) ) * xi + (var.lambda_/2 - param.a2) * xi .^ 2 - param.a3 * xi .^ 3;
xi_list = -sqrt(3):0.001:sqrt(3);
plot(xi_list,eqn(xi_list)); hold on;
yline(0)

fprintf("a0: %.5f \n a1: %.5f \n a2: %.5f \n a3: %.5f \n",param.a0,param.a1,param.a2,param.a3)
fprintf("- param.a3: %.5f \n (var.lambda_/2 - param.a2): %.5f \n - (param.a1 + var.gamma_): %.5f \n (var.mu_ - param.a0): %.5f \n",-param.a3,var.lambda_/2 - param.a2,- (param.a1 + var.gamma_),(var.mu_ - param.a0))