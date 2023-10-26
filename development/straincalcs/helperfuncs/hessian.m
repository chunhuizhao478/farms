function[out_gamma_r,xi_1,alpha_out1_func,alpha_out2_func] = hessian(param)
%input: param (struct); 
%output: gamma_r, xi_1, two functions (alpha_cr)
%%%Reproduce the instability plot%%%
%clear all; close all; clc
syms gamma_r alpha_ xi
%% Initialization %%
mu_o = param.mu_0;
lambda_o = mu_o;
xi_o = param.xi_0;
%% Parameters Definition %%
lambda_ = lambda_o;
mu_ = mu_o + alpha_ * xi_o * gamma_r;
gamma_ = alpha_ * gamma_r;
%% Equation 2 %%
expr2 = (2*mu_-gamma_*xi)^2+(2*mu_-gamma_*xi)*(3*lambda_-gamma_*xi)+(lambda_*gamma_*xi-gamma_^2)*(3-xi^2);
%%Solve for gamma_r
solve_gamma_r = expr2;
solve_gamma_r = subs(solve_gamma_r,alpha_,1);
solve_gamma_r = subs(solve_gamma_r,xi,xi_o);
out_gamma_r = solve(solve_gamma_r,gamma_r);
out_gamma_r = out_gamma_r(2);
%%Solve for alpha
solve_alpha2 = expr2;
solve_alpha2 = subs(solve_alpha2,gamma_r,out_gamma_r);
alpha_out2 = solve(solve_alpha2==0,alpha_);
alpha_out2 = alpha_out2(2);
%%Compute alpha_cr(xi=0)%%
alpha_cr_xi_0 = subs(alpha_out2,xi,0);
fprintf("alpha_cr at xi = 0: %d\n",eval(alpha_cr_xi_0))
%% Equation 1 %%
expr1 = 2 * mu_ - gamma_ * xi;
%%Solve for alpha
solve_alpha1 = expr1;
alpha_out1 = solve(solve_alpha1==0,alpha_);
alpha_out1 = subs(alpha_out1,gamma_r,out_gamma_r);
%%Compute xi_1%%
xi_1 = xi_o + sqrt(xi_o^2 + 2 * mu_o / lambda_o);
%% Plot figure %%
xi_list = linspace(xi_o,sqrt(3),100);
alpha_plot1 = subs(alpha_out1,xi,xi_list);
alpha_plot2 = subs(alpha_out2,xi,xi_list);
figure(1);
plot(xi_list,alpha_plot2); hold on
plot(xi_list,alpha_plot1);
xline(xi_o,'m-.',{'xi_o'}) 
xline(xi_1,'m-.',{'xi_1'}) 
xlim([-sqrt(3),sqrt(3)])
ylim([0,1])
title('Strain invariants ratio vs Damage')
xlabel('Strain invariants ratio')
ylabel('Damage')
legend('Eq15a','Eq15b')
%% Store the results
alpha_out1_func = matlabFunction(alpha_out1); %Create function w xi as variable
alpha_out2_func = matlabFunction(alpha_out2);

 