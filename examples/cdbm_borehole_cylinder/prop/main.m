clear all; close all; clc;
%% Input Parameters 
%%
%Borehole Breakouts Induced in Arkosic Sandstones and a Discrete Element Analysis
%Tablerock sandstone
% E = rho * Vs ^ 2 * ( 3 * Vp ^ 2 - 4 * Vs ^ 2 ) / ( Vp ^ 2 - Vs ^ 2 ) / 1e9;     %GPa
% nu = 0.5 * ( Vp ^ 2 - 2 * Vs ^ 2 ) / ( Vp ^ 2 - Vs ^ 2 );  %-
E = 67; %GPa
nu = 0.25; %-
rho = 2640; %kg/m^3
phi = 31;   %internal friction angle
dx = 0.001; %m
cd_constant = 10; %/s damage accumulation rate
%% Calculate Parameters
%calculate lambda, mu
% lambda = E * nu / ( ( 1 + nu )  * ( 1 - 2 * nu ));
% mu = E / ( 2  * ( 1 + nu ));
lambda = 32.04e9;
mu = 32.04e9;
%calculate 3d xi_o
xi_o = xiocalc3d(lambda,mu,phi);
%calculate shear wave speed
cs = sqrt(mu*1e9/rho);
%calculate pressure wave speed
cd = sqrt((lambda+2*mu)*1e9/rho);
%calculate other parameters
[param] = struct_param(lambda*1e9,mu*1e9,xi_o);
%% Output Parameters
fprintf("lambda (GPa): %.3f \n", lambda);
fprintf("mu (GPa): %.3f \n", mu);
fprintf("youngE (GPa): %.3f \n", E);
fprintf("nu (-): %.3f \n", nu);
fprintf("rho (kg/m^3): %.3f \n", rho);
fprintf("cs (m/s): %.3f \n", cs);
fprintf("cd (m/s): %.3f \n", cd);
fprintf("min time step (ns): %.3f \n", 0.1*dx/cd*1e9);
fprintf("xi_o (-): %.3f \n", xi_o);
fprintf("xi_1 (-): %.3f \n", param.xi_1);
fprintf("gamma_r (GPa): %.3f \n", param.gamma_r/1e9);
fprintf("a0 (GPa): %.3f \n", param.a0/1e9);
fprintf("a1 (GPa): %.3f \n", param.a1/1e9);
fprintf("a2 (GPa): %.3f \n", param.a2/1e9);
fprintf("a3 (GPa): %.3f \n", param.a3/1e9);
param.alpha_out1_func
param.alpha_out2_func