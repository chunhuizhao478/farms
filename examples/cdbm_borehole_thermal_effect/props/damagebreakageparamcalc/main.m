clear all; close all; clc;
%% Input Parameters 
%% Temp 150 degree
% Vs = 2435.88; %m/s
% Vp = 3981; %m/s
% rho = 2641.7; %kg/m^3
%% Temp 200 degree
Vs = 2413; %m/s
Vp = 3981; %m/s
rho = 2633.4; %kg/m^3
%% Temp 300 degree
% Vs = 2355.45; %m/s
% Vp = 3881; %m/s
% rho = 2625.15; %kg/m^3
%%
%
E = rho * Vs ^ 2 * ( 3 * Vp ^ 2 - 4 * Vs ^ 2 ) / ( Vp ^ 2 - Vs ^ 2 ) / 1e9;     %GPa
nu = 0.5 * ( Vp ^ 2 - 2 * Vs ^ 2 ) / ( Vp ^ 2 - Vs ^ 2 );  %-
phi = 46;   %internal friction angle
dx = 0.001; %m
cd_constant = 5e7; %/s damage accumulation rate
%% Calculate Parameters
%calculate lambda, mu
lambda = E * nu / ( ( 1 + nu )  * ( 1 - 2 * nu ));
mu = E / ( 2  * ( 1 + nu ));
% lambda = 15.6177;
% mu = 19.877;
%calculate 2d xi_o
xi_o = xiocalc2d(lambda,mu,phi);
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