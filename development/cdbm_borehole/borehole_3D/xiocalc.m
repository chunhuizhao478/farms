clear all; close all; clc;
E = 48.5;
nu = 0.22;
rho = 2640;
%calculate lambda, mu
lambda = E * nu / ( ( 1 + nu )  * ( 1 - 2 * nu ));
mu = E / ( 2  * ( 1 + nu ));
%calculate 3d xi_o
phi = 46;
r = 0.27;
xi_o = -sqrt(2) / sqrt(3*(lambda/mu)^2+4*lambda/mu+2-r*(3*lambda/mu+2)^2);
%calculate shear wave speed
cs = sqrt(mu*1e9/rho);
%calculate pressure wave speed
cd = sqrt((lambda+2*mu)*1e9/rho);