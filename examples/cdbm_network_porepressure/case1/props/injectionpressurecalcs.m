clear all; close all; clc;
%injection pressure
q = 10; %kg/s
rho_o = 10^3; %kg/m^3
c = 6.8e-2; %m^2/s
mu = 20e9; %Pa
nu = 0.25; %-
nu_u = 0.3; %-
k = 3e-16; %m^2
eta = 0.4e-3; %Pa s
%
lambda = 2 * mu * nu / ( 1 - 2 * nu );
lambda_u = 2 * mu * nu_u / ( 1 - 2 * nu_u );
%
x_min = -10^2;
x_max = 10^2;
y_min = -10^2;
y_max = 10^2;
n_points = 1000;
% Define the range for x and y
x = linspace(x_min, x_max, n_points);
y = linspace(y_min, y_max, n_points);

% Create the grid
[X, Y] = meshgrid(x, y);

% injection location
x_inject = 0;
y_inject = 0;

R = sqrt((X-x_inject).^2+(Y-y_inject).^2);

% pressure
% time
t = 1*24*60*60;

xi = R ./ sqrt(c .* t);

p = q .* eta ./ ( 4 .* pi .* rho_o .* R .* k ) .* erfc( 0.5 .* xi );

%surf
surf(X,Y,p); hold on;
shading interp;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('pressure');
zlim([0,20e6]);
title('3D Surface Plot of Pressure Analytical Solution For Fluid Injection');
colorbar;