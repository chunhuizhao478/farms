clear all; close all; clc;

%use struct_param_sw
struct_param_sw;

%plot alpha-xi
figure();
xi = linspace(-sqrt(2),sqrt(2),1000);
plot(xi,alpha_out1_func(xi));
hold on;
plot(xi,alpha_out2_func(xi));
ylim([0,1]);
xlim([-sqrt(2),sqrt(2)]);

%% Check correctness of forcing function implementation
%Given xi_given, alpha_given which makes it into granular state
xi_given = 0.83676;
I2_given = 0.000839612;
alpha_given = 0.872122;
B_given = 0;
dt = 0.0025;
gamma_given = 3.7150e10;
lambda_given = 3.204e10;
mu_given = 2.32e9;
%
alpha_cr = comp_alpha_cr(xi_given,param);
Prob = 1 / ( exp( (alpha_cr - alpha_given) / param.beta ) + 1 );
%B forcing func
Bforcingfunc = param.C_B * Prob * ( 1 - B_given ) * ( I2_given * (mu_given - gamma_given * xi_given + lambda_given/2 * xi_given ^ 2) - I2_given * (param.a0 + param.a1 * xi_given + param.a2 * xi_given ^ 2 + param.a3 * xi_given ^ 3));
%B_new
B_new = dt * Bforcingfunc;

%% Check energy curve vs xi
xi_vals = linspace(-sqrt(2),sqrt(2),1000);
Esolid = mu_given - gamma_given .* xi_given + lambda_given/2 * xi_given .^ 2;

Egranular = param.a0 + param.a1 .* xi_vals + param.a2 * xi_vals .^ 2 + param.a3 * xi_vals .^ 3;
Egranular_2 = 7.9677e+09 + (-2.2849e+10) .* xi_vals + 2.0269e+10 * xi_vals .^ 2 + (-5.1873e+09) * xi_vals .^ 3;

Esolid_bc = param.chi * (param.mu_0 - 0.0 .* xi_vals + param.lambda_0/2 * xi_vals .^ 2);
figure();
plot(xi_given,Esolid,'r*')
hold on;
plot(xi_vals,Egranular_2)
hold on;
plot(xi_vals,Esolid_bc)
hold on;
xline(xi_given)
legend("Esolid","Egranular","Esolidbc")
title("Energy Curve vs xi range")
ylabel("Energy")
xlabel("xi")

%% Check Esolid

xi_vals = linspace(-sqrt(2),sqrt(2),1000);
