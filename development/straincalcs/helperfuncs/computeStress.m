%if parameters can retrieve the desired stress
clear all; close all; clc;
%setup initial I1 I2 xi
%%given shear modulus and lambda
mu_o = 32.04e9;
lambda_o = 32.04e9;
%%given stress components
sts11 = -240e6;
sts12 = 70e6;
sts22 = -120e6;
%%compute young's modulus and poisson ratio
youngs_modulus_o = mu_o * ( 3 * lambda_o + 2 * mu_o ) / ( lambda_o + mu_o );
nu_o = lambda_o / ( 2 * ( lambda_o + mu_o ) );
%%compute strain components
%%in plane strain, strain_z = 0, don't forget stress_z = nu_o * ( stress_x
%%+ stress_y )
%%please check "plane_strain_plane_stress.pdf" in the reference folder!!!
eps11 = 1 / youngs_modulus_o * ( ( 1 + nu_o ) * sts11 - nu_o * ( sts11 + sts22 + nu_o * ( sts11 + sts22 ) ) );
eps22 = 1 / youngs_modulus_o * ( ( 1 + nu_o ) * sts22 - nu_o * ( sts11 + sts22 + nu_o * ( sts11 + sts22 ) ) );
eps12 = 1 / youngs_modulus_o * ( ( 1 + nu_o ) * sts12 );
%%compute I1 I2 xi
I1 = eps11 + eps22;
I2 = eps11 ^ 2 + eps22 ^ 2 + 2 * eps12 ^ 2;
xi = I1 / sqrt(I2);
%%%%%%%%%%%%%%%%%%%%
sigma_s_11 = lambda_o * I1 + 2 * mu_o * eps11;
sigma_s_22 = lambda_o * I1 + 2 * mu_o * eps22;
sigma_s_12 = 2 * mu_o * eps12;