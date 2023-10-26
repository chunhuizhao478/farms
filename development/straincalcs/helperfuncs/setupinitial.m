clear all; close all; clc;
%setup initial I1 I2 xi
%%given shear modulus and lambda
mu_o = 32.04e9;
lambda_o = 32.04e9;
%%given stress components
sts11 = -240e6;
sts12 = 81.6e6;
sts22 = -120e6;
%%compute young's modulus and poisson ratio
youngs_modulus_o = mu_o * ( 3 * lambda_o + 2 * mu_o ) / ( lambda_o + mu_o );
nu_o = lambda_o / ( 2 * ( lambda_o + mu_o ) );
%%compute strain components
eps11 = 1 / youngs_modulus_o * ( sts11 - nu_o * sts22 );
eps22 = 1 / youngs_modulus_o * ( sts22 - nu_o * sts11 );
eps12 = 1 / mu_o * sts12;
%%compute I1 I2 xi
I1 = eps11 + eps22;
I2 = eps11 ^ 2 + eps22 ^ 2 + 2 * eps12 ^ 2;
xi = I1 / sqrt(I2);
%%print results
fprintf("xi(computed): %.4f, xi_0: -0.8 \n",xi)
%%compute f for damage evolution
struct_param_sw;
B = 0;
alpha = 0;
if xi < param.xi_0
    %xi<xi_o
    f = ( 1 - B ) * ( param.C_1 * exp( alpha / param.C_2 ) * I2 * ( xi - param.xi_0 ) );
else
    %xi>xi_o
    f = ( 1 - B ) * ( param.C_d * I2 * ( xi - param.xi_0  ) );
end
%%print results
fprintf("f(computed): %.4f \n",f)
