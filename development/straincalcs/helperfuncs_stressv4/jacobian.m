%%Compute Jacobian of ComputeStrain function%%
clear all; close all; clc;
syms eps11p_pre eps22p_pre eps33p_pre eps12p_pre eps11p_inc eps22p_inc eps33p_inc eps12p_inc eps11t eps22t eps33t eps12t lambda gamma_damaged shear_modulus xi a0 a1 a2 a3 B m1 m2 C_g dt nu 
%define viscoelastic strain 
eps11p = eps11p_pre + eps11p_inc;
eps22p = eps22p_pre + eps22p_inc;
eps33p = eps33p_pre + eps33p_inc;
eps12p = eps12p_pre + eps12p_inc;
%define elastic strain
eps11e = eps11t - eps11p;
eps22e = eps22t - eps22p;
eps33e = eps33t - eps33p;
eps12e = eps12t - eps12p;
%represent I1 I2 xi using elastic strain
I1 = eps11e + eps22e + eps33e;
I2 = eps11e ^ 2 + eps22e ^ 2 + 2 * eps12e ^ 2 + eps33e ^ 2;
xi = I1 / sqrt(I2);
%Represent sigma (solid(s) + granular(b))
sigma11_s = ( lambda - gamma_damaged / xi ) * I1 ...
          + ( 2 * shear_modulus - gamma_damaged * xi ) * eps11e;
sigma11_b = ( 2 * a2 + a1 / xi + 3 * a3 * xi ) * I1 ...
          + ( 2 * a0 + a1 * xi - a3 * xi ^ 3 ) * eps11e;
sigma22_s = ( lambda - gamma_damaged / xi ) * I1 ...
          + ( 2 * shear_modulus - gamma_damaged * xi ) * eps22e;
sigma22_b = ( 2 * a2 + a1 / xi + 3 * a3 * xi ) * I1 ...
          + ( 2 * a0 + a1 * xi - a3 * xi ^ 3 ) * eps22e;
sigma33_s = ( lambda - gamma_damaged / xi ) * I1 ...
          + ( 2 * shear_modulus - gamma_damaged * xi ) * eps33e;
sigma33_b = ( 2 * a2 + a1 / xi + 3 * a3 * xi ) * I1 ...
          + ( 2 * a0 + a1 * xi - a3 * xi ^ 3 ) * eps33e;
sigma12_s = ( 2 * shear_modulus - gamma_damaged * xi ) * eps12e;
sigma12_b = ( 2 * a0 + a1 * xi - a3 * xi ^ 3 ) * eps12e;
%Represent total stress
sigma11_t = (1 - B) * sigma11_s + B * sigma11_b;
sigma22_t = (1 - B) * sigma22_s + B * sigma22_b;
sigma33_t = (1 - B) * sigma33_s + B * sigma33_b;
sigma12_t = (1 - B) * sigma12_s + B * sigma12_b;
%Represent deviatroic stress
sigma_d11 = sigma11_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
sigma_d22 = sigma22_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
sigma_d33 = sigma22_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
sigma_d12 = sigma12_t;
%Setup three equations, three unknowns
eqn1 = eps11p_inc - dt * C_g * power(B,m1) * power(sigma_d11,m2);
eqn2 = eps22p_inc - dt * C_g * power(B,m1) * power(sigma_d22,m2);
eqn3 = eps12p_inc - dt * C_g * power(B,m1) * power(sigma_d12,m2);
eqn4 = eps33p_inc - dt * C_g * power(B,m1) * power(sigma_d33,m2);
%jacobian
%
eqout11 = string(diff(eqn1,eps11p_inc));
eqout12 = string(diff(eqn1,eps22p_inc));
eqout13 = string(diff(eqn1,eps12p_inc));
eqout14 = string(diff(eqn1,eps33p_inc));
%
eqout21 = string(diff(eqn2,eps11p_inc));
eqout22 = string(diff(eqn2,eps22p_inc));
eqout23 = string(diff(eqn2,eps12p_inc));
eqout24 = string(diff(eqn2,eps33p_inc));
%
eqout31 = string(diff(eqn3,eps11p_inc));
eqout32 = string(diff(eqn3,eps22p_inc));
eqout33 = string(diff(eqn3,eps12p_inc));
eqout34 = string(diff(eqn3,eps33p_inc));
%
eqout41 = string(diff(eqn4,eps11p_inc));
eqout42 = string(diff(eqn4,eps22p_inc));
eqout43 = string(diff(eqn4,eps12p_inc));
eqout44 = string(diff(eqn4,eps33p_inc));
%
writematrix(eqout11,'./out_jacobian/eqout11.txt');
writematrix(eqout12,'./out_jacobian/eqout12.txt');
writematrix(eqout13,'./out_jacobian/eqout13.txt');
writematrix(eqout14,'./out_jacobian/eqout14.txt');
%
writematrix(eqout21,'./out_jacobian/eqout21.txt');
writematrix(eqout22,'./out_jacobian/eqout22.txt');
writematrix(eqout23,'./out_jacobian/eqout23.txt');
writematrix(eqout24,'./out_jacobian/eqout24.txt');
%
writematrix(eqout31,'./out_jacobian/eqout31.txt');
writematrix(eqout32,'./out_jacobian/eqout32.txt');
writematrix(eqout33,'./out_jacobian/eqout33.txt');
writematrix(eqout34,'./out_jacobian/eqout34.txt');
%
writematrix(eqout41,'./out_jacobian/eqout41.txt');
writematrix(eqout42,'./out_jacobian/eqout42.txt');
writematrix(eqout43,'./out_jacobian/eqout43.txt');
writematrix(eqout44,'./out_jacobian/eqout44.txt');