% clear all; close all; clc;
%1d_continuum_damage_breakage_model%
%%material_parameter%%
%Define all parameters in "param" struct
%% Prescribed values %%
param.P_c = 10 * 10 ^ 6;   %Pa  <confinement pressure> refer to "Lyak_BZ_JMPS14_splitstrain" Section 3.1
param.lambda_0 = 3.204e10;  %Pa  <initial 1st lame constant> refer to "Lyak_BZ_JMPS14_splitstrain" Section 3.1 | 10 ^ 10
param.mu_0 = 3.204e10;      %Pa  <initial 2nd lame constant> refer to "Lyak_BZ_JMPS14_splitstrain" Section 3.1 | 10 ^ 10
param.xi_0 = -0.8;         %--  <strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
param.xi_d = -0.9;         %--  <strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
param.xi_max = sqrt(2);    %--  <strain invariants ratio: maximum allowable value>: set boundary
param.xi_min = -sqrt(2);   %--  <strain invariants ratio: minimum allowable value>: set boundary
param.chi = 0.2;           %--  <ratio of two energy state: F_b/F_s = chi < 1>: ensure the energy transition from solid state to granular state.
param.C_d = 5000;           %1/s <coefficient gives positive damage   evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
param.C_1 = 300;           %1/s <coefficient of healing for damage evolution>: refer to "Lyakhovsky_2011_Hessian_Matrix" Section 3.4 %%1e-14 %%300
param.C_2 = 0.05;          %1/s <coefficient of healing for damage evolution>: refer to "Lyakhovsky_2011_Hessian_Matrix" Section 3.4 %%0.1 %%0.05
param.beta = 5e-3;         %--  <coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1 %%0.001%%
param.C_B = 30000; %1e-5   %1/s <coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
param.C_BH = 300000;  %1/s <coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14".
param.C_g  = 1e-5;         %1/(Pa*s) <material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
param.m1  = 10;            %--  <coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
param.m2  = 1;             %--  <coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
%% Compute gamma_r, xi_1 %%
[out_gamma_r,out_xi_1,alpha_out1_func,alpha_out2_func] = hessian(param);
param.gamma_r = double(out_gamma_r);
param.xi_1 = double(out_xi_1);
param.alpha_out1_func = alpha_out1_func;
param.alpha_out2_func = alpha_out2_func;
%% Compute a0 a1 a2 a3 %%
xi_critical = 0.5;
[a0_out,a1_out,a2_out,a3_out] = fb_coeff(param,xi_critical);
param.a0 = double(a0_out); %--  <coefficient of F_B energy>: refer to "fb_coeff" and "hessian"
param.a1 = double(a1_out); %--  <coefficient of F_B energy>: refer to "fb_coeff" and "hessian"
param.a2 = double(a2_out); %--  <coefficient of F_B energy>: refer to "fb_coeff" and "hessian"
param.a3 = double(a3_out); %--  <coefficient of F_B energy>: refer to "fb_coeff" and "hessian"
%%%%Reference%%%%
%"Lyakhovsky_Ben-Zion_P14": A Continuum Damage–Breakage Faulting Model and Solid-Granular Transitions VLADIMIR LYAKHOVSKY 1 and YEHUDA BEN-ZION 2
%"Lyak_BZ_JMPS14_splitstrain": Damage–breakage rheology model and solid-granular transition near brittle instability Vladimir Lyakhovskya,n, Yehuda Ben-Zionb
%"note_mar25": in 'Self Note' folder 