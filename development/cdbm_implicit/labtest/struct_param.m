function [param] = struct_param()
%1d_continuum_damage_breakage_model%
%%material_parameter%%
%Define all parameters in "param" struct
%% Prescribed values %%
param.P_c = 50 * 10 ^ 6;   %Pa  <confinement pressure> refer to "Lyak_BZ_JMPS14_splitstrain" Section 3.1
param.lambda_0 = 30e9; %Pa  <initial 1st lame constant> refer to "Lyak_BZ_JMPS14_splitstrain" Section 3.1 | 10 ^ 10 | 3.204e10
param.mu_0 = 30e9;     %Pa  <initial 2nd lame constant> refer to "Lyak_BZ_JMPS14_splitstrain" Section 3.1 | 10 ^ 10 | 3.204e10
param.xi_0 = -0.8;        %--  <strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
param.xi_d = -0.9;        %--  <strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
param.xi_max = sqrt(3);    %--  <strain invariants ratio: maximum allowable value>: set boundary
param.xi_min = -sqrt(3);   %--  <strain invariants ratio: minimum allowable value>: set boundary
param.chi = 1;          %--  <ratio of two energy state: F_b/F_s = chi < 1>: ensure the energy transition from solid state to granular state.
param.C_d = 10;          %1/s <coefficient gives positive damage   evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
param.C_1 = 300;           %1/s <coefficient of healing for damage evolution>: refer to "Lyakhovsky_2011_Hessian_Matrix" Section 3.4 %%1e-14 %%300
param.C_2 = 0.05;          %1/s <coefficient of healing for damage evolution>: refer to "Lyakhovsky_2011_Hessian_Matrix" Section 3.4 %%0.1 %%0.05
param.beta = 0.03;         %--  <coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1 %%0.001%%
param.C_B = 100 * param.C_d; %1e-5     %1/s <coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
param.C_BH = 1e4;  %1/s <coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14".
param.C_g  = 1e-8;         %1/(Pa*s) <material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
param.m1  = 10;            %--  <coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
param.m2  = 1;             %--  <coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
param.xi_given = 0.8;      %location where ensure the Fb < Fs
%% Compute gamma_r, xi_1 %%
[out_gamma_r,out_xi_1,alpha_out1_func,alpha_out2_func] = hessian(param);
param.gamma_r = double(out_gamma_r);
param.xi_1 = double(out_xi_1);
param.alpha_out1_func = alpha_out1_func;
param.alpha_out2_func = alpha_out2_func;
%% Compute a0 a1 a2 a3 %%
[a0_out,a1_out,a2_out,a3_out] = fb_coeff(param);
param.a0 = double(a0_out); %--  <coefficient of F_B energy>: refer to "fb_coeff" and "hessian"
param.a1 = double(a1_out); %--  <coefficient of F_B energy>: refer to "fb_coeff" and "hessian"
param.a2 = double(a2_out); %--  <coefficient of F_B energy>: refer to "fb_coeff" and "hessian"
param.a3 = double(a3_out); %--  <coefficient of F_B energy>: refer to "fb_coeff" and "hessian"
%%%%Reference%%%%
%"Lyakhovsky_Ben-Zion_P14": A Continuum Damage–Breakage Faulting Model and Solid-Granular Transitions VLADIMIR LYAKHOVSKY 1 and YEHUDA BEN-ZION 2
%"Lyak_BZ_JMPS14_splitstrain": Damage–breakage rheology model and solid-granular transition near brittle instability Vladimir Lyakhovskya,n, Yehuda Ben-Zionb
%"note_mar25": in 'Self Note' folder 

%% Plot the phase diagram
xi = -sqrt(3):0.01:sqrt(3);
figure();
plot(xi,param.alpha_out2_func(xi),'b-'); hold on;
plot(xi,param.alpha_out1_func(xi),'r-'); hold on;
xline(-0.8,'k-.'); hold on;
xline(-sqrt(3),'k-.'); hold on;
xline(sqrt(3),'k-.'); hold on;
xline(0.8248,'k-.'); hold on;
legend("Convexity loss : condition 1", "Convexity loss : condition 2","","","","",'Fontsize',15);
title("Phase Diagram : \fontsize{15}{0}\selectfont$\xi$ vs \fontsize{15}{0}\selectfont$\alpha$ ",Interpreter="latex");
ylabel("Damage Variable \fontsize{15}{0}\selectfont$\alpha$",Interpreter="latex",FontSize=25);
xlabel("Strain Invariant Ratio \fontsize{15}{0}\selectfont$\xi$",Interpreter="latex",FontSize=25);
ylim([0,1])
xlim([-1.8,1.8])
ax = gca;
ax.FontSize = 15;
%% Plot the energy diagram
xi = -sqrt(3):0.01:sqrt(3);
figure();
fb = param.a0 + param.a1 .* (xi) + param.a2 .* (xi) .^ 2 + param.a3 .* (xi) .^ 3;
fs = (param.mu_0 + 0.7 .* param.xi_0 .* param.gamma_r) - (0.7 .* param.gamma_r) .* (xi) + param.lambda_0/2 .* (xi) .^ 2;
plot(xi,fb,'b-'); hold on;
plot(xi,fs,'r-'); hold on;
plot(xi,0.2*fs+0.8*fb,'k-'); hold on;
legend("fb","fs","0.2*fs+0.8*fb",'Fontsize',15);
xlim([-1.8,1.8])
ax = gca;
ax.FontSize = 15;
xi_choice = 1.3;
fb_choice = param.a0 + param.a1 .* (xi_choice) + param.a2 .* (xi_choice) .^ 2 + param.a3 .* (xi_choice) .^ 3;
fs_choice = (param.mu_0 + 0.7 .* param.xi_0 .* param.gamma_r) - (0.7 .* param.gamma_r) .* (xi_choice) + param.lambda_0/2 .* (xi_choice) .^ 2;
ratio = fb_choice / (0.2 * fs_choice + 0.8 * fb_choice);
disp(ratio)
end