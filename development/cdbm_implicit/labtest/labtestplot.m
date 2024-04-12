clear all; close all; clc;
eps = readmatrix("./data_scec_plots/eps.txt");
sts = readmatrix("./data_scec_plots/sigma.txt");
alpha = readmatrix("./data_scec_plots/alpha.txt");
B = readmatrix("./data_scec_plots/B.txt");
xi = readmatrix("./data_scec_plots/xi.txt");

eps_extend = linspace(0.0045,0.008,10);
sts_extend = linspace(81.5143e6,81.5143e6,10);
alpha_extend = linspace(0.705743,0.705743,10);
B_extend = linspace(1.0,1.0,10);

eps_extend = reshape(eps_extend,[10,1]);
sts_extend = reshape(sts_extend,[10,1]);
alpha_extend = reshape(alpha_extend,[10,1]);
B_extend = reshape(B_extend,[1,10]);

eps = cat(1,eps,eps_extend);
sts = cat(1,sts,sts_extend);
alpha = cat(1,alpha,alpha_extend);
B = cat(2,B,B_extend);

subplot(2,1,1)
plot(eps,sts/1e6,'k.-',LineWidth=3); hold on;
xlim([0,8e-3])
ylim([80,420])
legend("stress in loading direction",FontSize=15,FontName='Times New Roman');
ylabel("total stress in loading direction (MPa)",FontSize=25,FontName='Times New Roman');
title("Stress/Damage/Breakage vs Strain Plot",FontSize=30,FontName='Times New Roman');
ax = gca;
ax.FontSize = 15; 

subplot(2,1,2)
plot(eps,alpha,'r-',LineWidth=3); hold on;
plot(eps,B,"b-",LineWidth=3); hold on;
xlim([0,8e-3])
xlim([4.07e-3,4.09e-3])
legend("\alpha","B",FontSize=15,FontName='Times New Roman');
xlabel("total strain",FontSize=25,FontName='Times New Roman');
ylabel("damage (\alpha), breakage variable (B)",FontSize=25,FontName='Times New Roman');
ax = gca;
ax.FontSize = 15; 

% subplot(3,1,3)
% plot(eps,xi,"m-");  hold on;
% xlim([0,1.2e-3])
