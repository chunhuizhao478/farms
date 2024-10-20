clear all; close all; clc;

%uguca results
datauguca = load("./TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx2880_s2.00_tf0.35_npc1-DataFiles/vel0.txt");
% datauguca = load("./TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx1440_s2.00_tf0.35_npc1-DataFiles/vel0.txt");
% datauguca = load("./TPV101_Nx/TPV101_Nx720_s2.00_tf0.35_npc1/TPV101_Nx720_s2.00_tf0.35_npc1-DataFiles/vel0.txt");
% ptrs = [310 361 412 335 387 285 437]; %100m (-5 0 5 -2.5 2.5 -7.5 7.5)
% ptrs = [622 721 820 671 771 570 872]; %50m
ptrs = [1240 1441 1641 1338 1541 1138 1741]; %25m
time_max = 4;
num_max = 40;
meshsize = "50m";
lenged1 = "uguca-2d-25m";
legend2 = "moose-2d-25m";
filename = "50m_921_damp05";
datauguca_time = linspace(0,time_max,num_max);
datauguca_ptr1 = 0.5*(datauguca(:,1240)+datauguca(:,1241));
datauguca_ptr2 = 0.5*(datauguca(:,1440)+datauguca(:,1441));
datauguca_ptr3 = 0.5*(datauguca(:,1641)+datauguca(:,1642));
datauguca_ptr4 = 0.5*(datauguca(:,1340)+datauguca(:,1341));
datauguca_ptr5 = 0.5*(datauguca(:,1540)+datauguca(:,1541));
datauguca_ptr6 = 0.5*(datauguca(:,1140)+datauguca(:,1141));
datauguca_ptr7 = 0.5*(datauguca(:,1740)+datauguca(:,1741));

%moose results
% ptr1_loc = 202 - 1 #-5012
% ptr2_loc = 403 - 1 #12.5
% ptr3_loc = 603 - 1 #5012
% ptr4_loc = 302 - 1 #-2512
% ptr5_loc = 503 - 1 #2512
% ptr6_loc = 102 - 1 #-7512
% ptr7_loc = 703 - 1 #7512
datatime = load("./files/"+filename+"/time.txt");
datasliprate_ptr1 = load("./files/"+filename+"/sliprate_neg5kmstrike0dot75dip.txt");
datasliprate_ptr2 = load("./files/"+filename+"/sliprate_0strike0dot75dip.txt");
datasliprate_ptr3 = load("./files/"+filename+"/sliprate_pos5kmstrike0dot75dip.txt");
datasliprate_ptr4 = load("./files/"+filename+"/sliprate_neg2dot5kmstrike0dot75dip.txt");
datasliprate_ptr5 = load("./files/"+filename+"/sliprate_pos2dot5kmstrike0dot75dip.txt");
datasliprate_ptr6 = load("./files/"+filename+"/sliprate_neg7dot5kmstrike0dot75dip.txt");
datasliprate_ptr7 = load("./files/"+filename+"/sliprate_pos7dot5kmstrike0dot75dip.txt");

%figures
%%sliprate
figure(1);
plot(datauguca_time(1:num_max),datauguca_ptr1(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr1(1:4:end),'b-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at -5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr1_'+meshsize+'.png')

figure(2);
plot(datauguca_time(1:num_max),datauguca_ptr2(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr2(1:4:end),'b-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at 0km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr2_'+meshsize+'.png')

figure(3);
plot(datauguca_time(1:num_max),datauguca_ptr3(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr3(1:4:end),'b-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at 5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr3_'+meshsize+'.png')

figure(4);
plot(datauguca_time(1:num_max),datauguca_ptr4(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr4(1:4:end),'b-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at -2.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr4_'+meshsize+'.png')

figure(5);
plot(datauguca_time(1:num_max),datauguca_ptr5(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr5(1:4:end),'b-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at 2.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr5_'+meshsize+'.png')

% figure(6);
% plot(datauguca_time(1:num_max),datauguca_ptr6(1:num_max)*2,'r-'); hold on;
% plot(datatime(1:4:end),datasliprate_ptr6(1:4:end),'b-')
% xlim([0,time_max])
% ylim([0,12.0])
% title("Slip Rate Time History at -7.5km","FontSize",30)
% xlabel("time (s)","FontSize",25)
% ylabel("slip rate (m/s)","FontSize",25)
% legend(lenged1,legend2,"FontSize",15,'Location','northwest');
% ax = gca;
% ax.FontSize = 15; 
% saveas(gcf,'./plots/sliprate_ptr6_'+meshsize+'.png')
% 
% figure(7);
% plot(datauguca_time(1:num_max),datauguca_ptr7(1:num_max)*2,'r-'); hold on;
% plot(datatime(1:4:end),datasliprate_ptr7(1:4:end),'b-')
% xlim([0,time_max])
% ylim([0,12.0])
% title("Slip Rate Time History at 7.5km","FontSize",30)
% xlabel("time (s)","FontSize",25)
% ylabel("slip rate (m/s)","FontSize",25)
% legend(lenged1,legend2,"FontSize",15,'Location','northwest');
% ax = gca;
% ax.FontSize = 15; 
% saveas(gcf,'./plots/sliprate_ptr7_'+meshsize+'.png')
