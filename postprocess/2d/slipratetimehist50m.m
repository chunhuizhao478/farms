clear all; close all; clc;

%uguca results
datauguca = load("./TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx2880_s2.00_tf0.35_npc1-DataFiles/vel0.txt");
ptrs = [1240 1441 1641 1338 1541 1138 1741 ]; %25m
time_max = 4;
num_max = 40;
meshsize = "50m";
lenged1 = "uguca-2d-25m";
legend2 = "moose-2d-50m";
filename = "50m_922_uguca_setup";
datauguca_time = linspace(0,time_max,num_max);
datauguca_ptr1 = 0.5*(datauguca(:,1240)+datauguca(:,1241));
datauguca_ptr2 = 0.5*(datauguca(:,1450)+datauguca(:,1451));
datauguca_ptr3 = 0.5*(datauguca(:,1641)+datauguca(:,1642));
datauguca_ptr4 = 0.5*(datauguca(:,1340)+datauguca(:,1341));
datauguca_ptr5 = 0.5*(datauguca(:,1540)+datauguca(:,1541));
datauguca_ptr6 = 0.5*(datauguca(:,1140)+datauguca(:,1141));
datauguca_ptr7 = 0.5*(datauguca(:,1741)+datauguca(:,1742));

datauguca_ptr8 = 0.5*(datauguca(:,1198)+datauguca(:,1199)); %6025
datauguca_ptr9 = 0.5*(datauguca(:,1681)+datauguca(:,1682));

% datauguca_ptr8 = 0.5*(datauguca(:,1259)+datauguca(:,1260)); %4525
% datauguca_ptr9 = 0.5*(datauguca(:,1621)+datauguca(:,1622));

% datauguca_ptr8 = 0.5*(datauguca(:,1247)+datauguca(:,1248)); %4825
% datauguca_ptr9 = 0.5*(datauguca(:,1634)+datauguca(:,1635));

% datauguca_ptr8 = 0.5*(datauguca(:,1245)+datauguca(:,1246)); %4875
% datauguca_ptr9 = 0.5*(datauguca(:,1636)+datauguca(:,1637));

% datauguca_ptr8 = 0.5*(datauguca(:,1243)+datauguca(:,1244)); %4925
% datauguca_ptr9 = 0.5*(datauguca(:,1637)+datauguca(:,1638));

% datauguca_ptr8 = 0.5*(datauguca(:,1241)+datauguca(:,1242)); %4975
% datauguca_ptr9 = 0.5*(datauguca(:,1639)+datauguca(:,1640));

%moose results
% ptr1_loc = 101 - 1 #-5025
% ptr2_loc = 202 - 1 #25
% ptr3_loc = 302 - 1 #5025
% ptr4_loc = 151 - 1 #-2525
% ptr5_loc = 252 - 1 #2525
% ptr6_loc = 51  - 1 #-7525
% ptr7_loc = 352 - 1 #7525
% ptr8_loc = 103 - 1 #-4925
% ptr9_loc = 300 - 1 #4925
datatime = load("./files/"+filename+"/time.txt");
datasliprate_ptr1 = load("./files/"+filename+"/sliprate_neg5kmstrike0dot75dip.txt");
datasliprate_ptr2 = load("./files/"+filename+"/sliprate_0strike0dot75dip.txt");
datasliprate_ptr3 = load("./files/"+filename+"/sliprate_pos5kmstrike0dot75dip.txt");
datasliprate_ptr4 = load("./files/"+filename+"/sliprate_neg2dot5kmstrike0dot75dip.txt");
datasliprate_ptr5 = load("./files/"+filename+"/sliprate_pos2dot5kmstrike0dot75dip.txt");
datasliprate_ptr6 = load("./files/"+filename+"/sliprate_neg7dot5kmstrike0dot75dip.txt");
datasliprate_ptr7 = load("./files/"+filename+"/sliprate_pos7dot5kmstrike0dot75dip.txt");

datasliprate_ptr8 = load("./files/"+filename+"/sliprate_neg6025mstrike0dot75dip.txt");
datasliprate_ptr9 = load("./files/"+filename+"/sliprate_pos6025mstrike0dot75dip.txt");


%figures
%%sliprate
figure(1);
plot(datauguca_time(1:num_max),datauguca_ptr1(1:num_max)*2,'r.-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr1(1:4:end),'b.-')
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
plot(datauguca_time(1:num_max),datauguca_ptr2(1:num_max)*2,'r.-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr2(1:4:end),'b.-')
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
plot(datauguca_time(1:num_max),datauguca_ptr3(1:num_max)*2,'r.-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr3(1:4:end),'b.-')
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
plot(datauguca_time(1:num_max),datauguca_ptr4(1:num_max)*2,'r.-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr4(1:4:end),'b.-')
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
plot(datauguca_time(1:num_max),datauguca_ptr5(1:num_max)*2,'r.-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr5(1:4:end),'b.-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at 2.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr5_'+meshsize+'.png')

figure(6);
plot(datauguca_time(1:num_max),datauguca_ptr6(1:num_max)*2,'r.-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr6(1:4:end),'b.-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at -7.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr6_'+meshsize+'.png')

figure(7);
plot(datauguca_time(1:num_max),datauguca_ptr7(1:num_max)*2,'r.-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr7(1:4:end),'b.-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at 7.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr7_'+meshsize+'.png')

% figure(8);
% plot(datauguca_time(1:num_max),datauguca_ptr8(1:num_max)*2,'r.-'); hold on;
% plot(datatime(1:4:end),datasliprate_ptr8(1:4:end),'b.-')
% xlim([0,time_max])
% ylim([0,12.0])
% title("Slip Rate Time History at -4.5km","FontSize",30)
% xlabel("time (s)","FontSize",25)
% ylabel("slip rate (m/s)","FontSize",25)
% legend(lenged1,legend2,"FontSize",15,'Location','northwest');
% ax = gca;
% ax.FontSize = 15; 
% saveas(gcf,'./plots/sliprate_ptr8_'+meshsize+'.png')
% 
% figure(9);
% plot(datauguca_time(1:num_max),datauguca_ptr9(1:num_max)*2,'r.-'); hold on;
% plot(datatime(1:4:end),datasliprate_ptr9(1:4:end),'b.-')
% xlim([0,time_max])
% ylim([0,12.0])
% title("Slip Rate Time History at 4.5km","FontSize",30)
% xlabel("time (s)","FontSize",25)
% ylabel("slip rate (m/s)","FontSize",25)
% legend(lenged1,legend2,"FontSize",15,'Location','northwest');
% ax = gca;
% ax.FontSize = 15; 
% saveas(gcf,'./plots/sliprate_ptr9_'+meshsize+'.png')