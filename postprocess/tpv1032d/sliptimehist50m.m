clear all; close all; clc;

%uguca results
datauguca = load("./TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx2880_s2.00_tf0.35_npc1-DataFiles/disp0.txt");
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

%moose results
datatime = load("./files/"+filename+"/time.txt");
datasliprate_ptr1 = load("./files/"+filename+"/slip_neg5kmstrike0dot75dip.txt");
datasliprate_ptr2 = load("./files/"+filename+"/slip_0strike0dot75dip.txt");
datasliprate_ptr3 = load("./files/"+filename+"/slip_pos5kmstrike0dot75dip.txt");
datasliprate_ptr4 = load("./files/"+filename+"/slip_neg2dot5kmstrike0dot75dip.txt");
datasliprate_ptr5 = load("./files/"+filename+"/slip_pos2dot5kmstrike0dot75dip.txt");
datasliprate_ptr6 = load("./files/"+filename+"/slip_neg7dot5kmstrike0dot75dip.txt");
datasliprate_ptr7 = load("./files/"+filename+"/slip_pos7dot5kmstrike0dot75dip.txt");

%figures
%%sliprate
figure(1);
plot(datauguca_time(1:num_max),datauguca_ptr1(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr1(1:4:end),'b-')
xlim([0,time_max])
ylim([0,7.0])
title("Slip Time History at -5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip (m)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/slip_ptr1_'+meshsize+'.png')

figure(2);
plot(datauguca_time(1:num_max),datauguca_ptr2(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr2(1:4:end),'b-')
xlim([0,time_max])
ylim([0,7.0])
title("Slip Time History at 0km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip (m)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/slip_ptr2_'+meshsize+'.png')

figure(3);
plot(datauguca_time(1:num_max),datauguca_ptr3(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr3(1:4:end),'b-')
xlim([0,time_max])
ylim([0,7.0])
title("Slip Time History at 5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip (m)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/slip_ptr3_'+meshsize+'.png')

figure(4);
plot(datauguca_time(1:num_max),datauguca_ptr4(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr4(1:4:end),'b-')
xlim([0,time_max])
ylim([0,7.0])
title("Slip Time History at -2.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip (m)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/slip_ptr4_'+meshsize+'.png')

figure(5);
plot(datauguca_time(1:num_max),datauguca_ptr5(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr5(1:4:end),'b-')
xlim([0,time_max])
ylim([0,7.0])
title("Slip Time History at 2.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip (m)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/slip_ptr5_'+meshsize+'.png')

figure(6);
plot(datauguca_time(1:num_max),datauguca_ptr6(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr6(1:4:end),'b-')
xlim([0,time_max])
ylim([0,7.0])
title("Slip Time History at -7.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip (m)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/slip_ptr6_'+meshsize+'.png')

figure(7);
plot(datauguca_time(1:num_max),datauguca_ptr7(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr7(1:4:end),'b-')
xlim([0,time_max])
ylim([0,7.0])
title("Slip Time History at 7.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip (m)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/slip_ptr7_'+meshsize+'.png')


