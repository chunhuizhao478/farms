clear all; close all; clc;

%uguca results
datauguca = load("./TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx2880_s2.00_tf0.35_npc1-DataFiles/vel0.txt");
datauguca = datauguca(:,1040:1841);
ptrs = [202 402 602 301 502 101 702 221 582]; %25m
time_max = 4;
num_max = 40;
meshsize = "25m";
lenged1 = "uguca-2d-25m";
legend2 = "moose-2d-25m";
filename = "25m_922_uguca_setup";
datauguca_time = linspace(0,time_max,num_max);
datauguca_ptr1 = datauguca(:,ptrs(1));
datauguca_ptr2 = datauguca(:,ptrs(2));
datauguca_ptr3 = datauguca(:,ptrs(3));
datauguca_ptr4 = datauguca(:,ptrs(4));
datauguca_ptr5 = datauguca(:,ptrs(5));
datauguca_ptr6 = datauguca(:,ptrs(6));
datauguca_ptr7 = datauguca(:,ptrs(7));
datauguca_ptr8 = datauguca(:,ptrs(8));
datauguca_ptr9 = datauguca(:,ptrs(9));

%moose results
% ptr1_loc = 201 - 1 #-5012
% ptr2_loc = 402 - 1 #12.5
% ptr3_loc = 602 - 1 #5012
% ptr4_loc = 301 - 1 #-2512
% ptr5_loc = 502 - 1 #2512
% ptr6_loc = 101 - 1 #-7512
% ptr7_loc = 702 - 1 #7512
% ptr8_loc = 221 - 1 #-4512
% ptr9_loc = 582 - 1 #4512
datatime = load("./files/"+filename+"/time.txt");
datasliprate_ptr1 = load("./files/"+filename+"/sliprate_neg5kmstrike0dot75dip.txt");
datasliprate_ptr2 = load("./files/"+filename+"/sliprate_0strike0dot75dip.txt");
datasliprate_ptr3 = load("./files/"+filename+"/sliprate_pos5kmstrike0dot75dip.txt");
datasliprate_ptr4 = load("./files/"+filename+"/sliprate_neg2dot5kmstrike0dot75dip.txt");
datasliprate_ptr5 = load("./files/"+filename+"/sliprate_pos2dot5kmstrike0dot75dip.txt");
datasliprate_ptr6 = load("./files/"+filename+"/sliprate_neg7dot5kmstrike0dot75dip.txt");
datasliprate_ptr7 = load("./files/"+filename+"/sliprate_pos7dot5kmstrike0dot75dip.txt");
datasliprate_ptr8 = load("./files/"+filename+"/sliprate_neg4dot5kmstrike0dot75dip.txt");
datasliprate_ptr9 = load("./files/"+filename+"/sliprate_pos4dot5kmstrike0dot75dip.txt");

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

figure(6);
plot(datauguca_time(1:num_max),datauguca_ptr6(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr6(1:4:end),'b-')
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
plot(datauguca_time(1:num_max),datauguca_ptr7(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr7(1:4:end),'b-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at 7.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr7_'+meshsize+'.png')

figure(8);
plot(datauguca_time(1:num_max),datauguca_ptr8(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr8(1:4:end),'b-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at -4.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr8_'+meshsize+'.png')

figure(9);
plot(datauguca_time(1:num_max),datauguca_ptr9(1:num_max)*2,'r-'); hold on;
plot(datatime(1:4:end),datasliprate_ptr9(1:4:end),'b-')
xlim([0,time_max])
ylim([0,12.0])
title("Slip Rate Time History at 4.5km","FontSize",30)
xlabel("time (s)","FontSize",25)
ylabel("slip rate (m/s)","FontSize",25)
legend(lenged1,legend2,"FontSize",15,'Location','northwest');
ax = gca;
ax.FontSize = 15; 
saveas(gcf,'./plots/sliprate_ptr9_'+meshsize+'.png')
