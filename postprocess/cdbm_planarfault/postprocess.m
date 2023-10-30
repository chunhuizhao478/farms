clear all; close all; clc;
%postprocess
%number of steps
num_steps = 21;
data_0 = readmatrix("outputs/tangent_jump_rate/tangent_jump_rate_0.txt");
%save data
arr_sliprate = zeros(size(data_0,1),num_steps);
arr_slip = zeros(size(data_0,1),num_steps);
%% xcoord
arr_xcoord = readmatrix("outputs/tangent_jump_rate/xcoord.txt");
%% sliprate
%loop over data files
for i = 1 : num_steps
    %read file
    data_i = readmatrix("outputs/tangent_jump_rate/tangent_jump_rate_"+string(i)+".txt");
    %save
    arr_sliprate(:,i) = data_i;
end
%add time
tplot=0:0.1:2;
colormap(cool);
cm=colormap(cool(length(tplot)));
%plot
figure(1);
for i = 1 : num_steps
    plot(arr_xcoord./10^3,arr_sliprate(:,i),'color',cm(i,:)); hold on;
end
xlabel("Distance along fault, x(km)")
ylabel("Slip rate (m/s)")
xlim([-10,10])
ylim([0,50])

cbh=colorbar;
cbh.Ticks =tplot/tplot(end) ;
cbh.TickLabels = num2cell(cbh.Ticks*tplot(end));

saveas(gcf, './outputs/pngs/SR');
saveas(gcf, './outputs/pngs/SR.png');
%% slip
%loop over data files
for i = 1 : num_steps
    %read file
    data_i = readmatrix("outputs/tangent_jump/tangent_jump_"+string(i)+".txt");
    %save
    arr_slip(:,i) = data_i;
end
%add time
tplot=0:0.1:2;
colormap(cool);
cm=colormap(cool(length(tplot)));
%plot
figure(2);
for i = 1 : num_steps
    plot(arr_xcoord./10^3,arr_slip(:,i),'color', cm(i,:)); hold on;
end
xlabel("Distance along fault, x(km)")
ylabel("Slip (m)")
xlim([-10,10])
ylim([0,10])

cbh=colorbar;
cbh.Ticks =tplot/tplot(end) ;
cbh.TickLabels = num2cell(cbh.Ticks*tplot(end));
colormap(cool)
saveas(gcf, './outputs/pngs/Slip');
saveas(gcf, './outputs/pngs/Slip.png');
%% front position vs time
%loop over data files
data_time_0 = readmatrix("outputs/front_time0.txt");
data_front_0 = readmatrix("outputs/front_xcoord0.txt");
data_time_1 = readmatrix("outputs/front_time1.txt");
data_front_1 = readmatrix("outputs/front_xcoord1.txt");
%
data_time_0 = data_time_0(1:21);
data_front_0 = data_front_0(1:21);
data_time_1 = data_time_1(1:21);
data_front_1 = data_front_1(1:21);
%plot
figure(3);
plot(data_front_0./10^3,data_time_0,'r-'); hold on;
plot(data_front_1./10^3,data_time_1,'r-');
xlabel("Distance along fault, x(km)")
ylabel("Time (s)")
xlim([-10,10]);
ylim([0,2]);
saveas(gcf, './outputs/pngs/SlipRateContour');
saveas(gcf, './outputs/pngs/SlipRateContour.png');
%
t=0:0.1:2;
xnn = 
ynn = t;