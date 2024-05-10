clear all; close all; clc;
load("vik.mat");
%postprocess
%number of steps
num_steps = 60;
data_0 = readmatrix("outputs/tangent_jump_rate/tangent_jump_rate_0.txt");
%save data
arr_sliprate = zeros(size(data_0,1),num_steps);
arr_slip = zeros(size(data_0,1),num_steps);
%save data refinement
arr_sliprate_refine = zeros(size(data_0,1),num_steps);
arr_slip_refine = zeros(size(data_0,1),num_steps);
%save traction
arr_traction = zeros(size(data_0,1),num_steps);
% xcoord
arr_xcoord = readmatrix("outputs/tangent_jump_rate/xcoord.txt");
%% traction x
for i = 1 : num_steps
    %read file
    data_i = readmatrix("outputs/traction_x/traction_x_"+string(i)+".txt");
    %save
    arr_traction(:,i) = data_i;
end
figure();
for i = 21 : num_steps
    plot(arr_xcoord./10^3,arr_traction(:,i),'r-','LineWidth',1.5); hold on;
end
ylim([10,130])
%% traction y
for i = 1 : num_steps
    %read file
    data_i = readmatrix("outputs/traction_y/traction_y_"+string(i)+".txt");
    %save
    arr_traction(:,i) = data_i;
end
figure();
for i = 21
    plot(arr_xcoord./10^3,arr_traction(:,i),'r-','LineWidth',1.5); hold on;
end
ylim([10,130])
%% sliprate
%loop over data files
for i = 1 : num_steps
    %read file
    data_i = readmatrix("outputs/tangent_jump_rate/tangent_jump_rate_"+string(i)+".txt");
    %save
    arr_sliprate(:,i) = data_i;
    %read file fine
    data_i_refine = readmatrix("../cdbm_planarfault/outputs/tangent_jump_rate/tangent_jump_rate_"+string(i)+".txt");
    %save fine
    arr_sliprate_refine(:,i) = data_i_refine;
end
%add time
tplot=0:0.1:2;
colormap(autumn);
cm=colormap(autumn(21));
%plot
figure(1);
for i = 1 : num_steps
    plot(arr_xcoord./10^3,arr_sliprate(:,i),'color',cm(i,:)); hold on;
end
% for i = 11
%     plot(arr_xcoord./10^3,arr_sliprate(:,i),'r-','LineWidth',1.5); hold on;
% %     plot(arr_xcoord./10^3,arr_sliprate_refine(:,i),'r-','LineWidth',1.5); hold on;
% end
xlabel("Distance along fault, x(km)",'FontSize',15)
ylabel("Slip rate (m/s)",'FontSize',15)
% legend("Refine","Uniform")
xlim([-10,10])
% ylim([0,50])

ax = gca;
ax.FontSize = 15; 

% cbh=colorbar;
% cbh.Ticks =tplot/tplot(end) ;
% cbh.TickLabels = num2cell(cbh.Ticks*tplot(end));
colormap(turbo)
saveas(gcf, './outputs/pngs/SR');
saveas(gcf, './outputs/pngs/SR.png');
%% slip
%loop over data files
for i = 1 : num_steps
    %read file
    data_i = readmatrix("outputs/tangent_jump/tangent_jump_"+string(i)+".txt");
    %save
    arr_slip(:,i) = data_i;
    %read file
    data_i_refine = readmatrix("../cdbm_planarfault/outputs/tangent_jump/tangent_jump_"+string(i)+".txt");
    %save
    arr_slip_refine(:,i) = data_i_refine;
end
%add time
tplot=0:0.1:2;
% colormap(turbo);
% cm=colormap(turbo(length(tplot)));
%plot
figure(2);
for i = 1 : num_steps
    plot(arr_xcoord./10^3,arr_slip(:,i),'r-','LineWidth',1.5); hold on;
%     plot(arr_xcoord./10^3,arr_slip_refine(:,i),'b-','LineWidth',1.5); hold on;
end
xlabel("Distance along fault, x(km)",'FontSize',15)
ylabel("Slip (m)",'FontSize',15)
% legend("Refine","Uniform")
xlim([-10,10])
ylim([0,10])

ax = gca;
ax.FontSize = 15; 

% cbh=colorbar;
% cbh.Ticks =tplot/tplot(end) ;
% cbh.TickLabels = num2cell(cbh.Ticks*tplot(end));
% colormap(turbo)
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
xnn = arr_xcoord./10^3;
ynn = t;
[X,Y]= ndgrid (xnn,ynn);
Z = arr_sliprate;
figure;
pcolor(X,Y,abs(Z));
shading interp;
colorbar;
caxis([1e-2 50]);
% set (gca, 'ColorScale', 'log');
xlabel ('Distance along Fault, x (km)','FontSize',15)
ylabel('Time (s)','FontSize',15)
set(gca, 'FontSize', 10)
set(gcf, 'PaperUnits', 'inches');
x_width=4 ;y_width=2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
%% plastic work
data = readmatrix("/Users/andyz/projects/farms/postprocess/cdbm_planarfault/outputs/plastic_work/plastic_work.txt");
time = 0:0.1:2;
figure();
plot(time,data,'-m', 'linewidth', 1);
xlabel('Time (s)')
ylabel('CDBM Plastic Work (N-m)')
set(gca, 'fontsize', 10)
 set(gcf, 'PaperUnits', 'inches');
x_width=4 ;y_width=3;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
xlim([0,2.5])
%% slip rate time hist
