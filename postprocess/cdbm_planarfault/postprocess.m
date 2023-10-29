clear all; close all; clc;
%postprocess
%number of steps
num_steps = 27;
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
%plot
figure(1);
for i = 1 : num_steps
    plot(arr_xcoord./10^3,arr_sliprate(:,i),'r-'); hold on;
end
xlabel("Distance along fault, x(km)", "FontSize", 15)
ylabel("Slip rate (m/s)", "FontSize", 15)
xlim([-10,10])
ylim([0,55])
ax = gca;
ax.FontSize = 15; 
%% slip
%loop over data files
for i = 1 : num_steps
    %read file
    data_i = readmatrix("outputs/tangent_jump/tangent_jump_"+string(i)+".txt");
    %save
    arr_slip(:,i) = data_i;
end
%plot
figure(2);
for i = 1 : num_steps
    plot(arr_xcoord./10^3,arr_slip(:,i),'r-'); hold on;
end
xlabel("Distance along fault, x(km)", "FontSize", 15)
ylabel("Slip (m)", "FontSize", 15)
xlim([-10,10])
ylim([0,15])
ax = gca;
ax.FontSize = 15; 