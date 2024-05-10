%paper revision, slip rate, traction x time hist plots, May 10 2024
close all; clear all; clc;
num_steps_linearelastic = 60;
time_linearelastic = 0.05:0.05:3;
locs = [320 300 280 260 240]; %2km 3km 4km 5km 6km 7km

sliprate_data_mat_linearelast = zeros(length(locs),num_steps_linearelastic);
traction_data_mat_linearelast = zeros(length(locs),num_steps_linearelastic);

%loop over data files
for i = 1 : num_steps_linearelastic
    %read file
    data_sliprate_i = readmatrix("outputs/tangent_jump_rate/tangent_jump_rate_"+string(i)+".txt");
    data_traction_i = readmatrix("outputs/traction_x/traction_x_"+string(i)+".txt");
    %save
    sliprate_data_mat_linearelast(:,i) = data_sliprate_i(locs);
    traction_data_mat_linearelast(:,i) = data_traction_i(locs);
end

time_linearelastic = time_linearelastic(2:2:54);
sliprate_data_mat_linearelast = sliprate_data_mat_linearelast(:,2:2:54);
traction_data_mat_linearelast = traction_data_mat_linearelast(:,2:2:54);

%
num_steps_damagmodel = 27;
time_damagemodel = 0.1:0.1:2.7;

sliprate_data_mat_damagemodel = zeros(length(locs),num_steps_damagmodel);
traction_data_mat_damagemodel = zeros(length(locs),num_steps_damagmodel);

%loop over data files
for i = 1 : num_steps_damagmodel
    %read file
    data_sliprate_i = readmatrix("../cdbm_planarfault_uniform/outputs/tangent_jump_rate/tangent_jump_rate_"+string(i)+".txt");
    data_traction_i = readmatrix("../cdbm_planarfault_uniform/outputs/traction_x/traction_x_"+string(i)+".txt");
    %save
    sliprate_data_mat_damagemodel(:,i) = data_sliprate_i(locs);
    traction_data_mat_damagemodel(:,i) = data_traction_i(locs);
end

figure(1);
for loc_i = 1 : length(locs)
    plot(time_linearelastic, traction_data_mat_linearelast(loc_i,:),'r-',LineWidth=1.5); hold on;
    plot(time_damagemodel, traction_data_mat_damagemodel(loc_i,:),'k-',LineWidth=1.5); hold on;
end
legend("Linear Elastic", "CDBM")
set(gca, 'YScale', 'log')
xlabel("Time (s)",FontSize=20)
ylabel("Shear stress (MPa)", FontSize=20)
ylim([8,75])
ax = gca;
ax.FontSize = 15; 
yline(81.24,'b-',label="shear strength",FontSize=15,LineStyle='--',LineWidth=1)
yline(12,'b-',label="residual strength",FontSize=15,LineStyle='--',LineWidth=1)
title("Shear stress time history at selected points"," ",FontSize=25)
