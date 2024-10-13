%exec_cycleloading
clear all; close all; clc;
%% create a folder
folderName = 'results_cycleloading/';
% Check if the folder exists
if exist(folderName, 'dir')
    % Remove the folder and its contents
    rmdir(folderName, 's'); % 's' option removes the folder and all subfolders/files
end

% Create a fresh folder
mkdir(folderName);
%% loop
C1_params  = [30 3 0.3];
CBH_params = [300 30 3];
applied_strains = [-0.005 -0.01 -0.02]; %compressive strain is negative
num_cycle_val = 100;
for num_strain_i = 1 : size(applied_strains,2)
    strain_i = applied_strains(num_strain_i);
    for num_C_i = 1 : size(C1_params, 2)
        C1_i = C1_params(num_C_i);
        CBH_i = CBH_params(num_C_i);
        folderName = 'results_cycleloading/';
        addname = 'strain_'+string(strain_i)+'_C1_'+string(C1_i)+'_CBH_'+string(CBH_i)+'/';
        fprintf("strain: %.4f, C1: %.4f, CBH: %.4f \n", strain_i, C1_i, CBH_i) 
        folderName = folderName + addname;
        main_cycleloading;
    end
end
disp("exec cycleloading completed!");