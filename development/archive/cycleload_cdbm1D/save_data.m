function [] = save_data(var,param,folderName)
    % Check if the folder exists
    if exist(folderName, 'dir')
        % Remove the folder and its contents
        rmdir(folderName, 's'); % 's' option removes the folder and all subfolders/files
    end
    
    % Create a fresh folder
    mkdir(folderName);
    %Store data
    writematrix(var.xi_list,folderName+"xi_"+"dt_"+string(var.dt)+".txt")
    writematrix(var.alpha_list,folderName+"alpha_"+"dt_"+string(var.dt)+".txt")
    writematrix(var.time_list,folderName+"time_"+"dt_"+string(var.dt)+".txt")
    writematrix([var.eps_x_list var.sigma_x_list],folderName+"stress_strain_"+"dt_"+string(var.dt)+".txt")
    writematrix(var.epsp_x_list,folderName+"viscostr_"+"dt_"+string(var.dt)+".txt")
end