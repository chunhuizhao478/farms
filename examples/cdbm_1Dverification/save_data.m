function [] = save_data(var,given_alpha)
    %Store data
    writematrix(var.xi_list,"xi_fh"+string(given_alpha)+"dt_"+string(var.dt)+".txt")
    writematrix(var.alpha_list,"alpha_fh"+string(given_alpha)+"dt_"+string(var.dt)+".txt")
    writematrix(var.time_list,"time_fh"+string(given_alpha)+"dt_"+string(var.dt)+".txt")
    writematrix([var.eps_x_list var.sigma_x_list],"stress_strain_fh"+string(given_alpha)+"dt_"+string(var.dt)+".txt")
    writematrix(var.epsp_x_list,"viscostr_fh"+string(given_alpha)+"dt_"+string(var.dt)+".txt")
end