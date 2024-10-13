%% 1D-Continuum-Damage-Breakage-Model %%
%% main code %%
% clear all; clc; close all;
format long;
%% import parameter %%
[param] = struct_param_cycleloading(C1_i,CBH_i);
%% import variable %%
[var] = struct_var();
%% initialize material variable %%
[var] = comp_modulus(param,var);
[var] = comp_young_nu(param,var,"init");
%% initialize strain tensor %%
[var] = init_str_eps(param,var);
[var] = comp_xi(var);
xi_pre = var.xi;
%% Init Strain %%
store_init_strain = var.eps(1,1);
store_init_stress = var.sigma(1,1);
%% define alpha threshold
given_total_strain = strain_i;
applied_total_strain = 0;
%% number of cycles
num_cycle = num_cycle_val;
cycle_index = 0;
var.cycle_index_list = [];
flag = 1;
%% while loop %%
while var.t < var.Tmax
      %% Set Time Step & strain increment %%
%       %{
      % if xi_pre > var.xi
      %     var.strain_rate = -1e-4;
      %     var.dt = 1e-2;
      %     var.depsx = var.strain_rate * var.dt;
      % else
      %     var.strain_rate = -1e-4;
      %     var.dt = 1e-1;
      %     var.depsx = var.strain_rate * var.dt;
      % end
      %% Set Time Step & strain increment %%
      if flag == 1
        % Loading phase: apply negative strain increment until given_total_strain is reached
        var.strain_rate = -1e-2;
        var.dt = 1e-2;
        var.depsx = var.strain_rate * var.dt;
        applied_total_strain = applied_total_strain + var.depsx;
        
        if applied_total_strain <= given_total_strain
            % Switch to unloading phase
            flag = -1;  % Change flag to indicate unloading phase
        end
      elseif flag == -1
        % Unloading phase: apply positive strain increment until strain (or stress) reaches 0
        var.strain_rate = 1e-2;
        var.dt = 1e-2;
        var.depsx = var.strain_rate * var.dt;
        applied_total_strain = applied_total_strain + var.depsx;
        
        if var.sigma(1,1)>=-1e7
            % Switch back to loading phase for next cycle
            flag = 1;  % Change flag to indicate loading phase
            cycle_index = cycle_index + 1;  % Increment cycle index
            var.cycle_index_list = [var.cycle_index_list cycle_index];
            applied_total_strain = 0;  % Reset total strain for next cycle
            disp(['Cycle: ', num2str(cycle_index), ' completed.']);
        end
      end
      % Exit loop if number of cycles is completed
      if cycle_index >= num_cycle
            disp('All cycles completed.');
            break;
      end
      %% Compute alpha & B %%
      if var.xi < param.xi_d
          case_index = '0';
      elseif var.xi >= param.xi_d && var.xi < param.xi_0
          case_index = '1';
      elseif var.xi >= param.xi_0 && var.xi <= param.xi_max
          case_index = '2';
      else
          error("xi is out of range !")
      end
      [tt,y] = ode45(@(t,y)comp_alphaB_odesys(y,param,var,case_index),[var.t,var.t+var.dt],[var.alpha_,var.B]);
      alpha_new_comp = y(end,1); B_new_comp = y(end,2);
      % Check whether alpha,B are within reasonable range
      [var] = check_alphaB(alpha_new_comp,B_new_comp,var);
      % Compute Prob
      alpha_cr = comp_alpha_cr(var.xi,param);
      Prob = 1 / ( exp( (alpha_cr - var.alpha_) / param.beta ) + 1 );
      var.Prob = Prob;
      % Update material variable %%
      [var] = comp_modulus(param,var);
      %% Compute strain %%
      xi_pre = var.xi;
      %[var] = comp_eps_sys_new(param, var);
      [var] = comp_eps_sys_fsolve(param,var);
      %[var] = comp_eps_sys_loop(param,var);
      %% Time Increment %%
      var.t = var.t + var.dt;
      %% Store results;
      [var] = store_res(var);
      %% Determine Which state / Printout variables %%
      alpha_cr_new = comp_alpha_cr(var.xi,param);
      if var.Prob == 0 %%var.alpha_ < alpha_cr_new 
          fprintf("simulaton time: %.4f, xi: %.10f, alpha: %.10f, B: %.10f, alpha_cr: %.10f eps_p %.10f, eps_e %.10f eps: %.10f, Prob: %d, state: %s \n", ...
                  var.t,var.xi,var.alpha_,var.B,alpha_cr_new,var.eps_p(1,1),var.eps_e(1,1),var.eps(1,1), var.Prob, "solid")
      elseif var.Prob > 0 %%(var.alpha_ >= alpha_cr_new && var.alpha_ <= 1)
          fprintf("simulaton time: %.4f, xi: %.10f, alpha: %.10f, B: %.10f, alpha_cr: %.10f eps_p %.10f  eps_e %.10f eps: %.10f, Prob: %d, state: %s \n", ...
                  var.t,var.xi,var.alpha_,var.B,alpha_cr_new,var.eps_p(1,1),var.eps_e(1,1),var.eps(1,1), var.Prob, "granular")
      else
          error("alpha is out of range !")
      end
end
%% Plot figures %%
plot_figure(var);
%% Save data %%
save_data(var,param,folderName);
disp("All done!")
