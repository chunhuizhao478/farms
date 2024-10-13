%% 1D-Continuum-Damage-Breakage-Model %%
%% main code %%
clear all; clc; close all;
format long;
%%save path
folder_name = 'results_monothlic/';
%% import parameter %%
[param] = struct_param();
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
%% while loop %%
while var.t < var.Tmax
      %% Set Time Step & strain increment %%
      var.strain_rate = -1e-2;
      var.dt = 1e-2;
      var.depsx = var.strain_rate * var.dt;
      % Exit loop if B = 1
      if var.B == 1.0
            disp('simulation completed.');
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
save_data(var,param,folder_name);
disp("All done!")
