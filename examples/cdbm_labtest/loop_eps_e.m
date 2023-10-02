function [var] = loop_eps_e(param,var)
%while loop for computing elastic strain/plastic strain/stress
%until the elastic does not change within a threshold value
%input: eps_new - new computed total strain
%       eps_p_pre - plastic strain not yet updated
%% 0.Initialize the difference between two iterations of elastic strain %%
var.eps_e = var.eps - var.eps_p;
eps_p_pre = var.eps_p;
eps_e_pre = var.eps_e;
diff = norm(eps_e_pre);
threshold = 1e-6;
%% 1.Start while loop/Check convergence %%
while diff > threshold
      %% 2.Compute I1,I2,xi %%
      [var] = comp_xi(var);
      %% 3.Compute stress %%
      [var] = comp_sigma(param,var);
      %% 4.Compute plastic strain %%
      [var] = comp_eps_p(eps_p_pre,param,var);
      %% 5.Recompute elastic strain %%
      var.eps_e = var.eps - var.eps_p;
      eps_e_new = var.eps_e;
      %% 6.Compute relative error
      diff = norm(eps_e_new-eps_e_pre)/norm(eps_e_pre);
      %% 7.Assign to eps_p_pre
      eps_e_pre = eps_e_new;
end
%% 8.Recompute all quantities based on converged eps_e %%
[var] = comp_xi(var);
[var] = comp_sigma(param,var);
end