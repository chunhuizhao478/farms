function [] = plot_figure(var)
%Plot figure
%% Strain invariant ratio %%
%%digi_sir = load("digitial_sir.txt");
figure(1);
tsize = size(var.xi_list);
%time = linspace(0,1,tsize(1));
%plot(time,var.xi_list, 'r-'); hold on
%plot(digi_sir(:,1),digi_sir(:,2), 'b-');
%------non-normalized---------%
plot(var.time_list,var.xi_list, 'r-.');
%-----------------------------%

title('time vs strain invariants ratio')
xlabel('time')
ylabel('strain invariants ratio')
legend("Our Simulation")
%% Damage Variable %%
%digi_damage = load("digitial_alpha.txt");
figure(2);
%plot(time,var.alpha_list,"r-"); holplotd on;
%plot(digi_damage(:,1),digi_damage(:,2),"b-")
plot(var.time_list,var.alpha_list, 'r-.');
title('time vs damage variable')
xlabel('time')
ylabel('damage variable')
legend("Our Simulation")
%% Stress strain plot %%
figure(3);
plot(-(var.eps_x_list),-(var.sigma_x_list),'r-*'); hold on;
plot(-(var.eps_x_list),-(var.sigmas_x_list),'b-'); hold on;
plot(-(var.eps_x_list),-(var.sigmab_x_list),'k-'); hold on;
title('strain vs stress')
xlabel('strain')
ylabel('stress')
legend("total","elastic",'granular')
ylim([0 23e6])
%% Strain invariants ratio vs. Damage variable plot %%
figure(4);
plot(var.xi_list,var.alpha_list,"r-.")
title('strain inariant ratio vs damage variable')
xlabel('strain inariant ratio')
ylabel('damage variable')
legend("Our Simulation")
%% Strain vs. Damage variable plot %%
figure(5);
plot(-(var.eps_x_list),var.alpha_list,"r-.")
title('strain vs damage variable')
xlabel('strain inariant ratio')
ylabel('damage variable')
legend("Our Simulation")
%% Plastic strain (nondecreasing function) %%
figure(6);
plot(var.time_list,-(var.epsp_x_list),"r-"); hold on;
plot(var.time_list,-(var.epse_x_list),"b-"); hold on;
plot(var.time_list,-(var.eps_x_list),"k-"); hold on;
title('time vs viscoelastic strain')
xlabel('time')
ylabel('viscoelastic strain')
legend("viscoelastic strain","elastic strain","total strain")
end