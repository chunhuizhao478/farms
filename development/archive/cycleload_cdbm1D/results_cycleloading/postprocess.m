% Post-processing
clear all; close all; clc;

% Plot params
plot_params.frontsize_titles = 20;
plot_params.frontsize_legend = 12;
plot_params.frontsize_ticklabels = 12;
plot_params.linewidth = 0.5;
plot_params.ploteveynsteps = 10;

%% stress-strain
% Load files
ss_1 = readmatrix("./strain_-0.02_C1_30_CBH_300/stress_strain_dt_0.01.txt");
ss_2 = readmatrix("./strain_-0.02_C1_3_CBH_30/stress_strain_dt_0.01.txt");
ss_3 = readmatrix("./strain_-0.02_C1_0.3_CBH_3/stress_strain_dt_0.01.txt");

%update params
plot_params.xlabel_str = "total strain $\epsilon_{xx}$";
plot_params.ylabel_str = "total stress $\sigma_{xx}$ (GPa)";
plot_params.plot_title = "stress-strain curve of cycle loading";
plot_params.legend_labels = {'C1 = 30 1/s, C_{BH} = 300 1/s', 'C1 = 3 1/s, C_{BH} = 30 1/s', 'C1 = 0.3 1/s, C_{BH} = 3 1/s'};
plot_params.line_styles = {'r-','b-','k-'};
plot_params.marker_styles = {'ro','bo','ko'};
plot_params.marker_facecolors = {'r','b','k'};
plot_params.framelimits = [0 0.1 0 1.8];

% Create data for plotting
x_data = {-ss_1(:,1), -ss_2(:,1), -ss_3(:,1)};
y_data = {-ss_1(:,2)/1e9, -ss_2(:,2)/1e9, -ss_3(:,2)/1e9};

% Plot stress-strain figure
plot_params.framelimits = [0 0.14 0 0.6];
postfunc_plotfigures(x_data, y_data, plot_params, './stress_strain_strains_large_strain_-0.02.png');
plot_params.framelimits = [0 0.03 0 0.6];
postfunc_plotfigures(x_data, y_data, plot_params, './stress_strain_strains_close_strain_-0.02.png');
% Generate GIF
% plot_params.framelimits = [0 0.14 0 0.6];
% postfunc_generategifs(x_data, y_data, plot_params,'./stress_strain_close_strain_-0.02.gif');
% plot_params.framelimits = [0 0.03 0 0.6];
% postfunc_generategifs(x_data, y_data, plot_params,'./stress_strain_large_strain_-0.02.gif');

%% stress-damage
% Load files
ss_1 = readmatrix("./strain_-0.02_C1_30_CBH_300/stress_strain_dt_0.01.txt");
ss_2 = readmatrix("./strain_-0.02_C1_3_CBH_30/stress_strain_dt_0.01.txt");
ss_3 = readmatrix("./strain_-0.02_C1_0.3_CBH_3/stress_strain_dt_0.01.txt");

d_1 = readmatrix("./strain_-0.02_C1_30_CBH_300/alpha_dt_0.01.txt");
d_2 = readmatrix("./strain_-0.02_C1_3_CBH_30/alpha_dt_0.01.txt");
d_3 = readmatrix("./strain_-0.02_C1_0.3_CBH_3/alpha_dt_0.01.txt");

% Create data for plotting
x_data = {d_1(:,1), d_2(:,1), d_3(:,1)};
y_data = {-ss_1(:,2)/1e9, -ss_2(:,2)/1e9, -ss_3(:,2)/1e9};

%update params
plot_params.xlabel_str = "damage $\alpha$";
plot_params.ylabel_str = "total stress $\sigma_{xx}$ (GPa)";
plot_params.plot_title = "stress-damage curve of cycle loading";
plot_params.legend_labels = {'C1 = 30 1/s, C_{BH} = 300 1/s', 'C1 = 3 1/s, C_{BH} = 30 1/s', 'C1 = 0.3 1/s, C_{BH} = 3 1/s'};
plot_params.line_styles = {'r-','b-','k-'};
plot_params.marker_styles = {'ro','bo','ko'};
plot_params.marker_facecolors = {'r','b','k'};
plot_params.framelimits = [0 0.6 0 0.6];

% Plot stress-strain figure
postfunc_plotfigures(x_data, y_data, plot_params, './stress_damage_strain_-0.02.png');
% Generate GIF
% postfunc_generategifs(x_data, y_data, plot_params,'./stress_damage_strain_-0.02.gif');

%% xi-damage
% Load files
xi_1 = readmatrix("./strain_-0.02_C1_30_CBH_300/xi_dt_0.01.txt");
xi_2 = readmatrix("./strain_-0.02_C1_3_CBH_30/xi_dt_0.01.txt");
xi_3 = readmatrix("./strain_-0.02_C1_0.3_CBH_3/xi_dt_0.01.txt");

d_1 = readmatrix("./strain_-0.02_C1_30_CBH_300/alpha_dt_0.01.txt");
d_2 = readmatrix("./strain_-0.02_C1_3_CBH_30/alpha_dt_0.01.txt");
d_3 = readmatrix("./strain_-0.02_C1_0.3_CBH_3/alpha_dt_0.01.txt");

% Create data for plotting
x_data = {xi_1, xi_2, xi_3};
y_data = {d_1(:,1), d_2(:,1), d_3(:,1)};

%update params
plot_params.xlabel_str = "strain invariant ratio $\xi$";
plot_params.ylabel_str = "damage $\alpha$";
plot_params.plot_title = "stress-damage curve of cycle loading";
plot_params.legend_labels = {'C1 = 30 1/s, C_{BH} = 300 1/s', 'C1 = 3 1/s, C_{BH} = 30 1/s', 'C1 = 0.3 1/s, C_{BH} = 3 1/s'};
plot_params.line_styles = {'r-','b-','k-'};
plot_params.marker_styles = {'ro','bo','ko'};
plot_params.marker_facecolors = {'r','b','k'};
plot_params.framelimits = [-sqrt(3) 0.6 0 0.6];

% Plot stress-strain figure
postfunc_plotfigures(x_data, y_data, plot_params, './xi_damage_strain_-0.02.png');
% Generate GIF
% postfunc_generategifs(x_data, y_data, plot_params,'./xi_damage_strain_-0.02.gif');