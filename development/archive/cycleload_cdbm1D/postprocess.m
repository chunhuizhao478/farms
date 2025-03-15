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
ss_1 = readmatrix("./results_sep29/stress_strain_fh-0.06dt_0.01_C1_30_CBH_300.txt");
ss_2 = readmatrix("./results_sep29/stress_strain_fh-0.06dt_0.01_C1_3_CBH_30.txt");
ss_3 = readmatrix("./results_sep29/stress_strain_fh-0.06dt_0.01_C1_0.3_CBH_3.txt");

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
plot_params.framelimits = [0 0.1 0 1.8];
postfunc_plotfigures(x_data, y_data, plot_params, 'figures_sep29/stress_strain_close.png');
plot_params.framelimits = [0 0.35 0 1.8];
postfunc_plotfigures(x_data, y_data, plot_params, 'figures_sep29/stress_strain_large.png');
% Generate GIF
plot_params.framelimits = [0 0.1 0 1.8];
postfunc_generategifs(x_data, y_data, plot_params,'figures_sep29/stress_strain_close.gif');
plot_params.framelimits = [0 0.35 0 1.8];
postfunc_generategifs(x_data, y_data, plot_params,'figures_sep29/stress_strain_large.gif');

%% stress-damage
% Load files
ss_1 = readmatrix("./results_sep29/stress_strain_fh-0.06dt_0.01_C1_30_CBH_300.txt");
ss_2 = readmatrix("./results_sep29/stress_strain_fh-0.06dt_0.01_C1_3_CBH_30.txt");
ss_3 = readmatrix("./results_sep29/stress_strain_fh-0.06dt_0.01_C1_0.3_CBH_3.txt");

d_1 = readmatrix("./results_sep29/alpha_fh-0.06dt_0.01_C1_30_CBH_300.txt");
d_2 = readmatrix("./results_sep29/alpha_fh-0.06dt_0.01_C1_3_CBH_30.txt");
d_3 = readmatrix("./results_sep29/alpha_fh-0.06dt_0.01_C1_0.3_CBH_3.txt");

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
plot_params.framelimits = [0 0.6 0 1.8];

% Plot stress-strain figure
postfunc_plotfigures(x_data, y_data, plot_params, 'figures_sep29/stress_damage.png');
% Generate GIF
postfunc_generategifs(x_data, y_data, plot_params,'figures_sep29/stress_damage.gif');

%% xi-damage
% Load files
xi_1 = readmatrix("./results_sep29/xi_fh-0.06dt_0.01_C1_30_CBH_300.txt");
xi_2 = readmatrix("./results_sep29/xi_fh-0.06dt_0.01_C1_3_CBH_30.txt");
xi_3 = readmatrix("./results_sep29/xi_fh-0.06dt_0.01_C1_0.3_CBH_3.txt");

d_1 = readmatrix("./results_sep29/alpha_fh-0.06dt_0.01_C1_30_CBH_300.txt");
d_2 = readmatrix("./results_sep29/alpha_fh-0.06dt_0.01_C1_3_CBH_30.txt");
d_3 = readmatrix("./results_sep29/alpha_fh-0.06dt_0.01_C1_0.3_CBH_3.txt");

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
plot_params.framelimits = [-sqrt(3) 0.5 0 0.6];

% Plot stress-strain figure
postfunc_plotfigures(x_data, y_data, plot_params, 'figures_sep29/xi_damage.png');
% Generate GIF
postfunc_generategifs(x_data, y_data, plot_params,'figures_sep29/xi_damage.gif');