clear all; close all; clc;
%Unified Parameter Choice For CBDM Complex Network Problem
%% Input Parameters
mu_s = 0.7;
mu_d = 0.5;
Dc = 0.1;
tau_S = 70e6;
sigma_N = 120e6;
G = 32.04e9;
%% Main Fault
tau_strength = mu_s * sigma_N;
mu = tau_S / sigma_N;
S = ( mu_s - mu ) / ( mu - mu_d );
L = G * Dc / ((mu_s - mu_d) * sigma_N);
%% Pintout
fprintf("mu_s (-): %.3f \n", mu_s);
fprintf("mu_d (-): %.3f \n", mu_d);
fprintf("D_c (-): %.3f \n", Dc);
fprintf("tau_S (MPa): %.3f \n", tau_S/1e6);
fprintf("sigma_N (MPa): %.3f \n", sigma_N/1e6);
fprintf("tau_strength (MPa): %.3f \n", tau_strength/1e6);
fprintf("pressure P needed to reach initial shear stress (MPa): %.3f \n", (sigma_N/1e6) - (tau_S/1e6) / mu_s);
fprintf("mu (-): %.3f \n", mu);
fprintf("S (-): %.3f \n", S);
fprintf("L (m): %.3f \n", L);