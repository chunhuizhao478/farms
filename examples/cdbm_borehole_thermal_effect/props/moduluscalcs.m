clear all; close all; clc;
%Vs,Vp
Vs = 2435.88;
Vp = 3981;
rho = 2641.7;
%
Ed = rho * Vs ^ 2 * ( 3 * Vp ^ 2 - 4 * Vs ^ 2 ) / ( Vp ^ 2 - Vs ^ 2 );
nud = 0.5 * ( Vp ^ 2 - 2 * Vs ^ 2 ) / ( Vp ^ 2 - Vs ^ 2 );
%
