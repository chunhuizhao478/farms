%check the slope of mean stress and deviatroic stress 2nd invariants

%if it is predicting that the strength is not too low...i.e. in the expected range for rocks. 
%Rocks have pressure sensitive strength. The proportionality between the strength and the pressure is the angle of internal friction. 
%Natural rocks have values for the tangent of the angle of internal friction to be between 0.5 and 0.9

% Meansts = [2.506e7 4.52e7 8.8e7 1.7525e8 3.503e8];
% devstsI2 = [1.3613e15 3.81e15 1.3829e16 5.44386e16 2.1723e17];
% 
% plot([0 0.2506],[0.06 0.368958],'b--');
% hold on;
% plot(Meansts/1e8,sqrt(devstsI2)/1e8,'b-*')
% ylabel("2nd invariant of deviatroic stress at peak stress",FontSize=20,FontName='Times New Roman')
% xlabel("mean stress at peak stress",FontSize=20,FontName='Times New Roman')
% title("I2 deviatroic stress vs mean stress at peak stress",FontSize=20,FontName='Times New Roman')

clear all; close all; clc;

term1  = [1.636e6 12.5e6 19.054e6 35.246e6 47.835e6 75.189e6 135.55e6 264.03e6 525.8e6];
term2  = [3.831e6 8.66e6 11.578e6 16.019e6 18.954e6  26.90e6  43.619e6 83.156e6 165e6];

% plot([0 23.56],[13.5 13.585],'b-.');
% hold on;
plot(term1/1e6,term2/1e6,'b-*',LineWidth=3)
ylabel("$\sqrt{J_2} (MPa)$",'Interpreter','latex',FontSize=25,FontName='Times New Roman')
xlabel("$I_1 (MPa)$" ,'Interpreter','latex',FontSize=25,FontName='Times New Roman')
title("Drucker–Prager Yield Criterion",FontSize=30,FontName='Times New Roman')
ax = gca;
ax.FontSize = 25; 
grid off;
%c  =  14.193MPa (16=A)
%phi = 39.84deg

A = 4.95;
B = 0.304;

syms phi c
phisol = vpasolve(B*sqrt(3)*(3-sin(phi))==2*sin(phi),phi);
eval(phisol(1)*180/pi)

csol = vpasolve(A*sqrt(3)*(3-sin(phisol(1)))==6*c*cos(phisol(1)),c);
eval(csol)

