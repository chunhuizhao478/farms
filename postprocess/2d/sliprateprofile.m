clear all; close all; clc;
%uguca results
% data0uguca = load("./TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx2880_s2.00_tf0.35_npc1-DataFiles/vel0.txt");

data1uguca = load("./TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx1440_s2.00_tf0.35_npc1-DataFiles/vel0.txt");

% data2uguca = load('faultEvo.txt');

%moose results
% data0 = load("./files/50m_922_uguca_setup/sliprate.txt");
% x0 = load("./files/50m_922_uguca_setup/interface_coordx.txt");
% 
data2 = load("./files/50m_922_uguca_setup/sliprate.txt");
x2 = load("./files/50m_922_uguca_setup/interface_coordx.txt");
% 
% data3 = load("./files/100m_920_dampcomp_0dot6/sliprate.txt");
% x3 = load("./files/100m_920_dampcomp_0dot6/interface_coordx.txt");
% 
% data4 = load("./files/100m_920_dampcomp_0dot8/sliprate.txt");
% x4 = load("./files/100m_920_dampcomp_0dot8/interface_coordx.txt");
% 
% data5 = load("./files/100m_920_dampcomp_1dot0/sliprate.txt");
% x5 = load("./files/100m_920_dampcomp_1dot0/interface_coordx.txt");

filename2 = 'sliprate.gif';
for i = 1:40
figure(1);
% plot(linspace(-36000,36000,1440*2)',data0uguca(i,:)*2,'-k');hold on;
plot(linspace(-35975,35975,720*2)',data1uguca(i,:)*2,'-k');hold on;
% plot(linspace(-35950,35950,360*2)',data2uguca(i,:)*2,'-k');hold on;
% plot(x0,data0(4*i-3,:),'-b'); hold on;
plot(x2,data2(4*i-3,:),'-b'); hold on;
% plot(x3,data3(4*i-3,:),'-k'); hold on;
% plot(x4,data4(4*i-3,:),'-m'); hold on;
% plot(x5,data5(4*i-3,:),'-b'); hold on;
hold off;
xlim([-15e3 15e3])
legend('uguca-50m','moose-50m','location','northeast',FontSize=10)
% legend('uguca-100m','MOOSE-100m-damp-1.0','Location','northoutside',FontSize=10)
% legend('uguca-100m','MOOSE-100m-damp-0.2','MOOSE-100m-damp-0.4','MOOSE-100m-damp-0.6','MOOSE-100m-damp-0.8','MOOSE-100m-damp-1.0',FontSize=10)
title("SCEC Benchmark TPV101","FontSize",30)
xlabel("Along Fault (m)","FontSize",25)
ylabel("Slip Rate (m/s)","FontSize",25)
ax = gca;
ax.FontSize = 15; 
drawnow; pause(0.1)
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
del =0.1;
if i == 1
imwrite(imind,cm,filename2,'gif','Loopcount',inf,'DelayTime',del);
else
imwrite(imind,cm,filename2,'gif','WriteMode','append','DelayTime',del);
end
end
