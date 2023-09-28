clear all; close all; clc;
%uguca results
data0uguca = load('./TPV103_Nx/TPV103_Nx2880_s2.00_tf0.10_npc5-DataFiles/vel0.txt');

% data1uguca = load('./TPV103_Nx/TPV103_Nx1440_s2.00_tf0.10_npc20-DataFiles/vel0.txt');

data2uguca = load('./TPV103_Nx/TPV103_Nx720_s2.00_tf0.10_npc5-DataFiles/vel0.txt');

%moose results
data0 = load("./files/100m_926/sliprate.txt");
x0 = load("./files/100m_926/interface_coordx.txt");

filename2 = 'sliprate.gif';
for i = 1:40
figure(1);
plot(linspace(-35987.5,35987.5,1440*2)',data0uguca(i,:)*2,'-k');hold on;
% plot(linspace(-35975,35975,720*2)',data1uguca(i,:)*2,'-k');hold on;
plot(linspace(-35950,35950,360*2)',data2uguca(i,:)*2,'-b');hold on;
% plot(x0,data0(4*i-3,:),'-b'); hold on;
% plot(x2,data2(4*i-3,:),'-b'); hold on;
% plot(x3,data3(4*i-3,:),'-k'); hold on;
% plot(x4,data4(4*i-3,:),'-m'); hold on;
% plot(x5,data5(4*i-3,:),'-b'); hold on;
hold off;
xlim([-15e3 15e3])
legend('uguca-25m','uguca-100m','location','northeast',FontSize=10)
% legend('uguca-100m','MOOSE-100m-damp-1.0','Location','northoutside',FontSize=10)
% legend('uguca-100m','MOOSE-100m-damp-0.2','MOOSE-100m-damp-0.4','MOOSE-100m-damp-0.6','MOOSE-100m-damp-0.8','MOOSE-100m-damp-1.0',FontSize=10)
title("SCEC Benchmark TPV103","FontSize",30)
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
