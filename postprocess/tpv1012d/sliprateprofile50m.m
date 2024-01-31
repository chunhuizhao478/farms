clear all; close all; clc;
%uguca results
data0uguca = load("./TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx2880_s2.00_tf0.35_npc1-DataFiles/vel0.txt");

%moose results
data2 = load("./files/50m_922_uguca_setup/sliprate.txt");
x2 = load("./files/50m_922_uguca_setup/interface_coordx.txt");

filename2 = 'sliprate50m.gif';
for i = 1:4:40
figure(1);

plot(linspace(-35987.5,35987.5,1440*2)',data0uguca(i,:)*2,'-k');hold on;
plot(x2,data2(4*i-3,:),'-b'); hold on;

% hold off;
xlim([-15e3 15e3])
legend('uguca-25m','moose-50m','location','northeast',FontSize=10)
title("SCEC Benchmark TPV101","FontSize",30)
xlabel("Along Fault (m)","FontSize",25)
ylabel("Slip Rate (m/s)","FontSize",25)
ylim([0,14])
ax = gca;
ax.FontSize = 15; 
end
% drawnow; pause(0.1)
% frame = getframe(1);
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% del =0.1;
% if i == 1
% imwrite(imind,cm,filename2,'gif','Loopcount',inf,'DelayTime',del);
% else
% imwrite(imind,cm,filename2,'gif','WriteMode','append','DelayTime',del);
% end
% end
