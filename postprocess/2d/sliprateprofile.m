% data1 = load("./TPV101_Nx/TPV101_Nx720_s2.00_tf0.35_npc1-DataFiles/vel1.txt");
% x1=load("./TPV101_Nx/TPV101_Nx720_s2.00_tf0.35_npc1.coord");
% x1=x1(:,1);
% x1=x1 - 35850;

data1 = load('faultEvo.txt');
x=-36e3/2:100:36e3/2;

data2 = load("./files/100m_919_fixnormal/sliprate.txt");
x2 = load("./files/100m_919_fixnormal/interface_coordx.txt");

% data3 = load("./files/100m_919_damp08/sliprate.txt");
% x3 = load("./files/100m_919_damp08/interface_coordx.txt");

% data4 = load("./files/100m_919_2/sliprate.txt");
% x4 = load("./files/100m_919_2/interface_coordx.txt");

filename2 = 'sliprate.gif';
for i = 1:30
figure(1); 
% plot(x1,data1(i,:)*2,'-r');hold on;
plot(linspace(-35950,35950,360*2)',data1(i,:)*2,'-r');hold on;
plot(x2,data2(4*i-3,:),'-b'); hold on;
% plot(x3,data3(4*i-3,:),'-k'); hold on;
% plot(x4,data4(4*i-3,:),'-k');
hold off;
xlim([-15e3 15e3])
legend('uguca-100m','MOOSE-100m-fixnormal')
title("SCEC Benchmark TPV101 ")
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
