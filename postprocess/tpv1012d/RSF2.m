data1 = load('faultEvo.txt');
x=-36e3/2:100:36e3/2

filename2 = 'sliprate.gif';
for i =  1:27
figure(1); 
plot(linspace(-18e3*2,18e3*2,360*2)',data1(i,:)*2,'-r');hold on;
 index = (i-1)*2
filename = join(['RSFData5/data2d_100m_traction_' num2str(index) '.txt']);
 data = readmatrix(filename);
  plot(data(1:3:end,3),smooth(data(1:3:end,2)),'-b');
%   filename = join(['RSFData3/data2d_100m_traction_' num2str(index) '.txt']);
%  data = readmatrix(filename);
%   plot(data(:,3),data(:,2),'--b');
  hold off;
  xlim([-15e3 15e3])

    legend('SBI','MOOSE_FARM','Interpreter',"Latex")
    title("SCEC Benchmark TPV104 ")
  drawnow; pause(0.1)


  frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      del =0.1;
      if i == 1;
        imwrite(imind,cm,filename2,'gif','Loopcount',inf,'DelayTime',del);
      else
        imwrite(imind,cm,filename2,'gif','WriteMode','append','DelayTime',del);
      end
    %count = count + 1;


end

