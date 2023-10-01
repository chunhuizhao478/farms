%% Make video

png_path = "./files/100m_918";

mp4_path = png_path;

% %% slip rate
writeObj = VideoWriter(mp4_path + "/sliprate.mp4","MPEG-4");
% writeObj = VideoWriter(mp4_path + "/sliprate" + current_fault + ".avi");
writeObj.FrameRate = 4;
writeObj.Quality = 100;
open(writeObj);
for i = 0 : 1 : 29
    filename = sprintf(png_path +  "/sliprate"+"%d.png", i);
    thisimage = imread(filename);
    drawnow;
    writeVideo(writeObj, thisimage);
end
close(writeObj);

% %% slip
% writeObj = VideoWriter(mp4_path + "/slip" + current_fault + ".mp4","MPEG-4");
% writeObj.FrameRate = 8;
% writeObj.Quality = 100;
% open(writeObj);
% for i = 0 : 300
%     filename = sprintf(png_path +  "/slip_" +"%d.png", i);
%     thisimage = imread(filename);
%     drawnow;
%     writeVideo(writeObj, thisimage);
% end
% close(writeObj);
% 
% %% tangent traction
% writeObj = VideoWriter(avi_path + "/tangtraction" + current_fault + ".avi");
% writeObj.FrameRate = 8;
% writeObj.Quality = 100;
% open(writeObj);
% for i = 0 : 171
%     filename = sprintf(png_path + "/tangtraction_" +"%d.png", i);
%     thisimage = imread(filename);
%     drawnow;
%     writeVideo(writeObj, thisimage);
% end
% close(writeObj);
% %% tangent traction global
% writeObj = VideoWriter(mp4_path + "/tangtraction" + current_fault + ".mp4","MPEG-4");
% writeObj.FrameRate = 8;
% writeObj.Quality = 100;
% open(writeObj);
% for i = 0 : 300
%     filename = sprintf(png_path + "/tangtraction_" +"%d.png", i);
%     thisimage = imread(filename);
%     drawnow;
%     writeVideo(writeObj, thisimage);
% end
% close(writeObj);
% %% normal jump rate
% writeObj = VideoWriter(avi_path + "/normaljumprate" + current_fault + ".avi");
% writeObj.FrameRate = 8;
% writeObj.Quality = 100;
% open(writeObj);
% for i = 0 : 700
%     filename = sprintf(png_path +  "/normaljumprate_" +"%d.png", i);
%     thisimage = imread(filename);
%     drawnow;
%     writeVideo(writeObj, thisimage);
% end
% close(writeObj);

% %% normal jump
% writeObj = VideoWriter(avi_path + "/normaljump" + current_fault + ".avi");
% writeObj.FrameRate = 8;
% writeObj.Quality = 100;
% open(writeObj);
% for i = 0 : 700
%     filename = sprintf(png_path +  "/normaljump_" +"%d.png", i);
%     thisimage = imread(filename);
%     drawnow;
%     writeVideo(writeObj, thisimage);
% end
% close(writeObj);

% %% normal traction
% writeObj = VideoWriter(mp4_path + "/normaltraction" + current_fault + ".mp4","MPEG-4");
% writeObj.FrameRate = 8;
% writeObj.Quality = 100;
% open(writeObj);
% for i = 0 : 300
%     filename = sprintf(png_path + "/normaltraction_" +"%d.png", i);
%     thisimage = imread(filename);
%     drawnow;
%     writeVideo(writeObj, thisimage);
% end
% close(writeObj);
% %% shear stress / strength
% writeObj = VideoWriter(mp4_path + "/tang_stress_strength" + current_fault + ".mp4","MPEG-4");
% writeObj.FrameRate = 8;
% writeObj.Quality = 100;
% open(writeObj);
% for i = 0 : 300
%     filename = sprintf(png_path + "/tang_stress_strength_" +"%d.png", i);
%     thisimage = imread(filename);
%     drawnow;
%     writeVideo(writeObj, thisimage);
% end
% close(writeObj);