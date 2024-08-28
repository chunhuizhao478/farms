%take out_jacobian files do post-process to match MOOSE c++ format
clear all; close all; clc;
%delete previous files
folderPath = './tangents_formatted';  % Specify the path to your folder
filePattern = fullfile(folderPath, '*');  % Create a pattern to match all files
delete(filePattern);  % Delete all files that match the pattern
%
nums = ['1111';'1122';'1133';'1123';'1113';'1112';'2222';'2233';'2223';'2213';'2212';'3333';'3323';'3313';'3312';'2323';'2313';'2312';'1313';'1312';'1212'];
for num = 1:size(nums,1)
    num_str = nums(num,:);
    num_names = append('./tangents/out_jacobian_',num_str,'.txt');
    disp(num_names)
    fileID = fopen(num_names,'r');
    formatSpec = '%c';
    v = fscanf(fileID,formatSpec);
    
    %replace string
    
    string_I2 = 'Ee00^2 + Ee01^2 + Ee02^2 + Ee10^2 + Ee11^2 + Ee12^2 + Ee20^2 + Ee21^2 + Ee22^2';
    string_I1 = 'Ee00 + Ee11 + Ee22';
    string_I2power0dot5 = '(I2)^(1/2)';
    string_I2power1dot5 = '(I2)^(3/2)';
    string_I2power2dot5 = '(I2)^(5/2)';
    string_I1power2dot0 = '(I1)^2';
    string_I1power3dot0 = '(I1)^3';
   

    v = replace(v,string_I2,'I2');
    v = replace(v,string_I1,'I1');
    v = replace(v,string_I2power0dot5,'pow(I2,0.5)');
    v = replace(v,string_I2power1dot5,'pow(I2,1.5)');
    v = replace(v,string_I2power2dot5,'pow(I2,2.5)');
    v = replace(v,string_I1power2dot0,'pow(I1,2.0)');
    v = replace(v,string_I1power3dot0,'pow(I1,3.0)');
    
    %replace label
    string_cg = 'C_g';
    string_a0 = 'a0';
    string_a1 = 'a1';
    string_a2 = 'a2';
    string_a3 = 'a3';
    
    v = replace(v,string_cg,'_C_g');
    v = replace(v,string_a0,'_a0');
    v = replace(v,string_a1,'_a1');
    v = replace(v,string_a2,'_a2');
    v = replace(v,string_a3,'_a3');
    
    %
    out_names = append('./tangents_formatted/out_jacobian_',num_str,'_formatted.txt');
    writematrix(v,out_names);
end

