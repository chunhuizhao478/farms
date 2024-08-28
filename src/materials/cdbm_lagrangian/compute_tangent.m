clear all; close all; clc;
%delete previous files
folderPath = './tangents';  % Specify the path to your folder
filePattern = fullfile(folderPath, '*');  % Create a pattern to match all files
delete(filePattern);  % Delete all files that match the pattern
%
syms Ee00 Ee01 Ee02 Ee10 Ee11 Ee12 Ee20 Ee21 Ee22 lambda damaged_modulus shear_modulus B_breakagevar a0 a1 a2 a3
Ee = [Ee00 Ee01 Ee02;
      Ee10 Ee11 Ee12;
      Ee20 Ee21 Ee22];
I = [1 0 0; 0 1 0; 0 0 1];
I1 = trace(Ee);
I2 = Ee00^2 + Ee01^2 + Ee02^2 + Ee10^2 + Ee11^2 + Ee12^2 + Ee20^2 + Ee21^2 + Ee22^2;
xi = I1 / sqrt(I2);
sigma_s = (lambda - damaged_modulus / xi) * I1 * I + (2*shear_modulus - damaged_modulus * xi) * Ee;
% sigma_b = (2 * a2 + a1 / xi + 3 * a3 * xi) * I1 * I + (2 * a0 + a1 * xi - a3 * xi^3) * Ee;
% sigma_total = (1-B_breakagevar) * sigma_s + B_breakagevar * sigma_b;

sigma_total = sigma_s;

T = zeros(3, 3, 3, 3);
for i = 1 : 3
    for j = 1 : 3
        for k = 1 : 3
            for l = 1 : 3
                Tijkl = string(diff(sigma_total(i,j),Ee(k,l)));
                writematrix(Tijkl,'./tangents/out_jacobian_'+string(i)+string(j)+string(k)+string(l)+'.txt');
            end
        end
    end
end