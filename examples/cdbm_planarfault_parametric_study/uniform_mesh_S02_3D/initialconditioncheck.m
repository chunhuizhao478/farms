clear all; close all; clc;
%check initial values
b22 = 0.926793;
b33 = 1.073206;
b23 = -0.169029;

zcoords = 0:10:15000;

pf = 1000 * 9.8 * zcoords;
sigma11 = -2670 * 9.8 * zcoords;

sigma22 = zeros(size(sigma11));
for i = 1 : numel(sigma22)
    if zcoords(i) <= 15600
        sigma22(i) = b22 .* (sigma11(i) + pf(i)) - pf(i);
    else
        sigma22(i) = sigma11(i);
    end
end


sigma33 = zeros(size(sigma11));
for i = 1 : numel(sigma33)
    if zcoords(i) <= 15600
        sigma33(i) = b33 .* (sigma11(i) + pf(i)) - pf(i);
    else
        sigma33(i) = sigma11(i);
    end
end

sigma23 = zeros(size(sigma11));
for i = 1 : numel(sigma23)
    if zcoords(i) <= 15600
        sigma23(i) = b23 .* (sigma11(i) + pf(i));
    else
        sigma23(i) = 0;
    end
end

sigma_xx = sigma22;
sigma_yy = sigma33;
sigma_zz = sigma11;
sigma_xy = sigma23;

figure();
plot(zcoords,sigma_xx,'r-'); hold on;
plot(zcoords,sigma_yy,'m-'); hold on;
plot(zcoords,sigma_zz,'b-'); hold on;
plot(zcoords,sigma_xy,'k-'); hold on;
legend("sigma_xx","sigma_yy","sigma_zz","sigma_xy");

fprintf("MAX sigma_xx (MPa): %.3f \n", max(abs(sigma_xx))/1e6);
fprintf("MAX sigma_yy (MPa): %.3f \n", max(abs(sigma_yy))/1e6);
fprintf("MAX sigma_zz (MPa): %.3f \n", max(abs(sigma_zz))/1e6);
fprintf("MAX sigma_xy (MPa): %.3f \n", max(sigma_xy)/1e6);