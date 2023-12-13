%take out_jacobian files do post-process to match MOOSE c++ format
clear all; close all; clc;
nums = ['11','12','13','14','21','22','23','24','31','32','33','34','41','42','43','44'];
for num = 1:length(nums)/2
    num_str = nums(2*num-1:2*num);
    num_names = append('./out_jacobian/eqout',num_str,'.txt');
    disp(num_names)
    fileID = fopen(num_names,'r');
    formatSpec = '%c';
    v = fscanf(fileID,formatSpec);
    
    %replace string
    
    string_minus_epse11 = 'eps11p_inc - eps11t + eps11p_pre';
    string_minus_epse22 = 'eps22p_inc - eps22t + eps22p_pre';
    string_minus_epse33 = 'eps33p_inc - eps33t + eps33p_pre';
    string_minus_epse12 = 'eps12p_inc - eps12t + eps12p_pre';
    string_minus_epse11_epse22_epse33 = 'eps11p_inc - eps22t - eps33t - eps11t + eps22p_inc + eps11p_pre + eps33p_inc + eps22p_pre + eps33p_pre';
    string_minus_2epse11 = '2*eps11p_inc - 2*eps11t + 2*eps11p_pre';
    string_minus_2epse22 = '2*eps22p_inc - 2*eps22t + 2*eps22p_pre';
    string_minus_2epse33 = '2*eps33p_inc - 2*eps33t + 2*eps33p_pre';
    string_minus_epseI2 = '(minus_eps11e)^2 + 2*(minus_eps12e)^2 + (minus_eps22e)^2 + (minus_eps33e)^2';
    string_minus_4epse12 = '4*eps12p_inc - 4*eps12t + 4*eps12p_pre';
    
    v = replace(v,string_minus_epse11,'minus_eps11e');
    v = replace(v,string_minus_epse22,'minus_eps22e');
    v = replace(v,string_minus_epse33,'minus_eps33e'); 
    v = replace(v,string_minus_epse12,'minus_eps12e');
    v = replace(v,string_minus_epse11_epse22_epse33,'minus_eps_tr');
    v = replace(v,string_minus_2epse11,'2*minus_eps11e');
    v = replace(v,string_minus_2epse22,'2*minus_eps22e');
    v = replace(v,string_minus_2epse33,'2*minus_eps33e');
    v = replace(v,string_minus_epseI2,'minus_epse_I2');
    v = replace(v,string_minus_4epse12,'4*minus_eps12e');
    
    %replace label
    string_cg = 'C_g';
    string_dt = 'dt';
    string_a0 = 'a0';
    string_a1 = 'a1';
    string_a2 = 'a2';
    string_a3 = 'a3';
    string_m2 = 'm2';
    string_m1 = 'm1';
    string_coeff = 'effec_sts_coeff';
    string_p = 'pressure';
    
    v = replace(v,string_cg,'_C_g');
    v = replace(v,string_dt,'_dt');
    v = replace(v,string_a0,'_a0');
    v = replace(v,string_a1,'_a1');
    v = replace(v,string_a2,'_a2');
    v = replace(v,string_a3,'_a3');
    v = replace(v,string_m2,'_m2');
    v = replace(v,string_m1,'_m1');
    v = replace(v,string_coeff,'_effec_sts_coeff');
    v = replace(v,string_p,'_pressure_neg');
    
    %replace power
    string_power_1 = '(minus_epse_I2)^(3/2)';
    string_power_2 = '(minus_epse_I2)^(1/2)';
    string_power_3 = '(minus_eps_tr)^3';
    string_power_4 = '(minus_epse_I2)^(5/2)';
    string_power_5 = '(minus_eps_tr)^2';
    string_power_6 = 'B^_m1';
    
    v = replace(v,string_power_1,'pow(minus_epse_I2,1.5)');
    v = replace(v,string_power_2,'pow(minus_epse_I2,0.5)');
    v = replace(v,string_power_3,'pow(minus_eps_tr,3.0)');
    v = replace(v,string_power_4,'pow(minus_epse_I2,2.5)');
    v = replace(v,string_power_5,'pow(minus_eps_tr,2.0)');
    v = replace(v,string_power_6,'pow(B,_m1)');
    
    %
    out_names = append('./out_format/eqout',num_str,'_formatted.txt');
    writematrix(v,out_names);
end

