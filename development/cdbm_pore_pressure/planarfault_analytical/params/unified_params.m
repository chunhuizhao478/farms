%Unified Parameter Choice For CBDM Complex Network Problem
%---------------------------------------------------------%
mu_s = 0.6;
mu_d = 0.5;
Dc = 0.1;
tau_S = 70e6;
sigma_N = 120e6;
G = 32.04e9;
%---------------------------------------------------------%
tau_strength = mu_s * sigma_N;
diff = tau_strength - tau_S;
mu = tau_S / sigma_N;
S = (mu_s - mu)/(mu - mu_d);
L = G * Dc / ( (mu_s - mu_d) * sigma_N );
%---------------------------------------------------------%
fprintf('tau_strength/1e6 = %.3f \n',tau_strength/1e6);
fprintf('diff/1e6 = %.3f \n',diff/1e6);
fprintf('mu = %.3f \n',mu);
fprintf('S = %.3f \n',S);
fprintf('L = %.3f \n',L);