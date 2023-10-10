clear all; close all; clc

xi = -1.4:0.01:1.4;
mu_0 = 32.04e9;
xi_0 = -0.8;
lambda = 32.04e9;
gamma_r = 3.7150e10;
a0 = 4.9526e9;
a1 = -1.8888e10;
a2 = 2.3960e10;
a3 = -1.0112e10;

%
alpha = 0;

gamma_damaged = alpha * gamma_r;
mu = mu_0 + alpha * xi_0 * gamma_r;

Fs = mu - gamma_damaged .* xi + 0.5 .* lambda .* xi .^ 2;
Fb = a0 + a1 .* xi + a2 .* xi .^ 2 + a3 .* xi .^ 3;

figure();
plot(xi,Fb/1e9,'b-',LineWidth=2); hold on;
plot(xi,Fs/1e9,'r-'); hold on;

%
alpha = 0.5;

gamma_damaged = alpha * gamma_r;
mu = mu_0 + alpha * xi_0 * gamma_r;

Fs = mu - gamma_damaged .* xi + 0.5 .* lambda .* xi .^ 2;

plot(xi,Fs/1e9,'r--'); hold on;

%
alpha = 0.75;

gamma_damaged = alpha * gamma_r;
mu = mu_0 + alpha * xi_0 * gamma_r;

Fs = mu - gamma_damaged .* xi + 0.5 .* lambda .* xi .^ 2;

plot(xi,Fs/1e9,'r-.'); hold on;

%
alpha = 1.0;

gamma_damaged = alpha * gamma_r;
mu = mu_0 + alpha * xi_0 * gamma_r;

Fs = mu - gamma_damaged .* xi + 0.5 .* lambda .* xi .^ 2;

plot(xi,Fs/1e9,'r.'); hold on;

title("Energy vs $\xi$ under various $\alpha$ values",Interpreter="latex",FontSize=20)
legend("Granular Energy $F_b$", "Solid Energy $F_s$, $\alpha$ = 0", "Solid Energy $F_s$, $\alpha$ = 0.5", "Solid Energy $F_s$, $\alpha$ = 0.75","Solid Energy $F_s$, $\alpha$ = 1",Interpreter="latex",FontSize=15);
xlabel("strain invariant ratio $\xi$",Interpreter="latex",FontSize=15)
ylabel("Energy $F_s$ ($F_b$) $GJ$",Interpreter="latex",FontSize=15)
grid on;
ax = gca;
ax.FontSize = 15;
xlim([-1.5 0]);
