function [var] = comp_xi(var)
%Compute strain invariant ratio, I1, I2
%input: var output: update in var

%%total strain -> for update damage variable
eps_kkt = sum(diag(var.eps));
eps_ijijt = sum(var.eps.^2,"all");

%var.I1t = eps_kkt;
%var.I2t = eps_ijijt;
%var.xit = eps_kkt / sqrt(eps_ijijt);

%%elastic strain
eps_kk = sum(diag(var.eps_e));
eps_ijij = sum(var.eps_e.^2,"all");

var.I1 = eps_kk;
var.I2 = eps_ijij;
var.xi = eps_kk / sqrt(eps_ijij);

end