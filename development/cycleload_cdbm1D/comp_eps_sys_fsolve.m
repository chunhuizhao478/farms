function [var] = comp_eps_sys_fsolve(param,var)
%Compute total strain / elastic strain / plastic strain
%% Get initial pressure
eps_p_pre = var.eps_p;
%% Assign prescribed strain
%eps11t = var.eps(1,1) + var.depsx;
eps11_inc = var.depsx;
%% Construct system of equations & initial value
%F = @(x)eps_func(x,eps11t,eps_p_pre,param,var);
%F = @(x)eps_func(x,eps11_inc,eps_p_pre,param,var);
%x0 = [0,0,0,0,0]; %use value from previous step
%% Set output options
%options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-30,'FunctionTolerance',1e-30,'StepTolerance',1e-30);
%% Solve the system of equations
%[x,fval,exitflag] = fsolve(F,x0,options);
%disp(fval)
%disp(exitflag)
%% Get individual variable
%eps22_inc = x(1);
%eps33_inc = x(2);
%eps_p11_inc = x(3);
%eps_p22_inc = x(4);
%eps_p33_inc = x(5);
%% Symbolic Solve
[eps22_inc,eps33_inc,eps_p11_inc,eps_p22_inc,eps_p33_inc] = eps_func_symbolic(eps11_inc,eps_p_pre,param,var);
%% Update total strain/ plastic strain/ elastic strain
var.eps   = [ var.eps(1,1)  + eps11_inc    0 0; 
              0 var.eps(2,2) + eps22_inc     0; 
              0 0 var.eps(3,3) + eps33_inc    ];
var.eps_p = [ var.eps_p(1,1) + eps_p11_inc 0 0; 
              0 var.eps_p(2,2) + eps_p22_inc 0; 
              0 0 var.eps_p(3,3) + eps_p33_inc];
var.eps_e = var.eps - var.eps_p;
%% Update I1, I2, xi & total stress
[var] = comp_xi(var);
[var] = comp_sigma(param,var);
end