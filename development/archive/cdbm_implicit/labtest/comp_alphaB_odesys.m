function [f_out] = comp_alphaB_odesys(y,param,var,case_index)
%%Solve the coupled system using ode solver (built-in time adaptivity)%%
%%y(1)=alpha_new y(2)=B_new

alpha_cr = comp_alpha_cr(var.xi,param);
Prob = 1 / ( exp( (alpha_cr - y(1)) / param.beta ) + 1 );

mu_new = param.mu_0 + y(1) * param.xi_0 * param.gamma_r;
gamma_new = y(1) * param.gamma_r;
lambda_new = param.lambda_0;

if case_index == '0'
f_out = [  ( param.C_1 * exp( y(1) / param.C_2 ) * var.I2 * ( var.xi - param.xi_0 ) );
           param.C_BH * var.I2 * ( var.xi - param.xi_0 )
           ];
         
end

if case_index == '1' 
f_out = [ ( param.C_1 * exp( y(1) / param.C_2 ) * var.I2 * ( var.xi - param.xi_0 ) ); 
           param.C_B * Prob * (1-var.B) * var.I2 * ( var.xi - param.xi_0) ];
           %( param.C_1 * exp( y(1) / param.C_2 ) * var.I2 * ( var.xi - param.xi_0 ) )  
           %param.C_B * Prob * ( var.I2 * (mu_new - gamma_new * var.xi + lambda_new/2 * var.xi ^ 2) - var.I2 * (param.a0 + param.a1 * var.xi + param.a2 * var.xi ^ 2 + param.a3 * var.xi ^ 3)) ]; %see ggw183.pdf
end

if case_index == '2' 
f_out = [ ( 1 - y(2) ) * ( param.C_d  * var.I2 * ( var.xi - param.xi_0 ) );
           param.C_B * Prob * (1-var.B) * var.I2 * ( var.xi - param.xi_0) ];
           %( param.C_d  * var.I2 * ( var.xi - param.xi_0 ) );
           %param.C_B * Prob * ( var.I2 * (mu_new - gamma_new * var.xi + lambda_new/2 * var.xi ^ 2) - var.I2 * (param.a0 + param.a1 * var.xi + param.a2 * var.xi ^ 2 + param.a3 * var.xi ^ 3)) ]; %see ggw183.pdf
end
end