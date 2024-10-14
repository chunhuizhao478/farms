function [alpha_cr_out] = comp_alpha_cr(xi_,param)
%Compute alpha_cr based on the current xi
%The two curves represent convexity loss
if xi_ <= param.xi_0 %no convexity loss
    alpha_cr_out = 1.0;
elseif xi_ > param.xi_0 && xi_ <= param.xi_1 %convexity loss 1 %refer to "hessian.m" for detailed compuation
    alpha_cr_out = param.alpha_out2_func(xi_);
elseif xi_ > param.xi_1 && xi_ <= sqrt(3) %convexity loss 2 %refer to "hessian.m" for detailed compuation
    alpha_cr_out = param.alpha_out1_func(xi_);
else %raise error if xi exceeds maximum value (sqrt(3))
    error("xi exceeds the maximum allowable range!!!")
end 
end