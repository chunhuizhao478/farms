function [var] = check_alphaB(alpha_new,B_new,var)
%Check whether damage and breakage parameters are within allowable range
%This seems necessary for eliminating numerical error.
%%convexity%%
%--------------------%
% alpha_size = size(alpha_new);
% B_size = size(B_new);
% if alpha_size(1) ~= 1 %multiple roots, convexity loss 
%     alpha_new = max(alpha_new);
% end
% if B_size(1) ~= 1 %multiple roots, convexity loss 
%     B_new = max(B_new);
% end
%--------------------%
%%alpha%%
%--------------------%
if alpha_new < 0
    alpha_new_out = 0;
elseif alpha_new > 1
    alpha_new_out = 1;
else
    alpha_new_out = alpha_new;
end
%--------------------%
%%B%%
%--------------------%
if B_new < 0
    B_new_out = 0;
elseif B_new > 1
    B_new_out = 1;
else
    B_new_out = B_new;
end
%--------------------%
var.B = B_new_out;
var.alpha_ = alpha_new_out;