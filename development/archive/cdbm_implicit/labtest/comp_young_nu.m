function [var] = comp_young_nu(param,var,case_index)
%Compute young's modulus and poisson ratio
%input: var output: update in var
if case_index == "init"
    lambda_cbdm = param.lambda_0;
    mu_cbdm = param.mu_0;
elseif case_index == "update"
    lambda_cbdm = var.lambda_ - var.gamma_ / var.xi;
    mu_cbdm = var.mu_ - ( var.gamma_ * var.xi ) / 2;
else
    error("Provide a valid case_index !");
end
    var.E = mu_cbdm * ( 3 * lambda_cbdm + 2 * mu_cbdm ) / ( lambda_cbdm + mu_cbdm );
    var.nu = lambda_cbdm / ( 2 * ( lambda_cbdm + mu_cbdm ));

%     var.E = var.mu_ * ( 3 * lambda_cbdm + 2 * var.mu_ ) / ( var.lambda_ + var.mu_ );
%     var.nu =  var.lambda_ / ( 2 * (  var.lambda_ + var.mu_ ));
end