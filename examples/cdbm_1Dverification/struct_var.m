function var = struct_var()
%1d_continuum_damage_breakage_model%
%%material_var%%
var.strain_rate = -1e-4;%-1e-4         %[prescribed] strain rate
var.t0 = 0;                            %[prescribed] initial time
var.dt = 1e-1;                         %[prescribed] time step
var.depsx = var.strain_rate * var.dt;  %[prescribed] strain increment

var.t = 0;                             %[variable] current time
var.Tmax = 1e5;                       %[variable] maximum time
var.alpha_ = 0;                        %[variable] damage parameter
var.B = 0;                             %[variable] breakage parameter
var.E = 0;                             %[variable] young's modulus
var.nu = 0;                            %[variable] poission's ratio

%%elastic strain
var.xi = 0;                            %[variable] strain invariants ratio
var.I1 = 0;                            %[variable] 1st invariant
var.I2 = 0;                            %[variable] 2nd invariant
%%total strain
%var.xit = 0;                           %[variable] strain invariants ratio
%var.I1t = 0;                           %[variable] 1st invariant
%var.I2t = 0;                           %[variable] 2nd invariant

var.lambda_ = 0;                       %[variable] 1st lame constant 
var.mu_ = 0;                           %[variable] 2nd lame constant (shear modulus)
var.gamma_ = 0;                        %[variable] damaged modulus

var.Prob = 0;                          %[variable] probability for material be in granular state                           

var.eps = zeros(3,3);                  %[tensor] total strain tensor
var.eps_e = zeros(3,3);                %[tensor] elastic strain tensor
var.eps_p = zeros(3,3);                %[tensor] plastic strain tensor
var.sigma_s = zeros(3,3);              %[tensor] stress in solid state
var.sigma_b = zeros(3,3);              %[tensor] stress in granular state
var.sigma = zeros(3,3);                %[tensor] total stress 

var.alpha_list = [];                   %[list] store damage variable
var.xi_list = [];                      %[list] store strain invariant ratio
var.time_list = []; 

var.sigma_x_list = [];                      %[list] store stress in x direction
var.sigma_y_list = [];
var.sigmad_x_list = [];

var.eps_x_list = [];                   %[list] store strain in x direction
var.eps_y_list = [];                   %[list] store time
var.epse_x_list = [];                  %[list] elastic strain in x direction
var.epse_y_list = [];                  %[list] elastic strain in y direction
var.epsp_x_list = [];                  %[list] plastic strain in x direction
var.epsp_y_list = [];                  %[list] plastic strain in y direction
var.sigmas_x_list = [];
var.sigmas_y_list = [];
var.sigmab_x_list = [];
var.sigmab_y_list = [];
var.B_list = [];

end