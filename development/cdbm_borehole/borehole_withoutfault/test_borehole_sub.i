##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.1
# For Main Fault, 
# mu = shear stress / normal stress = 70e6 / 120e6 = 0.583
# mu_s = 0.677
# S = ( mu_s - mu ) / ( mu - mu_d ) = ( 0.677 - 0.583 ) / ( 0.583 - 0.4 ) = 0.514
# Frictional Length Scale L = G Dc / ( ( mu_s - mu_d ) sigma_yy ) = 32.04e9 * 0.4 / (( 0.677 - 0.1) * 120e6) = 185m
# Use mesh size = 25m

# Diffusion Length Scale D = 5e5
# sqrt(5e5*185/3464) = 163. using 6~7, 25m mesh to resolve it

# Check CFL condition 0.1 * 25 / 6000 ~ 0.0001s
##########################################################

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../../meshgenerator/cdbm/borehole_wofaults/borehole_wofaults.msh'
  []
[]

[GlobalParams]
  
  ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  lambda_o = 3.204e10
    
  #initial shear modulus value (second lame constant) [Pa]
  shear_modulus_o = 3.204e10

  #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  xi_0 = -1.0

  #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  xi_d = -1.1

  #<strain invariants ratio: maximum allowable value>: set boundary
  #Xu_etal_P15-2D
  #may need a bit space, use 1.5 as boundary
  xi_max = 1.8

  #<strain invariants ratio: minimum allowable value>: set boundary
  #Xu_etal_P15-2D
  xi_min = -1.8

  #<coefficient gives positive damage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #under slow strain rate < low strain rate threshold
  # C_d_min = 10

  #if option 2, use Cd_constant
  Cd_constant = 1e4

  #power-law correction
  #index
  m = 0.9

  #low strain rate threshold
  # mechanical_strain_rate_threshold = -1e-4

  #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
  CdCb_multiplier = 10 

  #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
  CBCBH_multiplier = 10

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_1 = 300

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_2 = 0.05

  #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  beta_width = 0.03 #1e-3

  #critical point of three phases (strain invariants ratio vs damage)
  xi_1 = 0.73205

  ##Compute parameters in granular states
  #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  #check struct_param.m

  #--------------------------------------------------------------------------------#
  #Note: "computeAlphaCr" needs to change every time the related parameters changed
  #--------------------------------------------------------------------------------#

  #coefficients
  # chi = 0.75
  a0 = 3.92486e9
  a1 = -1.51257e10
  a2 = 1.93525e10
  a3 = -8.21570e9

  #diffusion coefficient #for structural stress coupling
  D = 0
  
[]

[Variables]
  [alpha_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  [B_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  #-high-order-dummy-#
  [alpha_sub_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [B_sub_dummy]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxVariables]
  [alpha_old]
    order = CONSTANT
    family = MONOMIAL
  []
  [B_old]
    order = CONSTANT
    family = MONOMIAL
  []
  #-high-order-dummy-#
  [alpha_old_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [B_old_dummy]
    order = FIRST
    family = MONOMIAL
  []
  #
  [xi_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [I2_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [mu_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [lambda_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [gamma_old]
      order = CONSTANT
      family = MONOMIAL
  []
  #checked
  [alpha_checked]
    order = CONSTANT
    family = MONOMIAL
  []
  [B_checked]
    order = CONSTANT
    family = MONOMIAL
  []
  #high-order-dummy
  [alpha_checked_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [B_checked_dummy]
    order = FIRST
    family = MONOMIAL
  []
  #grad_alpha
  [alpha_grad_x_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  [alpha_grad_y_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  #mechanical strain rate sub
  [mechanical_strain_rate_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  #Track Cd
  [track_Cd]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [./timederivative_alpha]
      type = TimeDerivative
      variable = alpha_sub
  []
  [./alpha_forcing_func]
      type = DamageVarForcingFuncDev
      option = 2
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      variable = alpha_sub
      healing = true
  []
  [./timederivative_B]
    type = TimeDerivative
    variable = B_sub
  []
  [./B_forcing_func]
      type = BreakageVarForcingFuncDevOld
      option = 2
      variable = B_sub
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      mu_old = mu_old
      gamma_old = gamma_old
      lambda_old = lambda_old
      healing = true
  []
  #-high-order-dummy-#
  [./timederivative_alpha_dummy]
    type = TimeDerivative
    variable = alpha_sub_dummy
  []
  [./alpha_forcing_func_dummy]
      type = DamageVarForcingFuncDev
      option = 2
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      variable = alpha_sub_dummy
      healing = true
  []
  [./timederivative_B_dummy]
    type = TimeDerivative
    variable = B_sub_dummy
  []
  [./B_forcing_func_dummy]
      type = BreakageVarForcingFuncDevOld
      option = 2
      variable = B_sub_dummy
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      mu_old = mu_old
      gamma_old = gamma_old
      lambda_old = lambda_old
      healing = true
  []
[]

[AuxKernels]
  [check_alpha]
    type = CheckAlphaB
    coupled = alpha_sub
    variable = alpha_checked
    execute_on = 'TIMESTEP_END'
  []
  [check_B]
    type = CheckAlphaB
    coupled = B_sub
    variable = B_checked
    execute_on = 'TIMESTEP_END'
  []
  #-high-order-dummy-#
  [check_alpha_dummy]
    type = CheckAlphaB
    coupled = alpha_sub_dummy
    variable = alpha_checked_dummy
    execute_on = 'TIMESTEP_END'
  []
  [check_B_dummy]
    type = CheckAlphaB
    coupled = B_sub_dummy
    variable = B_checked_dummy
    execute_on = 'TIMESTEP_END'
  []
  #
  [alpha_grad_x_sub]
    type = VariableGradientComponent
    variable = alpha_grad_x_sub
    component = x
    gradient_variable = alpha_checked_dummy
    execute_on = 'TIMESTEP_END'
 []
 [alpha_grad_y_sub]
   type = VariableGradientComponent
   variable = alpha_grad_y_sub
   component = y
   gradient_variable = alpha_checked_dummy
   execute_on = 'TIMESTEP_END'
  []
[]

#by default, subApp is using the same time step as mainApp
[Executioner]
  type = Transient
  [TimeIntegrator]
    type = ExplicitSSPRungeKutta
    order = 3
  []
[]