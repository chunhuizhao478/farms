##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.4
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

#sample test geometry
[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 800
    ny = 200
    xmin = -10000
    xmax = 10000
    ymin = -2500
    ymax = 2500
    elem_type = TRI3
  []
  [./new_block]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'y<0'
    block_id = 1
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
    split_interface = true
  []
[]
  
[GlobalParams]

  ##----continuum damage breakage model----##

  #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  xi_0 = -0.8

  #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  xi_d = -0.9

  #<strain invariants ratio: maximum allowable value>: set boundary
  #Xu_etal_P15-2D
  #may need a bit space, use 1.5 as boundary
  xi_max = 1.5

  #<strain invariants ratio: minimum allowable value>: set boundary
  #Xu_etal_P15-2D
  xi_min = -1.5

  #<coefficient gives positive damage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #under slow strain rate < low strain rate threshold
  C_d_min = 100

  #power-law correction
  #index
  m = 0.9

  #low strain rate threshold
  mechanical_strain_rate_threshold = 1e-4

  #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
  CdCb_multiplier = 100

  #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
  C_BH = 0

  #<coefficient of healing for damage evolution>: refer to "Lyakhovsky_2011_Hessian_Matrix" Section 3.4
  C_1 = 0

  #<coefficient of healing for damage evolution>: refer to "Lyakhovsky_2011_Hessian_Matrix" Section 3.4
  C_2 = 0.05

  #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  beta_width = 0.03 #1e-3

  #critical point of three phases (strain invariants ratio vs damage)
  xi_1 = 0.8248

  ##Compute parameters in granular states
  #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  #check struct_param.m

  #--------------------------------------------------------------------------------#
  #Note: "computeAlphaCr" needs to change every time the related parameters changed
  #--------------------------------------------------------------------------------#

  #coefficients
  a0 = 4.9526e9
  a1 = -1.8888e10
  a2 = 2.3960e10
  a3 = -1.0112e10

  #diffusion coefficient #self-defined value
  D = 1e3
  
[]

[Variables]
  [alpha_sub]
      order = FIRST
      family = LAGRANGE
  []
[]

[AuxVariables]
  [alpha_old]
      order = FIRST
      family = LAGRANGE
  []
  [B_old]
      order = CONSTANT
      family = MONOMIAL
  []
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
      order = FIRST
      family = LAGRANGE
  []
  [B_checked]
      order = CONSTANT
      family = MONOMIAL
  []
  #
  [B_sub]
    order = CONSTANT
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
      scale = 10
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      mechanical_strain_rate = mechanical_strain_rate_sub
      variable = alpha_sub
  []
[]

[AuxKernels]
  [compute_B]
    type = BreakageVarUpdateDev
    scale = 10
    variable = B_sub
    alpha_old = alpha_old
    B_old = B_old
    xi_old = xi_old
    I2_old = I2_old
    mu_old = mu_old
    gamma_old = gamma_old
    lambda_old = lambda_old
    mechanical_strain_rate = mechanical_strain_rate_sub
    execute_on = 'TIMESTEP_BEGIN'
  []
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
  [alpha_grad_x_sub]
    type = VariableGradientComponent
    variable = alpha_grad_x_sub
    component = x
    gradient_variable = alpha_checked
    execute_on = 'TIMESTEP_END'
 []
 [alpha_grad_y_sub]
   type = VariableGradientComponent
   variable = alpha_grad_y_sub
   component = y
   gradient_variable = alpha_checked
   execute_on = 'TIMESTEP_END'
  []
  #Track Cd
  [track_Cd_aux]
    type = TrackCdEvoAux
    variable = track_Cd
    mechanical_strain_rate = mechanical_strain_rate_sub
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