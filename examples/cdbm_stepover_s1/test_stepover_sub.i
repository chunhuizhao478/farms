[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../meshgenerator/cdbm/stepover/step_over_s.msh'
  []
  [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'x >= -1500 & x<= 3000 & y<=100 & y>=0'
      block_id = 1
  []
  [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'x >= -1500 & x<= 3000 & y<=0 & y>=-100'
      block_id = 2
  []
  [./new_block_3]
      type = ParsedSubdomainMeshGenerator
      input = new_block_2
      combinatorial_geometry = 'x >= 1500 & x<= 6000 & y<=-900 & y>=-1000'
      block_id = 3
  []
  [./new_block_4]
      type = ParsedSubdomainMeshGenerator
      input = new_block_3
      combinatorial_geometry = 'x >= 1500 & x<= 6000 & y<=-1000 & y>=-1100'
      block_id = 4
  []
  [./split]
      type = BreakMeshByBlockGenerator
      input = new_block_4
      split_interface = true
      block_pairs = '1 2; 3 4'
  []
  [./sidesets]
    input = split
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0'
    new_boundary = 'left right bottom top'
  []
[]

[GlobalParams]

  ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  lambda_o = 3.204e10
    
  #initial shear modulus value (second lame constant) [Pa]
  shear_modulus_o = 3.204e10

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
  C_d_min = 10

  #power-law correction
  #index
  m = 0.9

  #low strain rate threshold
  mechanical_strain_rate_threshold = 1e-4

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
  xi_1 = 0.8248

  ##Compute parameters in granular states
  #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  #check struct_param.m

  #--------------------------------------------------------------------------------#
  #Note: "computeAlphaCr" needs to change every time the related parameters changed
  #--------------------------------------------------------------------------------#

  #coefficients
  # a0 = 4.9526e9
  # a1 = -1.8888e10
  # a2 = 2.3960e10
  # a3 = -1.0112e10

  #chi = 0.95
  # a0 = 9.410e9
  # a1 = -2.474e10
  # a2 = 1.850e10
  # a3 = -2.832e9

  #ggw183.pdf
  a0 = 7.414e9
  a1 = -2.1344e10
  a2 = 1.9058e10
  a3 = -4.946e9

  #diffusion coefficient #for structural stress coupling
  D = 1e3
  
[]

[Variables]
  [alpha_sub]
    order = FIRST
    family = LAGRANGE
  []
  [B_sub]
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
    order = FIRST
    family = LAGRANGE
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
    order = FIRST
    family = LAGRANGE
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
      option = 1
      scale = 10
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      mechanical_strain_rate = mechanical_strain_rate_sub
      variable = alpha_sub
      healing = true
  []
  [./timederivative_B]
    type = TimeDerivative
    variable = B_sub
  []
  [./B_forcing_func]
      type = BreakageVarForcingFuncDev
      option = 1
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
      healing = false
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