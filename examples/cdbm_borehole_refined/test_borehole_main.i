##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.4 [Check!]
# For Main Fault, 
# mu = shear stress / normal stress = 70e6 / 120e6 = 0.583
# mu_s = 0.677

# under the stress state, 
# -135MPa 70MPa; 70MPa -120e6
# the principal stress is 198e6 

# S = ( mu_s - mu ) / ( mu - mu_d ) = ( 0.677 - 0.583 ) / ( 0.583 - 0.4 ) = 0.5
# Frictional Length Scale L = G Dc / ( ( mu_s - mu_d ) sigma_yy ) = 32.04e9 * 0.4 / (( 0.677 - 0.4) * 198e6) = 233.6m

# Use mesh size = 25m, resolved by more than 7 elements [Check!]

# Use 1% overstress

# Use dt = 2e-4
##########################################################

[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  '../../meshgenerator/cdbm/borehole_refined/borehole_wofaults.msh'
  []
[]

[GlobalParams]
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y'
  
  #damping ratio
  q = 0.2

  ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  lambda_o = 2.15e9
    
  #initial shear modulus value (second lame constant) [Pa]
  shear_modulus_o = 2.73e9

  #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  xi_0 = -0.976

  #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  xi_d = -0.976

  #<strain invariants ratio: maximum allowable value>: set boundary
  #Xu_etal_P15-2D
  #may need a bit space, use 1.5 as boundary
  xi_max = 1.8

  #<strain invariants ratio: minimum allowable value>: set boundary
  #Xu_etal_P15-2D
  xi_min = -1.8

  #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  C_g = 1e-10

  #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  m1 = 10

  #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
  m2 = 1

  ##Compute gamma_damaged_r, xi_1
  #Determine two parameters using convexity of Hessian matrix, positivity of eigenvalues
  #two equations [15a] = 0 [15b] = 0 solves gamma_damaged_r and xi_1 
  #check struct_param.m 

  #coefficient of damage solid modulus
  gamma_damaged_r = 2.99256e9

  #critical point of three phases (strain invariants ratio vs damage)
  xi_1 = 0.74231

  ##Compute parameters in granular states
  #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  #check struct_param.m

  #--------------------------------------------------------------------------------#
  #Note: "computeAlphaCr" needs to change every time the related parameters changed
  #--------------------------------------------------------------------------------#

  # #coefficients
  # chi = 0.75
  a0 = 5.15336e8
  a1 = -1.57828e9
  a2 = 1.44664e9
  a3 = -3.44462e8

  #diffusion coefficient #for structural stress coupling
  D = 0

[]

[AuxVariables]
  [./resid_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./resid_y]
      order = FIRST
      family = LAGRANGE
  []
  [./resid_slipweakening_x]
      order = FIRST
      family = LAGRANGE
  [../]
  [./resid_slipweakening_y]
      order = FIRST
      family = LAGRANGE
  [../]
  [./disp_slipweakening_x]
      order = FIRST
      family = LAGRANGE
  []
  [./disp_slipweakening_y]
      order = FIRST
      family = LAGRANGE
  []
  [./vel_slipweakening_x]
      order = FIRST
      family = LAGRANGE
  []
  [./vel_slipweakening_y]
      order = FIRST
      family = LAGRANGE
  []
  [./mu_s]
      order = CONSTANT
      family = MONOMIAL
  []
  [./mu_d]
      order = CONSTANT
      family = MONOMIAL
  []
  [./ini_shear_stress]
      order = FIRST
      family = LAGRANGE
  []
  [./tangent_jump_rate]
      order = CONSTANT
      family = MONOMIAL
  []
  [./tria_area_aux]
      order = CONSTANT
      family = MONOMIAL
  []
  [./nodal_area]
      order = FIRST
      family = LAGRANGE
  [../]
  #obtain parameters from MaterialRealAux, pass parameters to subApp
  # [./alpha_old]
  #   order = FIRST
  #   family = LAGRANGE
  # []
  [./B_old]
    order = FIRST
    family = LAGRANGE
  []
  [./xi_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [./I2_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [./mu_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [./lambda_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [./gamma_old]
      order = CONSTANT
      family = MONOMIAL
  []
  #updated alpha, B
  [./alpha_in]
    order = CONSTANT
    family = MONOMIAL
  []
  [./B_in]
    order = CONSTANT
    family = MONOMIAL
  []
  #high-order-dummy
  [./alpha_in_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [./B_in_dummy]
    order = FIRST
    family = MONOMIAL
  []
  #grad_alpha
  [./alpha_grad_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [./alpha_grad_y]
    order = CONSTANT
    family = MONOMIAL
  []
  #mechanical strain rate
  [./mechanical_strain_rate]
    order = CONSTANT
    family = MONOMIAL
  []
  #track Cd
  [./track_Cd]
    order = CONSTANT
    family = MONOMIAL
  []
  #pressure
  [pressure]
    order = FIRST
    family = LAGRANGE
  []
  #
  [./initial_stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [./initial_stress_xy]
    order = CONSTANT
    family = MONOMIAL
  []
  [./initial_stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  #
  [arctanyx]
    order = FIRST
    family = LAGRANGE
  []
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = SMALL
        add_variables = true
        planar_formulation = PLANE_STRAIN
        generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_xy strain_yy'
      [../]
    [../]
  [../]
[]

[AuxKernels]
  [Displacment_x]
    type = ProjectionAux
    variable = disp_slipweakening_x
    v = disp_x
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Displacement_y]
    type = ProjectionAux
    variable = disp_slipweakening_y
    v = disp_y
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Vel_x]
      type = CompVarRate
      variable = vel_slipweakening_x
      coupled = disp_x
      execute_on = 'TIMESTEP_BEGIN'
  []
  [Vel_y]
      type = CompVarRate
      variable = vel_slipweakening_y
      coupled = disp_y
      execute_on = 'TIMESTEP_BEGIN'
  []
  [arctanyx]
      type = CompArctanyx
      variable = arctanyx
      execute_on = 'INITIAL'
  []
  [get_xi_old]
      type = MaterialRealAux
      property = xi
      variable = xi_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [get_I2_old]
      type = MaterialRealAux
      property = I2
      variable = I2_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [get_mu_old]
      type = MaterialRealAux
      property = shear_modulus
      variable = mu_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [get_lambda_old]
      type = MaterialRealAux
      property = lambda
      variable = lambda_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [get_gamma_old]
      type = MaterialRealAux
      property = gamma_damaged
      variable = gamma_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    use_displaced_mesh = true
    variable = disp_x
  []
  [./inertia_y]
    type = InertialForce
    use_displaced_mesh = true
    variable = disp_y
  []
  [./Reactionx]
    type = StiffPropDamping
    variable = 'disp_x'
    component = '0'
  []
  [./Reactiony]
    type = StiffPropDamping
    variable = 'disp_y'
    component = '1'
  []
[]

[Materials]
  #damage breakage model
  [stress_medium]
    type = ComputeDamageBreakageStressv4
    alpha_in = alpha_in_dummy
    B_in = B_in_dummy
    alpha_grad_x = alpha_grad_x
    alpha_grad_y = alpha_grad_y
    output_properties = 'eps_p eps_e eps_total I1 sts_total'
    outputs = exodus
  []
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2650
  []
  #young's modulus: 6.67GPa
  #poisson's ratio: 0.22
  [elasticity]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 6.67e9
      poissons_ratio = 0.22
  []
  [./static_initial_stress_tensor]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor
      tensor_functions = 'func_stress_xx              func_stress_xy             func_initial_stress_00 
                          func_stress_xy              func_stress_yy             func_initial_stress_00
                          func_initial_stress_00      func_initial_stress_00     func_initial_stress_00'
  [../]
[]

[Functions]
  #pressure
  [func_pressure_x]
    type = InitialStressXYPressureLateralBC
  []
  [func_pressure_y]
    type = InitialStressXYPressureVerticalBC
  []
  [func_stress_xx]
    type = SolutionFunction
    solution = load_stress_xx
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  [../]
  [func_stress_xy]
    type = SolutionFunction
    solution = load_stress_xy
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  [../]
  [func_stress_yy]
    type = SolutionFunction
    solution = load_stress_yy
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  [../]
  [func_initial_stress_00]
    type = ConstantFunction
    value = 0
  []
[]

[Executioner]
  type = Transient
  dt = 3.0e-8
  end_time = 0.01
  num_steps = 400000
  # num_steps = 10
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

#for cluster run
[Outputs]
  exodus = true
  interval = 1000
  # interval = 1
  [sample_snapshots]
    type = Exodus
    interval = 20000
  []
  [snapshots]
    type = Exodus
    interval = 2000
    overwrite = true
  []
  [checkpoints]
    type = Checkpoint
    interval = 20000
    num_files = 2
  []
[]

[UserObjects]
  [load_stress_xx]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_xx_saved
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_xy]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_xy_saved
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_yy]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_yy_saved
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    timestep = LATEST
    force_preaux = true
  []
[]

[BCs]
  [Pressure_left_x]
    type = Pressure
    boundary = left
    function = func_pressure_x
    variable = disp_x
  []
  [Pressure_right_x]
    type = Pressure
    boundary = right
    function = func_pressure_x
    variable = disp_x
  []
  [Pressure_top_y]
    type = Pressure
    boundary = top
    function = func_pressure_y
    variable = disp_y
  []
  [Pressure_bottom_y]
    type = Pressure
    boundary = bottom
    function = func_pressure_y
    variable = disp_y
  []
[]

[MultiApps]
  [./sub_app]
      type = TransientMultiApp
      positions = '0 0 0'
      input_files = 'test_borehole_sub.i'
      execute_on = 'TIMESTEP_BEGIN'
  [../]
[]

[Transfers]
  [pull_resid]
      type = MultiAppCopyTransfer
      from_multi_app = sub_app
      source_variable = 'alpha_checked B_checked alpha_grad_x_sub alpha_grad_y_sub track_Cd alpha_checked_dummy B_checked_dummy'
      variable = 'alpha_in B_in alpha_grad_x alpha_grad_y track_Cd alpha_in_dummy B_in_dummy'
      execute_on = 'TIMESTEP_BEGIN'
  []
  #we actually don't need to pass alpha and B
  [push_disp]
      type = MultiAppCopyTransfer
      to_multi_app = sub_app
      source_variable = 'alpha_in B_in xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate alpha_in_dummy B_in_dummy'
      variable = 'alpha_old B_old xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate_sub alpha_old_dummy B_old_dummy'
      execute_on = 'TIMESTEP_BEGIN'
  []
[]