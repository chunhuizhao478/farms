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
      file =  '../../../meshgenerator/cdbm/borehole_wofaults/borehole_wofaults.msh'
  []
[]

[GlobalParams]
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y'
  
  #damping ratio
  q = 0.2
  
  #characteristic length (m)
  # Dc = 0.4

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
  gamma_damaged_r = 3.49061e10

  #critical point of three phases (strain invariants ratio vs damage)
  xi_1 = 0.73205

  ##Compute parameters in granular states
  #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  #check struct_param.m

  #--------------------------------------------------------------------------------#
  #Note: "computeAlphaCr" needs to change every time the related parameters changed
  #--------------------------------------------------------------------------------#

  # #coefficients
  # chi = 0.75
  a0 = 3.92486e9
  a1 = -1.51257e10
  a2 = 1.93525e10
  a3 = -8.21570e9

  #diffusion coefficient #for structural stress coupling
  D = 0

  #pressure analytical solution
  #Reference: Injection-induced seismicity: Poroelastic and earthquake nucleation effects (P. Segall1 and S. Lu2)
  # effec_sts_coeff = 0.31 #
  # flux_q = 5e2 #kg/s
  # density_rho_0 = 1e3 #kg/m^3
  # permeability_k = 3e-12 #m^2
  # viscosity_eta = 0.4e-3 #Pa s
  biotcoeff_alpha = 0.31 #-
  # undrained_nu_u = 0.3  #-
  # shear_modulus_mu = 32.04e9 #Pa
  # drained_nu = 0.25 #-
  
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
[]

# [Modules/TensorMechanics/CohesiveZoneMaster]
#   [./czm_ik]
#     boundary = 'czm'
#     strain = SMALL
#     generate_output='tangent_jump normal_jump'
#   [../]
# []


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
  #define shear strain material property (elastic) inside damage stress 
  #and compute its rate using "MaterialRateRealAux"
  [get_shear_strain_rate]
      type = MaterialRateRealAux
      property = principal_strain
      variable = mechanical_strain_rate
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  #fault length
  [fault_len]
      type = ConstantAux
      variable = nodal_area
      value = 25
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  #load steady solve stress as initial condition
  [initial_stress_xx]
    type = SolutionAux
    variable = initial_stress_xx
    solution = load_stress_xx
    execute_on = 'INITIAL'
  []
  [initial_stress_xy]
    type = SolutionAux
    variable = initial_stress_xy
    solution = load_stress_xy
    execute_on = 'INITIAL'
  []
  [initial_stress_yy]
    type = SolutionAux
    variable = initial_stress_yy
    solution = load_stress_yy
    execute_on = 'INITIAL'
  []
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_x
  []
  [./inertia_y]
    type = InertialForce
    use_displaced_mesh = false
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
    type = ComputeDamageBreakageStressBorehole
    alpha_in = alpha_in_dummy
    B_in = B_in_dummy
    alpha_grad_x = alpha_grad_x
    alpha_grad_y = alpha_grad_y
    pressure = pressure
    initial_stress_xx = initial_stress_xx
    initial_stress_xy = initial_stress_xy
    initial_stress_yy = initial_stress_yy
    output_properties = 'eps_p eps_e eps_total I1 sts_total'
    outputs = exodus
  []
  [density]
      type = GenericConstantMaterial
      prop_names = density
      prop_values = 2670
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
    type = InitialStressXYPressureBorehole_x_fast
  []
  [func_pressure_y]
    type = InitialStressXYPressureBorehole_y_fast
  []
  [func_stress_xx]
    type = SolutionFunction
    solution = load_stress_xx
    execute_on = 'INITIAL'
  [../]
  [func_stress_xy]
    type = SolutionFunction
    solution = load_stress_xy
    execute_on = 'INITIAL'
  [../]
  [func_stress_yy]
    type = SolutionFunction
    solution = load_stress_yy
    execute_on = 'INITIAL'
  [../]
  [func_initial_stress_00]
    type = ConstantFunction
    value = 0
  []
[]

[UserObjects]
  [load_stress_xx]
    type = SolutionUserObject
    mesh = ../borehole_staticsolve/test_borehole_main_out.e
    system_variables = stress_xx_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_xy]
    type = SolutionUserObject
    mesh = ../borehole_staticsolve/test_borehole_main_out.e
    system_variables = stress_xy_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_yy]
    type = SolutionUserObject
    mesh = ../borehole_staticsolve/test_borehole_main_out.e
    system_variables = stress_yy_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
[]

[Executioner]
  type = Transient
  dt = 2e-4
  end_time = 40.0
  # num_steps = 10
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

#for cluster run
[Outputs]
  exodus = true
  interval = 200
  [sample_snapshots]
    type = Exodus
    interval = 2000
  []
  [snapshots]
    type = Exodus
    interval = 2000
    overwrite = true
  []
  [checkpoints]
    type = Checkpoint
    interval = 2000
    num_files = 2
  []
[]

[BCs]
  [fix_top_x]
    type = DirichletBC
    boundary = top
    variable = disp_x
    value = 0
  []
  [fix_top_y]
    type = DirichletBC
    boundary = top
    variable = disp_y
    value = 0
  []
  [fix_bottom_x]
    type = DirichletBC
    boundary = bottom
    variable = disp_x
    value = 0
  []
  [fix_bottom_y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0
  []
  [fix_left_x]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0
  []
  [fix_left_y]
    type = DirichletBC
    boundary = left
    variable = disp_y
    value = 0
  []
  [fix_right_x]
    type = DirichletBC
    boundary = right
    variable = disp_x
    value = 0
  []
  [fix_right_y]
    type = DirichletBC
    boundary = right
    variable = disp_y
    value = 0
  []
  [./Pressure]
    [./borehole_x]
      boundary = borehole
      function = func_pressure_x
      displacements = 'disp_x disp_y'
    []
    [./borehole_y]
      boundary = borehole
      function = func_pressure_y
      displacements = 'disp_x disp_y'
    []
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