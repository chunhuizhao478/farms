[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  '../../../meshgenerator/cdbm/borehole_3d/borehole_example.msh'
  []
[]

[GlobalParams]
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y disp_z'
  
  #damping ratio
  # q = 0.2

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
  #
  # [./disp_cdbm_x]
  #     order = FIRST
  #     family = LAGRANGE
  # []
  # [./disp_cdbm_y]
  #     order = FIRST
  #     family = LAGRANGE
  # []
  # [./disp_cdbm_z]
  #     order = FIRST
  #     family = LAGRANGE
  # []
  #
  # [./vel_cdbm_x]
  #     order = FIRST
  #     family = LAGRANGE
  # []
  # [./vel_cdbm_y]
  #     order = FIRST
  #     family = LAGRANGE
  # []
  # [./vel_cdbm_z]
  #     order = FIRST
  #     family = LAGRANGE
  # []
  #
  # [./mu_s]
  #     order = CONSTANT
  #     family = MONOMIAL
  # []
  # [./mu_d]
  #     order = CONSTANT
  #     family = MONOMIAL
  # []
  # [./ini_shear_stress]
  #     order = FIRST
  #     family = LAGRANGE
  # []
  # [./tangent_jump_rate]
  #     order = CONSTANT
  #     family = MONOMIAL
  # []
  # [./tria_area_aux]
  #     order = CONSTANT
  #     family = MONOMIAL
  # []
  # [./nodal_area]
  #     order = FIRST
  #     family = LAGRANGE
  # [../]
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
  [./alpha_grad_z]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = SMALL
        add_variables = true
        # generate_output = 'stress_xx stress_yy stress_xy stress_xz stress_yz stress_zz strain_xx strain_xy strain_yy strain_xz strain_yz strain_zz'
      [../]
    [../]
  [../]
[]

[AuxKernels]
  #
  # [Displacment_x]
  #   type = ProjectionAux
  #   variable = disp_cdbm_x
  #   v = disp_x
  #   execute_on = 'TIMESTEP_BEGIN'
  # []
  # [Displacement_y]
  #   type = ProjectionAux
  #   variable = disp_cdbm_y
  #   v = disp_y
  #   execute_on = 'TIMESTEP_BEGIN'
  # []
  # [Displacement_z]
  #   type = ProjectionAux
  #   variable = disp_cdbm_z
  #   v = disp_z
  #   execute_on = 'TIMESTEP_BEGIN'
  # []
  #
  # [Vel_x]
  #   type = CompVarRate
  #   variable = vel_cdbm_x
  #   coupled = disp_x
  #   execute_on = 'TIMESTEP_BEGIN'
  # []
  # [Vel_y]
  #   type = CompVarRate
  #   variable = vel_cdbm_y
  #   coupled = disp_y
  #   execute_on = 'TIMESTEP_BEGIN'
  # []
  # [Vel_z]
  #   type = CompVarRate
  #   variable = vel_cdbm_z
  #   coupled = disp_z
  #   execute_on = 'TIMESTEP_BEGIN'
  # []
  #
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
  [./inertia_z]
    type = InertialForce
    use_displaced_mesh = true
    variable = disp_z
  []
  # [./Reactionx]
  #   type = StiffPropDamping
  #   variable = 'disp_x'
  #   component = '0'
  # []
  # [./Reactiony]
  #   type = StiffPropDamping
  #   variable = 'disp_y'
  #   component = '1'
  # []
  # [./Reactionz]
  #   type = StiffPropDamping
  #   variable = 'disp_z'
  #   component = '2'
  # []
[]

[Materials]
  #damage breakage model
  [stress_medium]
    type = ComputeDamageBreakageStress3D
    alpha_in = alpha_in_dummy
    B_in = B_in_dummy
    alpha_grad_x = alpha_grad_x
    alpha_grad_y = alpha_grad_y
    alpha_grad_z = alpha_grad_z
    # output_properties = 'eps_p eps_e eps_total I1 sts_total'
    outputs = nemesis
  []
  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2650
  []
  [./static_initial_stress_tensor]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor
      tensor_functions = 'func_stress_xx        func_stress_xy     func_stress_xz 
                          func_stress_xy        func_stress_yy     func_stress_yz
                          func_stress_xz        func_stress_yz     func_stress_zz'
  [../]
  [./static_initial_strain_tensor]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_strain_tensor
      tensor_functions = 'func_strain_xx        func_strain_xy     func_strain_xz 
                          func_strain_xy        func_strain_yy     func_strain_yz
                          func_strain_xz        func_strain_yz     func_strain_zz'
  [../]
[]

[Functions]
  #func_stress
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
  [func_stress_xz]
    type = SolutionFunction
    solution = load_stress_xz
    execute_on = 'INITIAL'
  [../]
  [func_stress_yz]
    type = SolutionFunction
    solution = load_stress_yz
    execute_on = 'INITIAL'
  [../]
  [func_stress_zz]
    type = SolutionFunction
    solution = load_stress_zz
    execute_on = 'INITIAL'
  [../]
  #func_strain
  [func_strain_xx]
    type = SolutionFunction
    solution = load_strain_xx
    execute_on = 'INITIAL'
  [../]
  [func_strain_xy]
    type = SolutionFunction
    solution = load_strain_xy
    execute_on = 'INITIAL'
  [../]
  [func_strain_yy]
    type = SolutionFunction
    solution = load_strain_yy
    execute_on = 'INITIAL'
  [../]
  [func_strain_xz]
    type = SolutionFunction
    solution = load_strain_xz
    execute_on = 'INITIAL'
  [../]
  [func_strain_yz]
    type = SolutionFunction
    solution = load_strain_yz
    execute_on = 'INITIAL'
  [../]
  [func_strain_zz]
    type = SolutionFunction
    solution = load_strain_zz
    execute_on = 'INITIAL'
  [../]
[]

[Executioner]
  type = Transient
  dt = 5.9e-8
  end_time = 1
  num_steps = 10000
  # num_steps = 5
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
    use_constant_mass = true
  []
[]

#for cluster run
[Outputs]
  interval = 100
  nemesis = true
  exodus = false
  show = 'disp_x disp_y disp_z alpha_in B_in mu_old'
[]

[UserObjects]
  #load_stress
  [load_stress_xx]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_xx_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_xy]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_xy_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_yy]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_yy_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_xz]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_xz_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_yz]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_yz_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_zz]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_zz_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  #load_strain
  [load_strain_xx]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = strain_xx_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_strain_xy]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = strain_xy_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_strain_yy]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = strain_yy_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_strain_xz]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = strain_xz_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_strain_yz]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = strain_yz_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
  []
  [load_strain_zz]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = strain_zz_saved
    execute_on = 'INITIAL'
    timestep = LATEST
    force_preaux = true
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
      source_variable = 'alpha_checked B_checked alpha_grad_x_sub alpha_grad_y_sub alpha_grad_z_sub alpha_checked_dummy B_checked_dummy'
      variable = 'alpha_in B_in alpha_grad_x alpha_grad_y alpha_grad_z alpha_in_dummy B_in_dummy'
      execute_on = 'TIMESTEP_BEGIN'
  []
  #we actually don't need to pass alpha and B
  [push_disp]
      type = MultiAppCopyTransfer
      to_multi_app = sub_app
      source_variable = 'alpha_in B_in xi_old I2_old mu_old lambda_old gamma_old alpha_in_dummy B_in_dummy'
      variable = 'alpha_old B_old xi_old I2_old mu_old lambda_old gamma_old alpha_old_dummy B_old_dummy'
      execute_on = 'TIMESTEP_BEGIN'
  []
[]