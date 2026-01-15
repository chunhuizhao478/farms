[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = '../../meshgenerator/tpv1013d/tpv1013d_100m.msh'
  []
  [./new_block_1]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'x >= -18000 & x <= 18000 & z >= -18000 & y < 0'
    block_id = 100
  []
  [./new_block_2]
    type = ParsedSubdomainMeshGenerator
    input = new_block_1
    combinatorial_geometry = 'x >= -18000 & x <= 18000 & z >= -18000 & y > 0'
    block_id = 200
  []
  [./split_1]
    type = BreakMeshByBlockGenerator
    input = new_block_2
    split_interface = true
    block_pairs = '100 200'
  []
  [./sidesets]
    input = split_1
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0
                0 0 -1
                0 0 1'
    new_boundary = 'left right bottom top back front'
  []
[]

[GlobalParams]
  #primary variables
  displacements = 'disp_x disp_y disp_z'
  q = 0.05
  ##rate-and-state coefficients
  f_o = 0.6
  rsf_a = 0.008
  rsf_b = 0.012
  rsf_L = 0.02
  delta_o = 1e-6

  ##initial normal traction (Pa)
  T2_o = 120e6

  ##initial strike shear traction (Pa)
  T1_o = 75e6

  ##initial dip shear traction (Pa)
  T3_o = 0

  ##initial sliprate (m/s)
  sliprate_strike_init = 1e-12

  ##initial state variable
  statevar_init = 1.606238999213454e9

  #element edge length (m)
  len = 200
[]

[AuxVariables]
  [resid_x]
    order = FIRST
    family = LAGRANGE
  []
  [resid_y]
    order = FIRST
    family = LAGRANGE
  []
  [resid_z]
    order = FIRST
    family = LAGRANGE
  []
  [resid_ratestate_x]
      order = FIRST
      family = LAGRANGE
  []
  [resid_ratestate_y]
      order = FIRST
      family = LAGRANGE
  []
  [resid_ratestate_z]
      order = FIRST
      family = LAGRANGE
  []
  [disp_ratestate_x]
      order = FIRST
      family = LAGRANGE
  []
  [disp_ratestate_y]
      order = FIRST
      family = LAGRANGE
  []
  [disp_ratestate_z]
    order = FIRST
    family = LAGRANGE
  []
  [vel_ratestate_x]
    order = FIRST
    family = LAGRANGE
  []
  [vel_ratestate_y]
      order = FIRST
      family = LAGRANGE
  []
  [vel_ratestate_z]
    order = FIRST
    family = LAGRANGE
  []
  [Ts_perturb]
    order = FIRST
    family = MONOMIAL
  []
  ###
  [jump_x_aux]
    order = FIRST
    family = MONOMIAL
  []
  [jump_x_rate_aux]
    order = FIRST
    family = MONOMIAL
  []
  [traction_x_aux]
    order = FIRST
    family = MONOMIAL
  []
  ###
  [statevar_aux]
    order = FIRST
    family = MONOMIAL
  []
[]

[Physics/SolidMechanics/CohesiveZone]
  [./czm_ik]
    boundary = 'Block100_Block200'
    strain = SMALL
    generate_output='traction_x traction_y traction_z jump_x jump_y jump_z'
  [../]
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        generate_output = 'stress_xx stress_yy stress_xy'
        extra_vector_tags = 'restore_tag'
      []
    []
  []
[]

[Problem]
  extra_tag_vectors = 'restore_tag'
[]

[Functions]
  [func_initial_strike_shear_stress]
    type = InitialStrikeShearStressPerturbRSF3D
  []
  #[func_initial_strike_shear_stress]
  #  type = ConstantFunction
  #  value = 0
  #[]
[]

[AuxKernels]
  [Displacment_x]
    type = ProjectionAux
    variable = disp_ratestate_x
    v = disp_x
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Displacement_y]
    type = ProjectionAux
    variable = disp_ratestate_y
    v = disp_y
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Displacement_z]
    type = ProjectionAux
    variable = disp_ratestate_z
    v = disp_z
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Vel_x]
    type = CompVarRate
    variable = vel_ratestate_x
    coupled = disp_x
    execute_on = 'TIMESTEP_END'
  []
  [Vel_y]
    type = CompVarRate
    variable = vel_ratestate_y
    coupled = disp_y
    execute_on = 'TIMESTEP_END'
  []
  [Vel_z]
    type = CompVarRate
    variable = vel_ratestate_z
    coupled = disp_z
    execute_on = 'TIMESTEP_END'
  []
  [Residual_x]
    type = ProjectionAux
    variable = resid_ratestate_x
    v = resid_x
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Residual_y]
    type = ProjectionAux
    variable = resid_ratestate_y
    v = resid_y
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Residual_z]
    type = ProjectionAux
    variable = resid_ratestate_z
    v = resid_z
    execute_on = 'TIMESTEP_BEGIN'
  []
  [restore_x]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_x'
    variable = 'resid_x'
  []
  [restore_y]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_y'
    variable = 'resid_y'
  []
  [restore_z]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_z'
    variable = 'resid_z'
  []
  ##
  [StrikeShearStress]
    type = FunctionAux
    variable = Ts_perturb
    function = func_initial_strike_shear_stress
    execute_on = 'TIMESTEP_BEGIN'
  []
  ##
  [get_jump_x_aux]
    type = MaterialRealAux
    property = jump_x
    variable = jump_x_aux
    boundary = 'Block100_Block200'
  []
  [get_jump_x_rate_aux]
    type = FDCompVarRate
    variable = jump_x_rate_aux
    coupled = jump_x
    execute_on = 'TIMESTEP_END'
    boundary = 'Block100_Block200'
  []
  [get_traction_x_aux]
    type = MaterialRealAux
    property = traction_x
    variable = traction_x_aux
    boundary = 'Block100_Block200'
  []
  #
  [get_statevar_aux]
    type = MaterialRealAux
    property = statevar
    variable = statevar_aux
    boundary = 'Block100_Block200'
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
  [./inertia_z]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_z
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
  [./Reactionz]
    type = StiffPropDamping
    variable = 'disp_z'
    component = '2'
  []
[]

[Materials]
  [elasticity]
      type = ComputeIsotropicElasticityTensor
      lambda = 32.04e9
      shear_modulus = 32.04e9
      use_displaced_mesh = false
  []
  [stress]
      type = ComputeLinearElasticStress
  []
  [density]
      type = GenericConstantMaterial
      prop_names = density
      prop_values = 2670
  []
  [./czm_mat]
      type = RateStateFrictionczm3d
      disp_x     = disp_ratestate_x
      disp_y     = disp_ratestate_y
      disp_z     = disp_ratestate_z
      vel_x      = vel_ratestate_x
      vel_y      = vel_ratestate_y
      vel_z      = vel_ratestate_z
      reaction_x = resid_ratestate_x
      reaction_y = resid_ratestate_y
      reaction_z = resid_ratestate_z
      Ts_perturb = Ts_perturb
      boundary = 'Block100_Block200'
  [../]
[]

[UserObjects]
  [recompute_residual_tag]
      type = ResidualEvaluationUserObject
      vector_tag = 'restore_tag'
      force_preaux = true
      execute_on = 'TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  dt = 0.0025
  end_time = 3.0 #for testing
  # num_steps = 10
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

[Outputs]
  exodus = true
  show = 'vel_ratestate_x vel_ratestate_y vel_ratestate_z disp_ratestate_x disp_ratestate_y disp_ratestate_z Ts_perturb statevar_aux'
  time_step_interval = 40
  [csv]
    type = CSV
    execute_on = 'timestep_end'
    time_step_interval = 20
  []
  [out]
    type = Checkpoint
    time_step_interval = 160
    num_files = 2
  []
[]

[BCs]
[]

[VectorPostprocessors]
[]
