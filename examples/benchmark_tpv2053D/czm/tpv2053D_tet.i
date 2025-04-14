# Verification of Benchmark Problem TPV14-2D from the SCEC Dynamic Rupture Validation exercises #
# Reference: #
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126. #
# [Note]: This serves as a test file, to run the full problem, please extend the domain size by modifying nx, ny, xmin, xmax, ymin, ymax

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = '../../../meshgenerator/tpv205/tpv2053d_200m.msh'
  []
  [./new_block_1]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'x >= -15000 & x <= 15000 & z >= -15000 & y < 0'
    block_id = 100
  []
  [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'x >= -15000 & x <= 15000 & z >= -15000 & y > 0'
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
  #damping ratio
  q = 0.1
  #characteristic length (m)
  Dc = 0.4
  #initial normal stress (Pa)
  T2_o = 120e6
  #initial shear stress along dip direction (Pa)
  T3_o = 0.0
  #dynamic friction coefficient
  mu_d = 0.525
  #element edge length (m)
  len = 200
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
  [./resid_z]
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
  [./resid_slipweakening_z]
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
  [./disp_slipweakening_z]
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
  [./vel_slipweakening_z]
    order = FIRST
    family = LAGRANGE
  []
  [./mu_s]
      order = CONSTANT
      family = MONOMIAL
  []
  [./ini_shear_stress]
    order = FIRST
    family = MONOMIAL
  []
  [./ini_normal_stress]
    order = FIRST
    family = MONOMIAL
  []
  #global quantities
  [./global_jump_x]
    order = FIRST
    family = MONOMIAL
  []
  [./global_jump_y]
    order = FIRST
    family = MONOMIAL
  []
  [./global_jump_z]
    order = FIRST
    family = MONOMIAL
  []
  [./global_traction_x]
    order = FIRST
    family = MONOMIAL
  []
  [./global_traction_y]
    order = FIRST
    family = MONOMIAL
  []
  [./global_traction_z]
    order = FIRST
    family = MONOMIAL
  []
  #local quantities
  [./local_shear_jump]
      order = FIRST
      family = MONOMIAL
  []
  [./local_shear_jump_rate]
      order = FIRST
      family = MONOMIAL
  []
  [./local_shear_traction]
    order = FIRST
    family = MONOMIAL
  []
  #
  [./local_normal_jump]
    order = FIRST
    family = MONOMIAL
  []
  [./local_normal_jump_rate]
    order = FIRST
    family = MONOMIAL
  []
  [./local_normal_traction]
    order = FIRST
    family = MONOMIAL
  []
  #
  [./local_dip_jump]
    order = FIRST
    family = MONOMIAL
  []
  [./local_dip_jump_rate]
    order = FIRST
    family = MONOMIAL
  []
  [./local_dip_traction]
    order = FIRST
    family = MONOMIAL
  []
  #
  [./normal_x]
    order = FIRST
    family = MONOMIAL
  []
  [./normal_y]
    order = FIRST
    family = MONOMIAL
  []
  [./normal_z]
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
  [func_static_friction_coeff_mus]
    type = InitialStaticFrictionCoeff
  []
  #the center is (-8000,0)
  #expand 1.5e3 on both side (-9500,0) and (-6500,0)
  [func_initial_strike_shear_stress]
    type = InitialShearStressTPV2053d
  []
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
  [Displacement_z]
    type = ProjectionAux
    variable = disp_slipweakening_z
    v = disp_z
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Vel_x]
    type = CompVarRate
    variable = vel_slipweakening_x
    coupled = disp_x
    execute_on = 'TIMESTEP_END'
  []
  [Vel_y]
    type = CompVarRate
    variable = vel_slipweakening_y
    coupled = disp_y
    execute_on = 'TIMESTEP_END'
  []
  [Vel_z]
    type = CompVarRate
    variable = vel_slipweakening_z
    coupled = disp_z
    execute_on = 'TIMESTEP_END'
  []
  [Residual_x]
    type = ProjectionAux
    variable = resid_slipweakening_x
    v = resid_x
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Residual_y]
    type = ProjectionAux
    variable = resid_slipweakening_y
    v = resid_y
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Residual_z]
    type = ProjectionAux
    variable = resid_slipweakening_z
    v = resid_z
    execute_on = 'TIMESTEP_BEGIN'
  []
  [restore_x]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_x'
    variable = 'resid_x'
    execute_on = 'TIMESTEP_END'
  []
  [restore_y]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_y'
    variable = 'resid_y'
    execute_on = 'TIMESTEP_END'
  []
  [restore_z]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_z'
    variable = 'resid_z'
    execute_on = 'TIMESTEP_END'
  []
  [StaticFricCoeff]
    type = FunctionAux
    variable = mu_s
    function = func_static_friction_coeff_mus
    execute_on = 'LINEAR TIMESTEP_BEGIN'
  []
  ##
  [StrikeShearStress]
    type = FunctionAux
    variable = ini_shear_stress
    function = func_initial_strike_shear_stress
    execute_on = 'LINEAR TIMESTEP_BEGIN'
  []
  ##
  [GetShearJump]
    type = FarmsMaterialRealAux3D
    material_property_name = 'local_shear_jump'
    variable = local_shear_jump
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetNormalJump]
    type = FarmsMaterialRealAux3D
    material_property_name = 'local_normal_jump'
    variable = local_normal_jump
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetDipJump]
    type = FarmsMaterialRealAux3D
    material_property_name = 'local_dip_jump'
    variable = local_dip_jump
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetShearJumpRate]
    type = FarmsMaterialRealAux3D
    material_property_name = 'local_shear_jump_rate'
    variable = local_shear_jump_rate
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetNormalJumpRate]
    type = FarmsMaterialRealAux3D
    material_property_name = 'local_normal_jump_rate'
    variable = local_normal_jump_rate
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetDipJumpRate]
    type = FarmsMaterialRealAux3D
    material_property_name = 'local_dip_jump_rate'
    variable = local_dip_jump_rate
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetShearTraction]
    type = FarmsMaterialRealAux3D
    material_property_name = 'local_shear_traction'
    variable = local_shear_traction
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetNormalTraction]
    type = FarmsMaterialRealAux3D
    material_property_name = 'local_normal_traction'
    variable = local_normal_traction
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetDipTraction]
    type = FarmsMaterialRealAux3D
    material_property_name = 'local_dip_traction'
    variable = local_dip_traction
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetNormalX]
    type = FarmsMaterialRealAux3D
    material_property_name = 'normal_x'
    variable = normal_x
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetNormalY]
    type = FarmsMaterialRealAux3D
    material_property_name = 'normal_y'
    variable = normal_y
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
  []
  [GetNormalZ]
    type = FarmsMaterialRealAux3D
    material_property_name = 'normal_z'
    variable = normal_z
    boundary = 'Block100_Block200'
    ini_normal_sts = ini_normal_stress
    ini_shear_sts = ini_shear_stress
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
      type = SlipWeakeningFrictionczm3d
      disp_slipweakening_x     = disp_slipweakening_x
      disp_slipweakening_y     = disp_slipweakening_y
      disp_slipweakening_z     = disp_slipweakening_z
      vel_slipweakening_x      = vel_slipweakening_x
      vel_slipweakening_y      = vel_slipweakening_y
      vel_slipweakening_z      = vel_slipweakening_z
      reaction_slipweakening_x = resid_slipweakening_x
      reaction_slipweakening_y = resid_slipweakening_y
      reaction_slipweakening_z = resid_slipweakening_z
      mu_s = mu_s
      ini_shear_sts = ini_shear_stress
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
  dt = 0.005
  end_time = 3.0
  # num_steps = 10
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    execute_on = 'timestep_end'
    time_step_interval = 20
  []
[]

[BCs]
[]

[VectorPostprocessors]
  [main_fault]
    type = SideValueSampler
    variable = 'local_shear_traction local_shear_jump local_shear_jump_rate local_dip_traction local_dip_jump local_dip_jump_rate' 
    boundary = 'Block100_Block200'
    sort_by = x
  []
[]