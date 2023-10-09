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
    type = FileMeshGenerator
    file =  '../../meshgenerator/cdbm/contmf/tria/contmfsmall.msh'
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
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y'
  
  #damping ratio
  q = 0.2
  
  #characteristic length (m)
  Dc = 0.4
  
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  []
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
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm_ik]
    boundary = 'Block0_Block1'
    strain = SMALL
    generate_output='tangent_jump'
  [../]
[]

[Problem]
  extra_tag_vectors = 'restore_tag'
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
  [StaticFricCoeff]
    type = FunctionAux
    variable = mu_s
    function = func_static_friction_coeff_mus
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [DynamicFricCoeff]
      type = FunctionAux
      variable = mu_d
      function = func_dynamic_friction_coeff_mud
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
  [StrikeShearStress]
    type = FunctionAux
    variable = ini_shear_stress
    function = func_initial_strike_shear_stress
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [TJump_rate]
    type = FDCompVarRate
    variable = tangent_jump_rate
    coupled = tangent_jump
    execute_on = 'TIMESTEP_BEGIN'
  []
  #fault length
  [fault_len]
      type = ConstantAux
      variable = nodal_area
      value = 25
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Kernels]
  [./stsdev_x]
    type = StressDivergenceTensors
    variable = disp_x
    displacements = 'disp_x disp_y'
    component = 0
    use_displaced_mesh = false
  []
  [./stsdev_y]
    type = StressDivergenceTensors
    variable = disp_y
    displacements = 'disp_x disp_y'
    component = 1
    use_displaced_mesh = false
  []
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
  [strain]
    type = ComputeIncrementalSmallStrain
    eigenstrain_names = eigenstrain_all
    displacements = 'disp_x disp_y'
  []
  [eigenstrain_all]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'func_initial_stress_xx           func_initial_strike_shear_stress 0  
                      func_initial_strike_shear_stress func_initial_stress_yy           0      
                      0 0 func_initial_stress_zz'
    eigenstrain_name = eigenstrain_all
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    lambda = 32.04e9
    shear_modulus = 32.04e9
    use_displaced_mesh = false
  []
  [./admissible]
    type = ComputeMultipleInelasticStress
    inelastic_models = cdp
    perform_finite_strain_rotations = false
  [../]
  [./cdp]
    type = CappedDruckerPragerStressUpdate
    DP_model = dp
    tensile_strength = ts
    compressive_strength = cs
    yield_function_tol = 1E-8
    tip_smoother = 4
    smoothing_tol = 1E-5
  [../]
  [density]
      type = GenericConstantMaterial
      prop_names = density
      prop_values = 2670
  []
  #SlipWeakeningMultifaults ONLY supports TRIA currently!
  [./czm_mat]
      type = SlipWeakeningMultifaults
      disp_slipweakening_x     = disp_slipweakening_x
      disp_slipweakening_y     = disp_slipweakening_y
      reaction_slipweakening_x = resid_slipweakening_x
      reaction_slipweakening_y = resid_slipweakening_y
      nodal_area = nodal_area
      mu_d = mu_d
      mu_s = mu_s
      tria_area = tria_area_aux
      boundary = 'Block0_Block1'
  [../]
  [./static_initial_stress_tensor_slipweakening]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor_slipweakening
      tensor_functions = 'func_initial_stress_xx         func_initial_strike_shear_stress           func_initial_stress_00 
      func_initial_strike_shear_stress         func_initial_stress_yy           func_initial_stress_00
                          func_initial_stress_00         func_initial_stress_00           func_initial_stress_00'
  [../]
  [./static_initial_stress_tensor]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor
      tensor_functions = 'func_initial_stress_xx             func_initial_stress_xy_const        func_initial_stress_00 
                          func_initial_stress_xy_const       func_initial_stress_yy              func_initial_stress_00
                          func_initial_stress_00             func_initial_stress_00              func_initial_stress_00'
  [../]
[]

[Functions]
  #mus constant value: 0.7
  [func_static_friction_coeff_mus]
      type = ConstantFunction
      value = 0.677
  []
  #mud constant value: 0.4
  [func_dynamic_friction_coeff_mud]
      type = ConstantFunction
      value = 0.1
  []
  #Note:restrict stress variation along the fault only
  #this function is used in czm only
  [func_initial_strike_shear_stress]
    type = InitialStressXYcontmfbfs
    # type = ConstantFunction
    # value = 70e6
  []
  #this function is used in medimum
  [func_initial_stress_xy_const]
    type = ConstantFunction
    value = 70e6
  []
  [./func_initial_stress_00]
    type = ConstantFunction
    value = 0.0
  []
  [./func_initial_stress_yy]
    type = ConstantFunction
    value = 120e6
  []
  #In problems with inelasticity, the sigma11 is important
  #This is different from pure tpv205 
  [./func_initial_stress_xx]
    type = ConstantFunction
    value = 135e6
  []
  [./func_initial_stress_zz]
    type = ConstantFunction
    value = 63.75e6
  []
[]

[UserObjects]
  [recompute_residual_tag]
      type = ResidualEvaluationUserObject
      vector_tag = 'restore_tag'
      force_preaux = true
      execute_on = 'TIMESTEP_END'
  []
  [./ts]
    type = TensorMechanicsHardeningConstant
    value = 1000
  [../]
  [./cs]
    type = TensorMechanicsHardeningConstant
    value = 1000
  [../]
  [./mc_coh]
    type = TensorMechanicsHardeningConstant
    value = 4.3e6
  [../]
  [./mc_phi]
    type = TensorMechanicsHardeningConstant
    value = 30.5
    convert_to_radians = true
  [../]
  [./mc_psi]
    type = TensorMechanicsHardeningConstant
    value = 0
  [../]
  [./dp]
    type = TensorMechanicsPlasticDruckerPrager
    mc_cohesion = mc_coh
    mc_friction_angle = mc_phi
    mc_dilation_angle = mc_psi
    mc_interpolation_scheme = outer_tip
    internal_constraint_tolerance = 1 # irrelevant here
    yield_function_tolerance = 1      # irrelevant here
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.0025
  end_time = 4.0
  num_steps = 10
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

[Outputs]
  exodus = true
  interval = 20
[]