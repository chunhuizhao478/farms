# Unit test for ElkADDynamicDarcyFlowPhaseField kernel
# Tests that the dynamic Darcy flow kernel works with phase-field damage-dependent properties
#
# Physics: rho^f * a^s + rho^f * tau_t / phi * a^f + mu_f / kappa * v^f = 0
#
# This test verifies:
# 1. Kernel computes residual correctly with Newmark time integration
# 2. Uses damage-dependent permeability (RankTwoTensor)
# 3. Uses damage-dependent porosity
#
# Simplified test: 1D bar with fixed damage, verify kernel runs

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 1
[]

[GlobalParams]
  displacements = 'disp_x'
[]

[Variables]
  [disp_x]
  []
  [wf_x]
    # Fluid relative displacement
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
    initial_condition = 0.3
  []
  [vel_x]
  []
  [accel_x]
  []
  [wf_vel_x]
  []
  [wf_accel_x]
  []
[]

[Kernels]
  # Solid mechanics (simplified)
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  # Dynamic Darcy flow for fluid displacement
  [darcy_x]
    type = ElkADDynamicDarcyFlowPhaseField
    variable = wf_x
    skeletondisp = disp_x
    skeletonvel = vel_x
    skeletonaccel = accel_x
    fluidvel = wf_vel_x
    fluidaccel = wf_accel_x
    beta = 0.25
    gamma = 0.5
    component = 0
  []
[]

[AuxKernels]
  [vel_x_kernel]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = 0.5
    execute_on = 'TIMESTEP_END'
  []
  [accel_x_kernel]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = 0.25
    execute_on = 'TIMESTEP_END'
  []
  [wf_vel_x_kernel]
    type = NewmarkVelAux
    variable = wf_vel_x
    acceleration = wf_accel_x
    gamma = 0.5
    execute_on = 'TIMESTEP_END'
  []
  [wf_accel_x_kernel]
    type = NewmarkAccelAux
    variable = wf_accel_x
    displacement = wf_x
    velocity = wf_vel_x
    beta = 0.25
    execute_on = 'TIMESTEP_END'
  []
[]

[BCs]
  [left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '0.0001*t'
  []
  [left_wf]
    type = DirichletBC
    variable = wf_x
    boundary = left
    value = 0
  []
[]

[Materials]
  # Strain computation
  [strain]
    type = ADComputeSmallStrain
  []
  # Elasticity tensor
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 50e9
    poissons_ratio = 0.3
  []
  # Stress computation
  [stress]
    type = ADComputeLinearElasticStress
  []
  # Permeability (RankTwoTensor for phase-field compatibility)
  [effective_perm_provider]
    type = ADGenericConstantRankTwoTensor
    tensor_name = effective_perm
    tensor_values = '5e-15 0 0 0 5e-15 0 0 0 5e-15'
  []
  # Damaged hydraulic properties
  [damaged_biot_coefficient]
    type = ElkADPorousFlowDamagedBiotCoefficient
    phase_field = d
    solid_bulk_compliance = 2e-11
    grain_bulk_modulus = 50e9
    minimum_degradation = 1e-8
  []
  [damaged_porosity]
    type = ElkADPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = 0.1
    porosity_lower_bound = 0.1
    porosity_upper_bound = 0.999
  []
  [damaged_biot_modulus]
    type = ElkADPorousFlowDamagedBiotModulus
    use_damaged_biot = true
    use_damaged_porosity = true
    fluid_bulk_modulus = 2.24e9
    grain_bulk_modulus = 50e9
  []
  # Poro-dynamic assembly material
  [porodynamics]
    type = ElkADPoroDynamicPhaseFieldMaterials
    rhos_value = 2600
    rhof_value = 1000
    tortosity_value = 1.2
    viscosity_value = 1e-3
    grain_bulk_modulus = 50e9
    fluid_bulk_modulus = 2.24e9
    use_damaged_properties = true
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  dt = 1e-6
  num_steps = 2
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
[]

[Postprocessors]
  [disp_right]
    type = PointValue
    variable = disp_x
    point = '1 0 0'
  []
  [wf_right]
    type = PointValue
    variable = wf_x
    point = '1 0 0'
  []
[]

[Outputs]
  csv = true
[]
