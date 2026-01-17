# Unit test for ElkADMassConservationNewmark with damaged Biot properties
# Tests that the mass conservation kernel uses damage-dependent Biot coefficient and modulus
#
# Physics: alpha(d) * M(d) * div(v^s) * test
#
# This test verifies:
# 1. Kernel uses damaged Biot coefficient (alpha)
# 2. Kernel uses damaged Biot modulus (M)
# 3. Pressure response varies with damage level
#
# Test setup: 1D compression with spatially varying damage

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = 0
  xmax = 1
[]

[GlobalParams]
  displacements = 'disp_x'
[]

[Variables]
  [disp_x]
  []
  [pressure]
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
  [vel_x]
  []
  [accel_x]
  []
  [biot_coeff_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [biot_modulus_aux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  # Spatially varying damage: d = x (0 at left, 1 at right)
  [damage_ic]
    type = FunctionIC
    variable = d
    function = 'x'
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
  [get_biot_coeff]
    type = ADMaterialRealAux
    variable = biot_coeff_aux
    property = biot_coefficient
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_biot_modulus]
    type = ADMaterialRealAux
    variable = biot_modulus_aux
    property = biot_modulus
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Kernels]
  # Solid mechanics
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  # Poromechanics coupling
  [poro_x]
    type = ElkADPoroMechanicsCoupling
    variable = disp_x
    porepressure = pressure
    component = 0
    multiply_biot_coefficient = true
  []
  # Mass conservation with Newmark (uses damaged biot coefficient and modulus)
  [mass_conservation]
    type = ElkADMassConservationNewmark
    variable = pressure
    displacement_x = disp_x
    velocity_x = vel_x
    acceleration_x = accel_x
    beta = 0.25
    gamma = 0.5
    multiply_biot_coefficient = true
  []
  # Pressure diffusion (for numerical stabilization)
  [pressure_diff]
    type = ADDiffusion
    variable = pressure
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
    function = '-0.001*t'  # Compression
  []
  [left_p]
    type = DirichletBC
    variable = pressure
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
  # Permeability (required by poro-dynamic material)
  [effective_perm_provider]
    type = ADGenericConstantRankTwoTensor
    tensor_name = effective_perm
    tensor_values = '5e-19 0 0 0 5e-19 0 0 0 5e-19'
  []
  # Damaged hydraulic properties
  [damaged_biot_coefficient]
    type = ElkADPorousFlowDamagedBiotCoefficient
    phase_field = d
    solid_bulk_compliance = 2e-11  # 1/K0 = 1/50e9
    grain_bulk_modulus = 50e9      # Ks = K0 => alpha = 1 - g
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
  # Poro-dynamic assembly (provides biot_coefficient and biot_modulus from damaged)
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
  num_steps = 3
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
[]

[Postprocessors]
  # Biot coefficient at different damage levels
  [biot_coeff_left]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.1 0 0'
  []
  [biot_coeff_right]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.9 0 0'
  []
  # Biot modulus at different damage levels
  [biot_modulus_left]
    type = PointValue
    variable = biot_modulus_aux
    point = '0.1 0 0'
  []
  [biot_modulus_right]
    type = PointValue
    variable = biot_modulus_aux
    point = '0.9 0 0'
  []
  # Pressure at different locations (should vary with damage)
  [pressure_left]
    type = PointValue
    variable = pressure
    point = '0.1 0 0'
  []
  [pressure_mid]
    type = PointValue
    variable = pressure
    point = '0.5 0 0'
  []
  [pressure_right]
    type = PointValue
    variable = pressure
    point = '0.9 0 0'
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
