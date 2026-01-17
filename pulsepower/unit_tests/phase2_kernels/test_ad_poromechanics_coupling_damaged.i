# Unit test for ElkADPoroMechanicsCoupling with damaged Biot coefficient
# Tests that the poromechanics coupling kernel uses the damage-dependent Biot coefficient
#
# Physics: -alpha(d) * p * grad(test)
#
# This test verifies:
# 1. Kernel correctly multiplies by damaged Biot coefficient
# 2. The coupling force is proportional to damage
#
# Test setup: 2D square with uniform pressure, measure displacement response
# At d=0: alpha=0, no coupling
# At d=0.5: alpha=0.75, partial coupling
# At d=1: alpha≈1, full coupling

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  ny = 3
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
  [pressure]
    family = LAGRANGE
    order = FIRST
    initial_condition = 1e6
  []
  [biot_coeff_aux]
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
  [get_biot_coeff]
    type = ADMaterialRealAux
    variable = biot_coeff_aux
    property = biot_coefficient
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
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
  []
  # Poromechanics coupling with damaged Biot coefficient
  [poro_x]
    type = ElkADPoroMechanicsCoupling
    variable = disp_x
    porepressure = pressure
    component = 0
    multiply_biot_coefficient = true
  []
  [poro_y]
    type = ElkADPoroMechanicsCoupling
    variable = disp_y
    porepressure = pressure
    component = 1
    multiply_biot_coefficient = true
  []
[]

[BCs]
  [left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
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
  # Poro-dynamic assembly (provides biot_coefficient from damaged)
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
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Postprocessors]
  # Biot coefficient at different damage levels
  [biot_coeff_left]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.1667 0.5 0'
  []
  [biot_coeff_mid]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.5 0.5 0'
  []
  [biot_coeff_right]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.8333 0.5 0'
  []
  # Displacement at right edge (should vary with damage)
  [disp_x_right]
    type = SideAverageValue
    variable = disp_x
    boundary = right
  []
  [disp_y_top]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
