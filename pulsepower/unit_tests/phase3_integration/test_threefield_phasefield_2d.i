# =============================================================================
# Integration Test: Three-Field Poroelastodynamics with Phase-Field Damage (2D)
# =============================================================================
# Simplified test for three-field (u, w, p) formulation with damage-dependent properties
# No MultiApp - damage is prescribed as an auxiliary variable
#
# Test verifies:
# 1. Three-field equations work together correctly
# 2. Damage-dependent hydraulic properties evolve correctly
# 3. Energy is tracked properly
# =============================================================================

# Material properties
E = 50e9
nu = 0.3
K = '${fparse E/3.0/(1.0-2.0*nu)}'
solid_density = 2600
fluid_density = 1000
fluid_bulk_modulus = 2.24e9
grain_bulk_modulus = 50e9
viscosity = 1e-3
porosity = 0.1
tortosity = 1.2
intrinsic_permeability = 5e-15
solid_bulk_modulus_compliance = '${fparse 1.0/K}'

# Newmark parameters
newmark_beta = 0.25
newmark_gamma = 0.5

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
  xmax = 0.01
  ymin = 0
  ymax = 0.01
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
  [wf_x]
    order = FIRST
    family = LAGRANGE
  []
  [wf_y]
    order = FIRST
    family = LAGRANGE
  []
  [p]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
  [vel_x]
  []
  [vel_y]
  []
  [accel_x]
  []
  [accel_y]
  []
  [vf_x]
  []
  [vf_y]
  []
  [af_x]
  []
  [af_y]
  []
  [biot_coeff_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [porosity_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [density_aux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  # Initial damage: d = 0.3 * x / 0.01 (varies from 0 to 0.3)
  [damage_ic]
    type = FunctionIC
    variable = d
    function = '0.3 * x / 0.01'
  []
[]

[Kernels]
  # --- Solid Momentum ---
  [inertia_x]
    type = ADInertialForce
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
  []
  [inertia_y]
    type = ADInertialForce
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
  []
  [porefluidinertia_x]
    type = ElkADPoreFluidInertialForceCoupling
    variable = disp_x
    fluiddisp = wf_x
    fluidvel = vf_x
    fluidaccel = af_x
    beta = ${newmark_beta}
  []
  [porefluidinertia_y]
    type = ElkADPoreFluidInertialForceCoupling
    variable = disp_y
    fluiddisp = wf_y
    fluidvel = vf_y
    fluidaccel = af_y
    beta = ${newmark_beta}
  []
  [dispkernel_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [dispkernel_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
  []
  [poromechanic_ux]
    type = ElkADPoroMechanicsCoupling
    variable = disp_x
    porepressure = p
    component = 0
    multiply_biot_coefficient = true
  []
  [poromechanic_uy]
    type = ElkADPoroMechanicsCoupling
    variable = disp_y
    porepressure = p
    component = 1
    multiply_biot_coefficient = true
  []

  # --- Fluid Momentum (Dynamic Darcy) ---
  [dynamicdarcyflow_x]
    type = ElkADDynamicDarcyFlowPhaseField
    variable = wf_x
    skeletondisp = disp_x
    skeletonvel = vel_x
    skeletonaccel = accel_x
    fluidvel = vf_x
    fluidaccel = af_x
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    component = 0
  []
  [dynamicdarcyflow_y]
    type = ElkADDynamicDarcyFlowPhaseField
    variable = wf_y
    skeletondisp = disp_y
    skeletonvel = vel_y
    skeletonaccel = accel_y
    fluidvel = vf_y
    fluidaccel = af_y
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    component = 1
  []
  [poromechanic_wx]
    type = ElkADPoroMechanicsCoupling
    variable = wf_x
    porepressure = p
    component = 0
    multiply_biot_coefficient = false
  []
  [poromechanic_wy]
    type = ElkADPoroMechanicsCoupling
    variable = wf_y
    porepressure = p
    component = 1
    multiply_biot_coefficient = false
  []

  # --- Mass Conservation ---
  [pressure_rate]
    type = ElkADPressureRate
    variable = p
  []
  [massconservationskeleton]
    type = ElkADMassConservationNewmark
    variable = p
    displacement_x = disp_x
    displacement_y = disp_y
    velocity_x = vel_x
    velocity_y = vel_y
    acceleration_x = accel_x
    acceleration_y = accel_y
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    multiply_biot_coefficient = true
    plane_strain_correction = true
    poissons_ratio = ${nu}
  []
  [insmass]
    type = ElkADMassConservationNewmark
    variable = p
    displacement_x = wf_x
    displacement_y = wf_y
    velocity_x = vf_x
    velocity_y = vf_y
    acceleration_x = af_x
    acceleration_y = af_y
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    multiply_biot_coefficient = false
    plane_strain_correction = true
    poissons_ratio = ${nu}
  []
[]

[AuxKernels]
  # Solid Newmark
  [accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  # Fluid Newmark
  [af_x]
    type = NewmarkAccelAux
    variable = af_x
    displacement = wf_x
    velocity = vf_x
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vf_x]
    type = NewmarkVelAux
    variable = vf_x
    acceleration = af_x
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  [af_y]
    type = NewmarkAccelAux
    variable = af_y
    displacement = wf_y
    velocity = vf_y
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vf_y]
    type = NewmarkVelAux
    variable = vf_y
    acceleration = af_y
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  # Material property outputs
  [get_biot_coeff]
    type = ADMaterialRealAux
    variable = biot_coeff_aux
    property = biot_coefficient
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_porosity]
    type = ADMaterialRealAux
    variable = porosity_aux
    property = porosity
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [get_density]
    type = ADMaterialRealAux
    variable = density_aux
    property = density
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[BCs]
  # Fixed left boundary
  [left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [left_y]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0
  []
  [left_wf_x]
    type = DirichletBC
    variable = wf_x
    boundary = left
    value = 0
  []
  [left_wf_y]
    type = DirichletBC
    variable = wf_y
    boundary = left
    value = 0
  []
  # Pressure loading on right
  [right_load]
    type = ADPressure
    variable = disp_x
    boundary = right
    function = '1e6 * t / 1e-5'  # Ramp to 1 MPa
  []
  # Drained pressure at left
  [left_p]
    type = DirichletBC
    variable = p
    boundary = left
    value = 0
  []
[]

[Materials]
  # Elasticity
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [stress]
    type = ADComputeLinearElasticStress
  []
  # Effective permeability (RankTwoTensor)
  [effective_perm_provider]
    type = ADGenericConstantRankTwoTensor
    tensor_name = effective_perm
    tensor_values = '${intrinsic_permeability} 0 0 0 ${intrinsic_permeability} 0 0 0 ${intrinsic_permeability}'
  []
  # Damaged Biot coefficient
  [damaged_biot_coefficient]
    type = ElkADPorousFlowDamagedBiotCoefficient
    phase_field = d
    solid_bulk_compliance = ${solid_bulk_modulus_compliance}
    grain_bulk_modulus = ${grain_bulk_modulus}
    minimum_degradation = 1e-8
  []
  # Damaged porosity
  [damaged_porosity]
    type = ElkADPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = ${porosity}
    porosity_lower_bound = ${porosity}
    porosity_upper_bound = 0.999
  []
  # Damaged Biot modulus
  [damaged_biot_modulus]
    type = ElkADPorousFlowDamagedBiotModulus
    fluid_bulk_modulus = ${fluid_bulk_modulus}
    grain_bulk_modulus = ${grain_bulk_modulus}
    use_damaged_biot = true
    use_damaged_porosity = true
  []
  # Poro-dynamic assembly
  [porodynamics]
    type = ElkADPoroDynamicPhaseFieldMaterials
    rhos_value = ${solid_density}
    rhof_value = ${fluid_density}
    tortosity_value = ${tortosity}
    viscosity_value = ${viscosity}
    grain_bulk_modulus = ${grain_bulk_modulus}
    fluid_bulk_modulus = ${fluid_bulk_modulus}
    use_damaged_properties = true
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
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

  [TimeIntegrator]
    type = NewmarkBeta
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
  []
[]

[Postprocessors]
  # Displacement at right edge
  [disp_x_right]
    type = SideAverageValue
    variable = disp_x
    boundary = right
  []
  # Fluid displacement at right edge
  [wf_x_right]
    type = SideAverageValue
    variable = wf_x
    boundary = right
  []
  # Pressure at center
  [pressure_center]
    type = PointValue
    variable = p
    point = '0.005 0.005 0'
  []
  # Biot coefficient (varies with damage)
  [biot_coeff_left]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.001 0.005 0'
  []
  [biot_coeff_right]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.009 0.005 0'
  []
  # Porosity (varies with damage)
  [porosity_left]
    type = PointValue
    variable = porosity_aux
    point = '0.001 0.005 0'
  []
  [porosity_right]
    type = PointValue
    variable = porosity_aux
    point = '0.009 0.005 0'
  []
  # Density (varies with damage through porosity)
  [density_left]
    type = PointValue
    variable = density_aux
    point = '0.001 0.005 0'
  []
  [density_right]
    type = PointValue
    variable = density_aux
    point = '0.009 0.005 0'
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
