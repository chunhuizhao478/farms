# DGRateStateFrictionMaterial Steady-State Test
# Test that the material correctly computes initial values
#
# This test verifies:
# 1. Initial state variable: θ = Dc/V_init = 0.004/1e-9 = 4e6 s
# 2. Initial friction coefficient from regularized law
# 3. Interface material property output via InterfaceMaterialRealAux

[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 4
    xmin = -100
    xmax = 0
    ymin = 0
    ymax = 200
    elem_type = QUAD4
    boundary_name_prefix = left
  []
  [left_block_id]
    type = SubdomainIDGenerator
    input = left_block
    subdomain_id = 1
  []
  [right_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 4
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 200
    elem_type = QUAD4
    boundary_name_prefix = right
  []
  [right_block_id]
    type = SubdomainIDGenerator
    input = right_block
    subdomain_id = 2
  []
  [stitch]
    type = StitchedMeshGenerator
    inputs = 'left_block_id right_block_id'
    stitch_boundaries_pairs = 'left_right right_left'
    clear_stitched_boundary_ids = false
  []
  [fault_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = stitch
    primary_block = 1
    paired_block = 2
    new_boundary = 'fault'
  []
[]

[Variables]
  [w]
    order = FIRST
    family = MONOMIAL
    block = '1 2'
  []
[]

[AuxVariables]
  [slip]
    order = CONSTANT
    family = MONOMIAL
  []
  [slip_rate]
    order = CONSTANT
    family = MONOMIAL
  []
  [state_variable]
    order = CONSTANT
    family = MONOMIAL
  []
  [fault_traction_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [friction_coeff]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [slip_aux]
    type = InterfaceMaterialRealAux
    variable = slip
    property = slip
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_aux]
    type = InterfaceMaterialRealAux
    variable = slip_rate
    property = slip_rate
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_aux]
    type = InterfaceMaterialRealAux
    variable = state_variable
    property = state_variable
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [traction_aux]
    type = InterfaceMaterialRealAux
    variable = fault_traction_aux
    property = fault_traction
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [friction_aux]
    type = InterfaceMaterialRealAux
    variable = friction_coeff
    property = friction_coefficient
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Kernels]
  [diffusion]
    type = Diffusion
    variable = w
    block = '1 2'
  []
[]

[DGKernels]
  [dg_elasticity]
    type = DGElasticityAntiplane
    variable = w
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
    block = '1 2'
  []
[]

[InterfaceKernels]
  [fault]
    type = DGAntiplaneFaultInterfaceKernel
    variable = w
    neighbor_var = w
    boundary = fault
    fault_traction = fault_traction
    dtraction_dslip = dtraction_dslip
    use_penalty = true
    penalty = 1e8
  []
[]

[Materials]
  [elastic]
    type = GenericConstantMaterial
    prop_names = 'shear_modulus'
    prop_values = '32.04e9'
    block = '1 2'
  []

  [rsf_material]
    type = DGRateStateFrictionMaterial
    displacement = w
    a = 0.015
    b = 0.020
    Dc = 0.004
    f0 = 0.6
    V0 = 1e-6
    sigma_n = 50e6
    shear_modulus = 32.04e9
    density = 2670
    initial_slip_rate = 1e-9
    boundary = fault
  []
[]

[BCs]
  [left_far]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_left
    function = '-1e-6'
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  [right_far]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_right
    function = '1e-6'
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  [top_left]
    type = DGElasticityNeumannBC
    variable = w
    boundary = left_top
    traction = 0
  []
  [top_right]
    type = DGElasticityNeumannBC
    variable = w
    boundary = right_top
    traction = 0
  []
  [bottom_left]
    type = DGElasticityNeumannBC
    variable = w
    boundary = left_bottom
    traction = 0
  []
  [bottom_right]
    type = DGElasticityNeumannBC
    variable = w
    boundary = right_bottom
    traction = 0
  []
[]

[Postprocessors]
  # Expected initial state: θ = Dc/V_init = 0.004/1e-9 = 4e6 s
  [state_avg]
    type = SideAverageValue
    variable = state_variable
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Expected initial slip rate: V_init = 1e-9 m/s
  [slip_rate_avg]
    type = SideAverageValue
    variable = slip_rate
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_avg]
    type = SideAverageValue
    variable = slip
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [traction_avg]
    type = SideAverageValue
    variable = fault_traction_aux
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [friction_avg]
    type = SideAverageValue
    variable = friction_coeff
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  dt = 1.0
  num_steps = 1
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
