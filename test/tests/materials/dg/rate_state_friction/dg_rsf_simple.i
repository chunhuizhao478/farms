# DGRateStateFrictionMaterial Simple Test
# Use a smaller domain and check if the solve works

[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
    xmin = -10
    xmax = 0
    ymin = 0
    ymax = 20
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
    ny = 2
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 20
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
    penalty = 1e12
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
    function = '0'
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  [right_far]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_right
    function = '0'
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
  [state_avg]
    type = SideAverageValue
    variable = state_variable
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
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
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  l_tol = 1e-4

  dt = 1e6
  num_steps = 3
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
