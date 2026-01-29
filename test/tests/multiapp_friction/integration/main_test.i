# MultiApp Integration Test - Main Elasticity App
#
# Simplified test to verify MultiApp coupling between elasticity and friction.
# Uses small mesh and few time steps for quick verification.

mu = 32.04e9
Vp = 1e-9

[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 5
    xmin = -1000
    xmax = 0
    ymin = 0
    ymax = 1000
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
    nx = 5
    ny = 5
    xmin = 0
    xmax = 1000
    ymin = 0
    ymax = 1000
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
  [slip_from_subapp]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  []
  [traction_to_subapp]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
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
    type = DGFaultSlipInterfaceKernel
    variable = w
    neighbor_var = w
    boundary = fault
    slip_prescribed = slip_from_subapp
    shear_modulus = shear_modulus
    penalty = 1e12
    epsilon = 1.0
    sigma = 6.0
  []
[]

[AuxKernels]
  [compute_traction]
    type = ElasticTractionAux
    variable = traction_to_subapp
    displacement = w
    shear_modulus = ${mu}
    normal_component = 0
    execute_on = 'TIMESTEP_END'
    boundary = fault
  []
[]

[Materials]
  [elastic]
    type = GenericConstantMaterial
    prop_names = 'shear_modulus'
    prop_values = '${mu}'
    block = '1 2'
  []
[]

[BCs]
  [left_far_field]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_left
    function = '-0.5 * ${Vp} * t'
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  [right_far_field]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_right
    function = '0.5 * ${Vp} * t'
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  [free_surface_left]
    type = DGElasticityNeumannBC
    variable = w
    boundary = left_bottom
    traction = 0
  []
  [free_surface_right]
    type = DGElasticityNeumannBC
    variable = w
    boundary = right_bottom
    traction = 0
  []
  [bottom_left]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_top
    function = '-0.5 * ${Vp} * t'
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  [bottom_right]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_top
    function = '0.5 * ${Vp} * t'
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
[]

[MultiApps]
  [friction_app]
    type = TransientMultiApp
    input_files = 'sub_test.i'
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[Transfers]
  [send_traction]
    type = MultiAppGeneralFieldNearestLocationTransfer
    to_multi_app = friction_app
    source_variable = traction_to_subapp
    variable = traction_received
    from_boundaries = 'fault'
  []
  [receive_slip]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = friction_app
    source_variable = slip
    variable = slip_from_subapp
    to_boundaries = 'fault'
  []
[]

[Postprocessors]
  [avg_slip]
    type = SideAverageValue
    variable = slip_from_subapp
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [avg_traction]
    type = SideAverageValue
    variable = traction_to_subapp
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
  dt = 1e6
  num_steps = 3
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
