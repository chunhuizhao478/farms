# DG Fault Interface Baseline Test
# Test with simplified DGFaultTractionMaterial to verify framework

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

  [fault_mat]
    type = DGFaultTractionMaterial
    displacement = w
    tau0 = 25e6
    dtau_ds = 0.0
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
  [w_left_avg]
    type = ElementAverageValue
    variable = w
    block = 1
  []
  [w_right_avg]
    type = ElementAverageValue
    variable = w
    block = 2
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
