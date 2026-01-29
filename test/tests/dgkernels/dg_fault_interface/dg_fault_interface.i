# DG Fault Interface Test
# Tests the fault interface kernel with prescribed traction
# Setup: Two blocks with a vertical fault at x=0.5
# The fault has a prescribed constant traction

[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 10
    xmin = 0
    xmax = 0.5
    ymin = 0
    ymax = 1
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
    nx = 5
    ny = 10
    xmin = 0.5
    xmax = 1.0
    ymin = 0
    ymax = 1
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
  []
[]

[Materials]
  [elastic]
    type = GenericConstantMaterial
    prop_names = 'shear_modulus'
    prop_values = '1.0'
    block = '1 2'
  []
  [fault_mat]
    type = DGFaultTractionMaterial
    displacement = w
    tau0 = 0.5  # Prescribed constant traction
    dtau_ds = 0.0  # No slip-traction coupling for this test
    boundary = fault
  []
[]

[BCs]
  [left]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_left
    value = 0
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  [right]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_right
    value = 1
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Postprocessors]
  # Check solution on left side of fault
  [w_left_of_fault]
    type = PointValue
    variable = w
    point = '0.45 0.5 0'
  []
  # Check solution on right side of fault
  [w_right_of_fault]
    type = PointValue
    variable = w
    point = '0.55 0.5 0'
  []
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
