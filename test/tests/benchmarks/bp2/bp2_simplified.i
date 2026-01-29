# BP2 Benchmark - Simplified Version
# 2D Antiplane Shear SEAS Problem with DG Formulation
#
# This is a simplified version of the SCEC BP2 benchmark for
# testing the DG SEAS framework. Full benchmark requires:
# - Larger domain (40 km fault depth)
# - Finer mesh (25-800 m resolution)
# - 1200 years simulation time
#
# This test uses:
# - Smaller domain (10 km x 10 km)
# - Coarse mesh for testing
# - Short simulation time

# Domain: [-5km, 5km] x [0, 10km]
# Fault at x = 0

[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 20
    xmin = -5000  # -5 km
    xmax = 0
    ymin = 0      # Free surface
    ymax = 10000  # 10 km depth
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
    nx = 10
    ny = 20
    xmin = 0
    xmax = 5000   # +5 km
    ymin = 0
    ymax = 10000
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
    epsilon = 1.0  # SIPG
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
    penalty = 1e10  # For stability
  []
[]

[Materials]
  # Bulk elastic properties (BP2 values)
  [elastic]
    type = GenericConstantMaterial
    prop_names = 'shear_modulus'
    prop_values = '32.04e9'  # Pa
    block = '1 2'
  []

  # Radiation damping
  [rad_damp]
    type = RadiationDampingMaterial
    shear_modulus = 32.04e9
    density = 2670
    block = '1 2'
  []

  # Fault friction - simplified with constant traction for now
  # Full rate-state would use DGRateStateFrictionMaterial
  [fault_mat]
    type = DGFaultTractionMaterial
    displacement = w
    tau0 = 25e6  # 25 MPa background shear stress
    dtau_ds = 0.0  # No slip-traction coupling for simplified test
    boundary = fault
  []
[]

[BCs]
  # Far-field: prescribed plate velocity loading
  # Left boundary: moving at -Vp/2
  [left_far_field]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_left
    function = '-0.5e-9*t'  # -Vp/2 * t, Vp = 1e-9 m/s
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  # Right boundary: moving at +Vp/2
  [right_far_field]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_right
    function = '0.5e-9*t'  # +Vp/2 * t
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  # Free surface at top (z = 0) - traction free
  [free_surface_left]
    type = DGElasticityNeumannBC
    variable = w
    boundary = left_top
    traction = 0
  []
  [free_surface_right]
    type = DGElasticityNeumannBC
    variable = w
    boundary = right_top
    traction = 0
  []
  # Bottom boundary - absorbing or fixed
  [bottom_left]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_bottom
    function = '-0.5e-9*t'
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
  [bottom_right]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_bottom
    function = '0.5e-9*t'
    epsilon = 1.0
    sigma = 6.0
    shear_modulus = shear_modulus
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  # Time stepping
  dt = 1e6  # 1e6 seconds ~ 11.6 days
  num_steps = 10  # Just 10 steps for testing
[]

[Postprocessors]
  # Monitor slip at mid-fault
  [w_left_fault]
    type = SideAverageValue
    variable = w
    boundary = fault
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Displacement at surface
  [w_surface_left]
    type = PointValue
    variable = w
    point = '-2500 0 0'
  []
  [w_surface_right]
    type = PointValue
    variable = w
    point = '2500 0 0'
  []
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
  [console]
    type = Console
    output_linear = false
  []
[]
