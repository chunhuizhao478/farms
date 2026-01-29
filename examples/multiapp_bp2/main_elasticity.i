# SCEC SEAS Benchmark Problem BP2-QD
# MultiApp Architecture - Main Application: Elasticity Solver
#
# This is the main app that solves the 2D antiplane elasticity problem.
# It uses a MultiApp to couple with the Friction SubApp.
#
# Mesh size: 800 m
#
# Reference: SEAS Benchmark Problem BP2-QD specification
# https://strike.scec.org/cvws/seas/

# Physical constants
mu = 32.04e9        # Shear modulus (Pa)

# Loading
Vp = 1e-9           # Plate rate (m/s)
Vinit = 1e-9        # Initial slip rate (m/s)

# Geometry (in meters)
L_x = 100000        # Half-width of domain (100 km)
L_z = 100000        # Total domain depth (100 km)
Wf = 40000          # Seismogenic depth (40 km) - fault ends here

# Pre-stress (for output only - SubApp handles friction balance)
tau0 = 26546122.0   # Initial shear stress (Pa)

# Mesh size
dz = 800            # Element size (m)
nz = ${fparse L_z / dz}  # Number of elements along depth = 125
nx = ${fparse L_x / dz}  # Number of elements each side = 125

# Following Tandem's mesh structure:
# - Fault only exists from y=0 to y=-Wf (seismogenic zone)
# - Below y=-Wf, both blocks meet but with NO slip allowed (continuous material)
# Strategy: Create two blocks spanning full depth, stitch at x=0, then create fault
# boundary only in the seismogenic zone (y > -Wf)
[Mesh]
  # Left block: x ∈ [-L_x, 0], y ∈ [-L_z, 0]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${nz}
    xmin = -${L_x}
    xmax = 0
    ymin = -${L_z}          # Deep boundary
    ymax = 0                # Free surface
    elem_type = QUAD4
    boundary_name_prefix = left
  []
  [left_block_id]
    type = SubdomainIDGenerator
    input = left_block
    subdomain_id = 1
  []

  # Right block: x ∈ [0, L_x], y ∈ [-L_z, 0]
  [right_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${nz}
    xmin = 0
    xmax = ${L_x}
    ymin = -${L_z}          # Deep boundary
    ymax = 0                # Free surface
    elem_type = QUAD4
    boundary_name_prefix = right
  []
  [right_block_id]
    type = SubdomainIDGenerator
    input = right_block
    subdomain_id = 2
  []

  # Stitch left and right blocks at x=0
  [stitch]
    type = StitchedMeshGenerator
    inputs = 'left_block_id right_block_id'
    stitch_boundaries_pairs = 'left_right right_left'
    clear_stitched_boundary_ids = false
  []

  # Create fault interface boundary between blocks 1 and 2 (full depth)
  [fault_full]
    type = SideSetsBetweenSubdomainsGenerator
    input = stitch
    primary_block = 1
    paired_block = 2
    new_boundary = 'fault_full'
  []

  # Restrict fault boundary to seismogenic zone (y > -Wf) using ParsedGenerateSideset
  # The 'fault' boundary will only include faces where y > -Wf
  [fault_sideset]
    type = ParsedGenerateSideset
    input = fault_full
    combinatorial_geometry = 'abs(x) < 0.01 & y > -${Wf}'
    new_sideset_name = 'fault'
    include_only_external_sides = false
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
  # Slip received from Friction SubApp
  [slip_from_subapp]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  []
  # Traction to send to SubApp (total traction = elastic + tau_pre)
  [traction_to_subapp]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = ${tau0}
  []
  # Slip rate received from SubApp (for output)
  [slip_rate_from_subapp]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = ${Vinit}
  []
  # State variable received from SubApp (for output)
  [state_from_subapp]
    order = CONSTANT
    family = MONOMIAL
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
    sigma = 200.0  # Increased for stability
    shear_modulus = shear_modulus
    block = '1 2'
  []
[]

[InterfaceKernels]
  # Enforce slip = prescribed_slip at fault interface
  [fault]
    type = DGFaultSlipInterfaceKernel
    variable = w
    neighbor_var = w
    boundary = fault
    slip_prescribed = slip_from_subapp
    shear_modulus = shear_modulus
    penalty = 1e14   # Increased for stability
    epsilon = 1.0
    sigma = 200.0    # Increased for stability
  []
[]

[AuxKernels]
  # Read traction from interface material (same approach as single app)
  # This ensures identical traction computation as single app
  [compute_traction]
    type = SEASTractionAux
    variable = traction_to_subapp
    traction_property = elastic_traction
    tau_pre = ${tau0}
    execute_on = 'TIMESTEP_END'
    boundary = fault
  []
[]

[Materials]
  # Bulk elastic properties
  [elastic]
    type = GenericConstantMaterial
    prop_names = 'shear_modulus'
    prop_values = '${mu}'
    block = '1 2'
  []
  # Elastic traction at fault interface
  # Use stiffness mode for correct quasi-static SEAS stress evolution
  [elastic_traction]
    type = DGElasticTractionMaterial
    displacement = w
    slip_prescribed = slip_from_subapp
    shear_modulus = ${mu}
    sigma = 200.0            # Must match DGFaultSlipInterfaceKernel
    penalty = 1e14           # Must match DGFaultSlipInterfaceKernel
    traction_scale = 1.0
    traction_name = elastic_traction
    boundary = fault
    # Use stiffness mode for quasi-static SEAS problems
    # This correctly captures stress evolution: τ = K * (Vp*t - slip)
    traction_mode = stiffness
    plate_rate = ${Vp}           # 1e-9 m/s
    fault_depth = ${Wf}          # 40 km
    domain_half_width = ${L_x}   # 100 km
  []
[]

[BCs]
  # Far-field Dirichlet BCs - plate motion drives loading
  [left_far_field]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_left
    function = far_field_left
    epsilon = 1.0
    sigma = 200.0
    shear_modulus = shear_modulus
  []
  [right_far_field]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_right
    function = far_field_right
    epsilon = 1.0
    sigma = 200.0
    shear_modulus = shear_modulus
  []

  # Free surface at y = 0 (top) - traction free
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

  # Deep boundary (y = -L_z) - far-field plate motion
  [deep_left]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_bottom
    function = far_field_left
    epsilon = 1.0
    sigma = 200.0
    shear_modulus = shear_modulus
  []
  [deep_right]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_bottom
    function = far_field_right
    epsilon = 1.0
    sigma = 200.0
    shear_modulus = shear_modulus
  []
[]

[Functions]
  # Far-field boundary conditions: ±Vp*t/2
  # Following Tandem's approach for antiplane shear
  [far_field_left]
    type = ParsedFunction
    expression = '-0.5 * ${Vp} * t'
  []
  [far_field_right]
    type = ParsedFunction
    expression = '0.5 * ${Vp} * t'
  []
[]

[MultiApps]
  [friction_app]
    type = TransientMultiApp
    input_files = 'sub_friction.i'
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

[Transfers]
  # Send traction to Friction SubApp (before SubApp solve)
  # Note: We transfer elastic traction only - SubApp adds tau_pre internally
  [send_traction]
    type = MultiAppGeneralFieldNearestLocationTransfer
    to_multi_app = friction_app
    source_variable = traction_to_subapp
    variable = traction_received
    from_boundaries = 'fault'
  []

  # Receive slip from Friction SubApp (after SubApp solve)
  [receive_slip]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = friction_app
    source_variable = slip
    variable = slip_from_subapp
    to_boundaries = 'fault'
  []

  # Receive slip rate for output and adaptive time stepping
  [receive_slip_rate]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = friction_app
    source_variable = slip_rate
    variable = slip_rate_from_subapp
    to_boundaries = 'fault'
  []

  # Receive state variable for output
  [receive_state]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = friction_app
    source_variable = state_variable
    variable = state_from_subapp
    to_boundaries = 'fault'
  []

  # Receive max slip rate for adaptive time stepping
  [receive_max_slip_rate]
    type = MultiAppPostprocessorTransfer
    from_multi_app = friction_app
    from_postprocessor = max_slip_rate
    to_postprocessor = max_slip_rate
    reduction_type = maximum
  []
[]

[Postprocessors]
  # Output at z = 0 km (free surface, y = 0) - VW region
  [slip_z0]
    type = PointValue
    variable = slip_from_subapp
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z0]
    type = PointValue
    variable = slip_rate_from_subapp
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z0]
    type = PointValue
    variable = traction_to_subapp
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z0]
    type = PointValue
    variable = state_from_subapp
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Output at z = 4.8 km (y = -4800 m) - VW region
  [slip_z4_8]
    type = PointValue
    variable = slip_from_subapp
    point = '0 -4800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z4_8]
    type = PointValue
    variable = slip_rate_from_subapp
    point = '0 -4800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z4_8]
    type = PointValue
    variable = traction_to_subapp
    point = '0 -4800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z4_8]
    type = PointValue
    variable = state_from_subapp
    point = '0 -4800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Output at z = 12 km (y = -12000 m) - VW region
  [slip_z12]
    type = PointValue
    variable = slip_from_subapp
    point = '0 -12000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z12]
    type = PointValue
    variable = slip_rate_from_subapp
    point = '0 -12000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z12]
    type = PointValue
    variable = traction_to_subapp
    point = '0 -12000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z12]
    type = PointValue
    variable = state_from_subapp
    point = '0 -12000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Output at z = 16.8 km (y = -16800 m) - transition zone
  [slip_z16_8]
    type = PointValue
    variable = slip_from_subapp
    point = '0 -16800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z16_8]
    type = PointValue
    variable = slip_rate_from_subapp
    point = '0 -16800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z16_8]
    type = PointValue
    variable = traction_to_subapp
    point = '0 -16800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z16_8]
    type = PointValue
    variable = state_from_subapp
    point = '0 -16800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Output at z = 24 km (y = -24000 m) - VS region
  [slip_z24]
    type = PointValue
    variable = slip_from_subapp
    point = '0 -24000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z24]
    type = PointValue
    variable = slip_rate_from_subapp
    point = '0 -24000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z24]
    type = PointValue
    variable = traction_to_subapp
    point = '0 -24000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z24]
    type = PointValue
    variable = state_from_subapp
    point = '0 -24000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Maximum slip rate for adaptive time stepping
  [max_slip_rate]
    type = Receiver
    default = ${Vinit}
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
  []

  # Time in years for convenience
  [time_years]
    type = FunctionValuePostprocessor
    function = 't / 31557600'  # seconds per year
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  # Use LU direct solver for DG - more robust than iterative methods
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  l_tol = 1e-8
  l_max_its = 100
  nl_max_its = 15

  # SEAS-specific adaptive time stepping based on slip rate
  [TimeStepper]
    type = SEASAdaptiveDT
    max_slip_rate_pp = max_slip_rate
    Dc = 0.004                  # Critical slip distance (m)
    C = 0.5                     # Safety factor
    dt_min = 1e-3               # Minimum dt = 1 ms
    dt_max = 1e6                # Maximum dt ~ 11.5 days
    V_seismic = 1e-3            # Seismic threshold (1 mm/s)
    dt_seismic = 0.01           # Target dt during seismic (10 ms)
    initial_dt = 100            # Initial dt = 100 s (solver-stable)
    growth_factor = 1.2         # Gradual growth during interseismic
  []

  # Run for 1 year initially for testing
  end_time = 3.1557e7           # 1 year in seconds
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [exodus]
    type = Exodus
    execute_on = 'INITIAL TIMESTEP_END'
    time_step_interval = 1000
  []
  [console]
    type = Console
    output_linear = false
  []
[]
