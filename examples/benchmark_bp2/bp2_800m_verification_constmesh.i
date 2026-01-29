# SCEC SEAS Benchmark Problem BP2-QD
# 2D Antiplane Shear with Rate-State Friction
# Staggered Solver Approach
#
# Mesh size: 800 m
#
# Reference: SEAS Benchmark Problem BP2-QD specification
# https://strike.scec.org/cvws/seas/
#
# Parameters from Table 1:
# ρ = 2670 kg/m³, cs = 3.464 km/s, μ = 32.04 GPa
# σn = 50 MPa, a0 = 0.010, amax = 0.025, b0 = 0.015
# Dc = 0.004 m, Vp = Vinit = 1e-9 m/s, V0 = 1e-6 m/s, f0 = 0.6
# H = 15 km, h = 3 km, Wf = 40 km

# Physical constants
mu = 32.04e9        # Shear modulus (Pa)
cs = 3464.0         # Shear wave speed (m/s)
# Note: ρ = 2670 kg/m³ (not used in quasi-dynamic formulation)

# Radiation damping: η = μ / (2 * cs)
eta = ${fparse mu / (2 * cs)}

# Rate-state parameters
a0 = 0.010          # VW region (near surface)
amax = 0.025        # VS region (deep)
b0 = 0.015          # Constant b
f0 = 0.6            # Reference friction coefficient
V0 = 1e-6           # Reference slip rate (m/s)
Dc = 0.004          # Critical slip distance (m)
sigma_n = 50e6      # Normal stress (Pa)

# Loading
Vp = 1e-9           # Plate rate (m/s)
Vinit = 1e-9        # Initial slip rate (m/s)

# Geometry (in meters)
H = 15000           # Depth extent of uniform VW region (15 km)
h = 3000            # Width of VW-VS transition zone (3 km)
Wf = 40000          # Width of rate-and-state fault (40 km)
L_x = 100000        # Half-width of domain (100 km)
L_z = 100000        # Total domain depth (100 km) - must extend below Wf for backslip

# Pre-stress calculation (from Eq. 11 in BP2 spec)
# τ⁰ = σn*amax*sinh⁻¹[Vinit/(2V0)*exp((f0 + b0*ln(V0/Vinit))/amax)] + η*Vinit
# τ⁰ ≈ 26.5461 MPa
tau0 = 26546122.0   # Initial shear stress (Pa)

# Mesh size
dz = 800            # Element size (m)
nz = ${fparse L_z / dz}  # Number of elements along depth = 125 (extended for backslip)
nx = ${fparse L_x / dz}  # Number of elements each side = 125

[Mesh]
  [left_block]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${nz}
    xmin = -${L_x}
    xmax = 0
    ymin = -${L_z}          # Deep boundary (100 km depth, negative y)
    ymax = 0                # Free surface at y = 0
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
    nx = ${nx}
    ny = ${nz}
    xmin = 0
    xmax = ${L_x}
    ymin = -${L_z}          # Deep boundary (100 km depth, negative y)
    ymax = 0                # Free surface at y = 0
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
  # Slip (accumulated displacement jump)
  [slip]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
  []
  # Slip rate (solved from traction balance)
  [slip_rate]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = ${Vinit}
  []
  # State variable - initialized using function
  [state_variable]
    order = CONSTANT
    family = MONOMIAL
  []
  # Elastic traction (from displacement field)
  [traction]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = ${tau0}
  []
  # Rate-state parameter a (depth-dependent)
  [a_param]
    order = CONSTANT
    family = MONOMIAL
  []
  # Friction stress (for output comparison with benchmark)
  # Computed from friction law: τ = σn * f(V, θ) + η * V
  [friction_stress]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  # Depth-dependent a(z) from Eq. 8
  # a = a0 for 0 <= z < H
  # a = a0 + (amax - a0)*(z - H)/h for H <= z < H + h
  # a = amax for H + h <= z < Wf
  # Note: y is negative (depth), so use abs(y) for depth z
  [a_func]
    type = ParsedFunction
    expression = 'if(abs(y) < ${H}, ${a0}, if(abs(y) < ${H} + ${h}, ${a0} + (${amax} - ${a0}) * (abs(y) - ${H}) / ${h}, ${amax}))'
  []

  # Initial state variable from Eq. 12 (using custom C++ function for accuracy)
  # θ(z,0) = (Dc/V0) * exp{(a/b) * ln[2V0/Vinit * sinh((τ⁰ - η*Vinit)/(a*σn))] - f0/b}
  [theta_func]
    type = BP2InitialStateFunction
    Dc = ${Dc}
    V0 = ${V0}
    Vinit = ${Vinit}
    b = ${b0}
    f0 = ${f0}
    tau0 = ${tau0}
    sigma_n = ${sigma_n}
    eta = ${eta}
    a0 = ${a0}
    amax = ${amax}
    H = ${H}
    h = ${h}
    depth_direction = 1  # y is the depth direction
  []

  # Far-field displacement functions
  # Standard convention: left side moves -z, right side moves +z
  [far_field_left]
    type = ParsedFunction
    expression = '-0.5 * ${Vp} * t'
  []
  [far_field_right]
    type = ParsedFunction
    expression = '0.5 * ${Vp} * t'
  []
[]

[ICs]
  [state_ic]
    type = FunctionIC
    variable = state_variable
    function = theta_func
  []
  [a_param_ic]
    type = FunctionIC
    variable = a_param
    function = a_func
  []
[]

[AuxKernels]
  # TIMESTEP_BEGIN: Update slip and state based on OLD values

  # 1. Update slip by integrating slip rate
  [update_slip]
    type = SEASSlipAux
    variable = slip
    slip_rate = slip_rate
    execute_on = 'TIMESTEP_BEGIN'
    boundary = fault
  []

  # 2. Update state variable using aging law
  [update_state]
    type = SEASStateAux
    variable = state_variable
    slip_rate = slip_rate
    Dc = ${Dc}
    execute_on = 'TIMESTEP_BEGIN'
    boundary = fault
  []

  # TIMESTEP_END: Compute traction, solve for new slip_rate

  # 3. Read elastic traction from interface material
  [read_traction]
    type = SEASTractionAux
    variable = traction
    traction_property = elastic_traction
    tau_pre = ${tau0}
    execute_on = 'TIMESTEP_END'
    boundary = fault
  []

  # 4. Solve for slip rate from traction balance (with variable a)
  #    Uses backslip loading: V = Vp for z > H+h (VS region)
  [solve_slip_rate]
    type = SEASSlipRateVarAAux
    variable = slip_rate
    traction = traction
    state_variable = state_variable
    a_var = a_param
    b = ${b0}
    Dc = ${Dc}
    f0 = ${f0}
    V0 = ${V0}
    sigma_n = ${sigma_n}
    eta = ${eta}
    # Backslip loading per BP2 Eq. 9: V(z,t) = Vp for z >= Wf
    # Combined with far-field BCs (±Vp*t/2) for correct loading
    use_backslip = true
    Vp = ${Vp}
    backslip_depth = ${Wf}  # 40 km per BP2 specification
    depth_direction = 1     # y is depth
    execute_on = 'TIMESTEP_END'
    boundary = fault
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
    sigma = 6.0    # Reduced (Tandem default)
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
    slip_prescribed = slip
    shear_modulus = shear_modulus
    penalty = 1e6    # Reduced penalty (Tandem uses ~3 for BR2)
    epsilon = 1.0
    sigma = 6.0      # Reduced sigma
  []
[]

[Materials]
  [elastic]
    type = GenericConstantMaterial
    prop_names = 'shear_modulus'
    prop_values = '${mu}'
    block = '1 2'
  []
  [elastic_traction]
    type = DGElasticTractionMaterial
    displacement = w
    slip_prescribed = slip
    shear_modulus = ${mu}
    sigma = 6.0
    penalty = 1e6
    traction_scale = 1.0
    traction_name = elastic_traction
    # Stiffness mode: τ = K × (Vp×t - slip) where K = μ/Wf
    # This is the standard quasi-static SEAS formulation
    traction_mode = stiffness
    fault_depth = ${Wf}       # 40 km
    plate_rate = ${Vp}        # 1e-9 m/s
    boundary = fault
  []
[]

[BCs]
  # Far-field Dirichlet BCs - plate motion drives loading
  [left_far_field]
    type = DGElasticityDirichletBC
    variable = w
    boundary = left_left
    function = far_field_left   # w = -Vp*t/2
    epsilon = 1.0
    sigma = 200.0
    shear_modulus = shear_modulus
  []
  [right_far_field]
    type = DGElasticityDirichletBC
    variable = w
    boundary = right_right
    function = far_field_right  # w = +Vp*t/2
    epsilon = 1.0
    sigma = 200.0
    shear_modulus = shear_modulus
  []

  # Free surface at y = 0 (ymax) - traction free
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

  # Deep boundary (y = -L_z, ymin) - consistent with far-field plate motion
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
    sigma = 100.0
    shear_modulus = shear_modulus
  []
[]

[Postprocessors]
  # Output at z = 0 km (free surface, y = 0) - VW region
  [slip_z0]
    type = PointValue
    variable = slip
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z0]
    type = PointValue
    variable = slip_rate
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z0]
    type = PointValue
    variable = traction
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z0]
    type = PointValue
    variable = state_variable
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Output at z = 4.8 km (y = -4800 m) - VW region
  [slip_z4_8]
    type = PointValue
    variable = slip
    point = '0 -4800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z4_8]
    type = PointValue
    variable = slip_rate
    point = '0 -4800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z4_8]
    type = PointValue
    variable = traction
    point = '0 -4800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z4_8]
    type = PointValue
    variable = state_variable
    point = '0 -4800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Output at z = 12 km (y = -12000 m) - VW region
  [slip_z12]
    type = PointValue
    variable = slip
    point = '0 -12000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z12]
    type = PointValue
    variable = slip_rate
    point = '0 -12000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z12]
    type = PointValue
    variable = traction
    point = '0 -12000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z12]
    type = PointValue
    variable = state_variable
    point = '0 -12000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Output at z = 16.8 km (y = -16800 m) - transition zone
  [slip_z16_8]
    type = PointValue
    variable = slip
    point = '0 -16800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z16_8]
    type = PointValue
    variable = slip_rate
    point = '0 -16800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z16_8]
    type = PointValue
    variable = traction
    point = '0 -16800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z16_8]
    type = PointValue
    variable = state_variable
    point = '0 -16800 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Output at z = 24 km (y = -24000 m) - VS region
  [slip_z24]
    type = PointValue
    variable = slip
    point = '0 -24000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [slip_rate_z24]
    type = PointValue
    variable = slip_rate
    point = '0 -24000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [shear_stress_z24]
    type = PointValue
    variable = traction
    point = '0 -24000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [state_z24]
    type = PointValue
    variable = state_variable
    point = '0 -24000 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # Maximum slip rate for adaptive time stepping
  [max_slip_rate]
    type = ElementExtremeValue
    variable = slip_rate
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
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
  # Uses dt = C * Dc / V_max for stability
  [TimeStepper]
    type = SEASAdaptiveDT
    max_slip_rate_pp = max_slip_rate
    Dc = 0.004                  # Critical slip distance (m)
    C = 0.1                     # Safety factor (very conservative)
    dt_min = 1e-3               # Minimum dt = 1 ms
    dt_max = 1e5                # Maximum dt ~ 1 day
    V_seismic = 1e-3            # Seismic threshold (1 mm/s)
    dt_seismic = 0.01           # Target dt during seismic (10 ms)
    initial_dt = 100            # Initial dt = 100 s
    growth_factor = 1.05        # Very conservative growth
  []

  # Run for 1 year initially for testing (31557600 seconds)
  end_time = 820.482e7  # 260 year in seconds for testing
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
