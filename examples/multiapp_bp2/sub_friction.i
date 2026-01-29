# SCEC SEAS Benchmark Problem BP2-QD
# MultiApp Architecture - Sub Application: Friction Solver
#
# This is the SubApp that solves rate-and-state friction on a 1D fault mesh.
# It receives traction from the Main App and sends slip back.
#
# The 1D mesh is along the fault depth direction (y in the main app).
# x coordinate in this SubApp corresponds to y (depth) in the main app.
#
# Equations solved:
#   State evolution (aging law): dtheta/dt = 1 - V*theta/Dc
#   Slip integration: ds/dt = V
#   Traction balance: tau = sigma_n * f(V,theta) + eta * V (solved for V)
#
# Reference: SCEC SEAS Benchmark Problem BP2-QD

# Physical constants
mu = 32.04e9        # Shear modulus (Pa)
cs = 3464.0         # Shear wave speed (m/s)

# Radiation damping: eta = mu / (2 * cs)
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
Vinit = 1e-9        # Initial slip rate (m/s)

# Geometry (in meters)
H = 15000           # Depth extent of uniform VW region (15 km)
h = 3000            # Width of VW-VS transition zone (3 km)
Wf = 40000          # Seismogenic depth (40 km) - fault ends here

# Pre-stress calculation (from Eq. 11 in BP2 spec)
tau0 = 26546122.0   # Initial shear stress (Pa)

# Mesh size
dz = 800            # Element size (m)
nz = ${fparse Wf / dz}  # Number of elements = 50 (only seismogenic zone)

# 2D Fault mesh positioned at x=0 to match MainApp fault location
# y coordinate matches MainApp depth direction
# y = 0 is free surface, y = -Wf is bottom of seismogenic zone
# Following Tandem: fault mesh only covers seismogenic zone (0 to Wf)
# Below Wf, the continuous elastic domain handles plate motion
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1             # Single element in x (fault is a line)
  ny = ${nz}         # Elements along depth (0 to -Wf)
  xmin = -0.1        # Thin strip around x=0
  xmax = 0.1
  ymin = -${Wf}      # Bottom of seismogenic zone (40 km depth)
  ymax = 0           # Free surface
[]

[Variables]
  # State variable theta (evolved implicitly via Kernel)
  [state_variable]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  # Traction received from MainApp (total traction = elastic + tau_pre)
  [traction_received]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${tau0}
  []
  # Slip rate (solved algebraically from traction balance)
  [slip_rate]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${Vinit}
  []
  # Accumulated slip
  [slip]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
  # Rate-state parameter a (depth-dependent)
  [a_param]
    order = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # Depth-dependent a(z) from Eq. 8
  # a = a0 for 0 <= z < H (near surface, VW)
  # a = a0 + (amax - a0)*(z - H)/h for H <= z < H + h (transition)
  # a = amax for H + h <= z < Wf (deep, VS)
  # Note: y is negative depth here, so use abs(y)
  [a_func]
    type = ParsedFunction
    expression = 'if(abs(y) < ${H}, ${a0}, if(abs(y) < ${H} + ${h}, ${a0} + (${amax} - ${a0}) * (abs(y) - ${H}) / ${h}, ${amax}))'
  []

  # Initial state variable from Eq. 12
  # theta(z,0) = (Dc/V0) * exp{(a/b) * ln[2V0/Vinit * sinh((tau0 - eta*Vinit)/(a*sigma_n))] - f0/b}
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
    depth_direction = 1  # y is the depth direction (matches MainApp)
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

[Kernels]
  # Time derivative of state variable
  [state_time]
    type = TimeDerivative
    variable = state_variable
  []
  # State evolution: dtheta/dt = 1 - V*theta/Dc
  [state_evolution]
    type = StateEvolutionKernel
    variable = state_variable
    slip_rate = slip_rate
    Dc = ${Dc}
    evolution_law = aging
  []
[]

[AuxKernels]
  # Solve for slip rate from traction balance (with variable a)
  # Following Tandem: fault only covers seismogenic zone (0 to Wf)
  # No backslip needed - the continuous domain below Wf handles plate motion
  # Execute at TIMESTEP_BEGIN before state evolution
  [solve_slip_rate]
    type = SEASSlipRateVarAAux
    variable = slip_rate
    traction = traction_received
    state_variable = state_variable
    a_var = a_param
    b = ${b0}
    Dc = ${Dc}
    f0 = ${f0}
    V0 = ${V0}
    sigma_n = ${sigma_n}
    eta = ${eta}
    # No backslip - fault only covers seismogenic zone (0 to Wf)
    use_backslip = false
    execute_on = 'TIMESTEP_BEGIN'
  []

  # Integrate slip: s += V * dt
  # Execute at TIMESTEP_END after state evolution
  # Uses slip_old_var coupled to self for proper state preservation in MultiApp
  [integrate_slip]
    type = SlipIntegrationAux
    variable = slip
    slip_rate = slip_rate
    slip_old_var = slip
    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  [avg_slip]
    type = ElementAverageValue
    variable = slip
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [avg_slip_rate]
    type = ElementAverageValue
    variable = slip_rate
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [avg_state]
    type = ElementAverageValue
    variable = state_variable
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [max_slip_rate]
    type = ElementExtremeValue
    variable = slip_rate
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [avg_traction]
    type = ElementAverageValue
    variable = traction_received
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Point values at z=0 (x=0 in SubApp)
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
  [state_z0]
    type = PointValue
    variable = state_variable
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [traction_z0]
    type = PointValue
    variable = traction_received
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_z0]
    type = PointValue
    variable = a_param
    point = '0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 20

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
  # Time stepping - inherited from MainApp via MultiApp
  # Note: This dt is just for initialization; actual dt comes from MainApp
[]

[Outputs]
  [csv]
    type = CSV
    file_base = sub_friction_out
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [console]
    type = Console
    output_linear = false
  []
[]
