# SEASAdaptiveDT Unit Test
# Test adaptive time stepping for SEAS simulations during aseismic phase
#
# This test verifies that the time stepper:
# 1. Uses initial_dt for the first step
# 2. Adapts dt based on max_slip_rate (dt = C * Dc / V)
# 3. Limits growth by growth_factor during aseismic phase
#
# With V = 1e-9 m/s (aseismic):
# dt_target = C * Dc / V = 0.5 * 0.004 / 1e-9 = 2e6 seconds
# But growth_factor = 1.2 limits: 100 -> 120 -> 144 -> 172.8 -> 207.4

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 2
[]

[Variables]
  [u]
    initial_condition = 0
  []
[]

[Kernels]
  [time]
    type = TimeDerivative
    variable = u
  []
  [source]
    type = BodyForce
    variable = u
    value = 1e-12
  []
[]

[Functions]
  # Aseismic slip rate: V = 1e-9 m/s
  # Expected dt = C * Dc / V = 0.5 * 0.004 / 1e-9 = 2e6 s
  # But limited by growth_factor = 1.2
  [slip_rate_function]
    type = ParsedFunction
    expression = '1e-9'
  []
[]

[Postprocessors]
  [max_slip_rate]
    type = FunctionValuePostprocessor
    function = slip_rate_function
    execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
  []
  [dt_value]
    type = TimestepSize
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_abs_tol = 1e-10

  [TimeStepper]
    type = SEASAdaptiveDT
    max_slip_rate_pp = max_slip_rate
    Dc = 0.004              # 4 mm (BP2 value)
    C = 0.5                 # Safety factor
    dt_min = 1.0            # 1 second minimum
    dt_max = 1e8            # Maximum
    V_seismic = 1e-3        # 1 mm/s threshold
    dt_seismic = 0.1        # 0.1 s during seismic
    initial_dt = 100.0      # Start with 100 s
    growth_factor = 1.2     # 20% growth limit
  []

  num_steps = 5
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
