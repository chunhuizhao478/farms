# SEASAdaptiveDT Seismic Phase Test
# Test adaptive time stepping during seismic phase (high slip rate)
#
# This test verifies that when slip rate exceeds V_seismic,
# the time stepper uses small time steps (dt_seismic or C*Dc/V)
#
# With V = 0.1 m/s (seismic, >> V_seismic = 1e-3):
# dt = min(C * Dc / V, dt_seismic) = min(0.5 * 0.004 / 0.1, 0.1) = min(0.02, 0.1) = 0.02 s

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
    value = 1e-6
  []
[]

[Functions]
  # High slip rate simulating seismic phase
  # V = 0.1 m/s >> V_seismic = 1e-3 m/s
  [slip_rate_function]
    type = ParsedFunction
    expression = '0.1'
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
    dt_min = 0.001          # 1 ms minimum
    dt_max = 1e8            # Maximum
    V_seismic = 1e-3        # 1 mm/s threshold
    dt_seismic = 0.1        # 0.1 s during seismic
    initial_dt = 100.0      # Start with 100 s
    growth_factor = 1.2
  []

  # Expected dt = 0.02 s during seismic phase
  num_steps = 5
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
