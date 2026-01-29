# Unit Test: SlipIntegrationAux
#
# Tests slip integration: S_new = S_old + dt * V
#
# For constant V:
#   S(t) = V * t
#
# Test case:
#   V = 1e-9 m/s (constant)
#   S0 = 0 m
#   After t = 1e8 s: S = 1e-9 * 1e8 = 0.1 m

V_const = 1e-9

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 1
[]

[Variables]
  [dummy]
    # Need at least one variable for MOOSE
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
[]

[AuxVariables]
  [slip_rate]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${V_const}
  []
  [slip]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
[]

[Kernels]
  [null]
    type = NullKernel
    variable = dummy
  []
[]

[AuxKernels]
  [integrate_slip]
    type = SlipIntegrationAux
    variable = slip
    slip_rate = slip_rate
    execute_on = 'TIMESTEP_END'
  []
[]

[Functions]
  [analytical_slip]
    type = ParsedFunction
    expression = '${V_const} * t'
  []
[]

[Postprocessors]
  [numerical_slip]
    type = ElementAverageValue
    variable = slip
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [analytical_slip_pp]
    type = FunctionValuePostprocessor
    function = analytical_slip
    point = '0.5 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [absolute_error]
    type = DifferencePostprocessor
    value1 = numerical_slip
    value2 = analytical_slip_pp
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  dt = 1e7
  end_time = 1e8
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
