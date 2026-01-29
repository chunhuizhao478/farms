# Unit Test: StateEvolutionKernel - Transient Evolution
#
# Tests the state evolution equation: dθ/dt = 1 - V*θ/Dc
#
# Analytical solution for constant V:
#   θ(t) = Dc/V + (θ0 - Dc/V) * exp(-V*t/Dc)
#
# Test case: Start away from steady state
#   V = 1e-6 m/s (constant)
#   Dc = 0.004 m
#   θ0 = 1e6 s (away from steady state Dc/V = 4000 s)
#
# After t = 0.01 s:
#   θ(0.01) = 4000 + (1e6 - 4000) * exp(-1e-6 * 0.01 / 0.004)
#           = 4000 + 996000 * exp(-2.5e-6)
#           ≈ 999997.5 s (very slow evolution)
#
# After t = 10000 s with larger V = 1e-3:
#   θ(10000) = 4 + (1e6 - 4) * exp(-1e-3 * 10000 / 0.004)
#            = 4 + 999996 * exp(-2500) ≈ 4 s (steady state)

Dc = 0.004
V_const = 1e-3
theta0 = 1000.0
theta_ss = ${fparse Dc / V_const}  # = 4 s

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 1
[]

[Variables]
  [state_variable]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${theta0}
  []
[]

[AuxVariables]
  [slip_rate]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${V_const}
  []
[]

[Kernels]
  [time_derivative]
    type = TimeDerivative
    variable = state_variable
  []
  [state_evolution]
    type = StateEvolutionKernel
    variable = state_variable
    slip_rate = slip_rate
    Dc = ${Dc}
    evolution_law = aging
  []
[]

[Functions]
  [analytical_theta]
    type = ParsedFunction
    expression = '${theta_ss} + (${theta0} - ${theta_ss}) * exp(-${V_const} * t / ${Dc})'
  []
[]

[Postprocessors]
  [numerical_state]
    type = ElementAverageValue
    variable = state_variable
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [analytical_state]
    type = FunctionValuePostprocessor
    function = analytical_theta
    point = '0.5 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [relative_error]
    type = RelativeDifferencePostprocessor
    value1 = numerical_state
    value2 = analytical_state
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
  dt = 1.0
  end_time = 100.0
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
