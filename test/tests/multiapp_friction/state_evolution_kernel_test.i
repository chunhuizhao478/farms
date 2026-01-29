# Unit Test: StateEvolutionKernel
#
# Tests the state evolution equation: dθ/dt = 1 - V*θ/Dc
#
# Analytical solution for constant V:
#   θ(t) = Dc/V + (θ0 - Dc/V) * exp(-V*t/Dc)
#
# At steady state: θ_ss = Dc/V
#
# Test case:
#   V = 1e-9 m/s (constant)
#   Dc = 0.004 m
#   θ0 = 4e6 s (= Dc/V, starting at steady state)
#
# Expected: θ should remain at 4e6 s

Dc = 0.004
V_const = 1e-9
theta0 = ${fparse Dc / V_const}

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
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

[Postprocessors]
  [avg_state]
    type = ElementAverageValue
    variable = state_variable
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [expected_state]
    type = FunctionValuePostprocessor
    function = '${theta0}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [error]
    type = DifferencePostprocessor
    value1 = avg_state
    value2 = expected_state
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
  dt = 1e8
  num_steps = 10
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
