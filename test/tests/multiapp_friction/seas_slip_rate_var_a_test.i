# Unit Test: SEASSlipRateVarAAux
#
# Tests slip rate computation with spatially-varying 'a' parameter.
# Uses the same logic as SEASSlipRateAux but reads 'a' from a coupled variable.
#
# At steady state with:
#   V = Vinit = 1e-9 m/s
#   θ = Dc/V0 = 4000 s
#   a = 0.025 (amax, VS region)
#
# The slip rate should remain approximately constant.

a_param = 0.025
b_param = 0.015
Dc = 0.004
f0 = 0.6
V0 = 1e-6
sigma_n = 50e6
eta = 4.634e6
tau_pre = 26546122.0
Vinit = 1e-9
theta0 = ${fparse Dc / V0}

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 1000
[]

[Variables]
  [dummy]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
[]

[AuxVariables]
  [traction]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${tau_pre}
  []
  [state_variable]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${theta0}
  []
  [a_var]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${a_param}
  []
  [slip_rate]
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
  [solve_slip_rate]
    type = SEASSlipRateVarAAux
    variable = slip_rate
    traction = traction
    state_variable = state_variable
    a_var = a_var
    b = ${b_param}
    Dc = ${Dc}
    f0 = ${f0}
    V0 = ${V0}
    sigma_n = ${sigma_n}
    eta = ${eta}
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Postprocessors]
  [computed_V]
    type = ElementAverageValue
    variable = slip_rate
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [traction_val]
    type = ElementAverageValue
    variable = traction
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [a_val]
    type = ElementAverageValue
    variable = a_var
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  num_steps = 1
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
