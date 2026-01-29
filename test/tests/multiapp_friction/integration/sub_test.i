# MultiApp Integration Test - Friction SubApp
#
# Simplified test to verify MultiApp coupling between elasticity and friction.

mu = 32.04e9
cs = 3464.0
eta = ${fparse mu / (2 * cs)}

a_param = 0.025
b_param = 0.015
f0 = 0.6
V0 = 1e-6
Dc = 0.004
sigma_n = 50e6
Vinit = 1e-9
tau_pre = 26546122.0
theta0 = ${fparse Dc / Vinit}

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 1000
[]

[Variables]
  [state_variable]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${theta0}
  []
[]

[AuxVariables]
  [traction_received]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
  [slip_rate]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${Vinit}
  []
  [slip]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
[]

[Kernels]
  [state_time]
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

[AuxKernels]
  [solve_slip_rate]
    type = FrictionSlipRateAux
    variable = slip_rate
    traction = traction_received
    state_variable = state_variable
    a = ${a_param}
    b = ${b_param}
    Dc = ${Dc}
    f0 = ${f0}
    V0 = ${V0}
    sigma_n = ${sigma_n}
    eta = ${eta}
    tau_pre = ${tau_pre}
    execute_on = 'TIMESTEP_BEGIN'
  []
  [integrate_slip]
    type = SlipIntegrationAux
    variable = slip
    slip_rate = slip_rate
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
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  dt = 1e6
  num_steps = 3
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
