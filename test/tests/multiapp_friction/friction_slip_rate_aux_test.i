# Unit Test: FrictionSlipRateAux
#
# Tests slip rate computation from traction balance:
#   τ = σn * f(V, θ) + η * V
#
# At steady state with:
#   V = V0 = 1e-6 m/s
#   θ = Dc/V0 = 4000 s
#
# The steady-state traction should be:
#   f_ss = a * asinh[V/(2*V0) * exp((f0 + b*ln(1))/a)]
#        = a * asinh[0.5 * exp(f0/a)]
#        = 0.025 * asinh[0.5 * exp(0.6/0.025)]
#        = 0.025 * asinh[0.5 * exp(24)]
#        ≈ 0.025 * 24 = 0.6 (approximately)
#
#   τ_ss = σn * f_ss + η * V
#        = 50e6 * 0.6 + 4.634e6 * 1e-6
#        ≈ 30e6 Pa
#
# Given this traction, the solver should return V ≈ 1e-6 m/s

a = 0.025
b = 0.015
Dc = 0.004
f0 = 0.6
V0 = 1e-6
sigma_n = 50e6
eta = 4.634e6

# At steady state: V = V0, θ = Dc/V0
V_expected = ${V0}
theta_ss = ${fparse Dc / V0}

# Pre-compute expected traction for steady state
# This is an approximation - we'll verify the solver finds V ≈ V0

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
  xmin = 0
  xmax = 1
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
    # Use steady-state traction based on expected V
    # τ = σn * f0 + η * V0 (simplified approximation)
    initial_condition = 30e6
  []
  [state_variable]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${theta_ss}
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
    type = FrictionSlipRateAux
    variable = slip_rate
    traction = traction
    state_variable = state_variable
    a = ${a}
    b = ${b}
    Dc = ${Dc}
    f0 = ${f0}
    V0 = ${V0}
    sigma_n = ${sigma_n}
    eta = ${eta}
    tau_pre = 0
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
  [state_val]
    type = ElementAverageValue
    variable = state_variable
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
