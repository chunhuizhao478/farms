# Minimal test: fixed dt=0.5 before t=3, dt=0.01 after t=3
# Using PiecewiseConstant function as timestep_limiting_function
# Expects output at sync_times: 1, 2, 3, 3.1, 3.2, 3.3, 3.4, 3.5

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 5
  []
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diff]
    type = ADDiffusion
    variable = u
  []
  [time]
    type = ADTimeDerivative
    variable = u
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = FunctionDirichletBC
    variable = u
    boundary = right
    function = 't'
  []
[]

[Functions]
  [dt_fn]
    type = PiecewiseConstant
    x = '0   3'
    y = '0.5 0.01'
    direction = LEFT_INCLUSIVE
  []
[]

[Postprocessors]
  [dt_pp]
    type = FunctionValuePostprocessor
    function = dt_fn
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  nl_rel_tol = 1e-8
  end_time = 3.5

  [TimeStepper]
    type = FunctionDT
    function = dt_fn
  []
[]

[Outputs]
  [exodus]
    type = Exodus
    sync_times = '1 2 3 3.1 3.2 3.3 3.4 3.5'
    sync_only = true
  []
  print_linear_residuals = false
[]
