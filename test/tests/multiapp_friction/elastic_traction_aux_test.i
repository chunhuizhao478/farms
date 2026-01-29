# Unit Test: ElasticTractionAux
#
# Tests elastic traction computation: τ = μ * ∂w/∂n
#
# For a linear displacement field w = x (in 2D):
#   ∂w/∂x = 1
#   τ = μ * 1 = μ
#
# Test case:
#   μ = 32.04e9 Pa
#   w = x (linear displacement)
#   Expected traction at any point: τ = 32.04e9 Pa

mu = 32.04e9

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[Variables]
  [w]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [traction]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [linear_disp]
    type = ParsedFunction
    expression = 'x'
  []
[]

[ICs]
  [w_ic]
    type = FunctionIC
    variable = w
    function = linear_disp
  []
[]

[Kernels]
  [diffusion]
    type = Diffusion
    variable = w
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = w
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = w
    boundary = right
    value = 1
  []
[]

[AuxKernels]
  [compute_traction]
    type = ElasticTractionAux
    variable = traction
    displacement = w
    shear_modulus = ${mu}
    normal_component = 0  # x-direction
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors]
  [avg_traction]
    type = ElementAverageValue
    variable = traction
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [expected_traction]
    type = FunctionValuePostprocessor
    function = '${mu}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [relative_error]
    type = RelativeDifferencePostprocessor
    value1 = avg_traction
    value2 = expected_traction
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
