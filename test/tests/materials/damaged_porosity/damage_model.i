# Regression test for ElkPorousFlowDamagedPorosity, DEFAULT (damage) update law.
#
# Locks in the existing damage-driven bounded-maximum porosity:
#   phi(d) = phi_0 + (1 - phi_0) * (1 - (1 - d)^2),   clamped to [lower, upper].
#
# The damage field d is imposed PIECEWISE-CONSTANT with breakpoints on the element
# boundaries (x = 1/3, 2/3) so each of the three elements carries a single, exact d
# value and the resulting porosity is element-uniform (hand-verifiable, no qp averaging).
#
# phi_0 = 0.008, lower = 0.008, upper = 0.999:
#   elem 0 (x in [0,1/3],   d = 1/6): g = (5/6)^2 = 0.69444  -> phi = 0.31111662
#   elem 1 (x in [1/3,2/3], d = 1/2): g = 0.25              -> phi = 0.752
#   elem 2 (x in [2/3,1],   d = 5/6): g = (1/6)^2 = 0.027779 -> phi = 0.97244334
# See EXPECTED_VALUES.md for the full derivation.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 3
    ny = 1
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 0.1
  []
[]

[Variables]
  # Trivial field so the Steady solve is non-empty.
  [u]
  []
[]

[AuxVariables]
  [d]
    family = MONOMIAL
    order = CONSTANT
  []
  [porosity_aux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  # Piecewise-constant damage; breakpoints coincide with element boundaries.
  [d_func]
    type = ParsedFunction
    expression = 'if(x < 0.333333333333, 0.166666666667, if(x < 0.666666666667, 0.5, 0.833333333333))'
  []
[]

[AuxKernels]
  [d_aux]
    type = FunctionAux
    variable = d
    function = d_func
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [porosity]
    type = MaterialRealAux
    variable = porosity_aux
    property = PorousFlow_porosity_qp_damaged
    execute_on = 'TIMESTEP_END'
  []
[]

[Kernels]
  [diff]
    type = Diffusion
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
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Materials]
  [porosity_damaged]
    type = ElkPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = 0.008
    porosity_lower_bound = 0.008
    porosity_upper_bound = 0.999
    # porosity_update_model omitted -> default 'damage'
  []
[]

[Postprocessors]
  [phi_elem0]
    type = PointValue
    variable = porosity_aux
    point = '0.1666667 0.05 0'
    execute_on = 'TIMESTEP_END FINAL'
  []
  [phi_elem1]
    type = PointValue
    variable = porosity_aux
    point = '0.5 0.05 0'
    execute_on = 'TIMESTEP_END FINAL'
  []
  [phi_elem2]
    type = PointValue
    variable = porosity_aux
    point = '0.8333333 0.05 0'
    execute_on = 'TIMESTEP_END FINAL'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
  execute_on = 'FINAL'
[]
