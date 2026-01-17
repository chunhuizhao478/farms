# Unit test for ElkADPorousFlowDamagedBiotCoefficient
# Tests that the damaged Biot coefficient is computed correctly:
#   alpha(d) = 1 - g(d) * K0 / Ks
# where g(d) = max((1-d)^2, eta)
#
# Test cases:
#   d = 0.0: g = 1.0, alpha = 1 - 1.0 * 50e9 / 50e9 = 0.0
#   d = 0.5: g = 0.25, alpha = 1 - 0.25 * 50e9 / 50e9 = 0.75
#   d = 1.0: g = eta, alpha = 1 - eta * 50e9 / 50e9 ≈ 1.0

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  ny = 1
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0.1
[]

[Variables]
  [dummy]
  []
[]

[AuxVariables]
  [d]
    family = LAGRANGE
    order = FIRST
  []
  [biot_coeff_aux]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  # Set damage values: d=0 at x=0, d=0.5 at x=0.5, d=1.0 at x=1.0
  [damage_ic]
    type = FunctionIC
    variable = d
    function = 'x'
  []
[]

[AuxKernels]
  [get_biot_coeff]
    type = ADMaterialRealAux
    variable = biot_coeff_aux
    property = biot_coefficient_damaged
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Materials]
  [damaged_biot_coefficient]
    type = ElkADPorousFlowDamagedBiotCoefficient
    phase_field = d
    solid_bulk_compliance = 2e-11  # 1/K0 = 1/50e9
    grain_bulk_modulus = 50e9      # Ks
    minimum_degradation = 1e-8
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = dummy
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = dummy
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = dummy
    boundary = right
    value = 1
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Postprocessors]
  [biot_coeff_left]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.1667 0.05 0'
  []
  [biot_coeff_mid]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.5 0.05 0'
  []
  [biot_coeff_right]
    type = PointValue
    variable = biot_coeff_aux
    point = '0.8333 0.05 0'
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
