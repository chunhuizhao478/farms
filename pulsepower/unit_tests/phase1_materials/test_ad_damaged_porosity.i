# Unit test for ElkADPorousFlowDamagedPorosity
# Tests that the damaged porosity is computed correctly:
#   phi(d) = phi_0 + (1 - phi_0) * (1 - g(d))
# where g(d) = (1-d)^2
#
# Test cases with phi_0 = 0.008:
#   d = 0.0: g = 1.0, phi = 0.008 + 0.992 * (1 - 1.0) = 0.008
#   d = 0.5: g = 0.25, phi = 0.008 + 0.992 * (1 - 0.25) = 0.008 + 0.744 = 0.752
#   d = 1.0: g = 0.0, phi = 0.008 + 0.992 * (1 - 0.0) = 1.0 (clamped to upper bound)

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
  [porosity_aux]
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
  [get_porosity]
    type = ADMaterialRealAux
    variable = porosity_aux
    property = PorousFlow_porosity_qp_damaged
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Materials]
  [damaged_porosity]
    type = ElkADPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = 0.008
    porosity_lower_bound = 0.008
    porosity_upper_bound = 0.999
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
  # Expected values:
  # At x=0.1667 (d≈0.1667): g=(1-0.1667)^2=0.694, phi=0.008+0.992*0.306=0.312
  # At x=0.5 (d=0.5): g=0.25, phi=0.008+0.992*0.75=0.752
  # At x=0.8333 (d≈0.8333): g=(1-0.8333)^2=0.0278, phi=0.008+0.992*0.972=0.972
  [porosity_left]
    type = PointValue
    variable = porosity_aux
    point = '0.1667 0.05 0'
  []
  [porosity_mid]
    type = PointValue
    variable = porosity_aux
    point = '0.5 0.05 0'
  []
  [porosity_right]
    type = PointValue
    variable = porosity_aux
    point = '0.8333 0.05 0'
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
