# Unit Test: ADPFFViscousResistance Kernel - Rate-Independent Limit
# Tests that eta = 0 recovers rate-independent phase field behavior
# Reference: Hofacker & Miehe (2012), IJNME 93:276-301
#
# With eta = 0, the viscous term vanishes and we should get
# the same results as the quasi-static phase field model.

# Test parameters - rate-independent limit
l_val = 0.01  # regularization length [m]
Gc_val = 100  # fracture toughness [J/m^2]
eta_val = 0   # ZERO viscosity - rate-independent limit

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
  [d]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [d_ic]
    type = ConstantIC
    variable = d
    value = 0.5
  []
[]

[AuxVariables]
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [psie_active_kernel]
    type = ConstantAux
    variable = psie_active
    value = 200
  []
[]

[Kernels]
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
  []
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
  # With eta = 0, this kernel should have no effect
  [viscous]
    type = ADPFFViscousResistance
    variable = d
    viscosity = eta
  []
[]

[Materials]
  [fracture_props]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc eta'
    prop_values = '${l_val} ${Gc_val} ${eta_val}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta_deg)+eta_deg
    phase_field = d
    parameter_names = 'p eta_deg'
    parameter_values = '2 1e-6'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
  []
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l + g*psie_active'
    coupled_variables = 'd psie_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    derivative_order = 1
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = d
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = d
    boundary = right
    value = 1
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  dt = 0.01
  num_steps = 5
  end_time = 0.05

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Postprocessors]
  [d_avg]
    type = ElementAverageValue
    variable = d
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [d_max]
    type = ElementExtremeValue
    variable = d
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [d_min]
    type = ElementExtremeValue
    variable = d
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
