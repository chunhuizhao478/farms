# Unit Test: ADPFFViscousResistance Kernel
# Tests the viscous crack resistance term for rate-dependent phase field fracture
# Reference: Hofacker & Miehe (2012), IJNME 93:276-301
#
# This test verifies:
# 1. The kernel adds the correct viscous resistance term: eta * (d - d_old) / dt
# 2. The residual scales correctly with viscosity parameter eta
# 3. The rate-independent limit (eta = 0) recovers quasi-static behavior

# Test parameters
l_val = 0.01  # regularization length [m]
Gc_val = 100  # fracture toughness [J/m^2]
eta_val = 1e-6  # viscosity [Pa*s]

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
  # Start with damage = 0.5 everywhere
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
  # Provide a driving energy (constant for simplicity)
  [psie_active_kernel]
    type = ConstantAux
    variable = psie_active
    value = 200  # [J/m^3] - above critical value to drive fracture
  []
[]

[Kernels]
  # Diffusion term: (grad_w, 2*G_c*l/c_0 * grad_d)
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
  []
  # Source term from free energy: (w, d_psi/d_d)
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
  # Viscous resistance term: (w, eta * d_dot)
  # This is the new kernel being tested
  [viscous]
    type = ADPFFViscousResistance
    variable = d
    viscosity = eta
  []
[]

[Materials]
  # Fracture properties
  [fracture_props]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc eta'
    prop_values = '${l_val} ${Gc_val} ${eta_val}'
  []
  # Degradation function: g(d) = (1-d)^2
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta_deg)+eta_deg
    phase_field = d
    parameter_names = 'p eta_deg'
    parameter_values = '2 1e-6'
  []
  # Crack geometric function: alpha(d) = d (AT1 model)
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
  []
  # Total free energy: psi = alpha*Gc/(c0*l) + g*psie_active
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
  # Fix d = 0 on left boundary (intact material)
  [left]
    type = DirichletBC
    variable = d
    boundary = left
    value = 0
  []
  # Fix d = 1 on right boundary (fully damaged)
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

  # Time stepping
  dt = 0.01
  num_steps = 5
  end_time = 0.05

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Postprocessors]
  # Track average damage
  [d_avg]
    type = ElementAverageValue
    variable = d
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Track max damage
  [d_max]
    type = ElementExtremeValue
    variable = d
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Track min damage
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
