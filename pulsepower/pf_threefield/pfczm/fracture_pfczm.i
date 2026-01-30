# =============================================================================
# PF-CZM Fracture Sub-Application for Three-Field Porodynamics
# Based on Wu's Unified Phase-Field Theory (JMPS 2017)
# =============================================================================
# Key changes from AT1:
#   1. Crack geometric function: alpha(d) = 2d - d^2 (vs d for AT1)
#   2. Normalization constant: c0 = pi (vs 8/3 for AT1)
#   3. Degradation function: Rational form with explicit failure strength
#   4. New parameters: psic, a2, a3, p_deg, eta
# =============================================================================

# Initial damage box 1 (horizontal notch)
bottom_left1 = '-0.0025 -3e-4 0'
top_right1 = '0.0025 3e-4 0'

# Initial damage box 2 (vertical notch)
bottom_left2 = '-3e-4 -0.0025 0'
top_right2 = '3e-4 0.0025 0'

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = '../mesh/fieldscale_test1_2d.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.01 0.01 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node = true
  []
  [./subdomain_id]
    type = SubdomainBoundingBoxGenerator
    bottom_left = ${bottom_left1}
    top_right = ${top_right1}
    location = INSIDE
    block_id = 1
    input = extranodeset1
  []
  [./subdomain_id2]
    type = SubdomainBoundingBoxGenerator
    bottom_left = ${bottom_left2}
    top_right = ${top_right2}
    location = INSIDE
    block_id = 1
    input = subdomain_id
  []
[]

[Variables]
  [d]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [mesh_size]
    order = CONSTANT
    family = MONOMIAL
  []
  [initial_damage_aux]
    family = LAGRANGE
    order = FIRST
  []
[]

[AuxKernels]
  [define_initial_damage_block1]
    type = ConstantAux
    variable = initial_damage_aux
    value = 0.9
    block = 1
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [define_initial_damage_block0]
    type = ConstantAux
    variable = initial_damage_aux
    value = 0
    block = '4 5'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Bounds]
  [irreversibility_first_step]
    type = VariableConstantIrreversibleBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
    bound_value = initial_damage_aux
  []
  [upper]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1
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
[]

[Materials]
  # =============================================================================
  # PF-CZM Material Properties
  # Note: xi and c0 are computed internally by CrackGeometricFunction
  # based on alpha(d) = 2d - d^2 (gives xi=2, c0=pi)
  # =============================================================================
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc psic'
    prop_values = '${l} ${Gc_const} ${psic}'
  []

  # =============================================================================
  # Crack Geometric Function: alpha(d) = 2d - d^2
  # For PF-CZM with xi = 2 (optimal for quasi-brittle materials)
  # This gives: c0 = pi, Du = pi*l/2 (finite support)
  # =============================================================================
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = '2*d - d^2'
    phase_field = d
  []

  # =============================================================================
  # Degradation Function: Rational form for PF-CZM
  # g(d) = (1-d)^p / [(1-d)^p + a1*d*P(d)] * (1-eta) + eta
  # where a1 = Gc/(psic*c0*l/xi) and P(d) = 1 + a2*d*(1 + a3*d)
  #
  # For linear softening: p=2, a2=-0.5, a3=0
  # For Cornelissen softening: p=2, a2=1.3868, a3=0.6567
  # =============================================================================
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+(Gc/psic*xi/c0/l)*d*(1+a2*d+a2*a3*d^2))*(1-eta)+eta
    phase_field = d
    material_property_names = 'Gc psic xi c0 l'
    parameter_names = 'p a2 a3 eta'
    parameter_values = '${p_deg} ${a2} ${a3} ${eta}'
  []

  # =============================================================================
  # Free Energy Functional
  # psi = psi_surface + psi_elastic
  #     = alpha(d)*Gc/(c0*l) + 0.5*g(d)*psie_active
  # Note: 0.5 factor is needed because psie_active already contains factor of 2
  # =============================================================================
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l+0.5*g*psie_active'
    coupled_variables = 'd psie_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    derivative_order = 1
  []

  # =============================================================================
  # Dissipated Energy Density for Energy Balance
  # Psi_d = Gc/c0 * (alpha/l + l*|grad(d)|^2)
  # =============================================================================
  [dissipated_energy]
    type = CrackDissipatedEnergyDensity
    phase_field = d
    alpha_name = alpha
    Gc_name = Gc
    l_name = l
    c0_name = c0
  []
[]

[AuxVariables]
  [dissipated_energy_density]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [dissipated_energy_density]
    type = ADMaterialRealAux
    variable = dissipated_energy_density
    property = dissipated_energy_density
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors]
  [dissipated_energy_dynamic]
    type = ADElementIntegralMaterialProperty
    mat_prop = dissipated_energy_density
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [dissipated_energy_first_step]
    type = FirstStepElementIntegralVariablePostprocessor
    variable = dissipated_energy_density
    execute_on = 'TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero -snes_type'
  petsc_options_value = 'gmres     hypre    boomeramg      True                       vinewtonrsls'

  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = false
  print_linear_residuals = false
[]
