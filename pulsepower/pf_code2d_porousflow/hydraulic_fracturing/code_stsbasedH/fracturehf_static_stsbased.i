# =============================================================================
# Stress-Based Driving Force Static Fracture Sub-Application
# Based on Miehe & Mauthe (CMAME 2016) + Wu's PF-CZM (JMPS 2017)
# =============================================================================
#
# Key Features:
#   - Uses stress-based driving force H instead of strain energy psie_active
#   - PF-CZM crack geometric function: α(d) = 2d - d²
#   - Phase field evolution (Miehe Eq. 13): η*ḋ = (1-d)*H - [d - l²Δd]
#
# Free Energy Functional:
#   ψ = α(d)*Gc/(c0*l) + 0.5*(1-d)² * (Gc/(c0*l)) * H
#
# This gives the correct (1-d)*H driving term in the evolution equation.
#
# Note: H is the stress-based driving force from main app (Miehe Eq. 56):
#   D = ζ * < Σ_a (<σ̃_eff^a>_+ / σ_c)² - 1 >_+
#   H = max(H_old, D)  (history variable for irreversibility)
#
# =============================================================================

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/square_with_hole_quicktest.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.01 0.01 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
[]

[Variables]
  [d]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [d_ic_domain]
    type = ConstantIC
    variable = d
    value = 0
    block = domain
  []
  [d_ic_main_fractures]
    type = ConstantIC
    variable = d
    value = 0.95
    block = main_fractures
  []
  [d_ic_branch_fractures]
    type = ConstantIC
    variable = d
    value = 0.95
    block = branch_fractures
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  # History variable H transferred from main app (stress-based driving force)
  [H_driving]
    order = CONSTANT
    family = MONOMIAL
  []
  [Gc_var]
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
  [define_initial_damage_block_domain]
    type = ConstantAux
    variable = initial_damage_aux
    value = 0
    block = domain
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [define_initial_damage_block_main_fractures]
    type = ConstantAux
    variable = initial_damage_aux
    value = 0.95
    block = main_fractures
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [define_initial_damage_block_sub_fractures]
    type = ConstantAux
    variable = initial_damage_aux
    value = 0.95
    block = branch_fractures
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
  # Stress-Based PF-CZM Material Properties
  # =============================================================================

  # Note: xi and c0 are computed internally by CrackGeometricFunction
  # based on alpha(d) = 2d - d^2 (gives xi=2, c0=pi)
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
  # Free Energy Functional with Stress-Based Driving Force (Miehe Eq. 13)
  #
  # Miehe's evolution equation:
  #   η*ḋ = (1-d)*H - [d - l²Δd]
  #
  # The ADPFFSource kernel computes: residual = -∂ψ/∂d
  #
  # With ψ_elastic = 0.5*(1-d)² * (Gc/(c0*l)) * H:
  #   ∂ψ_elastic/∂d = -(1-d) * (Gc/(c0*l)) * H
  #   -∂ψ_elastic/∂d = (1-d) * (Gc/(c0*l)) * H  ✓
  #
  # This gives the correct (1-d)*H driving term from Miehe's formulation.
  # =============================================================================
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l + 0.5*(1-d)^2 * (Gc/c0/l) * H_driving'
    coupled_variables = 'd H_driving'
    material_property_names = 'alpha(d) Gc c0 l'
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
  petsc_options_value = 'gmres     hypre  boomeramg True vinewtonrsls'

  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = false
  print_linear_residuals = false
[]
