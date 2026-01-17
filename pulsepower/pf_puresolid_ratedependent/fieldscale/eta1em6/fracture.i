# Rate-Dependent Phase Field Fracture - Phase Field Input
# Based on Hofacker & Miehe (2012) - IJNME 93:276-301
# This file handles the phase field (damage) evolution with viscous regularization
#
# Key modification from rate-independent case:
# Added ADPFFViscousResistance kernel implementing: eta * d_dot
# where d_dot = (d - d_old) / dt (backward Euler)
#
# The complete phase field equation becomes:
#   (G_c/l)[d - l^2*Laplacian(d)] + eta*d_dot = 2(1-d)*H

#initial damage box 1
bottom_left1 = '0 -4e-4 0'
top_right1 = '0.0925 4e-4 0'

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = '../../../pf_code2d_puresolid_pulseloadmaxmin/2dmeshfile/fieldscale_test1_2d_amr_fieldscale_debug_enlargedrefine.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '1.0 1.0 0'
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
  [define_initial_damage_block1]
    type = ConstantAux
    variable = initial_damage_aux
    value = 0.9
    block = 1
  []
  [define_initial_damage_block0]
    type = ConstantAux
    variable = initial_damage_aux
    value = 0
    block = '4 5'
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
  # Diffusion term: (grad_w, 2*G_c*l/c_0 * grad_d)
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
  []
  # Source term: (w, d_psi/d_d)
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
  # NEW: Viscous resistance term for rate-dependent fracture
  # Implements: (w, eta * d_dot) = (w, eta/dt * (d - d_old))
  # Reference: Hofacker & Miehe (2012), Eq. (46)
  [viscous]
    type = ADPFFViscousResistance
    variable = d
    viscosity = eta
  []
[]

[Materials]
  # Fracture properties including viscosity
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc eta'
    prop_values = '${l} ${Gc_const} ${eta_viscosity}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta_deg)+eta_deg
    phase_field = d
    parameter_names = 'p eta_deg'
    parameter_values = '2 1e-6'
  []
  [crack_geometric] #AT1 model
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
  []
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l+g*psie_active'
    coupled_variables = 'd psie_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    derivative_order = 1
  []
  # Dissipated energy density material
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
