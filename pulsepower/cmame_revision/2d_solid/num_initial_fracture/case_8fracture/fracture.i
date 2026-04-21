# Case: 8 initial fractures — 4 stripes at 0 deg, 45 deg, 90 deg, 135 deg
# (tips every 45 deg around the origin). Must match the mesh in elasticity.i.

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../../2d_mesh/2d_mesh.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.01 0.01 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
  # Stripe 1: 0 deg (along +x)
  [./subdomain_id_0deg]
    type = OrientedSubdomainBoundingBoxGenerator
    input = extranodeset1
    center = '0 0 0'
    length = 0.005
    width = 6e-4
    height = 1
    length_direction = '1 0 0'
    width_direction = '0 1 0'
    block_id = 1
  []
  # Stripe 2: 45 deg
  [./subdomain_id_45deg]
    type = OrientedSubdomainBoundingBoxGenerator
    input = subdomain_id_0deg
    center = '0 0 0'
    length = 0.005
    width = 6e-4
    height = 1
    length_direction = '0.7071068 0.7071068 0'
    width_direction = '-0.7071068 0.7071068 0'
    block_id = 1
  []
  # Stripe 3: 90 deg (along +y)
  [./subdomain_id_90deg]
    type = OrientedSubdomainBoundingBoxGenerator
    input = subdomain_id_45deg
    center = '0 0 0'
    length = 0.005
    width = 6e-4
    height = 1
    length_direction = '0 1 0'
    width_direction = '-1 0 0'
    block_id = 1
  []
  # Stripe 4: 135 deg
  [./subdomain_id_135deg]
    type = OrientedSubdomainBoundingBoxGenerator
    input = subdomain_id_90deg
    center = '0 0 0'
    length = 0.005
    width = 6e-4
    height = 1
    length_direction = '-0.7071068 0.7071068 0'
    width_direction = '-0.7071068 -0.7071068 0'
    block_id = 1
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
  [diff]
    type = ADPFFDiffusion #
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
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc'
    prop_values = '${l} ${Gc_const}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
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
  # Dissipated energy density material (pure solid):
  # computes Psi^d_density = Gc/c0 * ( alpha/l + l * |grad d|^2 )
  [dissipated_energy]
    type = CrackDissipatedEnergyDensity
    phase_field = d
    alpha_name = alpha
    Gc_name = Gc
    l_name = l
    c0_name = c0
    # property_name defaults to 'dissipated_energy_density'
  []
[]

[AuxVariables]
  # Aux to expose dissipated energy density for integration/output
  [dissipated_energy_density]
    order = CONSTANT
    family = MONOMIAL
  []
[]

# This will transfer back to elasticity.i for integration
[AuxKernels]
  [dissipated_energy_density]
    type = ADMaterialRealAux
    variable = dissipated_energy_density
    property = dissipated_energy_density
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# Integrate dissipated energy here (fracture app) to avoid transfer noise
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
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  # petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero -snes_type'
  petsc_options_value = 'gmres     hypre  boomeramg True vinewtonrsls'

  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = false
  # time_step_interval = 40
  print_linear_residuals = false
[]
