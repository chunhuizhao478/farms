[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/square_with_hole.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.01 0.01 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
[]

# d is now a Variable that evolves (damage evolution activated)
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

# Initial conditions for damage field
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

# AuxKernels for initial_damage_aux (used for bounds)
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

# Bounds for damage irreversibility: d can only increase
[Bounds]
  [irreversibility]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  []
  [upper]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1.0
  []
[]

# Phase-field damage evolution kernels
[Kernels]
  [diff]
    type = ADPFFDiffusion
    variable = d
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
    expression = 'alpha*Gc/c0/l+0.5*g*psie_active' #<-there is a 0.5 because 2 is taken within psie_active
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
  # Use vinewtonrsls for bounded variational inequality (damage irreversibility)
  petsc_options_iname = '-snes_type -ksp_type -pc_type -pc_hypre_type'
  petsc_options_value = 'vinewtonrsls gmres hypre boomeramg'

  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = false
  # time_step_interval = 40
  print_linear_residuals = false
[]
