# =============================================================================
# SPE-HFTC acrylic cube, AT1 phase-field sub-app (3D)
# Launched by elasticity_absorbing.i (cli_args pass Gc_const, l, jit, mesh_file); clone_parent_mesh = true.
# =============================================================================

# defaults for standalone parsing; overridden by the main app's cli_args
l = 6e-4
Gc_const = 2720.5714285714284
mesh_file = ../mesh/spe_cube_h0p3.msh
Lx = 51.11e-3
Ly = 51.92e-3

[Mesh]
  [msh]
    type = FileMeshGenerator
    file = ${mesh_file}
  []
  [corner]
    type = ExtraNodesetGenerator
    input = msh
    coord = '${fparse -Lx/2} ${fparse -Ly/2} 0'
    new_boundary = corner_ptr
    use_closest_node = true
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
  [dissipated_energy_density]
    order = CONSTANT
    family = MONOMIAL
  []
  # toughness field (constant here; the *_weibull.i / *_gaussian.i variants draw it at random);
  # copied to the main app for output
  [Gc_var]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [gc_var]
    type = ConstantIC
    variable = Gc_var
    value = ${Gc_const}
  []
[]

[AuxKernels]
  # homogeneous acrylic: no pre-seeded damage (lower bound of d is 0 everywhere)
  [define_initial_damage]
    type = ConstantAux
    variable = initial_damage_aux
    value = 0
  []
  # exposes the crack energy density for the energy budget (copied back to the main app)
  [dissipated_energy_density]
    type = ADMaterialRealAux
    variable = dissipated_energy_density
    property = dissipated_energy_density
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Bounds]
  [irreversibility]
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
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l'
    prop_values = '${l}'
  []
  [Gc_field]
    type = ADParsedMaterial
    property_name = Gc
    coupled_variables = 'Gc_var'
    expression = 'Gc_var'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 1e-6'
  []
  [crack_geometric]
    # AT1: alpha(d) = d, c0 = 8/3
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
  # crack energy density Gc/c0 * (alpha/l + l |grad d|^2)
  [dissipated_energy]
    type = CrackDissipatedEnergyDensity
    phase_field = d
    alpha_name = alpha
    Gc_name = Gc
    l_name = l
    c0_name = c0
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

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor -snes_type -ksp_gmres_restart'
  petsc_options_value = 'gmres hypre boomeramg 0.7 4 5 0.3 vinewtonrsls 100'

  automatic_scaling = true
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = false
  print_linear_residuals = false
  [perf]
    # the sub-app's own wall-time breakdown, same cadence as the main app's
    type = PerfGraphOutput
    execute_on = 'TIMESTEP_END FINAL'
    time_step_interval = 500
    level = 2
  []
[]
