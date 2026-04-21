[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 160
    ny = 160
    xmin = 0
    xmax = 40
    ymin = 0
    ymax = 40
    elem_type = QUAD4
  []
  [./damage_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 20 0'
    top_right = '4.0 20.25 0'
  [../]
  [./sidesets]
    input = damage_block
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0'
    new_boundary = 'left right bottom top'
  []
[]

# [Adaptivity]
#   max_h_level = 5
#   marker = 'combo'
#   cycles_per_step = 1
#   [Markers]
#       [./combo]
#         type = FarmsComboMarker
#         markers = 'damage_marker strain_energy_marker'
#         meshsize_marker = 'meshsize_marker'
#       [../]
#       [damage_marker]
#         type = ValueThresholdMarker
#         variable = d
#         refine = 0.01
#       []
#       [strain_energy_marker]
#         type = ValueThresholdMarker
#         variable = psie_active
#         refine = '${fparse 1.0*3/8*Gc_const/l}'
#       []   
#       # if mesh_size > dxmin, refine
#       # if mesh_size < dxmin/100, coarsen (which never happens)
#       # otherwise, do nothing
#       [meshsize_marker]
#         type = ValueThresholdMarker
#         variable = mesh_size
#         refine = '${dx_min}'
#         coarsen = '${fparse dx_min/100}'
#         third_state = DO_NOTHING
#       [] 
#   []
# []

[Variables]
  [d]
    family = LAGRANGE
    order = FIRST
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
    value = 1
    block = 1
  []
  [define_initial_damage_block0]
    type = ConstantAux
    variable = initial_damage_aux
    value = 0
    block = 0
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
  exodus = true
  # time_step_interval = 40
  print_linear_residuals = false
[]