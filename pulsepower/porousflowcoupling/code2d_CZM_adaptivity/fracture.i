bottom_left = '-0.0025 -2e-4 0'
top_right = '0.0025 2e-4 0'

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../2dmeshfile/fieldscale_test1_2d.msh'
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0.1 0.1 0'
    new_boundary = corner_ptr
    input = msh
    use_closest_node=true
  []
  [./subdomain_id]
    type = SubdomainBoundingBoxGenerator
    bottom_left = ${bottom_left}
    top_right = ${top_right}
    location = INSIDE
    block_id = 1
    input = extranodeset1
  []
[]

# [Adaptivity]
#   initial_marker = initial_marker
#   initial_steps = 2
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
#         refine = 0.25
#       []
#       [strain_energy_marker]
#         type = ValueThresholdMarker
#         variable = psie_active
#         refine = '${fparse 1.0*3/8*Gc/l}'
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
#       [initial_marker]
#         type = BoxMarker
#         bottom_left = ${bottom_left}
#         top_right = ${top_right}
#         outside = DO_NOTHING
#         inside = REFINE
#       []
#   []
# []

[Variables]
  [d]
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
  [a1_aux]
    family = MONOMIAL
    order = FIRST
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
  # [irreversibility]
  #   type = VariableOldValueBounds
  #   variable = bounds_dummy
  #   bounded_variable = d
  #   bound_type = lower
  # []
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
    normalization_constant = c_alpha
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
    prop_names =  'l Gc a1 a2 a3 p eta c_alpha'
    prop_values = '${l} ${Gc} ${a1} ${a2} ${a3} ${p} ${eta} ${c_alpha}'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = '2*d - d*d'
    phase_field = d
  []
  [degradation]
    type = RationalDegradationFunctionCZM
    property_name = g
    expression = (1-d)^p/((1-d)^p+a1*d*(1+a2*d+a2*a3*d^2))*(1-eta)+eta
    phase_field = d
    material_property_names = 'a1 a2 a3 p eta'
  []
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c_alpha/l+g*psie_active'
    coupled_variables = 'd psie_active'
    material_property_names = 'alpha(d) g(d) Gc c_alpha l'
    derivative_order = 1
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'

  # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero -snes_type'
  # petsc_options_value = 'gmres     hypre  boomeramg True vinewtonrsls'

  # automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = false
  # time_step_interval = 40
  print_linear_residuals = false
[]