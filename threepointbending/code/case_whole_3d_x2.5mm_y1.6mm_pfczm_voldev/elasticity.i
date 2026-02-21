E = 4.41e6
nu = 0.25
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'

l = 5e-5 # N * h, N: number of elements, h: element size

sigmat = 6.43e6

Gc = 3 #J/m^2
psic = '${fparse sigmat*sigmat/(2*E)}'

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc=${Gc};psic=${psic};l=${l}'
    execute_on = 'TIMESTEP_END'
    clone_parent_mesh = true
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = d
    source_variable = d
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = psie_active
    source_variable = psie_active
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

# Line support positions
left_support_x = 0.004
right_support_x = 0.024
loading_x = 0.014
z_center = 0.004  # extrude_z / 2

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../meshfile/mesh_whole_3d_x2.5mm_y1.6mm.msh'
  []
  [./elastic_region_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0.0035 0 0'
    top_right = '0.0045 0.001 1'
    block_id = 2
  []
  [./elastic_region_2]
    type = SubdomainBoundingBoxGenerator
    input = elastic_region_1
    bottom_left = '0.0235 0 0'
    top_right = '0.0245 0.001 1'
    block_id = 2
  []
  [./elastic_region_3]
    type = SubdomainBoundingBoxGenerator
    input = elastic_region_2
    bottom_left = '0.0135 0.007 0'
    top_right = '0.0145 0.008 1'
    block_id = 2
  []
  # Create nodesets for LINE supports using BoundingBoxNodeSetGenerator
  # Left support LINE: x = 0.004, y = 0, z = 0 to 0.008
  [./left_support_line]
    type = BoundingBoxNodeSetGenerator
    input = elastic_region_3
    bottom_left = '${fparse left_support_x - 1e-5} -1e-5 -1e-5'
    top_right = '${fparse left_support_x + 1e-5} 1e-5 0.009'
    new_boundary = left_support_line
  []
  # Right support LINE: x = 0.024, y = 0, z = 0 to 0.008
  [./right_support_line]
    type = BoundingBoxNodeSetGenerator
    input = left_support_line
    bottom_left = '${fparse right_support_x - 1e-5} -1e-5 -1e-5'
    top_right = '${fparse right_support_x + 1e-5} 1e-5 0.009'
    new_boundary = right_support_line
  []
  # Top loading LINE: x = 0.014, y = 0.008, z = 0 to 0.008
  [./top_loading_line]
    type = BoundingBoxNodeSetGenerator
    input = right_support_line
    bottom_left = '${fparse loading_x - 1e-5} 0.0079 -1e-5'
    top_right = '${fparse loading_x + 1e-5} 0.0081 0.009'
    new_boundary = top_loading_line
  []
  # Create nodesets for CENTER POINTS of support lines (z = z_center)
  # Left support center point
  [./left_support_center]
    type = ExtraNodesetGenerator
    input = top_loading_line
    coord = '${left_support_x} 0 ${z_center}'
    new_boundary = left_support_center
    use_closest_node = true
  []
  # Right support center point
  [./right_support_center]
    type = ExtraNodesetGenerator
    input = left_support_center
    coord = '${right_support_x} 0 ${z_center}'
    new_boundary = right_support_center
    use_closest_node = true
  []
  # Top loading center point
  [./top_loading_center]
    type = ExtraNodesetGenerator
    input = right_support_center
    coord = '${loading_x} 0.008 ${z_center}'
    new_boundary = top_loading_center
    use_closest_node = true
  []
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[AuxVariables]
  [d]
  []
[]

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
  []
  [solid_z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true
  []
[]

[Functions]
  [func_loading]
    type = ParsedFunction
    expression = '-1.667e-6 * t'
  []
  [dt_limit_fn]
    type = PiecewiseLinear
    x = '0     16.5   20'
    y = '1e10  0.001 0.001'
  []
[]

[Postprocessors]
  [dt_limit_pp]
    type = FunctionValuePostprocessor
    function = dt_limit_fn
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

# LINE SUPPORT AND LINE LOADING BOUNDARY CONDITIONS
[BCs]
  #=== TOP LOADING LINE ===
  [apply_load_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top_loading_line
    function = func_loading
  []
  [fix_load_center_x]
    type = DirichletBC
    variable = disp_x
    boundary = top_loading_center
    value = 0
  []
  [fix_load_center_z]
    type = DirichletBC
    variable = disp_z
    boundary = top_loading_center
    value = 0
  []

  #=== LEFT SUPPORT LINE ===
  [fix_left_support_line_z]
    type = DirichletBC
    variable = disp_z
    boundary = left_support_line
    value = 0
  []
  [fix_left_support_center_y]
    type = DirichletBC
    variable = disp_y
    boundary = left_support_center
    value = 0
  []

  #=== RIGHT SUPPORT LINE ===
  [fix_right_support_line_z]
    type = DirichletBC
    variable = disp_z
    boundary = right_support_line
    value = 0
  []
  [fix_right_support_center_y]
    type = DirichletBC
    variable = disp_y
    boundary = right_support_center
    value = 0
  []
[]

[Materials]
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G l Gc psic'
    prop_values = '${K} ${G} ${l} ${Gc} ${psic}'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
  []
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+(Gc/psic*xi/c0/l)*d*(1+a2*d+a2*a3*d^2))*(1-eta)+eta
    phase_field = d
    material_property_names = 'Gc psic xi c0 l '
    parameter_names = 'p a2 a3 eta '
    parameter_values = '2 -0.5 0 1e-6'
  []
  # Block 1: damage zone - large deformation
  [defgrad]
    type = ComputeDeformationGradient
    block = 1
  []
  [hencky]
    type = HenckyIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = VOLDEV
    output_properties = 'psie_active'
    outputs = exodus
    block = 1
  []
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky
    output_properties = 'stress'
    outputs = exodus
    block = 1
  []
  # Block 2: elastic zone - large deformation, no degradation
  [defgrad_elastic]
    type = ComputeDeformationGradient
    block = 2
  []
  [no_degradation]
    type = NoDegradation
    property_name = g_nodeg
    expression = 1
    phase_field = d
    block = 2
  []
  [hencky_elastic]
    type = HenckyIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g_nodeg
    block = 2
  []
  [stress_elastic]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky_elastic
    block = 2
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor -snes_type -ksp_gmres_restart'
  petsc_options_value = 'gmres hypre boomeramg 0.7 4 5 0.3 vinewtonrsls 100'

  line_search = basic

  automatic_scaling = true

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 50

  dt = 1
  end_time = 1000

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    optimal_iterations = 6
    iteration_window = 2
    growth_factor = 1.2
    cutback_factor = 0.5
    timestep_limiting_postprocessor = dt_limit_pp
  []

  fixed_point_max_its = 20
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8
[]

[Outputs]
  [./exodus]
    type = Exodus
    interval = 10
    show = 'd disp_x disp_y'
  [../]
  print_linear_residuals = false
[]
