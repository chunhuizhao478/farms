E = 40e6
nu = 0.25
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'

l = 5e-5 # N * h, N: number of elements, h: element size

sigmat = 6.43e6

Gc = '${fparse 8*l*sigmat*sigmat/(3*E)}'

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc=${Gc};l=${l}'
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
[]

# Line support positions
left_support_x = 0.004
right_support_x = 0.024
loading_x = 0.014
z_center = 0.004  # extrude_z / 2

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../meshfile_angles/mesh_whole_3d_x0.50mm_y1.6mm.msh'
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
    bottom_left = '0.0139 0.007 0'
    top_right = '0.0141 0.008 1'
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
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
  []
  [solid_z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
  []
[]

[Functions]
  [func_loading]
    type = ParsedFunction
    expression = '-1e-4 * t'
  []
  [dt_fn]
    type = PiecewiseConstant
    x = '0   14'
    y = '0.1 0.001'
    direction = LEFT_INCLUSIVE
  []
[]

# LINE SUPPORT AND LINE LOADING BOUNDARY CONDITIONS
# Boundary names from mesh:
#   top_loading_line: loading line at top center
#   top_loading_center: center point of loading line (for u_x=0, u_z=0)
#   left_support_line: left support line at bottom
#   right_support_line: right support line at bottom
#   left_support_center: center point of left support (for u_y=0)
#   right_support_center: center point of right support (for u_y=0)
[BCs]
  #=== TOP LOADING LINE ===
  # Apply displacement loading on the entire line
  [apply_load_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top_loading_line
    function = func_loading
  []
  # u_x = 0 only at center point of loading line
  [fix_load_center_x]
    type = DirichletBC
    variable = disp_x
    boundary = top_loading_center
    value = 0
  []
  # u_z = 0 only at center point of loading line
  [fix_load_center_z]
    type = DirichletBC
    variable = disp_z
    boundary = top_loading_center
    value = 0
  []

  #=== LEFT SUPPORT LINE ===
  # u_z = 0 on entire line (prevent out-of-plane motion)
  [fix_left_support_line_z]
    type = DirichletBC
    variable = disp_z
    boundary = left_support_line
    value = 0
  []
  # u_y = 0 only at center point (vertical constraint at middle)
  [fix_left_support_center_y]
    type = DirichletBC
    variable = disp_y
    boundary = left_support_center
    value = 0
  []

  #=== RIGHT SUPPORT LINE ===
  # u_z = 0 on entire line (prevent out-of-plane motion)
  [fix_right_support_line_z]
    type = DirichletBC
    variable = disp_z
    boundary = right_support_line
    value = 0
  []
  # u_y = 0 only at center point (vertical constraint at middle)
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
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 1e-6'
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = spectral
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
    block = 1
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
    block = 1
  []
  #elastic domain
  [./elastic_stress]
    type = ADComputeLinearElasticStress
    output_properties = 'stress'
    outputs = exodus
    block = 2
  [../]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
    block = 2
  [../]
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  #petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  #petsc_options_value = 'lu       superlu_dist                 '

  petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor -snes_type -ksp_gmres_restart'
  petsc_options_value = 'gmres hypre boomeramg 0.7 4 5 0.3 vinewtonrsls 100'

  line_search = basic

  automatic_scaling = true

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_max_its = 50
  end_time = 20

  [TimeStepper]
    type = FunctionDT
    function = dt_fn
  []

  fixed_point_max_its = 20
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8
[]

[Outputs]
  [./exodus]
    type = Exodus
    sync_times = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 14.1 14.2 14.3 14.4 14.5 14.6 14.7 14.8 14.9 15.0 15.1 15.2 15.3 15.4 15.5 15.6 15.7 15.8 15.9 16.0 16.1 16.2 16.3 16.4 16.5 16.6 16.7 16.8 16.9 17.0 17.1 17.2 17.3 17.4 17.5 17.6 17.7 17.8 17.9 18.0 18.1 18.2 18.3 18.4 18.5 18.6 18.7 18.8 18.9 19.0 19.1 19.2 19.3 19.4 19.5 19.6 19.7 19.8 19.9 20.0'
    sync_only = true
    show = 'd disp_x disp_y'
  [../]
  print_linear_residuals = false
[]
