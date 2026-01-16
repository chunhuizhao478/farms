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

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../meshfile/mesh_whole_3d_x2mm_y1.6mm.msh'
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
[]

#top_loading: id 4
#bottom_left_support: id 2
#bottom_right_support: id 3
[BCs]
  #top_loading
  [apply_load_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 4
    function = func_loading
  []
  [fix_load_x]
    type = DirichletBC
    variable = disp_x
    boundary = 4
    value = 0
  []
  [fix_load_z]
    type = DirichletBC
    variable = disp_z
    boundary = 4
    value = 0
  []
  #bottom_left_support
  [fix_bottom_left_support_y]
    type = DirichletBC
    variable = disp_y
    boundary = 2
    value = 0
  []
  [fix_bottom_left_support_z]
    type = DirichletBC
    variable = disp_z
    boundary = 2
    value = 0
  []
  #bottom_right_support
  [fix_bottom_right_support_y]
    type = DirichletBC
    variable = disp_y
    boundary = 3
    value = 0
  []
  [fix_bottom_right_support_z]
    type = DirichletBC
    variable = disp_z
    boundary = 3
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

  dt = 0.1
  end_time = 30

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    optimal_iterations = 6
    iteration_window = 2
    growth_factor = 1.2
    cutback_factor = 0.5
  []

  fixed_point_max_its = 20
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8
[]

[Outputs]
  exodus = true
  time_step_interval = 10
  print_linear_residuals = false
[]
