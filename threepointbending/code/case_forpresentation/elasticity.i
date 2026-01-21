E = 40e6
nu = 0.25
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../../meshfile/mesh_whole_3d_x1mm_y1.6mm_forpresentation.msh'
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
[]

[Materials]
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [strain]
    type = ADComputeSmallStrain
  []
  #elastic domain
  [./elastic_stress]
    type = ADComputeLinearElasticStress
    output_properties = 'stress'
    outputs = exodus
  [../]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
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
  num_steps = 1

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
  show = 'd stress_00'
  print_linear_residuals = false
[]
