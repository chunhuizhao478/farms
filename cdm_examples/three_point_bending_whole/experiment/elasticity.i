E = 40e6
nu = 0.25
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'

l = 1e-5 # N * h, N: number of elements, h: element size

sigmat = 6.43e6

Gc = '${fparse 8*l*sigmat*sigmat/(3*E)}' 

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc=${Gc};l=${l}'
    execute_on = 'TIMESTEP_END'
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
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../meshfile_refined/mesh_1.6up_1.5shrift.msh'
  []
  [./elastic_region_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0.003 0 0'
    top_right = '0.005 0.001 1'
    block_id = 1
  []
  [./elastic_region_2]
    type = SubdomainBoundingBoxGenerator
    input = elastic_region_1
    bottom_left = '0.023 0 0'
    top_right = '0.025 0.001 0'
    block_id = 1
  []
  [./elastic_region_3]
    type = SubdomainBoundingBoxGenerator
    input = elastic_region_2
    bottom_left = '0.013 0.007 0'
    top_right = '0.015 0.008 1'
    block_id = 1
  []
  [./extranodeset_0]
    type = ExtraNodesetGenerator
    coord = '0.004 0 0'
    new_boundary = support1
    input = elastic_region_3
    use_closest_node=true
  []
  [./extranodeset_1]
    type = ExtraNodesetGenerator
    coord = '0.024 0 0'
    new_boundary = support2
    input = extranodeset_0
    use_closest_node=true
  []
  [./extranodeset_2]
    type = ExtraNodesetGenerator
    coord = '0.0140 0.008 0'
    new_boundary = load
    input = extranodeset_1
    use_closest_node=true
  []
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
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
[]

[Functions]
  [func_loading]
    type = ParsedFunction
    expression = '-1e-4 * t'
  []
[]

[BCs]
  # [fix_support1_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = support1
  #   value = 0
  # []
  [fix_support1_y]
    type = DirichletBC
    variable = disp_y
    boundary = support1
    value = 0
  []
  # [fix_support2_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = support2
  #   value = 0
  # []
  [fix_support2_y]
    type = DirichletBC
    variable = disp_y
    boundary = support2
    value = 0
  []
  [fix_load]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = load
    function = func_loading
  []
  [fix_load_x]
    type = DirichletBC
    variable = disp_x
    boundary = load
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
    block = 0
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
    block = 0
  []
  #elastic domain
  [./elastic_stress]
    type = ADComputeLinearElasticStress
    output_properties = 'stress'
    outputs = exodus
    block = 1
  [../]      
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
    block = 1
  [../]
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '
  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  dt = 0.001
  end_time = 100

  fixed_point_max_its = 20
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-8
  fixed_point_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  time_step_interval = 10
  print_linear_residuals = false
[]
