[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/mesh.msh'
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
  [./extranodeset_1]
    type = ExtraNodesetGenerator
    coord = '0.004 0 0;
             0.024 0 0'
    new_boundary = support
    input = elastic_region_3
    use_closest_node=true
  []
  [./extranodeset_2]
    type = ExtraNodesetGenerator
    coord = '0.014 0.008 0'
    new_boundary = load
    input = extranodeset_1
    use_closest_node=true
  []
  displacements = 'disp_x disp_y'
[]

  [Variables]
    [disp_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE
    []
  []

  [AuxVariables]
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
    [./strength]
      order = CONSTANT
      family = MONOMIAL
    [../]
  []
  
  [Functions]
    [func_loading]
      type = ParsedFunction
      expression = '-0.0001 * t'
    []
  []

  [Kernels]
    [dispkernel_x]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
        use_displaced_mesh = true
    []
    [dispkernel_y]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
        use_displaced_mesh = true
    []
  []
  
  [BCs]
    # [fix_support_x]
    #   type = ADDirichletBC
    #   variable = disp_x
    #   boundary = support
    #   value = 0
    # []
    [fix_support_y]
      type = ADDirichletBC
      variable = disp_y
      boundary = support
      value = 0
    []
    [fix_load]
      type = ADFunctionDirichletBC
      variable = disp_y
      boundary = load
      function = func_loading
    []
    [fix_load_x]
      type = ADDirichletBC
      variable = disp_x
      boundary = load
      value = 0
    []
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ADComputeIsotropicElasticityTensor
      youngs_modulus = 1.5e9
      poissons_ratio = 0.25
    [../]
    [./elastic_stress]
      type = ADComputeSmearedCrackingStress
      cracking_stress = 6.43e6
      softening_models = abrupt_softening
      cracked_elasticity_type = FULL
      max_cracks = 1
      outputs = exodus
      block = '0'
    [../]
    [./elastic_stress_2]
      type = ADComputeLinearElasticStress
      block = '1'
      outputs = exodus
    [../]
    [strain]
        type = ADComputeFiniteStrain
        displacements = 'disp_x disp_y'
    []
    [./abrupt_softening]
      type = ADAbruptSoftening
    [../]
    [./exponential_softening]
      type = ADExponentialSoftening
    [../]
  []

  [Preconditioning]
    [smp]
      type = SMP
      full = true
    []
 []
  
  [Executioner]
    type = Transient
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    line_search = 'bt'
    automatic_scaling = true
    # num_steps = 1
    l_max_its = 100
    nl_max_its = 10
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    start_time = 0.0
    end_time = 100
    # dt = 0.005
    [./TimeIntegrator]
        type = ImplicitEuler
    [../]
    [TimeStepper]
      type = IterationAdaptiveDT
      optimal_iterations = 5
      linear_iteration_ratio = 20
      dt = 0.005
    []
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []