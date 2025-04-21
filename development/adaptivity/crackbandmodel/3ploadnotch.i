[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/meshwohole.msh'
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
    [./strength]
      order = CONSTANT
      family = MONOMIAL
      initial_condition = 8e6
    [../]
  []
  
  [Functions]
    [func_loading]
      type = ParsedFunction
      expression = '-1e-4 * t'
    []
  []

  [GlobalParams]
    displacements = 'disp_x disp_y'
  []
  
  [BCs]
    [fix_right_x]
      type = DirichletBC
      variable = disp_x
      boundary = right
      value = 0
    []
    [fix_right_y]
      type = DirichletBC
      variable = disp_y
      boundary = right
      value = 0
    []
    [fix_load_x]
      type = FunctionDirichletBC
      variable = disp_x
      boundary = left
      function = func_loading
    []
    [fix_left_y]
      type = DirichletBC
      variable = disp_y
      boundary = left
      value = 0
    []
  []

  [Kernels]
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
        use_displaced_mesh = true
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
        use_displaced_mesh = true
    []
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 20e9
      poissons_ratio = 0.2
    [../]
    [./damage_stress]
        type = FarmsCrackBandModel
        cracking_stress = strength
        hb = 4e-4 #m
        Gf = 90 #J/m^2
        output_properties = 'crack_damage stress elastic_strain crack_initiation_strain'
        outputs = exodus
    [../]    
    [strain]
      type = ComputeFiniteStrain
      displacements = 'disp_x disp_y'
    []
  []

  # [Adaptivity]
  #   max_h_level = 3
  #   marker = 'combo'
  #   # initial_steps = 2
  #   # initial_marker = initial_box
  #   [Indicators]
  #       [error]
  #         type = GradientJumpIndicator
  #           variable = crack_damage_aux
  #       []
  #   []
  #   [Markers]
  #       [./combo]
  #           type = ComboMarker
  #           markers = 'error_marker nonlocal_eqstrain_marker'
  #       [../]
  #       [./error_marker]
  #           type = ErrorFractionMarker
  #           indicator = error
  #           refine = 0.9
  #       [../]
  #       [./nonlocal_eqstrain_marker]
  #           type = ValueThresholdMarker
  #           variable = eqstrain_local
  #           refine = 5e-5
  #       []         
  #   []
  # []

  [Preconditioning]
    [smp]
      type = SMP
      full = true
    []
 []
  
  [Executioner]
    type = Transient
    solve_type = Newton
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    line_search = 'none'
    # num_steps = 1
    l_max_its = 200
    nl_max_its = 100
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    start_time = 0.0
    end_time = 100
    dt = 0.0025
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []