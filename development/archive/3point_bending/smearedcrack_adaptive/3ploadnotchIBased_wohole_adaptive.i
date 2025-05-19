[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/meshwohole.msh'
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
    [./strength]
      order = CONSTANT
      family = MONOMIAL
      initial_condition = 2.4e6
    [../]
    [crack_damage_aux]
      order = FIRST
      family = MONOMIAL
    []
    [eqstrain_nonlocal_aux]
      order = FIRST
      family = MONOMIAL
    []  
  []

  [AuxKernels]
    [get_damage]
      type = MaterialRealAux
      variable = crack_damage_aux
      property = crack_damage
      block = '0'
    []
    [get_eqstrain_nonlocal_aux]
      type = MaterialRealAux  
      variable = eqstrain_nonlocal_aux
      property = eqstrain_nonlocal
      execute_on = timestep_end
      block = 0
    []
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

  [Physics/SolidMechanics/QuasiStatic]
    [./all]
      strain = FINITE
      add_variables = true
      generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_zx'
    [../]
  []

  [UserObjects]
    [eqstrain_averaging]
        type = ElkRadialAverage
        length_scale = 5e-5
        prop_name = eqstrain_local
        radius = 5e-5
        weights = BAZANT
        execute_on = LINEAR
        block = 0
    []
  []

  [Adaptivity]
    max_h_level = 3
    marker = 'combo'
    # initial_steps = 2
    # initial_marker = initial_box
    [Indicators]
        [error]
          type = GradientJumpIndicator
            variable = crack_damage_aux
        []
    []
    [Markers]
        [./combo]
            type = ComboMarker
            markers = 'error_marker nonlocal_eqstrain_marker'
        [../]
        [./error_marker]
            type = ErrorFractionMarker
            indicator = error
            refine = 0.9
        [../]
        [./nonlocal_eqstrain_marker]
            type = ValueThresholdMarker
            variable = eqstrain_local
            refine = 5e-5
        []         
    []
  []
  
  [BCs]
    [fix_support1_x]
      type = DirichletBC
      variable = disp_x
      boundary = support1
      value = 0
    []
    [fix_support1_y]
      type = DirichletBC
      variable = disp_y
      boundary = support1
      value = 0
    []
    [fix_support2_x]
      type = DirichletBC
      variable = disp_x
      boundary = support2
      value = 0
    []
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
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 20e9
      poissons_ratio = 0.2
    [../]
    [./elastic_stress]
        type = FarmsComputeSmearedCrackingStress
        paramA = 0.99
        paramB = 500
        cracking_stress = strength
        output_properties = 'crack_damage stress eqstrain_local'
        outputs = exodus
        block = 0
    [../]
    [nonlocal_eqstrain]
      type = ElkNonlocalEqstrain
      average_UO = eqstrain_averaging
      output_properties = 'eqstrain_nonlocal'
      outputs = exodus
      block = 0
    []
    [./nondamage_stress]
      type = ComputeFiniteStrainElasticStress
      output_properties = 'stress'
      outputs = exodus
      block = 1
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
    solve_type = PJFNK #in smeared cracking w/o full jacobian, this is much more efficient than NEWTON
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    line_search = 'bt'
    # num_steps = 1
    l_max_its = 200
    nl_max_its = 10
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    start_time = 0.0
    end_time = 100
    dt = 0.001
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []